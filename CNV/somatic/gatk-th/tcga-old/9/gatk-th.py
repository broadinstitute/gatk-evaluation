#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd

import jax
import jax.numpy as jnp
import jax.scipy
import jax.ops

import numba
from numba import jit

import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from sklearn.utils.extmath import cartesian

import scipy.stats
import scipy.special

from collections import namedtuple

import tqdm

import emcee
import corner

#===============================================================================
# posterior fitting
#===============================================================================

eps_posterior_fitting = 1E-10

def fit_t_to_log2cr(modeled_segment, t_degrees_of_freedom):
    log2cr_10 = modeled_segment['LOG2_COPY_RATIO_POSTERIOR_10']
    log2cr_50 = modeled_segment['LOG2_COPY_RATIO_POSTERIOR_50']
    log2cr_90 = modeled_segment['LOG2_COPY_RATIO_POSTERIOR_90']
    
    if np.isnan(log2cr_50):
        return np.array([np.nan, np.nan])
    
    def loss(x):
        mu, sigma = x
        t_cdf = lambda x: scipy.stats.t._cdf((x - mu) / sigma, df=t_degrees_of_freedom) # use _cdf to avoid scipy.stats overhead
        return (0.1 - t_cdf(log2cr_10))**2 + \
               (0.5 - t_cdf(log2cr_50))**2 + \
               (0.9 - t_cdf(log2cr_90))**2
    
    mu, sigma = scipy.optimize.minimize(loss, x0=[log2cr_50, log2cr_90 - log2cr_10],
                                        method='L-BFGS-B',
                                        bounds=((None, None), (eps_posterior_fitting, None))).x
#     print(modeled_segment, mu, sigma)
    return np.array([mu, sigma])

def fit_beta_to_maf(modeled_segment):
    maf_10 = modeled_segment['MINOR_ALLELE_FRACTION_POSTERIOR_10']
    maf_50 = modeled_segment['MINOR_ALLELE_FRACTION_POSTERIOR_50']
    maf_90 = modeled_segment['MINOR_ALLELE_FRACTION_POSTERIOR_90']
    
    if np.isnan(maf_50):
        return np.array([np.nan, np.nan])
    
    def loss(x):
        a, b = x
        scaled_beta_cdf = lambda x: scipy.stats.beta._cdf(x / 0.5, a=a, b=b) # use _cdf to avoid scipy.stats overhead
        return (0.1 - scaled_beta_cdf(maf_10))**2 + \
               (0.5 - scaled_beta_cdf(maf_50))**2 + \
               (0.9 - scaled_beta_cdf(maf_90))**2
    
    a, b = scipy.optimize.minimize(loss, x0=[maf_50, maf_90 - maf_10],
                                   method='L-BFGS-B',
                                   bounds=((eps_posterior_fitting, None), (eps_posterior_fitting, None))).x
#     print(modeled_segment, a, b, min([max([0.5 * (a - 1) / (a + b - 2), 0.]), 0.5]))
    return np.array([a, b])

#===============================================================================
# logpdf utilities
#===============================================================================

# hypercube-to-simplex transform: https://arxiv.org/pdf/1010.3436.pdf

@numba.jit(nopython=True)
def to_simplex(z):
    m = z.size + 1
    x = np.zeros(m)
    x[0] = 1. - z[0]
    cum_prod = 1.
    for k in range(1, m - 1):
        cum_prod *= z[k - 1]
        x[k] = cum_prod * (1. - z[k]) 
    x[m - 1] = cum_prod * z[m - 2]
    return x

@numba.jit(nopython=True)
def to_hypercube(x):
    m = x.size
    z = np.zeros(m - 1)
    z[0] = 1. - x[0]
    cum_prod = 1.
    for k in range(1, m - 1):
        cum_prod *= z[k - 1]
        z[k] = 1. - x[k] / cum_prod
    return z
  
_norm_pdf_C = np.sqrt(2. * np.pi)
_norm_pdf_logC = np.log(_norm_pdf_C)
def norm_logpdf(x, loc=0., scale=1.):
    y = (x - loc) / scale
    return -y**2 / 2.0 - _norm_pdf_logC
  
def t_logpdf(x, df, loc=0., scale=1.):
    r = df * 1.
    y = (x - loc) / scale
    lPx = scipy.special.gammaln((r + 1) / 2) - scipy.special.gammaln(r / 2)
    lPx -= 0.5 * np.log(r * np.pi) + (r + 1) / 2 * np.log(1 + y**2 / r)
    return lPx
  
def t_logpdf_jax(x, df, loc=0., scale=1.):
    r = df * 1.
    y = (x - loc) / scale
    lPx = jax.scipy.special.gammaln((r + 1) / 2) - jax.scipy.special.gammaln(r / 2)
    lPx -= 0.5 * jnp.log(r * jnp.pi) + (r + 1) / 2 * jnp.log(1 + y**2 / r)
    return lPx
t_logpdf_jax_jit = jax.jit(t_logpdf_jax)
  
eps_beta_logpdf = 1E-3
def beta_logpdf(x, a, b, scale=1.):
    # up to factor independent of x
    x_bnd = np.minimum(scale - eps_beta_logpdf, np.maximum(eps_beta_logpdf, x)) / scale
    lPx = scipy.special.xlog1py(b - 1.0, -x_bnd) + scipy.special.xlogy(a - 1.0, x_bnd)
#     lPx -= scipy.special.betaln(a, b) + np.log(scale)
    return lPx

# slower than scipy implementation
# def beta_logpdf_jax(x, a, b, scale=1.):
#     # up to factor independent of x
#     x_bnd = jnp.minimum(scale - eps_beta_logpdf, jnp.maximum(eps_beta_logpdf, x)) / scale
#     lPx = jax.scipy.special.xlog1py(b - 1.0, -x_bnd) + jax.scipy.special.xlogy(a - 1.0, x_bnd)
#     return lPx
# beta_logpdf_jax_jit = jax.jit(beta_logpdf_jax)

def dirichlet_hypercube_logp(alpha, z):
    # up to factor involving Gammas independent of z
    alpha_tilde = np.sum(alpha) - np.cumsum(alpha)
    return np.sum(beta_logpdf(z, a=alpha_tilde[:-1], b=alpha[:-1]))
  
#===============================================================================
# discrete priors
#===============================================================================

GlobalDiscretePrior = namedtuple('GlobalDiscretePrior', ['allelic_copy_number_product_states_lij',
                                                         'unnorm_product_state_log_prior_li',
                                                         'product_state_log_prior_l',
                                                         'marginalization_allelic_copy_number_product_states_lij',
                                                         'marginalization_product_state_log_prior_l',
                                                         'marginalization_product_state_filter_l'])
  
def generate_label_ordered_product_states_and_log_prior(discrete_prior_config):
    num_overlapping_populations = discrete_prior_config.num_overlapping_populations
    num_alleles = discrete_prior_config.num_alleles
    allelic_copy_number_states = discrete_prior_config.allelic_copy_number_states
    normal_allelic_copy_number_state = discrete_prior_config.normal_allelic_copy_number_state
    allelic_copy_number_change_prior_prob = discrete_prior_config.allelic_copy_number_change_prior_prob
    hom_del_prior_prob = discrete_prior_config.hom_del_prior_prob

    num_allelic_copy_number_states = len(allelic_copy_number_states)

    # number of unordered product states = num_allelic_copy_number_states^(num_overlapping_populations * num_alleles)
    # each product state is num_overlapping_populations x num_alleles matrix of allelic copy number
    allelic_copy_number_product_states_lij = cartesian(np.tile(np.tile(allelic_copy_number_states, (num_alleles, 1)), (num_overlapping_populations, 1))).reshape(-1, num_overlapping_populations, num_alleles)
    print('Number of unordered product states:', len(allelic_copy_number_product_states_lij))

    # eliminate label switching between alleles by requiring the first allele to be the major allele in the population with the lowest index
    # out of those populations where the alleles are not identical (we assume that 1 <= num_alleles <= 2)
    def calculate_is_label_ordered(allelic_copy_number_product_state_ij):
        for i in range(num_overlapping_populations):
            if allelic_copy_number_product_state_ij[i][0] > allelic_copy_number_product_state_ij[i][1]:
                return True
            elif allelic_copy_number_product_state_ij[i][0] < allelic_copy_number_product_state_ij[i][1]:
                return False
        return True  # all alleles have equal copy number in all populations

    is_label_ordered_l = np.array([calculate_is_label_ordered(allelic_copy_number_product_state_ij) for allelic_copy_number_product_state_ij in allelic_copy_number_product_states_lij])
    label_ordered_allelic_copy_number_product_states_lij = allelic_copy_number_product_states_lij[is_label_ordered_l]
    print('Number of label-ordered product states:', len(label_ordered_allelic_copy_number_product_states_lij))

    unnorm_label_ordered_product_state_prior_li = (allelic_copy_number_change_prior_prob)**(np.array([np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 0, :] - normal_allelic_copy_number_state), axis=-1),
                                                                                                      np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 0, :] - label_ordered_allelic_copy_number_product_states_lij[:, 1, :]), axis=-1), 
                                                                                                      np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] - label_ordered_allelic_copy_number_product_states_lij[:, 2, :]), axis=-1)])).transpose()

    # remove states where clonal and subclonal populations have identical event (clonal states should be indicated by a zero subclonal cancer cell fraction with a normal subclone)
#     unnorm_label_ordered_product_state_prior_li[(np.all(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] == label_ordered_allelic_copy_number_product_states_lij[:, 2, :], axis=1)) * 
#                                                 (np.any(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] != normal_allelic_copy_number_state, axis=1))] = 0.

    # remove unphysical states where normal (clonal) deletions are reverted in clonal/subclonal (subclonal) population (except for clonal states)
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 1, :] != 0), axis=1)] = 0.
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 2, :] != 0), axis=1)] = 0.
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 1, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 2, :] != 0), axis=1) * 
                                                np.any(label_ordered_allelic_copy_number_product_states_lij[:, 2, :] != normal_allelic_copy_number_state, axis=1)] = 0.

    # heavily penalize hom dels
    unnorm_label_ordered_product_state_prior_li[np.all(label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == 0, axis=-1)] = hom_del_prior_prob
    unnorm_label_ordered_product_state_prior_li[np.all(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] == 0, axis=-1)] = hom_del_prior_prob
    unnorm_label_ordered_product_state_prior_li[np.all(label_ordered_allelic_copy_number_product_states_lij[:, 2, :] == 0, axis=-1)] = hom_del_prior_prob

    is_prior_nonzero_l = np.all(unnorm_label_ordered_product_state_prior_li > 0., axis=-1)
    label_ordered_allelic_copy_number_product_states_lij = label_ordered_allelic_copy_number_product_states_lij[is_prior_nonzero_l]
    unnorm_label_ordered_product_state_prior_li = unnorm_label_ordered_product_state_prior_li[is_prior_nonzero_l]

    unnorm_label_ordered_product_state_log_prior_li = np.log(unnorm_label_ordered_product_state_prior_li / (np.sum(np.product(unnorm_label_ordered_product_state_prior_li, axis=-1))))
    label_ordered_product_state_log_prior_l = np.sum(unnorm_label_ordered_product_state_log_prior_li, axis=-1)
    
    marginalization_product_state_log_prior_threshold = np.sort(label_ordered_product_state_log_prior_l)[::-1][np.minimum(discrete_prior_config.num_marginalization_product_states - 1, len(label_ordered_product_state_log_prior_l) - 1)]
    marginalization_product_state_filter_l = label_ordered_product_state_log_prior_l >= marginalization_product_state_log_prior_threshold
    marginalization_allelic_copy_number_product_states_lij = label_ordered_allelic_copy_number_product_states_lij[marginalization_product_state_filter_l]
    marginalization_product_state_log_prior_l = label_ordered_product_state_log_prior_l[marginalization_product_state_filter_l]
    print('Number of marginalization label-ordered product states:', len(marginalization_allelic_copy_number_product_states_lij))

    return GlobalDiscretePrior(
        allelic_copy_number_product_states_lij = label_ordered_allelic_copy_number_product_states_lij, 
        unnorm_product_state_log_prior_li = unnorm_label_ordered_product_state_log_prior_li, 
        product_state_log_prior_l = label_ordered_product_state_log_prior_l,
        marginalization_allelic_copy_number_product_states_lij = marginalization_allelic_copy_number_product_states_lij,
        marginalization_product_state_log_prior_l = marginalization_product_state_log_prior_l,
        marginalization_product_state_filter_l = marginalization_product_state_filter_l)
  
Data = namedtuple('Data', ['obs_log2cr_logp_k',
                           'obs_maf_logp_k',
                           'num_points_cr_k',
                           'num_points_maf_k',
                           'length_k',
                           'start_k',
                           'end_k',
                           'is_contig_end_k',
                           'is_maf_nan_k',
                           'allelic_copy_number_product_states_lij',
                           'marginalization_allelic_copy_number_product_states_lij'])
DiscretePrior = namedtuple('DiscretePrior', ['product_state_log_prior_kl',
                                             'marginalization_product_state_log_prior_kl'])

def prepare_data_and_discrete_prior(modeled_segments, global_discrete_prior, discrete_prior_config, likelihood_config):
    data_k = np.array([[fit_t_to_log2cr(modeled_segment, likelihood_config.t_degrees_of_freedom), fit_beta_to_maf(modeled_segment), 
                        modeled_segment['NUM_POINTS_COPY_RATIO'], modeled_segment['NUM_POINTS_ALLELE_FRACTION'], modeled_segment['END'] - modeled_segment['START'] + 1] 
                       for _, modeled_segment in modeled_segments.iterrows()])
    
    obs_log2cr_mu_sigma_k = np.vstack(data_k[:, 0])
    obs_log2cr_k = obs_log2cr_mu_sigma_k[:, 0]

    obs_maf_a_b_k = np.vstack(data_k[:, 1])
    obs_maf_k = np.array(list(map(lambda x: min([max([0.5 * (x[0] - 1) / (x[0] + x[1] - 2), 0.]), 0.5]), obs_maf_a_b_k)))
    
    length_k = data_k[:, 4].astype(int)
    end_k = np.cumsum(length_k) / np.sum(length_k)
    start_k = np.concatenate([[0.], end_k[:-1]])
    contig_ends = modeled_segments.groupby('CONTIG', sort=False)['END'].max().reset_index()
    
    product_state_prior_length_weights_ik = np.array([np.maximum(length_k / discrete_prior_config.normal_population_event_length_scale, 1.),
                                                      np.maximum(length_k / discrete_prior_config.tumor_population_event_length_scale, 1.),
                                                      np.maximum(length_k / discrete_prior_config.tumor_population_event_length_scale, 1.)])
    product_state_log_prior_kl = np.einsum('ik,li->kl', product_state_prior_length_weights_ik, global_discrete_prior.unnorm_product_state_log_prior_li)
    marginalization_product_state_log_prior_kl = product_state_log_prior_kl[:, global_discrete_prior.marginalization_product_state_filter_l]
    
    data = Data(
        obs_log2cr_logp_k = lambda x_k: t_logpdf_jax_jit(x_k, df=likelihood_config.t_degrees_of_freedom, loc=obs_log2cr_mu_sigma_k[:, 0], scale=obs_log2cr_mu_sigma_k[:, 1]),
        obs_maf_logp_k = lambda x_k: beta_logpdf(x_k, a=obs_maf_a_b_k[:, 0], b=obs_maf_a_b_k[:, 1], scale=0.5),
        num_points_cr_k = data_k[:, 2].astype(int),
        num_points_maf_k = data_k[:, 3].astype(int),
        length_k = length_k,
        start_k = start_k,
        end_k = end_k,
        is_contig_end_k = modeled_segments[['CONTIG', 'END']].apply(tuple, axis=1).isin(contig_ends.apply(tuple, axis=1)),
        is_maf_nan_k = np.all(np.isnan(obs_maf_a_b_k), axis=1),
        allelic_copy_number_product_states_lij = global_discrete_prior.allelic_copy_number_product_states_lij,
        marginalization_allelic_copy_number_product_states_lij = global_discrete_prior.marginalization_allelic_copy_number_product_states_lij)
    
    discrete_prior = DiscretePrior(
        product_state_log_prior_kl = product_state_log_prior_kl,
        marginalization_product_state_log_prior_kl = marginalization_product_state_log_prior_kl)
    
    return data, discrete_prior
  
#===============================================================================
# model
#===============================================================================

eps_model = 1E-10

Parameters = namedtuple('Parameters', ['subclonal_cancer_cell_fraction_s',
                                       'purity',
                                       'cr_norm'])
TransformedParameters = namedtuple('TransformedParameters', ['transformed_subclonal_cancer_cell_fraction_s',
                                                             'purity',
                                                             'cr_norm'])
def transformed_parameters_array_to_tuple(transformed_parameters_array):
    return TransformedParameters(
        transformed_subclonal_cancer_cell_fraction_s = transformed_parameters_array[:-2], 
        purity = transformed_parameters_array[-2], 
        cr_norm = transformed_parameters_array[-1])
  
def calculate_log2cr_maf_sl(parameters, data, use_marginalization_states):
    num_subclonal_populations = len(parameters.subclonal_cancer_cell_fraction_s)
    allelic_copy_number_product_states_lij = data.marginalization_allelic_copy_number_product_states_lij if use_marginalization_states else data.allelic_copy_number_product_states_lij
    # s = subclonal population, l = product state, j = allele
    overlapping_population_fraction_is = np.array([(1. - parameters.purity) * np.ones(num_subclonal_populations), 
                                                   parameters.purity * (1. - parameters.subclonal_cancer_cell_fraction_s),
                                                   parameters.purity * parameters.subclonal_cancer_cell_fraction_s])
    population_weighted_copy_number_slj = np.einsum('is,lij->slj', overlapping_population_fraction_is, allelic_copy_number_product_states_lij)
    total_cr_sl = np.sum(population_weighted_copy_number_slj, axis=-1)
    log2cr_sl = np.log2(total_cr_sl / (parameters.cr_norm + eps_model) + eps_model)
    maf_sl = np.min(population_weighted_copy_number_slj, axis=-1) / (total_cr_sl + eps_model)
    return log2cr_sl, maf_sl
  
def calculate_likelihood_logp_ksl(parameters, data, use_marginalization_states):
    log2cr_sl, maf_sl = calculate_log2cr_maf_sl(parameters, data, use_marginalization_states)
    
    likelihood_log2cr_logp_slk = data.obs_log2cr_logp_k(log2cr_sl[:, :, np.newaxis])
    likelihood_maf_logp_slk = data.obs_maf_logp_k(maf_sl[:, :, np.newaxis])
    likelihood_maf_logp_slk = np.array(likelihood_maf_logp_slk) # cannot use index assignment with jax
    likelihood_maf_logp_slk[:, :, data.is_maf_nan_k] = 0.
    
    likelihood_logp_slk = likelihood_log2cr_logp_slk + likelihood_maf_logp_slk
    likelihood_logp_ksl = np.transpose(likelihood_logp_slk, axes=(2, 0, 1))
    
    return likelihood_logp_ksl
  
def prior_logp(prior, transformed_parameters):
    logp = 0.
    logp += prior.continuous_prior.subclonal_cancer_cell_fraction_dirichlet_hypercube_logp(transformed_parameters.transformed_subclonal_cancer_cell_fraction_s)
    logp += prior.continuous_prior.purity_beta.logpdf(transformed_parameters.purity)
    logp += prior.continuous_prior.cr_norm_lognorm.logpdf(transformed_parameters.cr_norm)
    return logp
    
  
def logp(transformed_parameters_array, prior, data):
    transformed_parameters = transformed_parameters_array_to_tuple(transformed_parameters_array)
    
    if np.any(transformed_parameters.transformed_subclonal_cancer_cell_fraction_s < 0) or \
          np.any(transformed_parameters.transformed_subclonal_cancer_cell_fraction_s > 1) or \
          not 0 <= transformed_parameters.purity <= 1 or \
          transformed_parameters.cr_norm < 0:
        return -np.inf
      
    subclonal_cancer_cell_fraction_s = to_simplex(transformed_parameters.transformed_subclonal_cancer_cell_fraction_s)
    
    # enforce CCF sorting
    if not np.all(np.diff(subclonal_cancer_cell_fraction_s) <= 0):
        return -np.inf
      
    logp = 0.
    
    # prior
    logp += prior_logp(prior, transformed_parameters)
    if logp == -np.inf:
        return logp
      
    parameters = Parameters(
        subclonal_cancer_cell_fraction_s = subclonal_cancer_cell_fraction_s,
        purity = transformed_parameters.purity,
        cr_norm = transformed_parameters.cr_norm)
    
    # likelihood
    likelihood_logp_ksl = calculate_likelihood_logp_ksl(parameters, data, use_marginalization_states=True)
    
    # marginalize over subclonal population assignments and marginalization product states
    logp_ksl = prior.discrete_prior.marginalization_product_state_log_prior_kl[:, np.newaxis, :] + likelihood_logp_ksl
    logp += jnp.sum(jax.scipy.special.logsumexp(logp_ksl, axis=(-2, -1)))
    
    # L2 prior on equivalence between MAP ploidy and cr_norm
    flat_argmax_idx = np.argmax(np.reshape(logp_ksl, (logp_ksl.shape[0], -1)), axis=-1)
    map_s_k, map_l_k = np.unravel_index(flat_argmax_idx, logp_ksl.shape[-2:])
    map_copy_number_ijk = np.transpose(data.marginalization_allelic_copy_number_product_states_lij[map_l_k], axes=(1, 2, 0))
    map_z_sk = map_s_k
    subclonal_cancer_cell_fraction_k = parameters.subclonal_cancer_cell_fraction_s[map_z_sk]
    overlapping_population_fraction_ik = np.array([(1. - parameters.purity) * np.ones(len(map_z_sk)), 
                                                   parameters.purity * (1. - subclonal_cancer_cell_fraction_k),
                                                   parameters.purity * subclonal_cancer_cell_fraction_k])
    population_weighted_copy_number_jk = np.einsum('ik,ijk->jk', overlapping_population_fraction_ik, map_copy_number_ijk)
    total_cr_k = np.sum(population_weighted_copy_number_jk, axis=0)
    ploidy = np.sum(total_cr_k * data.num_points_cr_k) / np.sum(data.num_points_cr_k)
    logp += prior.continuous_prior.cr_norm_constraint_norm.logpdf(parameters.cr_norm - ploidy)
    
    return logp
  
DiscreteParameters = namedtuple('DiscreteParameters', ['copy_number_ijk',
                                                       'z_sk'])
def calculate_map_discrete_parameters(parameters, discrete_prior, data, use_marginalization_states):
    subclonal_cancer_cell_fraction_s, purity, cr_norm = parameters
    allelic_copy_number_product_states_lij = data.marginalization_allelic_copy_number_product_states_lij if use_marginalization_states else data.allelic_copy_number_product_states_lij
    product_state_log_prior_kl = discrete_prior.marginalization_product_state_log_prior_kl if use_marginalization_states else discrete_prior.product_state_log_prior_kl
    
    likelihood_logp_ksl = calculate_likelihood_logp_ksl(parameters, data, use_marginalization_states)
    logp_ksl = product_state_log_prior_kl[:, np.newaxis, :] + likelihood_logp_ksl
    
    flat_argmax_idx = np.argmax(np.reshape(logp_ksl, (logp_ksl.shape[0], -1)), axis=-1)
    map_s_k, map_l_k = np.unravel_index(flat_argmax_idx, logp_ksl.shape[-2:])
    
    map_copy_number_ijk = np.transpose(allelic_copy_number_product_states_lij[map_l_k], axes=(1, 2, 0))
    map_z_sk = map_s_k
    
    return DiscreteParameters(map_copy_number_ijk, map_z_sk)
  
#===============================================================================
# inference
#===============================================================================
  
def run_mcmc(inference_config, prior, data, logp):
    num_parameters, num_walkers, num_burn_in, num_samples, emcee_move = inference_config
    np.random.seed(0)
    
    # initialize walkers
    ensemble_position = np.random.rand(num_walkers, num_parameters)
    ensemble_position[:, -2] = prior.continuous_prior.purity_beta.rvs(num_walkers)
    ensemble_position[:, -1] = prior.continuous_prior.cr_norm_lognorm.rvs(num_walkers)
        
    # enforce CCF sorting
    for w, walker_position in enumerate(ensemble_position):
        transformed_subclonal_cancer_cell_fraction_s = walker_position[:-2]
        subclonal_cancer_cell_fraction_s = to_simplex(transformed_subclonal_cancer_cell_fraction_s)
        if not np.all(np.diff(subclonal_cancer_cell_fraction_s) <= 0):
            sorted_subclonal_cancer_cell_fraction_s = np.sort(subclonal_cancer_cell_fraction_s)[::-1]
            ensemble_position[w, :-2] = to_hypercube(sorted_subclonal_cancer_cell_fraction_s)

#     with Pool() as pool:
#           sampler = emcee.EnsembleSampler(num_walkers, num_parameters, logp, args=(prior, data), pool=pool, moves=emcee_move)
    sampler = emcee.EnsembleSampler(num_walkers, num_parameters, logp, args=(prior, data), moves=emcee_move)
    # burn-in sampling
    for ensemble_position, _, _ in tqdm.tqdm(sampler.sample(ensemble_position, iterations=num_burn_in), total=num_burn_in):
        pass
    # sampling   
    for ensemble_position, _, _ in tqdm.tqdm(sampler.sample(ensemble_position, iterations=num_samples), total=num_samples):
        pass
    
    acceptance_fraction = np.mean(sampler.acceptance_fraction)
    print('Acceptance fraction:', acceptance_fraction)
    
    transformed_parameter_samples = sampler.get_chain(flat=True)[num_burn_in:]
    logp_samples = sampler.get_log_prob(flat=True)[num_burn_in:]
    
    # postprocess samples
    subclonal_cancer_cell_fraction_s_samples = np.apply_along_axis(to_simplex, 1, transformed_parameter_samples[:, :-2])
    subclonal_cancer_cell_fraction_s_samples = np.sort(subclonal_cancer_cell_fraction_s_samples, axis=1)[:, ::-1]
    purity_samples = transformed_parameter_samples[:, -2][:, np.newaxis]
    cr_norm_samples = transformed_parameter_samples[:, -1][:, np.newaxis]
    parameter_samples = np.hstack([subclonal_cancer_cell_fraction_s_samples, purity_samples, cr_norm_samples])
    
    return sampler, parameter_samples, logp_samples
  
#===============================================================================
# vectorized methods
#===============================================================================
  
# we use p to index parameters and P to index non-transformed parameters

@numba.jit(nopython=True)
def to_simplex_w(z_wP):
    num_walkers, num_hypercube_dimensions = z_wP.shape
    num_simplex_dimensions = num_hypercube_dimensions + 1
    x_wp = np.zeros((num_walkers, num_simplex_dimensions))
    x_wp[:, 0] = 1. - z_wP[:, 0]
    cum_prod_w = np.ones(num_walkers)
    for k in range(1, num_hypercube_dimensions):
        cum_prod_w *= z_wP[:, k - 1]
        x_wp[:, k] = cum_prod_w * (1. - z_wP[:, k]) 
    x_wp[:, num_hypercube_dimensions] = cum_prod_w * z_wP[:, num_hypercube_dimensions - 1]
    return x_wp

@numba.jit(nopython=True)
def to_hypercube_w(x_wp):
    num_walkers, num_simplex_dimensions = x_wp.shape
    num_hypercube_dimensions = num_simplex_dimensions - 1
    z_wP = np.zeros((num_walkers, num_hypercube_dimensions))
    z_wP[:, 0] = 1. - x_wp[:, 0]
    cum_prod_w = np.ones(num_walkers)
    for k in range(1, num_hypercube_dimensions):
        cum_prod_w *= z_wP[:, k - 1]
        z_wP[:, k] = 1. - x_wp[:, k] / cum_prod_w
    return z_wP
  
def dirichlet_hypercube_logp_w(alpha, z_wP):
    # up to factor involving Gammas independent of z_wP
    alpha_tilde = np.sum(alpha) - np.cumsum(alpha)
    return np.sum(beta_logpdf(z_wP, a=alpha_tilde[:-1], b=alpha[:-1]), axis=1)
  
def calculate_log2cr_maf_wsl(ensemble_wp, data, use_marginalization_states):
    num_walkers, num_parameters = np.shape(ensemble_wp)
    num_subclonal_populations = num_parameters - 2
    subclonal_cancer_cell_fraction_ws = ensemble_wp[:, :-2]
    purity_ws = ensemble_wp[:, -2, np.newaxis]
    cr_norm_w = ensemble_wp[:, -1]
    allelic_copy_number_product_states_lij = data.marginalization_allelic_copy_number_product_states_lij if use_marginalization_states else data.allelic_copy_number_product_states_lij
    # s = subclonal population, l = product state, j = allele
    overlapping_population_fraction_iws = np.array([(1. - purity_ws) * np.ones((num_walkers, num_subclonal_populations)), 
                                                   purity_ws * (1. - subclonal_cancer_cell_fraction_ws),
                                                   purity_ws * subclonal_cancer_cell_fraction_ws])
    population_weighted_copy_number_wslj = np.einsum('iws,lij->wslj', overlapping_population_fraction_iws, allelic_copy_number_product_states_lij)
    total_cr_wsl = np.sum(population_weighted_copy_number_wslj, axis=-1)
    log2cr_wsl = np.log2(total_cr_wsl / (cr_norm_w[:, np.newaxis, np.newaxis] + eps_model) + eps_model)
    maf_wsl = np.min(population_weighted_copy_number_wslj, axis=-1) / (total_cr_wsl + eps_model)
    return log2cr_wsl, maf_wsl
  
def calculate_likelihood_logp_wksl(ensemble_wp, data, use_marginalization_states):
    log2cr_wsl, maf_wsl = calculate_log2cr_maf_wsl(ensemble_wp, data, use_marginalization_states)
    
    likelihood_log2cr_logp_wslk = data.obs_log2cr_logp_k(log2cr_wsl[:, :, :, np.newaxis])
    likelihood_maf_logp_wslk = data.obs_maf_logp_k(maf_wsl[:, :, :, np.newaxis])
    likelihood_maf_logp_wslk = np.array(likelihood_maf_logp_wslk) # cannot use index assignment with jax
    likelihood_maf_logp_wslk[:, :, :, data.is_maf_nan_k] = 0.
    
    likelihood_logp_wslk = likelihood_log2cr_logp_wslk + likelihood_maf_logp_wslk
    likelihood_logp_wksl = np.transpose(likelihood_logp_wslk, axes=(0, 3, 1, 2))
    
    return likelihood_logp_wksl
  
def prior_logp_w(prior, transformed_ensemble_wP):
    num_walkers, _ = np.shape(transformed_ensemble_wP)
    logp_w = np.zeros(num_walkers)
    logp_w += prior.continuous_prior.subclonal_cancer_cell_fraction_dirichlet_hypercube_logp(transformed_ensemble_wP[:, :-2])
    logp_w += prior.continuous_prior.purity_beta.logpdf(transformed_ensemble_wP[:, -2])
    logp_w += prior.continuous_prior.cr_norm_lognorm.logpdf(transformed_ensemble_wP[:, -1])
    return logp_w
  
def logp_w(transformed_ensemble_wP, prior, data):
    num_walkers, num_transformed_parameters = np.shape(transformed_ensemble_wP)
    num_parameters = num_transformed_parameters + 1
    
    is_in_transformed_boundary_w = \
          np.all(transformed_ensemble_wP[:, :-1] >= 0, axis=1) & np.all(transformed_ensemble_wP[:, :-1] <= 1, axis=1) & \
          (transformed_ensemble_wP[:, -1] > 0)
    
    subclonal_cancer_cell_fraction_ws = to_simplex_w(transformed_ensemble_wP[:, :-2])
    
    # enforce CCF sorting
    is_subclonal_cancer_cell_fraction_sorted_w = np.all(np.diff(subclonal_cancer_cell_fraction_ws, axis=1) <= 0, axis=1)
    
    is_allowed_w = is_in_transformed_boundary_w & is_subclonal_cancer_cell_fraction_sorted_w
    
    # prior
    logp_w = np.choose(is_allowed_w, [-np.inf, 0.])
    logp_w += prior_logp_w(prior, transformed_ensemble_wP)

    if np.all(logp_w == -np.inf):
        return logp_w
      
    ensemble_wp = np.zeros((num_walkers, num_parameters))
    ensemble_wp[:, :-2] = subclonal_cancer_cell_fraction_ws
    ensemble_wp[:, -2] = transformed_ensemble_wP[:, -2]
    ensemble_wp[:, -1] = transformed_ensemble_wP[:, -1]
    
    # likelihood
    likelihood_logp_wksl = calculate_likelihood_logp_wksl(ensemble_wp * is_allowed_w[:, np.newaxis], data, use_marginalization_states=True)
    
    # marginalize over subclonal population assignments and marginalization product states
    logp_wksl = prior.discrete_prior.marginalization_product_state_log_prior_kl[np.newaxis, :, np.newaxis, :] + likelihood_logp_wksl
    logp_w += jnp.sum(jax.scipy.special.logsumexp(logp_wksl, axis=(-2, -1)), axis=-1)
    
    # L2 prior on equivalence between MAP ploidy and cr_norm
    flat_argmax_idx = jnp.argmax(np.reshape(logp_wksl, (logp_wksl.shape[0], logp_wksl.shape[1], -1)), axis=-1)
    map_s_wk, map_l_wk = np.unravel_index(flat_argmax_idx, logp_wksl.shape[-2:])
    map_copy_number_wijk = np.transpose(data.marginalization_allelic_copy_number_product_states_lij[map_l_wk], axes=(0, 2, 3, 1))
    map_z_wsk = map_s_wk
    subclonal_cancer_cell_fraction_wk = np.take_along_axis(subclonal_cancer_cell_fraction_ws, map_z_wsk, axis=-1)
    overlapping_population_fraction_iwk = np.array([(1. - ensemble_wp[:, -2, np.newaxis]) * np.ones(map_z_wsk.shape), 
                                                   ensemble_wp[:, -2, np.newaxis] * (1. - subclonal_cancer_cell_fraction_wk),
                                                   ensemble_wp[:, -2, np.newaxis] * subclonal_cancer_cell_fraction_wk])
    population_weighted_copy_number_wjk = np.einsum('iwk,wijk->wjk', overlapping_population_fraction_iwk, map_copy_number_wijk)
    total_cr_wk = np.sum(population_weighted_copy_number_wjk, axis=1)
    ploidy_w = np.sum(total_cr_wk * data.num_points_cr_k, axis=-1) / np.sum(data.num_points_cr_k)
    logp_w += prior.continuous_prior.cr_norm_constraint_norm.logpdf(ensemble_wp[:, -1] - ploidy_w)
    
    return logp_w
  
def run_mcmc_vectorized(inference_config, prior, data, logp):
    num_parameters, num_walkers, num_burn_in, num_samples, emcee_move = inference_config
    np.random.seed(0)
    
    # initialize walkers
    transformed_ensemble_wP = np.random.rand(num_walkers, num_parameters)
    transformed_ensemble_wP[:, -2] = prior.continuous_prior.purity_beta.rvs(num_walkers)
    transformed_ensemble_wP[:, -1] = prior.continuous_prior.cr_norm_lognorm.rvs(num_walkers)
        
    # enforce CCF sorting
    for w, transformed_walker_P in enumerate(transformed_ensemble_wP):
        transformed_subclonal_cancer_cell_fraction_s = transformed_walker_P[:-2]
        subclonal_cancer_cell_fraction_s = to_simplex(transformed_subclonal_cancer_cell_fraction_s)
        if not np.all(np.diff(subclonal_cancer_cell_fraction_s) <= 0):
            sorted_subclonal_cancer_cell_fraction_s = np.sort(subclonal_cancer_cell_fraction_s)[::-1]
            transformed_ensemble_wP[w, :-2] = to_hypercube(sorted_subclonal_cancer_cell_fraction_s)

    sampler = emcee.EnsembleSampler(num_walkers, num_parameters, logp_w, args=(prior, data), moves=emcee_move, vectorize=True)
    # burn-in sampling
    for transformed_walker_P, _, _ in tqdm.tqdm(sampler.sample(transformed_ensemble_wP, iterations=num_burn_in), total=num_burn_in):
        pass
    # sampling   
    for transformed_walker_P, _, _ in tqdm.tqdm(sampler.sample(transformed_ensemble_wP, iterations=num_samples), total=num_samples):
        pass
    
    acceptance_fraction = np.mean(sampler.acceptance_fraction)
    print('Acceptance fraction:', acceptance_fraction)
    
    transformed_parameter_samples = sampler.get_chain(flat=True)[num_burn_in:]
    logp_samples = sampler.get_log_prob(flat=True)[num_burn_in:]
    
    # postprocess samples
    subclonal_cancer_cell_fraction_s_samples = np.apply_along_axis(to_simplex, 1, transformed_parameter_samples[:, :-2])
    subclonal_cancer_cell_fraction_s_samples = np.sort(subclonal_cancer_cell_fraction_s_samples, axis=1)[:, ::-1]
    purity_samples = transformed_parameter_samples[:, -2][:, np.newaxis]
    cr_norm_samples = transformed_parameter_samples[:, -1][:, np.newaxis]
    parameter_samples = np.hstack([subclonal_cancer_cell_fraction_s_samples, purity_samples, cr_norm_samples])
    
    return sampler, parameter_samples, logp_samples

#===============================================================================
# plotting
#===============================================================================

def plot_cr_maf_posteriors(modeled_segments, output_path, output_prefix, show=True):
    length_k = modeled_segments['END'].values - modeled_segments['START'].values + 1
    alpha_k = np.log10(length_k) / np.log10(np.max(length_k))
    colors_rgba = [np.array(matplotlib.colors.to_rgba(c)) for c in ['r', 'grey']]
    facecolor_k = np.take(colors_rgba, np.isnan(modeled_segments['MINOR_ALLELE_FRACTION_POSTERIOR_50']), axis=0)
    facecolor_k[:, 3] *= alpha_k
    
#     log2cr_min = np.floor(np.min(modeled_segments['LOG2_COPY_RATIO_POSTERIOR_10']))
    cr_min = 0.
    log2cr_max = np.ceil(np.max(modeled_segments['LOG2_COPY_RATIO_POSTERIOR_90']))
    maf_min = 0.
    maf_max = 0.5
    
    fig, ax = plt.subplots(1, figsize=(10, 5))
    
    def make_rectangle(modeled_segment):
        cr_10 = 2**modeled_segment['LOG2_COPY_RATIO_POSTERIOR_10']
        cr_90 = 2**modeled_segment['LOG2_COPY_RATIO_POSTERIOR_90']
        maf_10 = modeled_segment['MINOR_ALLELE_FRACTION_POSTERIOR_10']        
        maf_90 = modeled_segment['MINOR_ALLELE_FRACTION_POSTERIOR_90']        
        maf_10 = maf_min if np.isnan(maf_10) else maf_10
        maf_90 = maf_max if np.isnan(maf_90) else maf_90
        return Rectangle([cr_10, maf_10], cr_90 - cr_10, maf_90 - maf_10)
    
    rectangles = PatchCollection([make_rectangle(modeled_segment) for _, modeled_segment in modeled_segments.iterrows()],
                                 facecolor=facecolor_k,
                                 edgecolor='k')
    ax.add_collection(rectangles)

    plt.xlabel('copy ratio')
    plt.ylabel('minor-allele fraction')
    plt.xlim([cr_min, log2cr_max])
    plt.ylim([maf_min, maf_max])
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.scatter.png'))
    if show:
        plt.show()
    plt.close()
    
def plot_corner(parameter_samples, output_path, output_prefix, show=True):
    num_subclonal_populations = parameter_samples.shape[1] - 2

    fig = corner.corner(parameter_samples, 
                        labels=['$\pi_{}$'.format(i) for i in range(num_subclonal_populations)] + ['purity', 'cr_norm'],
                        quantiles=[0.1, 0.5, 0.9], show_titles=True)
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.corner.png'))
    if show:
        plt.show()
    plt.close()
    
def plot_copy_number_ijk_samples(copy_number_ijk_samples, data, allelic_copy_number_states, output_path, output_prefix, show=True):
    num_samples, num_populations, num_alleles, num_segments = np.shape(copy_number_ijk_samples)
    
    fig, axs = plt.subplots(nrows=num_populations, sharex=True, sharey=True, figsize=(18, 6))
    y_max = allelic_copy_number_states[-1] + 1
    plt.xlim([0, 1])
    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[2].set_xticks([])
    plt.ylim([0, y_max])
    axs[0].set_ylabel('normal allelic CN')
    axs[1].set_ylabel('clonal allelic CN')
    axs[2].set_ylabel('subclonal allelic CN')
    
    allele_rgba = [np.array(matplotlib.colors.to_rgba(c)) for c in ['b', 'g']]
    
    for k in range(num_segments):
        copy_number_ij_states, counts = np.unique(copy_number_ijk_samples[:, :, :, k], axis=0, return_counts=True)
        normalized_counts = counts / np.sum(counts)
        num_states = len(copy_number_ij_states)
        
        for i in range(num_populations):
            axs[i].vlines(data.end_k[data.is_contig_end_k], 0, y_max, color='grey', linestyle='dashed', alpha=0.5)
            for j in range(num_alleles):
                colors = np.tile(allele_rgba[j], (num_states, 1))
                colors[:, 3] *= normalized_counts
                lc = LineCollection([[[data.start_k[k], copy_number_ij_state[i][j]], [data.end_k[k], copy_number_ij_state[i][j]]] for copy_number_ij_state in copy_number_ij_states],
                                    color=colors, lw=4, alpha=0.5)
                axs[i].add_collection(lc)
    
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.acn.png'))
    if show:
        plt.show()
    plt.close()
    
def plot_subclonal_diagram(parameters, discrete_parameters, data, normal_allelic_copy_number_state, output_path, output_prefix, show=True):
    copy_number_ijk, z_sk = discrete_parameters
    _, _, num_segments = np.shape(copy_number_ijk)
    subclonal_cancer_cell_fraction_k = parameters.subclonal_cancer_cell_fraction_s[z_sk]
    is_normal_k = np.all(copy_number_ijk == normal_allelic_copy_number_state, axis=(0, 1))
    is_clonal_k = np.all(copy_number_ijk[1, :, :] == copy_number_ijk[2, :, :], axis=0)
    
    fig, ax = plt.subplots(figsize=(18, 3))
    plt.xlim([0, 1])
    ax.set_xticks([])
    ax.set_ylim([0, 1])
    ax.set_ylabel('subclonal CCF')
    
    ax.vlines(data.end_k[data.is_contig_end_k], 0, 4, color='grey', linestyle='dashed', alpha=0.5)
    
    lc = LineCollection([[[data.start_k[k], 0], [data.end_k[k], 0]] for k in np.nonzero(is_normal_k)[0]],
                        color='b', lw=4, alpha=0.5, label='normal')
    ax.add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 0], [data.end_k[k], 0]] for k in np.nonzero(~is_normal_k & is_clonal_k)[0]],
                        color='g', lw=4, alpha=0.5, label='clonal')
    ax.add_collection(lc)
    lc = LineCollection([[[data.start_k[k], subclonal_cancer_cell_fraction_k[k]], [data.end_k[k], subclonal_cancer_cell_fraction_k[k]]] for k in np.nonzero(~is_normal_k & ~is_clonal_k)[0]],
                        color='r', lw=4, alpha=0.5)
    ax.add_collection(lc)
    
    ax.legend(loc='lower right', bbox_to_anchor= (1.08, 0.))
    
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.ccf.png'))
    if show:   
        plt.show()
    plt.close()
  
def plot_fit(parameters, discrete_parameters, data, modeled_segments, output_path, output_prefix, show=True):
    subclonal_cancer_cell_fraction_s, purity, cr_norm = parameters
    copy_number_ijk, z_sk = discrete_parameters
    
    _, _, num_segments = np.shape(copy_number_ijk)
    
    subclonal_cancer_cell_fraction_k = subclonal_cancer_cell_fraction_s[z_sk]
    overlapping_population_fraction_ik = np.array([(1. - purity) * np.ones(num_segments), 
                                                   purity * (1. - subclonal_cancer_cell_fraction_k),
                                                   purity * subclonal_cancer_cell_fraction_k])
    population_weighted_copy_number_jk = np.einsum('ik,ijk->jk', overlapping_population_fraction_ik, copy_number_ijk)
    total_cr_k = np.sum(population_weighted_copy_number_jk, axis=0)
    log2cr_k = np.log2(total_cr_k / (cr_norm + eps_model) + eps_model)
    maf_k = np.min(population_weighted_copy_number_jk, axis=0) / (total_cr_k + eps_model)
    
    figs, axs = plt.subplots(nrows=2, sharex=True, figsize=(18, 6))
    plt.xlim([0, 1])
    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[0].set_ylim([0, 3])
    axs[1].set_ylim([0, 0.5])
    axs[0].set_ylabel('copy ratio')
    axs[1].set_ylabel('minor-allele fraction')
    
    axs[0].vlines(data.end_k[data.is_contig_end_k], 0, 4, color='grey', linestyle='dashed', alpha=0.5)
    axs[1].vlines(data.end_k[data.is_contig_end_k], 0, 0.5, color='grey', linestyle='dashed', alpha=0.5)
    
    lc = LineCollection([[[data.start_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_10']], [data.end_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_10']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.5)
    axs[0].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_50']], [data.end_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_50']]] for k in range(num_segments)],
                        color='r', lw=4, alpha=0.5)
    axs[0].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_90']], [data.end_k[k], 2**modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_90']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.5)
    axs[0].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 2**log2cr_k[k]], [data.end_k[k], 2**log2cr_k[k]]] for k in range(num_segments)],
                        color='b', lw=4, alpha=0.5)
    axs[0].add_collection(lc)
    
    lc = LineCollection([[[data.start_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_10']], [data.end_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_10']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.5)
    axs[1].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_50']], [data.end_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_50']]] for k in range(num_segments)],
                        color='r', lw=4, alpha=0.5, label='data')
    axs[1].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_90']], [data.end_k[k], modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_90']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.5)
    axs[1].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], maf_k[k]], [data.end_k[k], maf_k[k]]] for k in range(num_segments)],
                        color='b', lw=4, alpha=0.5, label='model')
    axs[1].add_collection(lc)
    
    axs[1].legend(loc='lower right', bbox_to_anchor= (1.08, 0.))
    
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.fit.png'))
    if show:
        plt.show()    
    plt.close()
    
# num_overlapping_populations and num_alleles are essentially fixed, due to assumptions made in discrete prior
DiscretePriorConfig = namedtuple('DiscretePriorConfig', ['num_overlapping_populations',
                                                         'num_alleles',
                                                         'allelic_copy_number_states',
                                                         'normal_allelic_copy_number_state',
                                                         'allelic_copy_number_change_prior_prob',
                                                         'hom_del_prior_prob',
                                                         'num_marginalization_product_states',
                                                         'normal_population_event_length_scale',
                                                         'tumor_population_event_length_scale'])
discrete_prior_config = DiscretePriorConfig(
    num_overlapping_populations = 3,
    num_alleles = 2,
    allelic_copy_number_states = np.arange(6 + 1),
    normal_allelic_copy_number_state = 1,
    allelic_copy_number_change_prior_prob = 0.99,
    hom_del_prior_prob = 1E-6,
    num_marginalization_product_states = 200,
    normal_population_event_length_scale = 1E4,
    tumor_population_event_length_scale = 1E7)

global_discrete_prior = generate_label_ordered_product_states_and_log_prior(discrete_prior_config)

ContinuousPriorConfig = namedtuple('ContinuousPriorConfig', ['num_subclonal_populations',
                                                             'subclonal_cancer_cell_fraction_alpha',
                                                             'purity_a',
                                                             'purity_b',
                                                             'cr_norm_s',
                                                             'cr_norm_scale',
                                                             'cr_norm_constraint_scale'])
continuous_prior_config = ContinuousPriorConfig(
    num_subclonal_populations = 4,
    subclonal_cancer_cell_fraction_alpha = 1E-2,
    purity_a = 1.,
    purity_b = 10.,
    cr_norm_s = 0.1,
    cr_norm_scale = 2.,
    cr_norm_constraint_scale = np.sqrt(1E-3))

ContinuousPrior = namedtuple('ContinuousPrior', ['subclonal_cancer_cell_fraction_dirichlet_hypercube_logp',
                                                 'purity_beta',
                                                 'cr_norm_lognorm',
                                                 'cr_norm_constraint_norm'])

continuous_prior = ContinuousPrior(
    subclonal_cancer_cell_fraction_dirichlet_hypercube_logp = \
      lambda transformed_subclonal_cancer_cell_fraction_s: dirichlet_hypercube_logp(
        continuous_prior_config.subclonal_cancer_cell_fraction_alpha * np.ones(continuous_prior_config.num_subclonal_populations), 
        transformed_subclonal_cancer_cell_fraction_s),
    purity_beta = scipy.stats.beta(
        a=continuous_prior_config.purity_a, 
        b=continuous_prior_config.purity_b),
    cr_norm_lognorm = scipy.stats.lognorm(
        s=continuous_prior_config.cr_norm_s, 
        scale=continuous_prior_config.cr_norm_scale),
    cr_norm_constraint_norm = scipy.stats.norm(
        scale=continuous_prior_config.cr_norm_constraint_scale))

Prior = namedtuple('Prior', ['continuous_prior',
                             'discrete_prior'])

LikelihoodConfig = namedtuple('LikelihoodConfig', ['t_degrees_of_freedom'])
likelihood_config = LikelihoodConfig(t_degrees_of_freedom = 5)

InferenceConfig = namedtuple('InferenceConfig', ['num_parameters',
                                                 'num_walkers',
                                                 'num_burn_in',
                                                 'num_samples',
                                                 'emcee_move'])
inference_config = InferenceConfig(
    num_parameters = (continuous_prior_config.num_subclonal_populations - 1) + 2,
    num_walkers = 128,
    num_burn_in = 200,
    num_samples = 50,
#    emcee_move = emcee.moves.WalkMove()
#    emcee_move = emcee.moves.StretchMove()
    emcee_move = emcee.moves.KDEMove()
#    emcee_move = emcee.moves.DEMove()
)

GlobalConfig = namedtuple('GlobalConfig', ['discrete_prior_config',
                                           'global_discrete_prior',
                                           'continuous_prior',
                                           'likelihood_config',
                                           'inference_config'])

global_config = GlobalConfig(
    discrete_prior_config = discrete_prior_config,
    global_discrete_prior = global_discrete_prior,
    continuous_prior = continuous_prior,
    likelihood_config = likelihood_config,
    inference_config = inference_config
)

continuous_prior_vectorized = ContinuousPrior(
    subclonal_cancer_cell_fraction_dirichlet_hypercube_logp = \
      lambda transformed_subclonal_cancer_cell_fraction_ws: dirichlet_hypercube_logp_w(
        continuous_prior_config.subclonal_cancer_cell_fraction_alpha * np.ones(continuous_prior_config.num_subclonal_populations), 
        transformed_subclonal_cancer_cell_fraction_ws),
    purity_beta = scipy.stats.beta(
        a=continuous_prior_config.purity_a, 
        b=continuous_prior_config.purity_b),
    cr_norm_lognorm = scipy.stats.lognorm(
        s=continuous_prior_config.cr_norm_s, 
        scale=continuous_prior_config.cr_norm_scale),
    cr_norm_constraint_norm = scipy.stats.norm(
        scale=continuous_prior_config.cr_norm_constraint_scale))


global_config_vectorized = GlobalConfig(
    discrete_prior_config = discrete_prior_config,
    global_discrete_prior = global_discrete_prior,
    likelihood_config = likelihood_config,
    continuous_prior = continuous_prior_vectorized,
    inference_config = inference_config
)

def run_th(modeled_segments_path, output_prefix, output_path, global_config, vectorize=True, show=True):
    discrete_prior_config, global_discrete_prior, continuous_prior, likelihood_config, inference_config = global_config
  
    # load modeled segments=====================================================
    modeled_segments = pd.read_csv(modeled_segments_path, sep='\t', comment='@', dtype={'CONTIG': str})

    print('Output prefix:', output_prefix)
    num_segments = len(modeled_segments)
    print('Number of segments:', num_segments)
    
    # output CR-MAF posteriors scatter plot=================================
    plot_cr_maf_posteriors(modeled_segments, output_path, output_prefix, show)
    
    # fit posteriors and apply segment-length weights to discrete priors========
    data, discrete_prior = prepare_data_and_discrete_prior(modeled_segments, global_discrete_prior, discrete_prior_config, likelihood_config)
    prior = Prior(continuous_prior, discrete_prior)
    
    # perform inference=========================================================
    
    sampler, parameter_samples, logp_samples = run_mcmc_vectorized(inference_config, prior, data, logp_w) if vectorize else run_mcmc(inference_config, prior, data, logp)
    
    # output corner plot========================================================
    
    plot_corner(parameter_samples, output_path, output_prefix, show)
    
    # get MAP discrete parameters===============================================
    
    map_sample_index = np.argmax(logp_samples)
    map_parameters = Parameters(
        parameter_samples[map_sample_index, :-2],
        parameter_samples[map_sample_index, -2],
        parameter_samples[map_sample_index, -1])
    map_discrete_parameters = calculate_map_discrete_parameters(map_parameters, discrete_prior, data, use_marginalization_states=False)
    
    print('MAP subclonal_cancer_cell_fraction_s:', map_parameters.subclonal_cancer_cell_fraction_s)
    print('MAP purity:', map_parameters.purity)
    print('MAP cr_norm:', map_parameters.cr_norm)
    print('MAP ploidy:', (map_parameters.cr_norm - 2 * (1 - map_parameters.purity)) / (map_parameters.purity + eps_model))
    
    # output MAP CN plot========================================================
    
    plot_copy_number_ijk_samples(np.array([map_discrete_parameters.copy_number_ijk]), data, discrete_prior_config.allelic_copy_number_states, 
                                 output_path, output_prefix, show)
    
    # output MAP subclonal plot=================================================
    
    plot_subclonal_diagram(map_parameters, map_discrete_parameters, data, discrete_prior_config.normal_allelic_copy_number_state,
                           output_path, output_prefix, show)
    
    # output MAP log2CR-MAF fit plot============================================
    
    plot_fit(map_parameters, map_discrete_parameters, data, modeled_segments, output_path, output_prefix, show)
    
    # output posterior table====================================================
    
    # save results==============================================================
    
    np.save(os.path.join(output_path, output_prefix + '.samples.npy'), parameter_samples)
    np.save(os.path.join(output_path, output_prefix + '.logp.npy'), logp_samples)
    
def main():
    parser = argparse.ArgumentParser(
        description='Model tumor heterogeneity given ModelSegments result.',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--output_path',
                        type=str,
                        help='Path to output directory.')

    parser.add_argument('--output_prefix',
                        type=str,
                        help='Output prefix.')

    parser.add_argument('--modeled_segments_path',
                        type=str,
                        help='Path to ModelSegments file.')
                        
    args = parser.parse_args()

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    run_th(args.modeled_segments_path, args.output_prefix, args.output_path, global_config_vectorized, vectorize=True, show=False)

if __name__ == '__main__':
    main()
