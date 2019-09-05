#!/usr/bin/env python

import os
import argparse
import re
import dill as pickle
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

from collections import namedtuple, OrderedDict

import tqdm

import emcee
import corner

#===============================================================================
# posterior fitting
#===============================================================================

eps_cr = 1E-3
eps_posterior_fitting = 1E-10

def fit_t_to_log2cr(modeled_segment, t_degrees_of_freedom):
    log2cr_10 = np.log2(2**modeled_segment['LOG2_COPY_RATIO_POSTERIOR_10'] + eps_cr)
    log2cr_50 = np.log2(2**modeled_segment['LOG2_COPY_RATIO_POSTERIOR_50'] + eps_cr)
    log2cr_90 = np.log2(2**modeled_segment['LOG2_COPY_RATIO_POSTERIOR_90'] + eps_cr)
    
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
    return np.array([a, b])

#===============================================================================
# logpdf utilities
#===============================================================================

def t_logpdf_jax(x, df, loc=0., scale=1.):
    r = df * 1.
    y = (x - loc) / scale
    lPx = jax.scipy.special.gammaln((r + 1) / 2) - jax.scipy.special.gammaln(r / 2)
    lPx -= 0.5 * jnp.log(r * jnp.pi) + (r + 1) / 2 * jnp.log(1 + y**2 / r)
    return lPx
t_logpdf_jax_jit = jax.jit(t_logpdf_jax)
  
eps_beta_logpdf = 1E-3
def beta_logpdf_jax(x, a, b, scale=1.):
    # up to factor independent of x
    x_bnd = jnp.minimum(scale - eps_beta_logpdf, jnp.maximum(eps_beta_logpdf, x)) / scale
    lPx = jax.scipy.special.xlog1py(b - 1.0, -x_bnd) + jax.scipy.special.xlogy(a - 1.0, x_bnd)
    return lPx
beta_logpdf_jax_jit = jax.jit(beta_logpdf_jax)

#===============================================================================
# discrete priors and data
#===============================================================================

def read_sequence_dictionary(modeled_segments_path):
    sequence_dictionary = OrderedDict()
    with open(modeled_segments_path, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                break
            contig_match = re.search("^@SQ.*SN:(.*)[\t]*.*$", line, re.M)
            length_match = re.search("^@SQ.*LN:(.*)[\t]*.*$", line, re.M)
            if contig_match is None or length_match is None:
                continue
            contig = contig_match.groups()[0].split('\t')[0]
            length = int(length_match.groups()[0].split('\t')[0])
            sequence_dictionary[contig] = length
    return sequence_dictionary

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
    copy_number_event_prior_penalty = discrete_prior_config.copy_number_event_prior_penalty
    allelic_copy_number_change_prior_penalty = discrete_prior_config.allelic_copy_number_change_prior_penalty
    hom_del_prior_penalty = discrete_prior_config.hom_del_prior_penalty

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
    
    unnorm_label_ordered_product_state_prior_li = np.ones((len(label_ordered_allelic_copy_number_product_states_lij), num_overlapping_populations))
    
    # penalize non-normal states
    is_non_normal_l = np.any(label_ordered_allelic_copy_number_product_states_lij != normal_allelic_copy_number_state, axis=(1, 2))
    unnorm_label_ordered_product_state_prior_li *= (1. - copy_number_event_prior_penalty)**is_non_normal_l[:, np.newaxis]
    
    # penalize non-germline states
#    is_non_germline_l = ~np.all((label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == label_ordered_allelic_copy_number_product_states_lij[:, 1, :]) *
#                                (label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == label_ordered_allelic_copy_number_product_states_lij[:, 2, :]), axis=-1)
#    unnorm_label_ordered_product_state_prior_li *= (1. - copy_number_event_prior_penalty)**is_non_germline_l[:, np.newaxis]
    
    # penalize subclonal non-clonal states
    is_subclonal_l = ~np.all(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] == label_ordered_allelic_copy_number_product_states_lij[:, 2, :], axis=-1)
    unnorm_label_ordered_product_state_prior_li *= (1. - copy_number_event_prior_penalty)**is_subclonal_l[:, np.newaxis]

    # penalize copy-number changes from population to population
    delta_copy_number_li = (np.array([np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 0, :] - normal_allelic_copy_number_state), axis=-1),
                                      np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 0, :] - label_ordered_allelic_copy_number_product_states_lij[:, 1, :]), axis=-1), 
                                      np.sum(np.abs(label_ordered_allelic_copy_number_product_states_lij[:, 1, :] - label_ordered_allelic_copy_number_product_states_lij[:, 2, :]), axis=-1)])).transpose()
    unnorm_label_ordered_product_state_prior_li *= (1. - allelic_copy_number_change_prior_penalty)**delta_copy_number_li
    
    # remove unphysical states where normal (clonal) deletions are reverted in clonal/subclonal (subclonal) population
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 1, :] != 0), axis=1)] = 0.
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 0, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 2, :] != 0), axis=1)] = 0.
    unnorm_label_ordered_product_state_prior_li[np.any((label_ordered_allelic_copy_number_product_states_lij[:, 1, :] == 0) * (label_ordered_allelic_copy_number_product_states_lij[:, 2, :] != 0), axis=1)] = 0.

    # heavily penalize hom dels
    unnorm_label_ordered_product_state_prior_li[np.any(np.all(label_ordered_allelic_copy_number_product_states_lij == 0, axis=-1), axis=-1)] *= (1. - hom_del_prior_penalty)

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

Data = namedtuple('Data', ['modeled_segments',
                           'sequence_dictionary',
                           'obs_log2cr_logp_k',
                           'obs_maf_logp_k',
                           'num_points_cr_k',
                           'num_points_maf_k',
                           'length_k',
                           'start_k',
                           'end_k',
                           'contig_ends',
                           'is_maf_nan_k',
                           'allelic_copy_number_product_states_lij',
                           'marginalization_allelic_copy_number_product_states_lij'])
DiscretePrior = namedtuple('DiscretePrior', ['product_state_log_prior_kl',
                                             'marginalization_product_state_log_prior_kl'])

def prepare_data_and_discrete_prior(modeled_segments, sequence_dictionary,
                                    global_discrete_prior, discrete_prior_config, likelihood_config):
    def translate_genomic_coordinate_to_plotting_coordinate(contig, position, 
                                                            total_genomic_length,
                                                            contig_starts_dictionary):
        return contig_starts_dictionary[contig] + position / total_genomic_length
        
    sequence_dictionary_for_plotting = OrderedDict((contig, sequence_dictionary[contig])
                                                   for contig in modeled_segments['CONTIG'].unique())
    contig_genomic_lengths = np.fromiter(sequence_dictionary_for_plotting.values(), dtype=int)
    total_genomic_length = np.sum(contig_genomic_lengths)
    contig_ends = np.cumsum(contig_genomic_lengths) / total_genomic_length
    contig_starts = np.concatenate([[0.], contig_ends[:-1]])
    contig_starts_dictionary = OrderedDict((contig, contig_starts[i])
                                           for i, contig in enumerate(sequence_dictionary_for_plotting.keys()))

    data_k = np.array([[fit_t_to_log2cr(modeled_segment, likelihood_config.t_degrees_of_freedom), fit_beta_to_maf(modeled_segment)]
                       for _, modeled_segment in modeled_segments.iterrows()])

    obs_log2cr_mu_sigma_k = np.vstack(data_k[:, 0])
    obs_maf_a_b_k = np.vstack(data_k[:, 1])
    length_k = (modeled_segments['END'] - modeled_segments['START'] + 1).values.astype(int)
    start_k = np.array([translate_genomic_coordinate_to_plotting_coordinate(
                            contig, start, total_genomic_length, contig_starts_dictionary)
                        for contig, start in zip(modeled_segments['CONTIG'], modeled_segments['START'])])
    end_k = np.array([translate_genomic_coordinate_to_plotting_coordinate(
                            contig, end, total_genomic_length, contig_starts_dictionary)
                        for contig, end in zip(modeled_segments['CONTIG'], modeled_segments['END'])])

    product_state_prior_length_weights_ik = np.array([np.maximum(length_k / discrete_prior_config.normal_population_event_length_scale, 1.),
                                                      np.maximum(length_k / discrete_prior_config.tumor_population_event_length_scale, 1.),
                                                      np.maximum(length_k / discrete_prior_config.tumor_population_event_length_scale, 1.)])
    product_state_log_prior_kl = np.einsum('ik,li->kl', product_state_prior_length_weights_ik, global_discrete_prior.unnorm_product_state_log_prior_li)
    marginalization_product_state_log_prior_kl = product_state_log_prior_kl[:, global_discrete_prior.marginalization_product_state_filter_l]
    
    data = Data(
        modeled_segments = modeled_segments,
        sequence_dictionary = sequence_dictionary,
        obs_log2cr_logp_k = lambda x_k: t_logpdf_jax_jit(x_k, df=likelihood_config.t_degrees_of_freedom, loc=obs_log2cr_mu_sigma_k[:, 0], scale=obs_log2cr_mu_sigma_k[:, 1]),
        obs_maf_logp_k = lambda x_k: beta_logpdf_jax_jit(x_k, a=obs_maf_a_b_k[:, 0], b=obs_maf_a_b_k[:, 1], scale=0.5),
        num_points_cr_k = modeled_segments['NUM_POINTS_COPY_RATIO'].values.astype(int),
        num_points_maf_k = modeled_segments['NUM_POINTS_ALLELE_FRACTION'].values.astype(int),
        length_k = length_k,
        start_k = start_k,
        end_k = end_k,
        contig_ends = contig_ends,
        is_maf_nan_k = np.all(np.isnan(obs_maf_a_b_k), axis=1),
        allelic_copy_number_product_states_lij = global_discrete_prior.allelic_copy_number_product_states_lij,
        marginalization_allelic_copy_number_product_states_lij = global_discrete_prior.marginalization_allelic_copy_number_product_states_lij)
    
    discrete_prior = DiscretePrior(
        product_state_log_prior_kl = product_state_log_prior_kl,
        marginalization_product_state_log_prior_kl = marginalization_product_state_log_prior_kl)
    
    return data, discrete_prior
  
#===============================================================================
# model/inference (vectorized methods)
#===============================================================================

eps_model = 1E-10

# hypercube-to-simplex transform: https://arxiv.org/pdf/1010.3436.pdf
# we use p to index parameters and P to index non-transformed (hypercube) parameters

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

def calculate_log2cr_maf_wsl(ensemble_wp, data, use_marginalization_states):
    num_walkers, num_parameters = np.shape(ensemble_wp)
    num_subclonal_populations = (num_parameters - 2) // 2
    subclonal_cancer_cell_fraction_weight_ws = ensemble_wp[:, :num_subclonal_populations]
    subclonal_cancer_cell_fraction_ws = ensemble_wp[:, num_subclonal_populations:-2]
    purity_ws = ensemble_wp[:, -2, np.newaxis]
    cr_norm_w = ensemble_wp[:, -1]
    allelic_copy_number_product_states_lij = data.marginalization_allelic_copy_number_product_states_lij if use_marginalization_states else data.allelic_copy_number_product_states_lij
    # s = subclonal population, l = product state, j = allele
    overlapping_population_fraction_iws = np.array([(1. - purity_ws) * np.ones((num_walkers, num_subclonal_populations)), 
                                                   purity_ws * (1. - subclonal_cancer_cell_fraction_ws),
                                                   purity_ws * subclonal_cancer_cell_fraction_ws])
    population_weighted_copy_number_wslj = np.einsum('iws,lij->wslj', overlapping_population_fraction_iws, allelic_copy_number_product_states_lij)
    total_cr_wsl = np.sum(population_weighted_copy_number_wslj, axis=-1)
    log2cr_wsl = np.log2(total_cr_wsl / (cr_norm_w[:, np.newaxis, np.newaxis] + eps_model) + eps_cr)
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

def prior_logp_w(prior, ensemble_wp):
    num_walkers, num_parameters = np.shape(ensemble_wp)
    num_subclonal_populations = (num_parameters - 2) // 2
    logp_w = np.zeros(num_walkers)
    logp_w += prior.continuous_prior.subclonal_cancer_cell_fraction_weight_dirichlet.logpdf(ensemble_wp[:, :num_subclonal_populations].transpose())
    logp_w += prior.continuous_prior.subclonal_cancer_cell_fraction_beta.logpdf(ensemble_wp[:, num_subclonal_populations:-2]).sum(axis=1)
    logp_w += prior.continuous_prior.purity_beta.logpdf(ensemble_wp[:, -2])
    logp_w += prior.continuous_prior.cr_norm_lognorm.logpdf(ensemble_wp[:, -1])
    return logp_w
  
def logp_w(transformed_ensemble_wP, prior, data, tumor_ploidy_samples=[]):
    num_walkers, num_transformed_parameters = np.shape(transformed_ensemble_wP)
    num_parameters = num_transformed_parameters + 1
    num_subclonal_populations = (num_parameters - 2) // 2
    
    is_in_transformed_boundary_w = \
          np.all(transformed_ensemble_wP[:, :-1] >= 0, axis=1) & np.all(transformed_ensemble_wP[:, :-1] <= 1, axis=1) & \
          (transformed_ensemble_wP[:, -1] > 0)
    
    subclonal_cancer_cell_fraction_weight_ws = to_simplex_w(transformed_ensemble_wP[:, :num_subclonal_populations - 1])
    subclonal_cancer_cell_fraction_ws = transformed_ensemble_wP[:, num_subclonal_populations - 1:-2]
    
#    # enforce CCF weight sorting
##    is_subclonal_cancer_cell_fraction_weight_sorted_w = np.all(np.diff(subclonal_cancer_cell_fraction_weight_ws, axis=1) <= 0, axis=1)
#    subclonal_cancer_cell_fraction_weight_ws = np.sort(subclonal_cancer_cell_fraction_weight_ws, axis=-1)[:, ::-1]
    # enforce CCF sorting
#    is_subclonal_cancer_cell_fraction_weight_sorted_w = np.all(np.diff(subclonal_cancer_cell_fraction_weight_ws, axis=1) <= 0, axis=1)
    subclonal_cancer_cell_fraction_ws = np.sort(subclonal_cancer_cell_fraction_ws, axis=-1)[:, ::-1]
    # TODO sort weights and fractions accordingly
    
    is_allowed_w = is_in_transformed_boundary_w #& is_subclonal_cancer_cell_fraction_sorted_w
    
    ensemble_wp = np.zeros((num_walkers, num_parameters))
    ensemble_wp[:, :num_subclonal_populations] = subclonal_cancer_cell_fraction_weight_ws
    ensemble_wp[:, num_subclonal_populations:-2] = transformed_ensemble_wP[:, num_subclonal_populations - 1:-2]
    ensemble_wp[:, -2] = transformed_ensemble_wP[:, -2]
    ensemble_wp[:, -1] = transformed_ensemble_wP[:, -1]
    # modify CCF weights of disallowed walkers to ensure valid arguments for scipy.stats.dirichlet.logpdf
    ensemble_wp[~is_allowed_w, :num_subclonal_populations] = eps_model
    ensemble_wp[~is_allowed_w, 0] = 1.
    
    # prior
    logp_w = np.choose(is_allowed_w, [-np.inf, 0.])
    logp_w += prior_logp_w(prior, ensemble_wp)

    if np.all(logp_w == -np.inf):
        return logp_w
      
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
    ploidy_w = np.sum(total_cr_wk * data.length_k, axis=-1) / np.sum(data.length_k)
    logp_w += prior.continuous_prior.cr_norm_constraint_norm.logpdf(ensemble_wp[:, -1] - ploidy_w)
    
    normal_ploidy_wj = np.einsum('wjk,k->wj', map_copy_number_wijk[:, 0, :, :], data.length_k) / np.sum(data.length_k)
    logp_w += np.sum(scipy.stats.lognorm(s=0.01, scale=1.).logpdf(normal_ploidy_wj), axis=-1)
    tumor_ploidy_wj = np.einsum('wjk,k->wj', (1. - subclonal_cancer_cell_fraction_wk[:, np.newaxis, :]) * map_copy_number_wijk[:, 1, :, :]
                                                 + subclonal_cancer_cell_fraction_wk[:, np.newaxis, :] * map_copy_number_wijk[:, 2, :, :], 
                                data.length_k) / np.sum(data.length_k)
    tumor_ploidy_samples.extend(np.sum(tumor_ploidy_wj, axis=-1))
    
    return logp_w
    
def run_mcmc_vectorized(inference_config, prior, data, logp_w):
    num_parameters, num_walkers, num_burn_in, num_samples, emcee_move = inference_config
    num_transformed_parameters = num_parameters - 1
    num_subclonal_populations = (num_parameters - 2) // 2
    np.random.seed(0)
    tumor_ploidy_samples = []
    
    # initialize walkers
    transformed_ensemble_wP = np.random.rand(num_walkers, num_transformed_parameters)
    transformed_ensemble_wP[:, -2] = prior.continuous_prior.purity_beta.rvs(num_walkers)
    transformed_ensemble_wP[:, -1] = prior.continuous_prior.cr_norm_lognorm.rvs(num_walkers)
        
    # enforce CCF sorting
    subclonal_cancer_cell_fraction_weight_ws = to_simplex_w(transformed_ensemble_wP[:, :num_subclonal_populations - 1])
    sorted_subclonal_cancer_cell_fraction_weight_ws = np.sort(subclonal_cancer_cell_fraction_weight_ws, axis=-1)[:, ::-1]
    transformed_ensemble_wP[:, :num_subclonal_populations - 1] = to_hypercube_w(sorted_subclonal_cancer_cell_fraction_weight_ws)

    sampler = emcee.EnsembleSampler(num_walkers, num_transformed_parameters, logp_w, args=(prior, data, tumor_ploidy_samples), moves=emcee_move, vectorize=True)
    # burn-in sampling
    for _, _, _ in tqdm.tqdm(sampler.sample(transformed_ensemble_wP, iterations=num_burn_in), total=num_burn_in):
        pass
    sampler.reset()
    tumor_ploidy_samples.clear()
    # sampling
    for _, _, _ in tqdm.tqdm(sampler.sample(transformed_ensemble_wP, iterations=num_samples), total=num_samples):
        pass
    
    acceptance_fraction = np.mean(sampler.acceptance_fraction)
    print('Acceptance fraction:', acceptance_fraction)
    
    transformed_parameter_samples = sampler.get_chain(flat=True)
    logp_samples = sampler.get_log_prob(flat=True)
    
    # postprocess samples
    subclonal_cancer_cell_fraction_weight_s_samples = to_simplex_w(transformed_parameter_samples[:, :num_subclonal_populations - 1])
    subclonal_cancer_cell_fraction_weight_s_samples = np.sort(subclonal_cancer_cell_fraction_weight_s_samples, axis=1)[:, ::-1]
    subclonal_cancer_cell_fraction_s_samples = transformed_parameter_samples[:, num_subclonal_populations - 1:-2]
    purity_samples = transformed_parameter_samples[:, -2][:, np.newaxis]
    cr_norm_samples = transformed_parameter_samples[:, -1][:, np.newaxis]
    parameter_samples = np.hstack([subclonal_cancer_cell_fraction_weight_s_samples, subclonal_cancer_cell_fraction_s_samples, purity_samples, cr_norm_samples])
    tumor_ploidy_samples = np.array(tumor_ploidy_samples[num_walkers:]) # remove initial ensemble sample
    
    return sampler, parameter_samples, logp_samples, tumor_ploidy_samples
    
Parameters = namedtuple('Parameters', ['subclonal_cancer_cell_fraction_weight_s',
                                       'subclonal_cancer_cell_fraction_s',
                                       'purity',
                                       'cr_norm'])
    
DiscreteParameters = namedtuple('DiscreteParameters', ['copy_number_ijk',
                                                       'z_sk'])
def calculate_map_discrete_parameters(parameters, discrete_prior, data, use_marginalization_states):
    subclonal_cancer_cell_fraction_weight_s, subclonal_cancer_cell_fraction_s, purity, cr_norm = parameters
    allelic_copy_number_product_states_lij = data.marginalization_allelic_copy_number_product_states_lij if use_marginalization_states else data.allelic_copy_number_product_states_lij
    product_state_log_prior_kl = discrete_prior.marginalization_product_state_log_prior_kl if use_marginalization_states else discrete_prior.product_state_log_prior_kl
    
    ensemble_wp = np.array([list(subclonal_cancer_cell_fraction_weight_s) + list(subclonal_cancer_cell_fraction_s) + [purity] + [cr_norm]])
    likelihood_logp_ksl = calculate_likelihood_logp_wksl(ensemble_wp, data, use_marginalization_states)[0]
    logp_ksl = product_state_log_prior_kl[:, np.newaxis, :] + likelihood_logp_ksl
    
    flat_argmax_idx = np.argmax(np.reshape(logp_ksl, (logp_ksl.shape[0], -1)), axis=-1)
    map_s_k, map_l_k = np.unravel_index(flat_argmax_idx, logp_ksl.shape[-2:])
    
    map_copy_number_ijk = np.transpose(allelic_copy_number_product_states_lij[map_l_k], axes=(1, 2, 0))
    map_z_sk = map_s_k
    
    return DiscreteParameters(map_copy_number_ijk, map_z_sk)

#===============================================================================
# plotting
#===============================================================================

def plot_corner(parameter_samples, output_path, output_prefix, show=True):
    num_subclonal_populations = (parameter_samples.shape[1] - 2) // 2

    fig = corner.corner(parameter_samples, 
                        labels=['$\pi_{}$'.format(i) for i in range(num_subclonal_populations)] + \
                               ['$\phi_{}$'.format(i) for i in range(num_subclonal_populations)] + \
                               ['purity', 'cr_norm'],
                        quantiles=[0.1, 0.5, 0.9], show_titles=True)
    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.th.corner.png'))
    if show:
        plt.show()
    plt.close()
    
def big_end(data, k, min_length=0.005):
    return max(data.end_k[k], data.start_k[k] + min_length)

def plot_fit(axs, parameters, discrete_parameters, data):
    def make_cr_rectangle(data, k):
        start = data.start_k[k]
        end = big_end(data, k)
        cr_10 = 2**data.modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_10']
        cr_90 = 2**data.modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_90']
        return Rectangle([start, cr_10], end - start, cr_90 - cr_10)
        
    def make_maf_rectangle(data, k):
        start = data.start_k[k]
        end = big_end(data, k)
        maf_10 = data.modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_10']
        maf_90 = data.modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_90']
        return Rectangle([start, maf_10], end - start, maf_90 - maf_10)

    assert len(axs) == 2
        
    subclonal_cancer_cell_fraction_weight_s, subclonal_cancer_cell_fraction_s, purity, cr_norm = parameters
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
    
    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[0].set_ylim([0, 3])
    axs[1].set_ylim([0, 0.5])
    axs[0].set_ylabel('copy ratio')
    axs[1].set_ylabel('minor-allele fraction')
    
    axs[0].vlines(data.contig_ends, 0, 4, color='grey', linestyle='dashed', alpha=0.5)
    axs[1].vlines(data.contig_ends, 0, 0.5, color='grey', linestyle='dashed', alpha=0.5)
    
    pc = PatchCollection([make_cr_rectangle(data, k) for k in range(num_segments)],
                         color='r', alpha=0.25)
    axs[0].add_collection(pc)
    lc = LineCollection([[[data.start_k[k], 2**data.modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_50']], [big_end(data, k), 2**data.modeled_segments.iloc[k]['LOG2_COPY_RATIO_POSTERIOR_50']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.75)
    axs[0].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 2**log2cr_k[k]], [big_end(data, k), 2**log2cr_k[k]]] for k in range(num_segments)],
                        color='b', lw=2, alpha=0.75)
    axs[0].add_collection(lc)
    
    pc = PatchCollection([make_maf_rectangle(data, k) for k in range(num_segments)],
                         color='r', alpha=0.25)
    axs[1].add_collection(pc)
    lc = LineCollection([[[data.start_k[k], data.modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_50']], [big_end(data, k), data.modeled_segments.iloc[k]['MINOR_ALLELE_FRACTION_POSTERIOR_50']]] for k in range(num_segments)],
                        color='r', lw=2, alpha=0.75, label='data')
    axs[1].add_collection(lc)
    lc = LineCollection([[[data.start_k[k], maf_k[k]], [big_end(data, k), maf_k[k]]] for k in range(num_segments)],
                        color='b', lw=2, alpha=0.75, label='model')
    axs[1].add_collection(lc)
    
    axs[1].legend(loc='lower right', bbox_to_anchor= (1.08, 0.))
    
def plot_copy_number_ijk_samples(axs, copy_number_ijk_samples, data, allelic_copy_number_states):
    assert len(axs) == 3
    
    num_samples, num_populations, num_alleles, num_segments = np.shape(copy_number_ijk_samples)
    
    y_max = allelic_copy_number_states[-1] + 1
    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[2].set_xticks([])
    axs[0].set_ylabel('normal allelic CN')
    axs[1].set_ylabel('clonal allelic CN')
    axs[2].set_ylabel('subclonal allelic CN')
    axs[0].set_ylim([0, y_max])
    axs[1].set_ylim([0, y_max])
    axs[2].set_ylim([0, y_max])
    
    allele_rgba = [np.array(matplotlib.colors.to_rgba(c)) for c in ['b', 'g']]
    
    for k in range(num_segments):
        copy_number_ij_states, counts = np.unique(copy_number_ijk_samples[:, :, :, k], axis=0, return_counts=True)
        normalized_counts = counts / np.sum(counts)
        num_states = len(copy_number_ij_states)
        
        for i in range(num_populations):
            axs[i].vlines(data.contig_ends, 0, y_max, color='grey', linestyle='dashed', alpha=0.5)
            for j in range(num_alleles):
                colors = np.tile(allele_rgba[j], (num_states, 1))
                colors[:, 3] *= normalized_counts
                lc = LineCollection([[[data.start_k[k], copy_number_ij_state[i][j]], [big_end(data, k), copy_number_ij_state[i][j]]] for copy_number_ij_state in copy_number_ij_states],
                                    color=colors, lw=4, alpha=0.5)
                axs[i].add_collection(lc)
    
def plot_subclonal_diagram(ax, parameters, discrete_parameters, data, normal_allelic_copy_number_state):
    copy_number_ijk, z_sk = discrete_parameters
    _, _, num_segments = np.shape(copy_number_ijk)
    subclonal_cancer_cell_fraction_k = parameters.subclonal_cancer_cell_fraction_s[z_sk]
    is_normal_k = np.all(copy_number_ijk == normal_allelic_copy_number_state, axis=(0, 1))
    is_clonal_k = np.all(copy_number_ijk[1, :, :] == copy_number_ijk[2, :, :], axis=0)
    
    ax.set_xticks([])
    ax.set_ylim([0, 1])
    ax.set_ylabel('subclonal CCF')
    
    ax.vlines(data.contig_ends, 0, 4, color='grey', linestyle='dashed', alpha=0.5)
    
    lc = LineCollection([[[data.start_k[k], 0], [big_end(data, k), 0]] for k in np.nonzero(is_normal_k)[0]],
                        color='b', lw=4, alpha=0.5, label='normal')
    ax.add_collection(lc)
    lc = LineCollection([[[data.start_k[k], 0], [big_end(data, k), 0]] for k in np.nonzero(~is_normal_k & is_clonal_k)[0]],
                        color='g', lw=4, alpha=0.5, label='clonal')
    ax.add_collection(lc)
    lc = LineCollection([[[data.start_k[k], subclonal_cancer_cell_fraction_k[k]], [big_end(data, k), subclonal_cancer_cell_fraction_k[k]]] for k in np.nonzero(~is_normal_k & ~is_clonal_k)[0]],
                        color='r', lw=4, alpha=0.5)
    ax.add_collection(lc)
    
    ax.legend(loc='lower right', bbox_to_anchor= (1.08, 0.))

def plot_all_genomic_results(parameters, discrete_parameters, copy_number_ijk_samples,
                             data, discrete_prior_config,
                             output_path, output_prefix, show=True):
    figs, axs = plt.subplots(nrows=6, sharex=True, figsize=(18, 6 * 3))
    plt.xlim([0, 1])

    plot_fit(axs[:2], parameters, discrete_parameters, data)
    plot_copy_number_ijk_samples(axs[2:5], copy_number_ijk_samples, data, discrete_prior_config.allelic_copy_number_states)
    plot_subclonal_diagram(axs[-1], parameters, discrete_parameters, data, discrete_prior_config.normal_allelic_copy_number_state)

    plt.suptitle(output_prefix)
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    plt.savefig(os.path.join(output_path, output_prefix + '.th.png'))
    if show:
        plt.show()
    plt.close()
    
#===============================================================================
# configs
#===============================================================================
    
# num_overlapping_populations and num_alleles are essentially fixed, due to assumptions made in discrete prior
DiscretePriorConfig = namedtuple('DiscretePriorConfig', ['num_overlapping_populations',
                                                         'num_alleles',
                                                         'allelic_copy_number_states',
                                                         'normal_allelic_copy_number_state',
                                                         'copy_number_event_prior_penalty',
                                                         'allelic_copy_number_change_prior_penalty',
                                                         'hom_del_prior_penalty',
                                                         'num_marginalization_product_states',
                                                         'normal_population_event_length_scale',
                                                         'tumor_population_event_length_scale'])
discrete_prior_config = DiscretePriorConfig(
    num_overlapping_populations = 3,
    num_alleles = 2,
    allelic_copy_number_states = np.arange(6 + 1),
    normal_allelic_copy_number_state = 1,
    copy_number_event_prior_penalty = 0.001,
    allelic_copy_number_change_prior_penalty = 0.001,
    hom_del_prior_penalty = 0.,
    num_marginalization_product_states = 250,
    normal_population_event_length_scale = 1E4,
    tumor_population_event_length_scale = 1E7)

global_discrete_prior = generate_label_ordered_product_states_and_log_prior(discrete_prior_config)

ContinuousPriorConfig = namedtuple('ContinuousPriorConfig', ['num_subclonal_populations',
                                                             'subclonal_cancer_cell_fraction_weight_alpha',
                                                             'subclonal_cancer_cell_fraction_a',
                                                             'subclonal_cancer_cell_fraction_b',
                                                             'purity_a',
                                                             'purity_b',
                                                             'cr_norm_s',
                                                             'cr_norm_scale',
                                                             'cr_norm_constraint_scale'])
continuous_prior_config = ContinuousPriorConfig(
    num_subclonal_populations = 4,
    subclonal_cancer_cell_fraction_weight_alpha = 1E-3,
    subclonal_cancer_cell_fraction_a = 1.,
    subclonal_cancer_cell_fraction_b = 10.,
    purity_a = 1.,
    purity_b = 10.,
    cr_norm_s = 0.1,
    cr_norm_scale = 2.,
    cr_norm_constraint_scale = np.sqrt(1E-3))

ContinuousPrior = namedtuple('ContinuousPrior', ['subclonal_cancer_cell_fraction_weight_dirichlet',
                                                 'subclonal_cancer_cell_fraction_beta',
                                                 'purity_beta',
                                                 'cr_norm_lognorm',
                                                 'cr_norm_constraint_norm'])
continuous_prior = ContinuousPrior(
    subclonal_cancer_cell_fraction_weight_dirichlet = scipy.stats.dirichlet(
        alpha=continuous_prior_config.subclonal_cancer_cell_fraction_weight_alpha * np.ones(continuous_prior_config.num_subclonal_populations)),
    subclonal_cancer_cell_fraction_beta = scipy.stats.beta(
        a=continuous_prior_config.subclonal_cancer_cell_fraction_a, 
        b=continuous_prior_config.subclonal_cancer_cell_fraction_b),
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
likelihood_config = LikelihoodConfig(t_degrees_of_freedom = 10)

InferenceConfig = namedtuple('InferenceConfig', ['num_parameters',
                                                 'num_walkers',
                                                 'num_burn_in',
                                                 'num_samples',
                                                 'emcee_move'])
inference_config = InferenceConfig(
    num_parameters = 2 * continuous_prior_config.num_subclonal_populations + 2,
    num_walkers = 128,
    num_burn_in = 200,
    num_samples = 50,
#    emcee_move = emcee.moves.WalkMove()
#    emcee_move = emcee.moves.StretchMove()
    emcee_move = emcee.moves.KDEMove()
#    emcee_move = emcee.moves.DEMove()
)

GlobalConfig = namedtuple('GlobalConfig', ['discrete_prior_config',
                                           'continuous_prior_config',
                                           'global_discrete_prior',
                                           'continuous_prior',
                                           'likelihood_config',
                                           'inference_config'])
global_config = GlobalConfig(
    discrete_prior_config = discrete_prior_config,
    continuous_prior_config = continuous_prior_config,
    global_discrete_prior = global_discrete_prior,
    likelihood_config = likelihood_config,
    continuous_prior = continuous_prior,
    inference_config = inference_config
)

#===============================================================================
# main
#===============================================================================

def run_th(modeled_segments_path, output_prefix, output_path, global_config, show=True):
    discrete_prior_config, continuous_prior_config, global_discrete_prior, continuous_prior, likelihood_config, inference_config = global_config
  
    # load modeled segments and sequence dictionary=============================
    modeled_segments = pd.read_csv(modeled_segments_path, sep='\t', comment='@', dtype={'CONTIG': str})
    sequence_dictionary = read_sequence_dictionary(modeled_segments_path)

    print('Output prefix:', output_prefix)
    num_segments = len(modeled_segments)
    print('Number of segments:', num_segments)
    
    # fit posteriors and apply segment-length weights to discrete priors========
    data, discrete_prior = prepare_data_and_discrete_prior(modeled_segments, sequence_dictionary,
                                                           global_discrete_prior, discrete_prior_config, likelihood_config)
    prior = Prior(continuous_prior, discrete_prior)
    with open(os.path.join(output_path, output_prefix + '.th.data.pkl'), 'wb') as f:
        pickle.dump(data, f)
    with open(os.path.join(output_path, output_prefix + '.th.prior.pkl'), 'wb') as f:
        pickle.dump(prior, f)
    
    # perform inference=========================================================
    
    sampler, parameter_samples, logp_samples, tumor_ploidy_samples = run_mcmc_vectorized(inference_config, prior, data, logp_w)

    np.save(os.path.join(output_path, output_prefix + '.th.samples.parameters.npy'), parameter_samples)
    np.save(os.path.join(output_path, output_prefix + '.th.samples.logp.npy'), logp_samples)
    np.save(os.path.join(output_path, output_prefix + '.th.samples.ploidy.npy'), tumor_ploidy_samples)
    
    # output corner plot========================================================
    
    plot_corner(parameter_samples, output_path, output_prefix, show)
    
    # calculate and output MAP parameters=======================================
    
    map_sample_index = np.argmax(logp_samples)
    map_parameters = Parameters(
        parameter_samples[map_sample_index, :continuous_prior_config.num_subclonal_populations],
        parameter_samples[map_sample_index, continuous_prior_config.num_subclonal_populations:-2],
        parameter_samples[map_sample_index, -2],
        parameter_samples[map_sample_index, -1])
    map_discrete_parameters = calculate_map_discrete_parameters(map_parameters, discrete_prior, data, use_marginalization_states=False)
#    map_tumor_ploidy = tumor_ploidy_samples[map_sample_index]
    
    print('MAP subclonal_cancer_cell_fraction_weight_s:', map_parameters.subclonal_cancer_cell_fraction_weight_s)
    print('MAP subclonal_cancer_cell_fraction_s:', map_parameters.subclonal_cancer_cell_fraction_s) 
    print('MAP purity:', map_parameters.purity)
    print('MAP cr_norm:', map_parameters.cr_norm)
#    print('MAP ploidy:', map_tumor_ploidy)

    map_result = tuple([map_parameters, map_discrete_parameters])
    with open(os.path.join(output_path, output_prefix + '.th.map.pkl'), 'wb') as f:
        pickle.dump(map_result, f)
    
    # output genomic plots======================================================

    plot_all_genomic_results(map_parameters, map_discrete_parameters, np.array([map_discrete_parameters.copy_number_ijk]),
                             data, discrete_prior_config, output_path, output_prefix, show)
    
    # output posterior table====================================================

    # TODO
    
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

    run_th(args.modeled_segments_path, args.output_prefix, args.output_path, global_config, show=False)

if __name__ == '__main__':
    main()
