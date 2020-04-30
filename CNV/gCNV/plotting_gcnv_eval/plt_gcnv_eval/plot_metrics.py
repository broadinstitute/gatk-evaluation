import os
import matplotlib
import numpy as np
import pandas as pd
import math
import vcf
from typing import List
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import io_plt
import constants


def _get_plot_full_path(output_dir: str, filename: str):
    if not os.path.exists(os.path.join(output_dir, constants.PLOT_OUTPUT_DIR_NAME)):
        os.mkdir(os.path.join(output_dir, constants.PLOT_OUTPUT_DIR_NAME))
    return os.path.join(output_dir, constants.PLOT_OUTPUT_DIR_NAME, filename)


def plot_number_of_events_distribution(output_dir: str, gcnv_segment_vcfs: list):
    variant_num_list = []
    for gcnv_vcf in gcnv_segment_vcfs:
        variant_num_list.append(_count_num_variants(gcnv_vcf))
    fig = plt.figure(figsize=(10, 6))
    plt.title("Event number distribution")
    plt.hist(variant_num_list, bins=100, edgecolor='black', linewidth=1)
    plt.xlabel("Number of segments")
    fig.savefig(_get_plot_full_path(output_dir, "event_number_distribution.png"), dpi=120)
    io_plt.log("Mean number of segments in each sample: %d" % (np.mean(np.array(variant_num_list))))


def _count_num_variants(vcf_file):
    num_vars = 0
    with open(vcf_file, 'r') as vcf_reader:
        for line in vcf_reader.readlines():
            if not line.startswith("#"):
                num_vars += 1
    return num_vars


def plot_quality_metric_distribution(output_dir: str, gcnv_segment_vcfs: list):
    qs_values = []
    qa_values = []

    for segment_vcf in gcnv_segment_vcfs:
        vcf_reader = vcf.Reader(open(segment_vcf, 'r'))
        sample_name = vcf_reader.samples[0]
        for record in vcf_reader:
            attributes = {'QS': int(record.genotype(sample_name)['QS']),
                          'QA': int(record.genotype(sample_name)['QA'])}
            qs_values.append(attributes['QS'])
            qa_values.append(attributes['QA'])

    fig = plt.figure(figsize=(10, 6))
    plt.title("QS score distribution")
    plt.hist(qs_values, bins=100, edgecolor='black', linewidth=1)
    plt.xlabel("QS score")
    fig.savefig(_get_plot_full_path(output_dir, "QS_distribution.png"), dpi=120)

    fig = plt.figure(figsize=(10, 6))
    plt.title("QA score distribution")
    plt.hist(qa_values, bins=100, edgecolor='black', linewidth=1)
    plt.xlabel("QA score")
    fig.savefig(_get_plot_full_path(output_dir, "QA_distribution.png"), dpi=120)


def plot_event_size_distribution(output_dir: str, gcnv_segment_vcfs: list):
    num_targets_events = []
    num_bases_events = []

    for segment_vcf in gcnv_segment_vcfs:
        vcf_reader = vcf.Reader(open(segment_vcf, 'r'))
        sample_name = vcf_reader.samples[0]
        for record in vcf_reader:
            num_targets_events.append(int(record.genotype(sample_name)['NP']))
            num_bases_events.append(int(record.INFO['END']) - int(record.POS))

    fig = plt.figure(figsize=(10, 6))
    plt.title("Event size distribution")
    plt.hist(num_targets_events, bins=100, edgecolor='black', linewidth=1)
    plt.xlabel("Number of exons")
    fig.savefig(_get_plot_full_path(output_dir, "event_size_distribution_exons.png"), dpi=120)

    fig = plt.figure(figsize=(10, 6))
    plt.title("Event size distribution")
    plt.hist(num_bases_events, bins=100, edgecolor='black', linewidth=1)
    plt.xlabel("Number of base pairs")
    fig.savefig(_get_plot_full_path(output_dir, "event_size_distribution_bp.png"), dpi=120)


def plot_ard_components(output_dir: str, models_dir: str, num_shards: int):
    num_plots_per_row = 4
    num_rows = math.ceil(num_shards/num_plots_per_row)
    fig = plt.figure(figsize=(num_plots_per_row * 5, num_rows * 5))
    for i in range(num_shards):
        ard_file = models_dir + "/shard-" + str(i) + "/mu_ard_u_log__.tsv"
        ard_df = pd.read_csv(open(ard_file, 'r'), comment="@", sep='\t', header=None)
        ard_vector = np.transpose(ard_df.as_matrix())[0]
        plt.subplot(num_rows, num_plots_per_row, i + 1)
        plt.plot(range(ard_vector.size), np.sort(ard_vector))
        plt.plot(range(ard_vector.size), np.zeros(ard_vector.size), 'r--')
        plt.ylim(-2, 2)
        plt.title("Shard " + str(i) + " ARD mean values")
    plt.tight_layout()
    fig.savefig(_get_plot_full_path(output_dir, "ard_components.png"), dpi=120)


def plot_precision_recall_curve(output_dir: str, precision_values: List[float], recall_values: List[float]):
    fig = plt.figure(figsize=(10, 10), dpi=120)
    plt.title("Precision-Recall curve")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.plot(recall_values, precision_values, 'r-')
    fig.savefig(_get_plot_full_path(output_dir, "precision_recall_curve.png"), dpi=120)

