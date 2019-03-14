#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import os
from typing import List

summary_table_file = 'cohorts_and_cases_summary.tsv'
entity_ids_file_suffix = '-entity_ids.txt'
read_count_paths_file_suffix = '-read_count_paths.txt'
training_tag = '-training'
case_tag = '-case'

def determine_cohorts_and_cases_and_write_results(output_dir: str,
                                                  entity_ids_file: List[str],
                                                  read_count_paths_file: List[str],
                                                  clustering_table_file: str,
                                                  training_blacklist_file: str,
                                                  number_of_training_samples_per_model: int):
    os.makedirs(output_dir, exist_ok=True)
    assert number_of_training_samples_per_model > 0, "Number of training samples per model must be positive."

    #load all files
    entity_ids = np.loadtxt(entity_ids_file, dtype=bytes).astype(str)
    read_count_paths = np.loadtxt(read_count_paths_file, dtype=bytes).astype(str)
    assert len(entity_ids) == len(read_count_paths), "Number of entity IDs and number of read count paths must be equal."
    clustering_table_df = pd.read_csv(clustering_table_file, sep='\t')
    cluster_names = clustering_table_df.columns[1:]
    training_blacklist = np.loadtxt(training_blacklist_file, dtype=bytes).astype(str)

    #construct summary table
    summary_df = pd.DataFrame(columns=['SAMPLE_NAME',
                                       'READ_COUNT_PATH',
                                       'CLUSTER_MEMBERSHIP',
                                       'IN_TRAINING_BLACKLIST',
                                       'IS_RESCUE_SAMPLE',
                                       'IS_TRAINING_SAMPLE'])
    summary_df['SAMPLE_NAME'] = entity_ids
    summary_df['READ_COUNT_PATH'] = read_count_paths
    summary_df['CLUSTER_MEMBERSHIP'] = clustering_table_df[cluster_names].idxmax(axis=1)        # use MAP responsibility
    summary_df['IN_TRAINING_BLACKLIST'] = summary_df['SAMPLE_NAME'].isin(training_blacklist)
    number_of_blacklisted_samples_per_cluster = summary_df.groupby('CLUSTER_MEMBERSHIP')['IN_TRAINING_BLACKLIST'].sum()
    number_of_available_samples_per_cluster = summary_df['CLUSTER_MEMBERSHIP'].value_counts() - number_of_blacklisted_samples_per_cluster
    rescue_clusters = number_of_available_samples_per_cluster[number_of_available_samples_per_cluster < number_of_training_samples_per_model].keys()
    summary_df['IS_RESCUE_SAMPLE'] = summary_df['CLUSTER_MEMBERSHIP'].isin(rescue_clusters)
    training_samples = summary_df[~summary_df['IN_TRAINING_BLACKLIST'] & ~summary_df['IS_RESCUE_SAMPLE']] \
                                 .groupby('CLUSTER_MEMBERSHIP') \
                                 .apply(lambda x: x.sample(number_of_training_samples_per_model))['SAMPLE_NAME']
    summary_df['IS_TRAINING_SAMPLE'] = summary_df['SAMPLE_NAME'].isin(training_samples)

    #sanity checks
    assert not any(summary_df['IN_TRAINING_BLACKLIST'] & summary_df['IS_TRAINING_SAMPLE'])
    assert not any(summary_df['IS_RESCUE_SAMPLE'] & summary_df['IS_TRAINING_SAMPLE'])

    #output entity id and read count path files for both training and case samples for each cluster
    for cluster_name, cluster_groupby in summary_df.groupby('CLUSTER_MEMBERSHIP'):
        if len(cluster_groupby) > 0:
            cluster_training_samples = cluster_groupby[cluster_groupby['IS_TRAINING_SAMPLE']]
            cluster_case_samples = cluster_groupby[~cluster_groupby['IS_TRAINING_SAMPLE']]

            np.savetxt(os.path.join(output_dir, cluster_name + training_tag + entity_ids_file_suffix),
                       cluster_training_samples['SAMPLE_NAME'], fmt='%s')
            np.savetxt(os.path.join(output_dir, cluster_name + case_tag + entity_ids_file_suffix),
                       cluster_case_samples['SAMPLE_NAME'], fmt='%s')

            np.savetxt(os.path.join(output_dir, cluster_name + training_tag + read_count_paths_file_suffix),
                       cluster_training_samples['READ_COUNT_PATH'], fmt='%s')
            np.savetxt(os.path.join(output_dir, cluster_name + case_tag + read_count_paths_file_suffix),
                       cluster_case_samples['READ_COUNT_PATH'], fmt='%s')

    #output summary table
    summary_df['IN_TRAINING_BLACKLIST'] = summary_df['IN_TRAINING_BLACKLIST'].astype(int)
    summary_df['IS_RESCUE_SAMPLE'] = summary_df['IS_RESCUE_SAMPLE'].astype(int)
    summary_df['IS_TRAINING_SAMPLE'] = summary_df['IS_TRAINING_SAMPLE'].astype(int)
    summary_df.to_csv(os.path.join(output_dir, summary_table_file), sep='\t', index=False)

def test_determine_cohorts_and_cases_and_write_results():
    #make fake data and files for testing
    np.random.seed(0)

    test_output_dir = 'output'
    test_entity_ids_file = 'test_files/entity_ids.txt'
    test_read_count_paths_file = 'test_files/read_count_paths.txt'
    test_clustering_table_file = 'test_files/clustering_table.tsv'
    test_training_blacklist_file = 'test_files/training_blacklist.txt'
    test_number_of_training_samples_per_model = 300

    os.makedirs('test_files', exist_ok=True)

    test_num_samples = 10000
    test_entity_ids = np.array(['sample_{:d}'.format(i)
                                for i in range(1, test_num_samples + 1)])
    np.savetxt(test_entity_ids_file, test_entity_ids, fmt='%s')

    test_read_count_paths = ['gs://fake-bucket/sample_{:d}.counts.tsv'.format(i)
                             for i in range(1, test_num_samples + 1)]
    np.savetxt(test_read_count_paths_file, test_read_count_paths, fmt='%s')

    test_num_clusters = 30
    test_dirichlet_alpha = 0.001
    test_around_num_digits = 3
    test_cluster_names = ['CLUSTER_{:d}'.format(i)
                          for i in range(1, test_num_clusters + 1)]
    test_clustering_df = pd.DataFrame(columns=['SAMPLE_NAME'] + test_cluster_names)
    test_unnorm_responsibilities = np.around(np.random.dirichlet(test_dirichlet_alpha * np.ones(test_num_clusters),
                                                                 size=test_num_samples),
                                             test_around_num_digits)
    test_responsibilities = test_unnorm_responsibilities / np.sum(test_unnorm_responsibilities, axis=1)[:, np.newaxis]
    test_clustering_df['SAMPLE_NAME'] = test_entity_ids
    test_clustering_df[test_cluster_names] = test_responsibilities
    test_clustering_df.to_csv(test_clustering_table_file, sep='\t', index=False)

    test_training_blacklist_fraction = 0.02
    test_training_blacklist = test_entity_ids[1 == np.random.choice([0, 1],
                                                                    p=[1 - test_training_blacklist_fraction, test_training_blacklist_fraction],
                                                                    size=test_num_samples)]
    np.savetxt(test_training_blacklist_file, test_training_blacklist, fmt='%s')

    #run on test files
    determine_cohorts_and_cases_and_write_results(test_output_dir,
                                                  test_entity_ids_file,
                                                  test_read_count_paths_file,
                                                  test_clustering_table_file,
                                                  test_training_blacklist_file,
                                                  test_number_of_training_samples_per_model)

def main():
    parser = argparse.ArgumentParser(
        description='For each cluster with a sufficient number of non-blacklisted training samples, outputs the following files to the output directory:\n'
                    + '\n'
                    + '1) A file containing entity IDs       for training samples (CLUSTER_NAME{0}{1})\n'.format(training_tag, entity_ids_file_suffix)
                    + '2) A file containing read count paths for training samples (CLUSTER_NAME{0}{1})\n'.format(training_tag, read_count_paths_file_suffix)
                    + '3) A file containing entity IDs       for case samples     (CLUSTER_NAME{0}{1})\n'.format(case_tag, entity_ids_file_suffix)
                    + '4) A file containing read count paths for case samples     (CLUSTER_NAME{0}{1})\n'.format(case_tag, read_count_paths_file_suffix)
                    + '\n'
                    + 'A table summarizing all results (cohorts_and_cases_summary.tsv) is also output.',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--output_dir',
                        type=str,
                        help='Output directory.')

    parser.add_argument('--entity_ids_file',
                        type=str,
                        help='File containing entity IDs corresponding to read count paths.')

    parser.add_argument('--read_count_paths_file',
                        type=str,
                        help='File containing read count paths (output of GATK CollectReadCounts) corresponding to entity IDs.')

    parser.add_argument('--clustering_table_file',
                        type=str,
                        help='Table of clustered samples (output of ClusterSamplesFromCoverage workflow).')

    parser.add_argument('--training_blacklist_file',
                        type=str,
                        help='File containing blacklist of entity IDs for samples not to be used for training models.')

    parser.add_argument('--number_of_training_samples_per_model',
                        type=int,
                        help='Number of samples used to train each model; samples in clusters with a number of ' +
                             'non-blacklisted samples less than this will be indicated by a IS_RESCUE tag in the summary table.')

    args = parser.parse_args()

    output_dir = args.output_dir
    entity_ids_file = args.entity_ids_file
    read_count_paths_file = args.read_count_paths_file
    clustering_table_file = args.clustering_table_file
    training_blacklist_file = args.training_blacklist_file
    number_of_training_samples_per_model = args.number_of_training_samples_per_model

    determine_cohorts_and_cases_and_write_results(output_dir,
                                                  entity_ids_file,
                                                  read_count_paths_file,
                                                  clustering_table_file,
                                                  training_blacklist_file,
                                                  number_of_training_samples_per_model)

if __name__ == '__main__':
    #test_determine_cohorts_and_cases_and_write_results()
    main()

