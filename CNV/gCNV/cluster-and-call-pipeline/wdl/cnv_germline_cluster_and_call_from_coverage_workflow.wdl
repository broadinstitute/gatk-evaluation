import "cluster_samples_from_coverage_workflow.wdl" as ClusterSamplesFromCoverageWorkflow
import "cnv_germline_cohort_from_coverage_workflow.wdl" as CNVGermlineCohortFromCoverageWorkflow
import "cnv_germline_case_scattered_from_coverage_workflow.wdl" as CNVGermlineCaseScatteredFromCoverageWorkflow

workflow CNVGermlineClusterAndCallFromCoverageWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    Array[String]+ entity_ids
    Array[String]+ read_count_paths
    Array[String]+ read_count_paths_for_clustering
    Int maximum_number_of_clusters
    File training_blacklist
    Int number_of_training_samples_per_model
    File intervals
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    Int num_intervals_per_scatter
    Int num_samples_per_scatter_block
    Int ref_copy_number_autosomal_contigs
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts

    #############################################################
    #### optional shared arguments for CNVGermline workflows ####
    #############################################################
    File? blacklist_intervals
    Float? ploidy_mapping_error_rate
    Float? ploidy_sample_psi_scale
    Float? gcnv_p_alt
    Float? gcnv_cnv_coherence_length
    Int? gcnv_max_copy_number
    Float? gcnv_mapping_error_rate
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    String? gcnv_copy_number_posterior_expectation_mode
    Int? gcnv_active_class_padding_hybrid_mode
    Array[String]? allosomal_contigs

    call ClusterSamplesFromCoverageWorkflow.ClusterSamplesFromCoverageWorkflow {
        input:
            entity_ids = entity_ids,
            read_count_files = read_count_paths_for_clustering,
            maximum_number_of_clusters = maximum_number_of_clusters,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    call DetermineCohortsAndCases {
        input:
            entity_ids = entity_ids,
            read_count_paths = read_count_paths,
            clustering_table = ClusterSamplesFromCoverageWorkflow.clustering_table,
            training_blacklist = training_blacklist,
            number_of_training_samples_per_model = number_of_training_samples_per_model,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    #process globs

    #call cohort mode

    #call case mode

    output {
        File clustering_table = ClusterSamplesFromCoverageWorkflow.clustering_table
        File sample_clusters_plot = ClusterSamplesFromCoverageWorkflow.sample_clusters_plot
    }
}

task DetermineCohortsAndCases {
    Array[String]+ entity_ids
    Array[String]+ read_count_paths
    File clustering_table
    File training_blacklist
    Int number_of_training_samples_per_model

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 3]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    command <<<
        set -e
        python <<EOF
            import argparse
            import os
            import pandas as pd
            import numpy as np
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
                                    help='File containing list of entity IDs corresponding to read count paths.')

                parser.add_argument('--read_count_paths_file',
                                    type=str,
                                    help='File containing list of read count paths (output of GATK CollectReadCounts) corresponding to entity IDs.')

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
                main()
        EOF \
            --output_dir . \
            --entity_ids_file ${write_lines(entity_ids)} \
            --read_count_paths_file ${write_lines(read_count_paths)} \
            --clustering_table_file ${clustering_table} \
            --training_blacklist_file ${training_blacklist} \
            --number_of_training_samples_per_model ${number_of_training_samples_per_model}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        #update this
    }
}