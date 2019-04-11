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
    Array[String] training_blacklist_entity_ids
    Int number_of_training_samples_per_model
    File intervals
    File contig_ploidy_priors
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    Int num_intervals_per_scatter
    Int num_samples_per_scatter_block
    Int ref_copy_number_autosomal_contigs
    String gatk_docker

    ############################
    #### optional arguments ####
    ############################
    File? clustering_prior_table        # currently not used
    Int? mem_gb_for_cluster_samples

    #############################################################
    #### shared optional arguments for CNVGermline workflows ####
    #############################################################

    ##################################
    #### optional basic arguments ####
    ##################################
    File? blacklist_intervals
    # If true, AnnotateIntervals will be run to create GC annotations and explicit
    # GC correction will be performed by the model generated by
    Boolean? do_explicit_gc_correction
    File? gatk4_jar_override
    Int? preemptible_attempts

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? bin_length

    ##################################################
    #### optional arguments for AnnotateIntervals ####
    ##################################################
    File? mappability_track_bed
    File? mappability_track_bed_idx
    File? segmental_duplication_track_bed
    File? segmental_duplication_track_bed_idx
    Int? feature_query_lookahead
    Int? mem_gb_for_annotate_intervals

    #################################################
    #### optional arguments for FilterIntervals ####
    ################################################
    File? blacklist_intervals_for_filter_intervals
    Float? minimum_gc_content
    Float? maximum_gc_content
    Float? minimum_mappability
    Float? maximum_mappability
    Float? minimum_segmental_duplication_content
    Float? maximum_segmental_duplication_content
    Int? low_count_filter_count_threshold
    Float? low_count_filter_percentage_of_samples
    Float? extreme_count_filter_minimum_percentile
    Float? extreme_count_filter_maximum_percentile
    Float? extreme_count_filter_percentage_of_samples
    Int? mem_gb_for_filter_intervals

    ##############################################################
    #### optional arguments for DetermineGermlineContigPloidy ####
    ##############################################################
    Float? ploidy_mean_bias_standard_deviation
    Float? ploidy_mapping_error_rate
    Float? ploidy_global_psi_scale
    Float? ploidy_sample_psi_scale
    Int? mem_gb_for_determine_germline_contig_ploidy
    Int? cpu_for_determine_germline_contig_ploidy
    Int? disk_for_determine_germline_contig_ploidy

    ##################################################
    #### optional arguments for GermlineCNVCaller ####
    ##################################################
    Float? gcnv_p_alt
    Float? gcnv_p_active
    Float? gcnv_cnv_coherence_length
    Float? gcnv_class_coherence_length
    Int? gcnv_max_copy_number
    Int? mem_gb_for_germline_cnv_caller
    Int? cpu_for_germline_cnv_caller
    Int? disk_for_germline_cnv_caller

    # optional arguments for germline CNV denoising model
    Int? gcnv_max_bias_factors
    Float? gcnv_mapping_error_rate
    Float? gcnv_interval_psi_scale
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    Float? gcnv_log_mean_bias_standard_deviation
    Float? gcnv_init_ard_rel_unexplained_variance
    Int? gcnv_num_gc_bins
    Float? gcnv_gc_curve_standard_deviation
    String? gcnv_copy_number_posterior_expectation_mode
    Boolean? gcnv_enable_bias_factors
    Int? gcnv_active_class_padding_hybrid_mode

    # optional arguments for Hybrid ADVI
    Float? gcnv_learning_rate
    Float? gcnv_adamax_beta_1
    Float? gcnv_adamax_beta_2
    Int? gcnv_log_emission_samples_per_round
    Float? gcnv_log_emission_sampling_median_rel_error
    Int? gcnv_log_emission_sampling_rounds
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_max_training_epochs
    Float? gcnv_initial_temperature
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_convergence_snr_averaging_window
    Float? gcnv_convergence_snr_trigger_threshold
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_max_calling_iters
    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_caller_internal_admixing_rate
    Float? gcnv_caller_external_admixing_rate
    Boolean? gcnv_disable_annealing

    ############################################################
    #### optional arguments for PostprocessGermlineCNVCalls ####
    ############################################################
    Array[String]? allosomal_contigs

    call ClusterSamplesFromCoverageWorkflow.ClusterSamplesFromCoverageWorkflow {
        input:
            entity_ids = entity_ids,
            read_count_files = read_count_paths_for_clustering,
            maximum_number_of_clusters = maximum_number_of_clusters,
            clustering_prior_table = clustering_prior_table,
            gatk_docker = gatk_docker,
            mem_gb_for_cluster_samples = mem_gb_for_cluster_samples,
            preemptible_attempts = preemptible_attempts
    }

    call DetermineCohortsAndCases {
        input:
            entity_ids = entity_ids,
            read_count_paths = read_count_paths,
            clustering_table = ClusterSamplesFromCoverageWorkflow.clustering_table,
            training_blacklist_entity_ids = training_blacklist_entity_ids,
            number_of_training_samples_per_model = number_of_training_samples_per_model,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    #process globs
    Array[File] training_entity_ids_per_cluster_ = DetermineCohortsAndCases.training_entity_ids_per_cluster
    Array[File] training_read_count_paths_per_cluster_ = DetermineCohortsAndCases.training_read_count_paths_per_cluster
    Array[File] case_entity_ids_per_cluster_ = DetermineCohortsAndCases.case_entity_ids_per_cluster
    Array[File] case_read_count_paths_per_cluster_ = DetermineCohortsAndCases.case_read_count_paths_per_cluster
    Int num_clusters = length(training_entity_ids_per_cluster_)

    scatter (cluster_index in range(num_clusters)) {
        #call cohort mode
        Array[String] training_entity_ids = read_lines(training_entity_ids_per_cluster_[cluster_index])
        if (len(training_entity_ids) > 0) {
            call CNVGermlineCohortFromCoverageWorkflow.CNVGermlineCohortFromCoverageWorkflow as CohortTraining {
                input:
                    intervals = intervals,
                    blacklist_intervals = blacklist_intervals,
                    entity_ids = training_entity_ids,
                    read_count_files = read_lines(training_read_count_paths_per_cluster_[cluster_index]),
                    cohort_entity_id = basename(training_entity_ids_per_cluster_[cluster_index], "-training-entity_ids.txt"),
                    contig_ploidy_priors = contig_ploidy_priors,
                    num_intervals_per_scatter = num_intervals_per_scatter,
                    ref_fasta_dict = ref_fasta_dict,
                    ref_fasta_fai = ref_fasta_fai,
                    ref_fasta = ref_fasta,
                    gatk_docker = gatk_docker,
                    do_explicit_gc_correction = do_explicit_gc_correction,
                    gatk4_jar_override = gatk4_jar_override,
                    preemptible_attempts = preemptible_attempts,
                    padding = padding,
                    bin_length = bin_length,
                    mappability_track_bed = mappability_track_bed,
                    mappability_track_bed_idx = mappability_track_bed_idx,
                    segmental_duplication_track_bed = segmental_duplication_track_bed,
                    segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
                    feature_query_lookahead = feature_query_lookahead,
                    mem_gb_for_annotate_intervals = mem_gb_for_annotate_intervals,
                    blacklist_intervals_for_filter_intervals = blacklist_intervals_for_filter_intervals,
                    minimum_gc_content = minimum_gc_content,
                    maximum_gc_content = maximum_gc_content,
                    minimum_mappability = minimum_mappability,
                    maximum_mappability = maximum_mappability,
                    minimum_segmental_duplication_content = minimum_segmental_duplication_content,
                    maximum_segmental_duplication_content = maximum_segmental_duplication_content,
                    low_count_filter_count_threshold = low_count_filter_count_threshold,
                    low_count_filter_percentage_of_samples = low_count_filter_percentage_of_samples,
                    extreme_count_filter_minimum_percentile = extreme_count_filter_minimum_percentile,
                    extreme_count_filter_maximum_percentile = extreme_count_filter_maximum_percentile,
                    extreme_count_filter_percentage_of_samples = extreme_count_filter_percentage_of_samples,
                    mem_gb_for_filter_intervals = mem_gb_for_filter_intervals,
                    ploidy_mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
                    ploidy_mapping_error_rate = ploidy_mapping_error_rate,
                    ploidy_global_psi_scale = ploidy_global_psi_scale,
                    ploidy_sample_psi_scale = ploidy_sample_psi_scale,
                    mem_gb_for_determine_germline_contig_ploidy = mem_gb_for_determine_germline_contig_ploidy,
                    cpu_for_determine_germline_contig_ploidy = cpu_for_determine_germline_contig_ploidy,
                    disk_for_determine_germline_contig_ploidy = disk_for_determine_germline_contig_ploidy,
                    gcnv_p_alt = gcnv_p_alt,
                    gcnv_p_active = gcnv_p_active,
                    gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
                    gcnv_class_coherence_length = gcnv_class_coherence_length,
                    gcnv_max_copy_number = gcnv_max_copy_number,
                    mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
                    cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
                    disk_for_germline_cnv_caller = disk_for_germline_cnv_caller,
                    gcnv_max_bias_factors = gcnv_max_bias_factors,
                    gcnv_mapping_error_rate = gcnv_mapping_error_rate,
                    gcnv_interval_psi_scale = gcnv_interval_psi_scale,
                    gcnv_sample_psi_scale = gcnv_sample_psi_scale,
                    gcnv_depth_correction_tau = gcnv_depth_correction_tau,
                    gcnv_log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
                    gcnv_init_ard_rel_unexplained_variance = gcnv_init_ard_rel_unexplained_variance,
                    gcnv_num_gc_bins = gcnv_num_gc_bins,
                    gcnv_gc_curve_standard_deviation = gcnv_gc_curve_standard_deviation,
                    gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                    gcnv_enable_bias_factors = gcnv_enable_bias_factors,
                    gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                    gcnv_learning_rate = gcnv_learning_rate,
                    gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
                    gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
                    gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                    gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                    gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                    gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                    gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                    gcnv_min_training_epochs = gcnv_min_training_epochs,
                    gcnv_max_training_epochs = gcnv_max_training_epochs,
                    gcnv_initial_temperature = gcnv_initial_temperature,
                    gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                    gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                    gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                    gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                    gcnv_max_calling_iters = gcnv_max_calling_iters,
                    gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                    gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                    gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                    gcnv_disable_annealing = gcnv_disable_annealing,
                    ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                    allosomal_contigs = allosomal_contigs
            }
        }

        #call scattered-case mode
        Array[String] case_entity_ids = read_lines(case_entity_ids_per_cluster_[cluster_index])
        if (len(case_entity_ids) > 0) {
            call CNVGermlineCaseScatteredFromCoverageWorkflow.CNVGermlineCaseScatteredFromCoverageWorkflow as ScatteredCase {
                input:
                    blacklist_intervals = blacklist_intervals,
                    filtered_intervals = CohortTraining.filtered_intervals[cluster_index],
                    entity_ids = read_lines(case_entity_ids_per_cluster_[cluster_index]),
                    read_count_files = read_lines(case_read_count_paths_per_cluster_[cluster_index]),
                    contig_ploidy_model_tar = CohortTraining.contig_ploidy_model_tar[cluster_index],
                    gcnv_model_tars = CohortTraining.gcnv_model_tars[cluster_index],
                    num_intervals_per_scatter = num_intervals_per_scatter,
                    ref_fasta_dict = ref_fasta_dict,
                    ref_fasta_fai = ref_fasta_fai,
                    ref_fasta = ref_fasta,
                    gatk_docker = gatk_docker,
                    num_samples_per_scatter_block = num_samples_per_scatter_block,
                    gatk4_jar_override = gatk4_jar_override,
                    preemptible_attempts = preemptible_attempts,
                    ploidy_mapping_error_rate = ploidy_mapping_error_rate,
                    ploidy_sample_psi_scale = ploidy_sample_psi_scale,
                    mem_gb_for_determine_germline_contig_ploidy = mem_gb_for_determine_germline_contig_ploidy,
                    cpu_for_determine_germline_contig_ploidy = cpu_for_determine_germline_contig_ploidy,
                    disk_for_determine_germline_contig_ploidy = disk_for_determine_germline_contig_ploidy,
                    gcnv_p_alt = gcnv_p_alt,
                    gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
                    gcnv_max_copy_number = gcnv_max_copy_number,
                    mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
                    cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
                    disk_for_germline_cnv_caller = disk_for_germline_cnv_caller,
                    gcnv_mapping_error_rate = gcnv_mapping_error_rate,
                    gcnv_sample_psi_scale = gcnv_sample_psi_scale,
                    gcnv_depth_correction_tau = gcnv_depth_correction_tau,
                    gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                    gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                    gcnv_learning_rate = gcnv_learning_rate,
                    gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
                    gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
                    gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                    gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                    gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                    gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                    gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                    gcnv_min_training_epochs = gcnv_min_training_epochs,
                    gcnv_max_training_epochs = gcnv_max_training_epochs,
                    gcnv_initial_temperature = gcnv_initial_temperature,
                    gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                    gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                    gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                    gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                    gcnv_max_calling_iters = gcnv_max_calling_iters,
                    gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                    gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                    gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                    gcnv_disable_annealing = gcnv_disable_annealing,
                    ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                    allosomal_contigs = allosomal_contigs
            }
        }
    }

    output {
        File clustering_table = ClusterSamplesFromCoverageWorkflow.clustering_table
        File sample_clusters_plot = ClusterSamplesFromCoverageWorkflow.sample_clusters_plot

        Array[File] training_entity_ids_per_cluster = training_entity_ids_per_cluster_
        Array[File] training_read_count_paths_per_cluster = training_read_count_paths_per_cluster_
        Array[File] case_entity_ids_per_cluster = case_entity_ids_per_cluster_
        Array[File] case_read_count_paths_per_cluster = case_read_count_paths_per_cluster_
        File summary_table = DetermineCohortsAndCases.summary_table

        File preprocessed_intervals = CohortTraining.preprocessed_intervals[0]
        File? annotated_intervals = CohortTraining.annotated_intervals[0]
        Array[File] training_filtered_intervals_per_cluster = CohortTraining.filtered_intervals
        Array[File] training_contig_ploidy_model_tar_per_cluster = CohortTraining.contig_ploidy_model_tar
        Array[File] training_contig_ploidy_calls_tar_per_cluster = CohortTraining.contig_ploidy_calls_tar
        Array[Array[File]] training_gcnv_model_tars_per_cluster = CohortTraining.gcnv_model_tars
        Array[Array[Array[File]]] training_gcnv_calls_tars_per_cluster = CohortTraining.gcnv_calls_tars
        Array[Array[File]] training_gcnv_tracking_tars_per_cluster = CohortTraining.gcnv_tracking_tars
        Array[Array[File]] training_genotyped_intervals_vcfs_per_cluster = CohortTraining.genotyped_intervals_vcfs
        Array[Array[File]] training_genotyped_segments_vcfs_per_cluster = CohortTraining.genotyped_segments_vcfs

        Array[Array[File]] scattered_case_contig_ploidy_calls_tars_per_cluster = ScatteredCase.contig_ploidy_calls_tars
        Array[Array[Array[Array[File]]]] scattered_case_gcnv_calls_tars_per_cluster = ScatteredCase.gcnv_calls_tars
        Array[Array[Array[File]]] scattered_case_gcnv_tracking_tars_per_cluster = ScatteredCase.gcnv_tracking_tars
        Array[Array[Array[File]]] scattered_case_genotyped_intervals_vcfs_per_cluster = ScatteredCase.genotyped_intervals_vcfs
        Array[Array[Array[File]]] scattered_case_genotyped_segments_vcfs_per_cluster = ScatteredCase.genotyped_segments_vcfs
    }
}

task DetermineCohortsAndCases {
    Array[String]+ entity_ids
    Array[String]+ read_count_paths
    File clustering_table
    Array[String] training_blacklist_entity_ids
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
        python - --output_dir . \
                 --entity_ids_file ${write_lines(entity_ids)} \
                 --read_count_paths_file ${write_lines(read_count_paths)} \
                 --clustering_table_file ${clustering_table} \
                 --training_blacklist_file ${write_lines(training_blacklist_entity_ids)} \
                 --number_of_training_samples_per_model ${number_of_training_samples_per_model} \
                 <<-'EOF'
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
                            + '1) A file containing entity IDs       for training samples (<CLUSTER_NAME>{0}{1})\n'.format(training_tag, entity_ids_file_suffix)
                            + '2) A file containing read count paths for training samples (<CLUSTER_NAME>{0}{1})\n'.format(training_tag, read_count_paths_file_suffix)
                            + '3) A file containing entity IDs       for case samples     (<CLUSTER_NAME>{0}{1})\n'.format(case_tag, entity_ids_file_suffix)
                            + '4) A file containing read count paths for case samples     (<CLUSTER_NAME>{0}{1})\n'.format(case_tag, read_count_paths_file_suffix)
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
        EOF
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] training_entity_ids_per_cluster = glob('*-training-entity_ids.txt')
        Array[File] training_read_count_paths_per_cluster = glob('*-training-read_count_paths.txt')
        Array[File] case_entity_ids_per_cluster = glob('*-case-entity_ids.txt')
        Array[File] case_read_count_paths_per_cluster = glob('*-case-read_count_paths.txt')
        File summary_table = 'cohorts_and_cases_summary.tsv'
    }
}