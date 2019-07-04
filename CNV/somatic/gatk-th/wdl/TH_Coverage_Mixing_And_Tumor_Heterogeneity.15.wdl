import "https://api.firecloud.org/ga4gh/v1/tools/slee:TH_Coverage_Mixing/versions/13/plain-WDL/descriptor" as CoverageMixingWorkflow
import "https://api.firecloud.org/ga4gh/v1/tools/slee:TH_Tumor_Heterogeneity_From_Coverage/versions/3/plain-WDL/descriptor" as THFromCoverageWorkflow

workflow CoverageMixingAndTHWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    File common_sites
    File intervals
    File? blacklist_intervals
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File read_count_pon
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    Array[Int] mixing_tumor_percentages
    File python_script
    File environment_yml
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts
    # Use as a last resort to increase the disk given to every task in case of ill behaving data
    Int? emergency_extra_disk

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? bin_length
    Int? mem_gb_for_preprocess_intervals

    ##############################################
    #### optional arguments for CollectCounts ####
    ##############################################
    String? collect_counts_format
    Int? mem_gb_for_collect_counts

    #####################################################
    #### optional arguments for CollectAllelicCounts ####
    #####################################################
    String? minimum_base_quality
    Int? mem_gb_for_collect_allelic_counts

    ##################################################
    #### optional arguments for DenoiseReadCounts ####
    ##################################################
    Int? number_of_eigensamples
    Int? mem_gb_for_denoise_read_counts

    ##############################################
    #### optional arguments for ModelSegments ####
    ##############################################
    Int? max_num_segments_per_chromosome
    Int? min_total_allele_count
    Int? min_total_allele_count_normal
    Float? genotyping_homozygous_log_ratio_threshold
    Float? genotyping_base_error_rate
    Float? kernel_variance_copy_ratio
    Float? kernel_variance_allele_fraction
    Float? kernel_scaling_allele_fraction
    Int? kernel_approximation_dimension
    Array[Int]+? window_sizes = [8, 16, 32, 64, 128, 256]
    Float? num_changepoints_penalty_factor
    Float? minor_allele_fraction_prior_alpha
    Int? num_samples_copy_ratio
    Int? num_burn_in_copy_ratio
    Int? num_samples_allele_fraction
    Int? num_burn_in_allele_fraction
    Float? smoothing_threshold_copy_ratio
    Float? smoothing_threshold_allele_fraction
    Int? max_num_smoothing_iterations
    Int? num_smoothing_iterations_per_fit
    Int? mem_gb_for_model_segments

    ######################################################
    #### optional arguments for CallCopyRatioSegments ####
    ######################################################
    Float? neutral_segment_copy_ratio_lower_bound
    Float? neutral_segment_copy_ratio_upper_bound
    Float? outlier_neutral_segment_copy_ratio_z_score_threshold
    Float? calling_copy_ratio_z_score_threshold
    Int? mem_gb_for_call_copy_ratio_segments

    #########################################
    #### optional arguments for plotting ####
    #########################################
    Int? minimum_contig_length
    Int? mem_gb_for_plotting

    call CoverageMixingWorkflow.CoverageMixingWorkflow as CoverageMixingWorkflow {
        input:
            common_sites = common_sites,
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta = ref_fasta,
            mixing_tumor_percentages = mixing_tumor_percentages,
            gatk_docker = gatk_docker,
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            emergency_extra_disk = emergency_extra_disk,
            padding = padding,
            bin_length = bin_length,
            mem_gb_for_preprocess_intervals = mem_gb_for_preprocess_intervals,
            collect_counts_format = collect_counts_format,
            mem_gb_for_collect_counts = mem_gb_for_collect_counts,
            minimum_base_quality = minimum_base_quality,
            mem_gb_for_collect_allelic_counts = mem_gb_for_collect_allelic_counts
    }

    scatter (mixing_index in range(length(CoverageMixingWorkflow.mixed_entity_ids))) {
        call THFromCoverageWorkflow.THFromCoverageWorkflow as THFromCoverageWorkflow {
            input:
                entity_id = CoverageMixingWorkflow.mixed_entity_ids[mixing_index],
                tumor_read_counts = CoverageMixingWorkflow.mixed_read_counts[mixing_index],
                tumor_allelic_counts = CoverageMixingWorkflow.mixed_allelic_counts[mixing_index],
                normal_allelic_counts = CoverageMixingWorkflow.normal_allelic_counts,
                read_count_pon = read_count_pon,
                ref_fasta_dict = ref_fasta_dict,
                python_script = python_script,
                environment_yml = environment_yml,
                gatk_docker = gatk_docker,
                gatk4_jar_override = gatk4_jar_override,
                preemptible_attempts = preemptible_attempts,
                emergency_extra_disk = emergency_extra_disk,
                number_of_eigensamples = number_of_eigensamples,
                mem_gb_for_denoise_read_counts = mem_gb_for_denoise_read_counts,
                max_num_segments_per_chromosome = max_num_segments_per_chromosome,
                min_total_allele_count = min_total_allele_count,
                min_total_allele_count_normal = min_total_allele_count_normal,
                genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
                genotyping_base_error_rate = genotyping_base_error_rate,
                kernel_variance_copy_ratio = kernel_variance_copy_ratio,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
                kernel_approximation_dimension = kernel_approximation_dimension,
                window_sizes = window_sizes,
                num_changepoints_penalty_factor = num_changepoints_penalty_factor,
                minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
                num_samples_copy_ratio = num_samples_copy_ratio,
                num_burn_in_copy_ratio = num_burn_in_copy_ratio,
                num_samples_allele_fraction = num_samples_allele_fraction,
                num_burn_in_allele_fraction = num_burn_in_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                max_num_smoothing_iterations = max_num_smoothing_iterations,
                num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
                mem_gb_for_model_segments = mem_gb_for_model_segments,
                neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
                neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
                outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
                calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
                mem_gb_for_call_copy_ratio_segments = mem_gb_for_call_copy_ratio_segments,
                minimum_contig_length = minimum_contig_length,
                mem_gb_for_plotting = mem_gb_for_plotting
        }
    }
}
