version 1.0

import "cnv_validation.wdl" as cnv_validation

workflow MultiCNVValidation {
    input {
        File ice_intervals
        File wgs_intervals
        File common_sites
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai

        File ice_read_count_pon
        File wgs_read_count_pon

        String gatk_docker

        File? gatk4_jar_override

        # In case you are using a branch that is making changes to the evaluation tools (e.g. CombineSegmentBreakpoints)
        File? gatk4_jar_override_evaluation
        Int? wgs_bin_length
        Int wxs_bin_length = 0

        String? group_id
        String group_id_final = select_first([group_id, "input_segs"])

        String eval_docker

        ### Validation parameters
        Float? num_changepoints_penalty_factor_normal
        Float? kernel_variance_allele_fraction
        Float? smoothing_threshold_allele_fraction
        Float? smoothing_threshold_copy_ratio
        Float? calling_copy_ratio_z_score_threshold

        ### parameter-fu
        ### WGS concordance
        # columns of interest in the ground truth files
#        Array[String]? wgs_columns_of_interest
#        Array[String]? wgs_columns_of_interest_seg_calls
#        Array[String] wgs_columns_of_interest_or_default = select_first([wgs_columns_of_interest, [
#            "final_total_cn", "final_major_cn", "final_minor_cn", "consensus_total_cn", "consensus_major_cn", "consensus_minor_cn",
#            "star", "level", "absolute_broad_major_cn", "absolute_broad_minor_cn", "battenberg_major_cn", "battenberg_minor_cn", "sclust_major_cn", "sclust_minor_cn",
#            "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90"
#        ]])
#        Array[String] wgs_columns_of_interest_seg_calls_or_default = select_first([wgs_columns_of_interest_seg_calls, [
#            "final_total_cn", "final_major_cn", "final_minor_cn", "consensus_total_cn", "consensus_major_cn", "consensus_minor_cn",
#            "star", "level", "absolute_broad_major_cn", "absolute_broad_minor_cn", "battenberg_major_cn", "battenberg_minor_cn", "sclust_major_cn", "sclust_minor_cn",
#            "CALL", "MEAN_LOG2_COPY_RATIO"
#        ]])
#        Array[File] wgs_tumor_bam_files
#        Array[File] wgs_tumor_bam_indices
#        Array[File] wgs_normal_bam_files
#        Array[File] wgs_normal_bam_indices
#        Array[File] wgs_gt_seg_files
#        Int num_wgs_bam_files = length(wgs_tumor_bam_files)

        ### purity files
        Array[String]? purity_columns_of_interest
        Array[String]? purity_columns_of_interest_seg_calls
        Array[String] purity_columns_of_interest_or_default = select_first([purity_columns_of_interest, [
            "cn", "NM_id", "gene_sym", "strand", "width",
            "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90",
            "MINOR_ALLELE_FRACTION_POSTERIOR_50", "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_90"
        ]])
        Array[String] purity_columns_of_interest_seg_calls_or_default = select_first([purity_columns_of_interest_seg_calls, [
            "cn", "NM_id", "gene_sym", "strand", "width",
            "CALL", "MEAN_LOG2_COPY_RATIO"
        ]])
        Array[File] purity_tumor_bam_files
        Array[File] purity_tumor_bam_indices
        Array[File] purity_normal_bam_files
        Array[File] purity_normal_bam_indices
        Array[File] purity_gt_seg_files

        Int num_purity_bam_files = length(purity_tumor_bam_files)

        ### Clinical Sensitivity
        Array[String]? clinical_columns_of_interest
        Array[String]? clinical_columns_of_interest_seg_calls
        Array[String] clinical_columns_of_interest_or_default = select_first([clinical_columns_of_interest, [
            "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90", "Segment_Call", "Segment_Mean"
        ]])
        Array[String] clinical_columns_of_interest_seg_calls_or_default = select_first([clinical_columns_of_interest_seg_calls, [
            "CALL", "Segment_Call", "Segment_Mean", "MEAN_LOG2_COPY_RATIO"
        ]])
        Array[File] clinical_tumor_bams
        Array[File] clinical_tumor_bais
        Int num_clinical_bams = length(clinical_tumor_bams)
        # These cannot be released publicly.
        Array[File] clinical_seg_gts
        #####

#        File centromere_track

        # SM-74P4M and SM-74NF5
        Array[Int] reproducibility_indexes = [5, 10]
        Int index1 = reproducibility_indexes[0]
        Int index2 = reproducibility_indexes[1]

        ### CallCopyRatioSegments parameters
        Float? neutral_segment_copy_ratio_lower_bound
        Float? neutral_segment_copy_ratio_upper_bound
    }

    # Run the purity, clinical, and CR concordance (CRSP) evaluations.  CR Concordance gets run with purity.
    scatter (i in range(num_purity_bam_files)) {
        call cnv_validation.CNVValidation as cnvValidationPurity {
            input:
                intervals = ice_intervals,
                common_sites = common_sites,
                tumor_bam = purity_tumor_bam_files[i],
                tumor_bam_idx = purity_tumor_bam_indices[i],
                normal_bam = purity_normal_bam_files[i],
                normal_bam_idx = purity_normal_bam_indices[i],
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                read_count_pon = ice_read_count_pon,
                gatk_docker = gatk_docker,
                gatk4_jar_override =  gatk4_jar_override,
                gatk4_jar_override_evaluation = gatk4_jar_override_evaluation,
                bin_length = wxs_bin_length,
                columns_of_interest = purity_columns_of_interest_or_default,
                columns_of_interest_seg_calls = purity_columns_of_interest_seg_calls_or_default,
                gt_seg_file = purity_gt_seg_files[i],
                num_changepoints_penalty_factor_normal = num_changepoints_penalty_factor_normal,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
                calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
                neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
                neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound
        }
    }

    call PurityValidation {
        input:
            combined_purity_series_segs = cnvValidationPurity.combined_seg_cr_calls_file,
            group_id = group_id_final,
            eval_docker = eval_docker
    }



}

task PurityValidation {
    input {
        Array[File] combined_purity_series_segs
        String group_id
        String eval_docker
    }
    command <<<
        set -e

        python /root/plot_purity_series_hcc1143.py -O ${group_id}/purity/  ${sep=" " combined_purity_series_segs}
        echo "Doing tar..."
        tar zcvf ${group_id}_purity.tar.gz ${group_id}/purity/
    >>>

    runtime {
        docker: "${eval_docker}"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File final_purity_validation_tar_gz = "${group_id}_purity.tar.gz"
    }
}

