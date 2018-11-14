import "cnv_validation.wdl" as cnv_validation

import "https://api.firecloud.org/ga4gh/v1/tools/mkanaszn:mkn_CNV_validation_CallModeledSegments/versions/1/plain-WDL/descriptor" as cnv_validation

workflow MultiCNVValidation {
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
    
    ### Python files for validation
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py

    ### Validation parameters
    Float? num_changepoints_penalty_factor_normal
    Float? kernel_variance_allele_fraction
    Float? smoothing_threshold_allele_fraction
    Float? smoothing_threshold_copy_ratio
    Float? calling_copy_ratio_z_score_threshold
    
    ### CallModeledSegments parameters
    Boolean? load_copy_ratio
    Boolean? load_allele_fraction
    Float? normal_minor_allele_fraction_threshold
    Float? copy_ratio_peak_min_relative_height
    Float? copy_ratio_kernel_density_bandwidth
    Float? min_fraction_of_points_in_normal_allele_fraction_region
    Float? min_weight_first_cr_peak_cr_data_only
    Float? max_phred_score_normal
    Int? n_inference_iterations
    Float? inference_total_grad_norm_constraint
    Int? n_extra_gaussians_mixture_model
    Int? max_n_peaks_in_copy_ratio
    Int? mem_gb_for_call_modeled_segments    

    ### parameter-fu
    ### WGS concordance
    # columns of interest in the ground truth files
    Array[String]? wgs_columns_of_interest
    Array[String]? wgs_columns_of_interest_seg_calls
    Array[String] wgs_columns_of_interest_or_default = select_first([wgs_columns_of_interest, [
        "final_total_cn", "final_major_cn", "final_minor_cn", "consensus_total_cn", "consensus_major_cn", "consensus_minor_cn",
        "star", "level", "absolute_broad_major_cn", "absolute_broad_minor_cn", "battenberg_major_cn", "battenberg_minor_cn", "sclust_major_cn", "sclust_minor_cn",
        "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90"
    ]])
    Array[String] wgs_columns_of_interest_seg_calls_or_default = select_first([wgs_columns_of_interest_seg_calls, [
        "final_total_cn", "final_major_cn", "final_minor_cn", "consensus_total_cn", "consensus_major_cn", "consensus_minor_cn",
        "star", "level", "absolute_broad_major_cn", "absolute_broad_minor_cn", "battenberg_major_cn", "battenberg_minor_cn", "sclust_major_cn", "sclust_minor_cn",
        "CALL", "LOG2_COPY_RATIO_POSTERIOR_50"
    ]])
    Array[File] wgs_tumor_bam_files
    Array[File] wgs_tumor_bam_indices
    Array[File] wgs_normal_bam_files
    Array[File] wgs_normal_bam_indices
    Array[File] wgs_gt_seg_files
    Int num_wgs_bam_files = length(wgs_tumor_bam_files)

    ### purity files
    Array[String]? purity_columns_of_interest
    Array[String]? purity_columns_of_interest_seg_calls
    Array[String] purity_columns_of_interest_or_default = select_first([purity_columns_of_interest, [
        "cn", "NM_id", "gene_sym", "strand", "width",
        "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90"
    ]])
    Array[String] purity_columns_of_interest_seg_calls_or_default = select_first([purity_columns_of_interest_seg_calls, [
        "cn", "NM_id", "gene_sym", "strand", "width",
        "CALL", "LOG2_COPY_RATIO_POSTERIOR_50"
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
        "CALL", "Segment_Call", "Segment_Mean", "LOG2_COPY_RATIO_POSTERIOR_50"
    ]])
    Array[File] clinical_tumor_bams
    Array[File] clinical_tumor_bais
    Int num_clinical_bams = length(clinical_tumor_bams)
    # These cannot be released publicly.
    Array[File] clinical_seg_gts
    #####
    
    File centromere_track

    # SM-74P4M and SM-74NF5
    Array[Int] reproducibility_indexes = [5, 10]
    Int index1 = reproducibility_indexes[0]
    Int index2 = reproducibility_indexes[1]

    ### Run WGS and combine the blacklists in the evaluation.
    scatter (i in range(num_wgs_bam_files)) {
        call cnv_validation.CNVValidation as cnvValidationWGS {
            input:
                intervals = wgs_intervals,
                common_sites = common_sites,
                tumor_bam = wgs_tumor_bam_files[i],
                tumor_bam_idx = wgs_tumor_bam_indices[i],
                normal_bam = wgs_normal_bam_files[i],
                normal_bam_idx = wgs_normal_bam_indices[i],
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                read_count_pon = wgs_read_count_pon,
                gatk_docker = gatk_docker,
                gatk4_jar_override =  gatk4_jar_override,
                gatk4_jar_override_evaluation = gatk4_jar_override_evaluation,
                bin_length = wgs_bin_length,
                columns_of_interest = wgs_columns_of_interest_or_default,
                columns_of_interest_seg_calls = wgs_columns_of_interest_seg_calls_or_default,
                gt_seg_file = wgs_gt_seg_files[i],
                num_changepoints_penalty_factor_normal = num_changepoints_penalty_factor_normal,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
            	load_copy_ratio = load_copy_ratio,
                load_allele_fraction = load_allele_fraction,
                normal_minor_allele_fraction_threshold = normal_minor_allele_fraction_threshold,
                copy_ratio_peak_min_relative_height = copy_ratio_peak_min_relative_height,
                copy_ratio_kernel_density_bandwidth = copy_ratio_kernel_density_bandwidth,
                min_fraction_of_points_in_normal_allele_fraction_region = min_fraction_of_points_in_normal_allele_fraction_region,
                min_weight_first_cr_peak_cr_data_only = min_weight_first_cr_peak_cr_data_only,
                max_phred_score_normal = max_phred_score_normal,
                n_inference_iterations = n_inference_iterations,
                inference_total_grad_norm_constraint = inference_total_grad_norm_constraint,
                n_extra_gaussians_mixture_model = n_extra_gaussians_mixture_model,
                max_n_peaks_in_copy_ratio = max_n_peaks_in_copy_ratio,
                mem_gb_for_call_modeled_segments = mem_gb_for_call_modeled_segments
        }

        call CombineTracks {
            input:
                combined_seg = cnvValidationWGS.combined_seg_file,
                matched_normal_called_seg = select_first([cnvValidationWGS.called_copy_ratio_segments_normal, "null"]),
                gatk4_jar_override  = gatk4_jar_override_evaluation,
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                gatk_docker = gatk_docker,
                centromere_track = centromere_track
        }
    }

    call BpConcordanceValidation {
        input:
            combined_segs = CombineTracks.combined_segs_with_tracks,
            gt_seg_files = wgs_gt_seg_files,
            group_id = group_id_final,
            eval_docker = eval_docker,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
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
            	load_copy_ratio = load_copy_ratio,
                load_allele_fraction = load_allele_fraction,
                normal_minor_allele_fraction_threshold = normal_minor_allele_fraction_threshold,
                copy_ratio_peak_min_relative_height = copy_ratio_peak_min_relative_height,
                copy_ratio_kernel_density_bandwidth = copy_ratio_kernel_density_bandwidth,
                min_fraction_of_points_in_normal_allele_fraction_region = min_fraction_of_points_in_normal_allele_fraction_region,
                min_weight_first_cr_peak_cr_data_only = min_weight_first_cr_peak_cr_data_only,
                max_phred_score_normal = max_phred_score_normal,
                n_inference_iterations = n_inference_iterations,
                inference_total_grad_norm_constraint = inference_total_grad_norm_constraint,
                n_extra_gaussians_mixture_model = n_extra_gaussians_mixture_model,
                max_n_peaks_in_copy_ratio = max_n_peaks_in_copy_ratio,
                mem_gb_for_call_modeled_segments = mem_gb_for_call_modeled_segments
        }
    }

    scatter (i in range(num_clinical_bams)) {
        call cnv_validation.CNVValidation as cnvValidationClinical {
            input:
                intervals = ice_intervals,
                common_sites = common_sites,
                tumor_bam = clinical_tumor_bams[i],
                tumor_bam_idx = clinical_tumor_bais[i],
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                read_count_pon = ice_read_count_pon,
                gatk_docker = gatk_docker,
                gatk4_jar_override =  gatk4_jar_override,
                gatk4_jar_override_evaluation = gatk4_jar_override_evaluation,
                bin_length = wxs_bin_length,
                columns_of_interest = clinical_columns_of_interest_or_default,
                columns_of_interest_seg_calls = clinical_columns_of_interest_seg_calls_or_default,
                gt_seg_file = clinical_seg_gts[i],
                num_changepoints_penalty_factor_normal = num_changepoints_penalty_factor_normal,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
            	load_copy_ratio = load_copy_ratio,
                load_allele_fraction = load_allele_fraction,
                normal_minor_allele_fraction_threshold = normal_minor_allele_fraction_threshold,
                copy_ratio_peak_min_relative_height = copy_ratio_peak_min_relative_height,
                copy_ratio_kernel_density_bandwidth = copy_ratio_kernel_density_bandwidth,
                min_fraction_of_points_in_normal_allele_fraction_region = min_fraction_of_points_in_normal_allele_fraction_region,
                min_weight_first_cr_peak_cr_data_only = min_weight_first_cr_peak_cr_data_only,
                max_phred_score_normal = max_phred_score_normal,
                n_inference_iterations = n_inference_iterations,
                inference_total_grad_norm_constraint = inference_total_grad_norm_constraint,
                n_extra_gaussians_mixture_model = n_extra_gaussians_mixture_model,
                max_n_peaks_in_copy_ratio = max_n_peaks_in_copy_ratio,
                mem_gb_for_call_modeled_segments = mem_gb_for_call_modeled_segments
        }

        call ClinicalSensitivityPrep {
            input:
                entity_id = basename(clinical_tumor_bams[i]),
                called_merged_seg_file = cnvValidationClinical.combined_seg_cr_calls_file,
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                intervals = ice_intervals,
                gatk4_jar_override = gatk4_jar_override_evaluation,
                eval_docker = gatk_docker
        }
    }

    call PurityValidation {
        input:
            combined_purity_series_segs = cnvValidationPurity.combined_seg_cr_calls_file,
            group_id = group_id_final,
            eval_docker = eval_docker,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
       }
       
    call ReproducibilityValidationPrep {
        input:
            called_segs_1 = cnvValidationPurity.combined_seg_cr_calls_file[index1],
            called_segs_2 = cnvValidationPurity.combined_seg_cr_calls_file[index2],
            group_id = group_id_final,
            targets_file = cnvValidationPurity.denoised_copy_ratios_tumor[index1],
            gatk4_jar_override  = gatk4_jar_override,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk_docker = gatk_docker,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
    }
    
    call ReproducibilityValidation {
        input:
            called_segs_1 = cnvValidationPurity.combined_seg_cr_calls_file[index1],
            called_segs_2 = cnvValidationPurity.combined_seg_cr_calls_file[index2],
            group_id = group_id_final,
            targets_file = cnvValidationPurity.denoised_copy_ratios_tumor[index1],
            reproducibility_targets = ReproducibilityValidationPrep.reproducibility_targets,
            gatk4_jar_override  = gatk4_jar_override_evaluation,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            eval_docker = eval_docker,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
    }

    call ClinicalSensitivity {
        input:
            merged_seg_and_gt_files = ClinicalSensitivityPrep.clinical_sensitivity_seg_file,
            group_id = group_id_final,
            eval_docker = eval_docker,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
    }

    call CreateHtmlReport {
        input:
            eval_docker = eval_docker,
            gatk_docker = gatk_docker,
            clinical_sensitivity_tar_gz = ClinicalSensitivity.final_clinical_sensitivity_tar_gz,
            bp_concordance_tar_gz = BpConcordanceValidation.final_bp_concordance_validation_outputs_tar_gz,
            reproducibility_tar_gz = ReproducibilityValidation.final_reproducibility_validation_tar_gz,
            purity_tar_gz = PurityValidation.final_purity_validation_tar_gz,
            group_id = group_id_final,
            cli_utils_py = cli_utils_py,
            clopper_pearson_py = clopper_pearson_py,
            plot_bp_concordance_pcawg_pilot_py = plot_bp_concordance_pcawg_pilot_py,
            plot_purity_series_hcc1143_py = plot_purity_series_hcc1143_py,
            run_plot_reproducibility_py = run_plot_reproducibility_py,
            plot_clinical_sensitivity_py = plot_clinical_sensitivity_py,
            run_html_report_py = run_html_report_py
    }
}

task CombineTracks {
    File combined_seg
    File matched_normal_called_seg
    File? gatk4_jar_override
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    String gatk_docker
    File centromere_track

    String output_name = basename(combined_seg)

    command <<<
    set -e
    echo "need to add --columns-of-interest POSSIBLE_GERMLINE when introducing germline tagging"

    echo "======= GERMLINE TAGGING"
    java -jar "/root/gatk.jar" TagGermlineEvents \
            --segments ${combined_seg} --called-matched-normal-seg-file ${matched_normal_called_seg} \
            -O ${output_name}.combined.germline_tagged.seg -R ${ref_fasta}

    echo "======= Centromeres "
    java -jar "/root/gatk.jar" CombineSegmentBreakpoints \
            --segments ${output_name}.combined.germline_tagged.seg --segments ${centromere_track}  \
            --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_10 \
            --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50 --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_90 \
            --columns-of-interest absolute_broad_major_cn --columns-of-interest absolute_broad_minor_cn --columns-of-interest battenberg_major_cn \
            --columns-of-interest battenberg_minor_cn --columns-of-interest consensus_major_cn --columns-of-interest consensus_minor_cn \
            --columns-of-interest consensus_total_cn --columns-of-interest final_major_cn --columns-of-interest final_minor_cn \
            --columns-of-interest final_total_cn --columns-of-interest level --columns-of-interest sclust_major_cn --columns-of-interest sclust_minor_cn --columns-of-interest star \
             --columns-of-interest type \
            --columns-of-interest POSSIBLE_GERMLINE \
            -O ${output_name}.final.seg -R ${ref_fasta}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.0.7.0"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
        bootDiskSizeInGb: "40"
    }

    output {
        File combined_segs_with_tracks = "${output_name}.final.seg"
    }
}

task BpConcordanceValidation {

    # This parameter is only optional to fix a coercion error.  Please treat as required
    Array[File?] combined_segs
    Array[File] gt_seg_files
    String group_id
    String eval_docker
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py

    command <<<
        set -e
        pushd .
        cd /root/
        # SUPER HACK:  We get the barcode from the ground truth filename
        python ${plot_bp_concordance_pcawg_pilot_py} -I ${sep=" -I " combined_segs} \
         -G ${sep=" -G " gt_seg_files} \
         -O ${group_id}/bp_concordance/
        tar zcvf ${group_id}_wgs_concordance.tar.gz ${group_id}/bp_concordance/
        popd
        cp /root/${group_id}_wgs_concordance.tar.gz .
    >>>

    runtime {
        docker: "${eval_docker}"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File final_bp_concordance_validation_outputs_tar_gz = "${group_id}_wgs_concordance.tar.gz"
    }
}

task PurityValidation {
    Array[File] combined_purity_series_segs
    String group_id
    String eval_docker
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py
	
    command <<<
        set -e

        python ${plot_purity_series_hcc1143_py} -O ${group_id}/purity/  ${sep=" " combined_purity_series_segs}
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

task ReproducibilityValidationPrep {
    File called_segs_1
    File called_segs_2
    String group_id
    String gatk_docker
    File targets_file
    File? gatk4_jar_override
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py

    # This should be a optional, but cromwell 30 croaks.
    Float ploidy = 3.7

    String base_targets_file = basename(targets_file)
    String sample1_name = basename(called_segs_1)
    String sample2_name = basename(called_segs_2)
    Boolean is_cr = false
    command <<<
        set -e
        # TODO: We need runtime parameters

        # Changing extension to work with CombineSegmentBreakpoints
        cp ${targets_file} targets_file.seg

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
            --segments ${called_segs_1} --segments ${called_segs_2}  \
			--columns-of-interest CALL --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50 \
            -O reproducibility.tsv.seg -R ${ref_fasta}

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
            --segments reproducibility.tsv.seg --segments targets_file.seg  \
			--columns-of-interest CALL_1 --columns-of-interest CALL_2 --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50_1 --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50_2 --columns-of-interest LOG2_COPY_RATIO \
            -O reproducibility_targets_tmp.tsv.seg -R ${ref_fasta}

        egrep -v "^\@" reproducibility_targets_tmp.tsv.seg > ${group_id}_reproducibility_targets.seg
        
    >>>
    
    runtime {
        docker: "${gatk_docker}"
        memory: "4 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }
    
    output {
       File reproducibility_targets = "${group_id}_reproducibility_targets.seg"
    }
}

task ReproducibilityValidation {
    File called_segs_1
    File called_segs_2
    File reproducibility_targets
    String group_id
    String eval_docker
    File targets_file
    File? gatk4_jar_override
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py
    
    # This should be a optional, but cromwell 30 croaks.
    Float ploidy = 3.7
    String sample1_name = basename(called_segs_1)
    String sample2_name = basename(called_segs_2)
    Boolean is_cr = false
    command <<<
        set -e
        
        echo "Plotting...."
        python ${run_plot_reproducibility_py} \
            ${reproducibility_targets} \
            ${sample1_name} \
            ${sample2_name} \
            ${group_id}/reproducibility/ \
            ${ploidy} \
            ${true='--cr' false='' is_cr}

        tar zcvf ${group_id}_reproducibility.tar.gz ${group_id}/reproducibility/
    >>>
	
    runtime {
        docker: "${eval_docker}"
        memory: "4 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File final_reproducibility_validation_tar_gz = "${group_id}_reproducibility.tar.gz"
    }
}

task ClinicalSensitivity {
    Array[File] merged_seg_and_gt_files
    String group_id
    String eval_docker
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py

    command <<<
        set -e
        mkdir -p ${group_id}
        python ${plot_clinical_sensitivity_py} -I ${sep=" -I " merged_seg_and_gt_files} \
            -O ${group_id}/clinical/

        tar zcvf ${group_id}_clinical_sensitivity.tar.gz ${group_id}/clinical/
    >>>
    runtime {
        docker: "${eval_docker}"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File final_clinical_sensitivity_tar_gz = "${group_id}_clinical_sensitivity.tar.gz"
    }
}

task ClinicalSensitivityPrep {
    String entity_id
    File called_merged_seg_file
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File intervals
    File? gatk4_jar_override
    String eval_docker

    command <<<

        # Need to convert the interval file.  FRAGILE
        egrep "^\@" ${intervals} > ${entity_id}_tmp_intervals.seg
        echo -e "CONTIG\tSTART\tEND\tstrand\tNAME" >> ${entity_id}_tmp_intervals.seg
        egrep -v "^\@" ${intervals} >> ${entity_id}_tmp_intervals.seg

        java -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
            --segments ${entity_id}_tmp_intervals.seg --segments ${called_merged_seg_file}  \
			--columns-of-interest CALL \
			--columns-of-interest Segment_Call \
			--columns-of-interest Segment_Mean \
			--columns-of-interest NAME \
            -O ${entity_id}.clinical_sensitivity_by_target.seg -R ${ref_fasta}

    >>>
    runtime {
        docker: "${eval_docker}"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File clinical_sensitivity_seg_file = "${entity_id}.clinical_sensitivity_by_target.seg"
        File combine_input_file = "${entity_id}_tmp_intervals.seg"
    }
}

task CreateHtmlReport {
    String eval_docker
    String gatk_docker
    File clinical_sensitivity_tar_gz
    File bp_concordance_tar_gz
    File reproducibility_tar_gz
    File purity_tar_gz
    String group_id
    File cli_utils_py
    File clopper_pearson_py
    File plot_bp_concordance_pcawg_pilot_py
    File plot_purity_series_hcc1143_py
    File run_plot_reproducibility_py
    File plot_clinical_sensitivity_py
    File run_html_report_py

    command <<<
        set -e
        pushd .
        mkdir /root/report
        cd /root/report
        tar zxvf ${clinical_sensitivity_tar_gz}
        tar zxvf ${bp_concordance_tar_gz}
        tar zxvf ${reproducibility_tar_gz}
        tar zxvf ${purity_tar_gz}
        cp /root/style.css .

        cd /root/report/

        python ${run_html_report_py} --eval-docker ${eval_docker} --gatk-docker ${gatk_docker} \
        --purity-dir ${group_id}/purity/ \
        --clinical-dir ${group_id}/clinical/ \
        --reproducibility-dir ${group_id}/reproducibility/ \
        --bp-concordance-dir ${group_id}/bp_concordance/ \
        --html_template /root/aggregate_template.html \
        /root/report/

        cd ..
        tar zcvf GATK_CNV_Report_${group_id}.tar.gz report/

        popd
        cp /root/GATK_CNV_Report_${group_id}.tar.gz .
    >>>

    runtime {
        docker: "${eval_docker}"
        memory: "1 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
      File report_tar_gz = "GATK_CNV_Report_${group_id}.tar.gz"
    }
}
