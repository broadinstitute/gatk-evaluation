version 1.0

import "cnv_validation.wdl" as cnv_validation

workflow PurityCNVValidation {
    input {
        File ice_intervals
        File common_sites
        File? blacklist_intervals

        Array[File] purity_tumor_bam_files
        Array[File] purity_tumor_bam_indices
        File purity_normal_bam_file
        File purity_normal_bam_index

        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
        File ice_read_count_pon

        String gatk_docker
        File? gatk4_jar_override
        File? gatk4_jar_override_evaluation

        Int wxs_bin_length = 0

        Array[String]? purity_columns_of_interest
        Array[String]? purity_columns_of_interest_seg_calls
        Array[String] purity_columns_of_interest_or_default = select_first([purity_columns_of_interest, [
                                                                                                        "cn", "NM_id", "gene_sym", "strand", "width",
                                                                                                        "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_90"
                                                                                                        ]])
        Array[String] purity_columns_of_interest_seg_calls_or_default = select_first([purity_columns_of_interest_seg_calls, [
                                                                                                                            "cn", "NM_id", "gene_sym", "strand", "width",
                                                                                                                            "CALL", "MEAN_LOG2_COPY_RATIO"
                                                                                                                            ]])
        File purity_gt_seg_file

        Float? num_changepoints_penalty_factor_normal
        Float? kernel_variance_allele_fraction
        Float? smoothing_threshold_allele_fraction
        Float? smoothing_threshold_copy_ratio
        Float? calling_copy_ratio_z_score_threshold

        String eval_docker
        String? group_id
        String group_id_final = select_first([group_id, "input_segs"])

        Array[Int] reproducibility_indexes = [5, 10]
        Int index1 = reproducibility_indexes[0]
        Int index2 = reproducibility_indexes[1]

        Int num_purity_bam_files = length(purity_tumor_bam_files)
    }

    scatter (i in range(num_purity_bam_files)) {
        call cnv_validation.CNVValidation as cnvValidationPurity {
            input:
                intervals = ice_intervals,
                blacklist_intervals = blacklist_intervals,
                common_sites = common_sites,
                tumor_bam = purity_tumor_bam_files[i],
                tumor_bam_idx = purity_tumor_bam_indices[i],
                normal_bam = purity_normal_bam_file,
                normal_bam_idx = purity_normal_bam_index,
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
                gt_seg_file = purity_gt_seg_file,
                num_changepoints_penalty_factor_normal = num_changepoints_penalty_factor_normal,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
                calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold
        }
    }

    call PurityValidation {
        input:
            combined_purity_series_segs = cnvValidationPurity.combined_seg_cr_calls_file,
            group_id = group_id_final,
            eval_docker = eval_docker
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
            gatk_docker = gatk_docker
    }

    call ReproducibilityValidation {
        input:
            called_segs_1 = cnvValidationPurity.combined_seg_cr_calls_file[index1],
            called_segs_2 = cnvValidationPurity.combined_seg_cr_calls_file[index2],
            reproducibility_targets = ReproducibilityValidationPrep.reproducibility_targets,
            group_id = group_id_final,
            targets_file = cnvValidationPurity.denoised_copy_ratios_tumor[index1],
            gatk4_jar_override  = gatk4_jar_override_evaluation,
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

task ReproducibilityValidationPrep {
    input {
        File called_segs_1
        File called_segs_2
        String group_id
        String gatk_docker
        File targets_file
        File? gatk4_jar_override
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
    }
    command <<<
        set -e
        # TODO: We need runtime parameters

        # Changing extension to work with CombineSegmentBreakpoints
        cp ${targets_file} targets_file.seg

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
        --segments ${called_segs_1} --segments ${called_segs_2}  \
        --columns-of-interest CALL --columns-of-interest MEAN_LOG2_COPY_RATIO \
        -O reproducibility.tsv.seg -R ${ref_fasta}

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
        --segments reproducibility.tsv.seg --segments targets_file.seg  \
        --columns-of-interest CALL_1 --columns-of-interest CALL_2 --columns-of-interest MEAN_LOG2_COPY_RATIO_1 \
        --columns-of-interest MEAN_LOG2_COPY_RATIO_2 --columns-of-interest LOG2_COPY_RATIO \
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
    input {
        File called_segs_1
        File called_segs_2
        File reproducibility_targets
        String group_id
        String eval_docker
        File targets_file
        File? gatk4_jar_override
    }
    # This should be a optional, but cromwell 30 croaks.
    Float ploidy = 3.7

    String sample1_name = basename(called_segs_1)
    String sample2_name = basename(called_segs_2)
    Boolean is_cr = false
    command <<<
        set -e

        echo "Plotting...."
        python /root/run_plot_reproducibility.py \
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