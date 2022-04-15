version 1.0

import "gatk_wdls/somatic/cnv_somatic_pair_workflow.wdl" as CNVSomaticPairWorkflow

workflow CNVValidation {
    input {
        ### CNV parameters
        File intervals
        File? blacklist_intervals
        File common_sites
        File tumor_bam
        File tumor_bam_idx
        File? normal_bam
        File? normal_bam_idx
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
        File read_count_pon
        String gatk_docker

        File? gatk4_jar_override
        Int? bin_length

        File? gatk4_jar_override_evaluation

        ### Validation parameters
        Array[String] columns_of_interest
        Array[String] columns_of_interest_seg_calls
        File gt_seg_file
        Float? num_changepoints_penalty_factor_normal
        Float? kernel_variance_allele_fraction
        Float? smoothing_threshold_allele_fraction
        Float? smoothing_threshold_copy_ratio
        Float? calling_copy_ratio_z_score_threshold
        #########

        Float? neutral_segment_copy_ratio_lower_bound
        Float? neutral_segment_copy_ratio_upper_bound
    }
    call CNVSomaticPairWorkflow.CNVSomaticPairWorkflow as cnvPair {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            common_sites = common_sites,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            read_count_pon = read_count_pon,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            bin_length = bin_length,
            kernel_variance_allele_fraction = kernel_variance_allele_fraction,
            smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
            smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
            mem_gb_for_model_segments = 31,
            calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
            neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
            neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound
    }

    call FixGtSegFile {
        input:
            seg_file = gt_seg_file
    }

    call CombineSegmentBreakpoints {
        input:
            seg_files = [cnvPair.modeled_segments_tumor, FixGtSegFile.cleaned_seg_file],
            labels = [basename(tumor_bam), "gt_seg_file"],
            columns_of_interest = columns_of_interest,
            gatk4_jar_override = gatk4_jar_override_evaluation,
            entity_id = basename(tumor_bam),
            gatk_docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai
    }

    call CombineSegmentBreakpoints as CombineSegmentBreakpointsCRCalls {
        input:
            seg_files = [cnvPair.called_copy_ratio_segments_tumor, FixGtSegFile.cleaned_seg_file],
            labels = [basename(tumor_bam), "gt_seg_file"],
            columns_of_interest = columns_of_interest_seg_calls,
            gatk4_jar_override = gatk4_jar_override_evaluation,
            entity_id = basename(tumor_bam) + ".calls",
            gatk_docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai
    }

    output {
        File combined_seg_file = CombineSegmentBreakpoints.combined_seg_file
        File combined_seg_cr_calls_file = CombineSegmentBreakpointsCRCalls.combined_seg_file
        File denoised_copy_ratios_tumor = cnvPair.denoised_copy_ratios_tumor
        File? called_copy_ratio_segments_normal = cnvPair.called_copy_ratio_segments_normal
    }

}

task FixGtSegFile {
    input {
        File seg_file
    }
    String base_seg_name = basename(seg_file)

    command <<<
        set -e

        cp ~{seg_file} ~{base_seg_name}.fixed.seg

    >>>

    runtime {
        docker: "ubuntu:16.04"
        preemptible: 2
    }

    output {
        File cleaned_seg_file = "${base_seg_name}.fixed.seg"
    }
}

task CombineSegmentBreakpoints {
    input {
        Array[File]+ seg_files
        Array[String]+ labels
        Array[String]+ columns_of_interest

        File? gatk4_jar_override
        String entity_id

        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai

        # Runtime parameters
        String gatk_docker
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? mem
    }
    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command <<<
        set -e
        # Use GATK Jar override if specified
        GATK_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx~{command_mem}m -jar $GATK_JAR CombineSegmentBreakpoints \
            --segments ~{sep=" --segments " seg_files} --labels ~{sep=" --labels " labels} \
            --columns-of-interest ~{sep=" --columns-of-interest " columns_of_interest} -O ~{entity_id}.seg -R ~{ref_fasta}
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File combined_seg_file = "${entity_id}.seg"
    }
}