#Workflow for running gCNV SFARI evaluation

import "cnv_common_tasks.wdl" as CNVTasks

workflow GermlineCNVSFARIEvaluation {
    ##################################
    #### required basic arguments ####
    ##################################
    File segmental_duplication_track_bed
    File segmental_duplication_track_bed_idx
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker

    File sfari_calculate_metrics_script
    Array[File]+ genotyped_segments_vcfs
    File filtered_intervals
    File truth_bed
    File sample_id_map
    Int threshold_QS
    String gcnv_evaluation_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts

    ##################################################
    #### optional arguments for AnnotateIntervals ####
    ##################################################
    File? mappability_track_bed
    File? mappability_track_bed_idx
    Int? feature_query_lookahead
    Int? mem_gb_for_annotate_intervals

    call CNVTasks.AnnotateIntervals {
        input:
            intervals = filtered_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            mappability_track_bed = mappability_track_bed,
            mappability_track_bed_idx = mappability_track_bed_idx,
            segmental_duplication_track_bed = segmental_duplication_track_bed,
            segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
            feature_query_lookahead = feature_query_lookahead,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_annotate_intervals,
            preemptible_attempts = preemptible_attempts
    }

    call CalculateMetrics {
        input:
            sfari_calculate_metrics_script = sfari_calculate_metrics_script,
            genotyped_segments_vcfs = genotyped_segments_vcfs,
            filtered_annotated_intervals = AnnotateIntervals.annotated_intervals,
            truth_bed = truth_bed,
            sample_id_map = sample_id_map,
            threshold_QS = threshold_QS,
            gcnv_evaluation_docker = gcnv_evaluation_docker
    }

    output {
        File metrics_tar = CalculateMetrics.metrics_tar
    }
}

task CalculateMetrics {
    File sfari_calculate_metrics_script
    
    Array[File]+ genotyped_segments_vcfs
    File filtered_annotated_intervals
    File truth_bed
    File sample_id_map
    Int threshold_QS

    #Runtime parameters
    String gcnv_evaluation_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int num_cores = 1

    String dollar = "$" #WDL workaround for using array[@], see https://github.com/broadinstitute/cromwell/issues/1819
    
    command <<<
        set -e
        
        # copy calls to calls directory
        mkdir calls
        genotyped_segments_vcfs_array=(${sep=" " genotyped_segments_vcfs})
        for genotyped_segments_vcf in ${dollar}{genotyped_segments_vcfs_array[@]}; do
            cp $genotyped_segments_vcf calls
        done

        # run evaluation script
        Rscript ${sfari_calculate_metrics_script} \
            calls \
            ${filtered_annotated_intervals} \
            ${truth_bed} \
            ${sample_id_map} \
            ${threshold_QS} \
            ${num_cores} \
            output
            
        tar czf metrics.tar.gz -C output .
    >>>

    runtime {
        docker: "${gcnv_evaluation_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File metrics_tar = "metrics.tar.gz"
    }
}
