# Evaluate performance of GATK gCNV against a callset produced by Talkowski lab SV pipeline
####

workflow EvaluateGermlineCNVCalls {

    ##################################
    #### required basic arguments ####
    ##################################
    Array[File]+ genotyped_segments_vcfs
    File truth_bed_sample_ids
    File padded_intervals
    String gcnv_evaluation_docker
    File blacklisted_intervals_truth
    File ref_fasta_dict

    ##################################
    #### optional basic arguments ####
    ##################################
    Int? preemptible_attempts
    String gcnv_eval_script = "/root/plot_evaluation_metrics.py"

    call EvaluateCalls {
        input:
            genotyped_segments_vcfs = genotyped_segments_vcfs,
            truth_bed_sample_ids = truth_bed_sample_ids,
            padded_intervals = padded_intervals,
            gcnv_evaluation_docker = gcnv_evaluation_docker,
            blacklisted_intervals_truth = blacklisted_intervals_truth,
            ref_fasta_dict = ref_fasta_dict,
            gcnv_eval_script = gcnv_eval_script,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File confusion_matrix = EvaluateCalls.confusion_matrix
        File confusion_matrix_bounded_filters = EvaluateCalls.confusion_matrix_bounded_filters
        File area_under_roc = EvaluateCalls.area_under_roc
        Array[File] metrics_plots = EvaluateCalls.metrics_plots
    }
}

task EvaluateCalls {
    Array[File]+ genotyped_segments_vcfs
    File truth_bed_sample_ids
    File padded_intervals
    String gcnv_eval_script
    File blacklisted_intervals_truth
    File ref_fasta_dict
    Array[File]+ gcnv_model_tars
    
    Array[String] callset_filter_names = ['QS', 'QA']
    Array[Float] callset_filter_max_values = [3100, 200]
    Array[Int] callset_filter_num_bins = [1000, 50]
    String attribute_for_roc_creation = 'QS'

    #Runtime parameters
    String gcnv_evaluation_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    Int num_model_tars = length(gcnv_model_tars)

    String dollar = "$" #WDL workaround for using array[@], see https://github.com/broadinstitute/cromwell/issues/1819
    command <<<
        set -e
        mkdir -p out/plots

        # untar models to shard-0, shard-1, etc directories 
        gcnv_model_tar_array=(${sep=" " gcnv_model_tars})
        for index in ${dollar}{!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${dollar}{gcnv_model_tar_array[$index]}
            mkdir shard-$index
            tar xzf $gcnv_model_tar -C shard-$index
        done


        python ${gcnv_eval_script} \
          --output_dir out \
          --ref_dict ${ref_fasta_dict} \
          --gcnv_segment_vcfs ${sep=' ' genotyped_segments_vcfs} \
          --sorted_truth_calls_bed ${truth_bed_sample_ids} \
          --padded_intervals ${padded_intervals} \
          --blacklisted_intervals_truth ${blacklisted_intervals_truth} \
          --callset_filter_names ${sep=' ' callset_filter_names} \
          --callset_filter_max_values ${sep=' ' callset_filter_max_values} \
          --callset_filter_num_bins ${sep=' ' callset_filter_num_bins} \
          --attribute_for_roc_creation ${attribute_for_roc_creation} \
          --gcnv_models_directory . \
          --num_model_shards ${num_model_tars}
    >>>

    runtime {
        docker: "${gcnv_evaluation_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File confusion_matrix = "out/confusion_matrix.tsv"
        File confusion_matrix_bounded_filters = "out/confusion_matrix_bounded_filters.tsv"
        File area_under_roc = "out/area_under_roc.tsv"
        Array[File] metrics_plots = glob("out/plots/*")
    }
}