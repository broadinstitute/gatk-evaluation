# Workflow to cluster a set of samples
#
#############

import "https://api.firecloud.org/ga4gh/v1/tools/broad-firecloud-dsde:GATK_Germline_CNV_Validation_Common_Tasks/versions/16/plain-WDL/descriptor" as CNVTasks

workflow ClusterSamples {

    ##################################
    #### required basic arguments ####
    ##################################
    File intervals_for_clustering
    Array[String]+ normal_bams
    Array[String]+ normal_bais
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    File clustering_script
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts

    ##############################################
    #### optional arguments for CollectCounts ####
    ##############################################
    String? collect_counts_format
    Int? mem_gb_for_collect_counts

    ###############################################
    #### optional arguments for ClusterSamples ####
    ###############################################
    File? clustering_prior_table
    Int? max_num_clusters
    Int? mem_gb_for_cluster_samples

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = intervals_for_clustering,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts
        }
    }

    call PerformClustering {
        input:
            read_count_files = CollectCounts.counts,
            entity_ids = CollectCounts.entity_id,
            clustering_prior_table = clustering_prior_table,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_cluster_samples,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File clustering_table = PerformClustering.clustering_table
        File sample_clusters_plot = PerformClustering.sample_clusters_plot
    }
}

task PerformClustering {
    Array[File] read_count_files
    Array[String] entity_ids
    File? clustering_prior_table
    File clustering_script
    Int? max_num_clusters

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
        python ${clustering_script} \
            --output_dir . \
            --max_clusters ${default=10 max_num_clusters} \
            --read_count_files ${sep=" " read_count_files} \
            --entity_ids ${sep=" " entity_ids}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File clustering_table = "clustering_table.tsv"
        File sample_clusters_plot = "sample_clusters.png"
    }
}