#Workflow for running CODEX pipeline on a single contig

workflow CODEXSingleContig {
    File bam_and_bai_list_file
    File sample_name_list_file
    File target_bed_file
    String contig
    String cohort_name
    
    #optional parameters for coverage collection
    Int? min_mapping_quality
    
    #optional parameters for running CODEX
    Int? min_median_coverage
    Int? max_median_coverage
    Int? min_length
    Int? max_length
    Int? min_GC
    Int? max_GC
    Float? min_mappability
    Int? max_num_latent_factors
    
    File codex_coverage_R_script
    File codex_R_script
    String codex_docker
    
    Array[Array[File]] bam_and_bais = read_tsv(bam_and_bai_list_file)
    Array[String] sample_names = read_lines(sample_name_list_file)

    Array[Pair[Array[File], String]] samples = zip(bam_and_bais, sample_names)
    
    scatter (sample in samples) {
        call CollectCoverage {
            input:
                bam = sample.left[0],
                bam_idx = sample.left[1],
                target_bed_file = target_bed_file,
                sample_name = sample.right,
                contig = contig,
                min_mapping_quality = min_mapping_quality,
                codex_coverage_R_script = codex_coverage_R_script,
                codex_docker = codex_docker
        }
    }
    
    call RunCODEX {
        input:
            coverage_files = CollectCoverage.coverage_file,
            contig = contig,
            cohort_name = cohort_name,
            min_median_coverage = min_median_coverage,
            max_median_coverage = max_median_coverage,
            min_length = min_length,
            max_length = max_length,
            min_GC = min_GC,
            max_GC = max_GC,
            min_mappability = min_mappability,
            max_num_latent_factors = max_num_latent_factors,
            codex_R_script = codex_R_script,
            codex_docker = codex_docker
    }

    output {
        Array[File] coverage_files = CollectCoverage.coverage_file
        File qcmat_file = RunCODEX.qcmat_file
        File choiceofK_file = RunCODEX.choiceofK_file
        File segment_file = RunCODEX.segment_file
    }
}

task CollectCoverage {
    File bam
    File bam_idx
    File target_bed_file
    String sample_name
    String contig
    
    Int? min_mapping_quality
    
    File codex_coverage_R_script

    #Runtime parameters
    String codex_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        set -e
        mv ${bam_idx} ${bam}.bai
        Rscript ${codex_coverage_R_script} \
            ${bam} \
            ${target_bed_file} \
            ${sample_name} \
            ${contig} \
            ${default="20" min_mapping_quality}
    >>>

    runtime {
        docker: "${codex_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File coverage_file = "${sample_name}_${contig}_coverage.tsv"
    }
}

task RunCODEX {
    Array[File] coverage_files
    String contig
    String cohort_name
    String? output_dir

    Int? min_median_coverage
    Int? max_median_coverage
    Int? min_length
    Int? max_length
    Int? min_GC
    Int? max_GC
    Float? min_mappability
    Int? max_num_latent_factors
    
    File codex_R_script

    #Runtime parameters
    String codex_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000

    String output_dir_ = select_first([output_dir, "output"])
 
    command <<<
        set -e
        mkdir ${output_dir_}
        Rscript ${codex_R_script} \
            ${write_lines(coverage_files)} \
            ${contig} \
            ${cohort_name} \
            ${output_dir_} \
            ${default="20" min_median_coverage} \
            ${default="4000" max_median_coverage} \
            ${default="20" min_length} \
            ${default="2000" max_length} \
            ${default="20" min_GC} \
            ${default="80" max_GC} \
            ${default="0.9" min_mappability} \
            ${default="9" max_num_latent_factors}
    >>>

    runtime {
        docker: "${codex_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File qcmat_file = "${output_dir_}/${cohort_name}_${contig}_qcmat.txt"
        File choiceofK_file = "${output_dir_}/${cohort_name}_${contig}_choiceofK.pdf"
        File segment_file = glob("${output_dir_}/${cohort_name}_${contig}_*_CODEX_integer.seg")[0]
    }
}
