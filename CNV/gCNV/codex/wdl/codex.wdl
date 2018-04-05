#Workflow for running CODEX pipeline

import "codex-single-contig.wdl" as CODEXSingleContig

workflow CODEX {
    File bam_and_bai_list_file
    File sample_name_list_file
    File target_bed_file
    File contig_list_file
    String cohort_name
    File codex_vcf_header_file
    
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
    File codex_vcf_converter_python_script
    String codex_docker
    
    Array[String] contigs = read_lines(contig_list_file)

    scatter (contig in contigs) {
        call CODEXSingleContig.CODEXSingleContig as CODEXSingleContig {
            input:
                bam_and_bai_list_file = bam_and_bai_list_file,
                sample_name_list_file = sample_name_list_file,
                target_bed_file = target_bed_file,
                contig = contig,
                cohort_name = cohort_name, 
                min_mapping_quality = min_mapping_quality,
                min_median_coverage = min_median_coverage,
                max_median_coverage = max_median_coverage,
                min_length = min_length,
                max_length = max_length,
                min_GC = min_GC,
                max_GC = max_GC,
                min_mappability = min_mappability,
                max_num_latent_factors = max_num_latent_factors,
                codex_coverage_R_script = codex_coverage_R_script,
                codex_R_script = codex_R_script,
                codex_docker = codex_docker
        }
    }

    call ConvertSegmentFilesToVCF {
        input:
            segment_files = CODEXSingleContig.segment_file,
            sample_name_list_file = sample_name_list_file,
            cohort_name = cohort_name,
            codex_vcf_header_file = codex_vcf_header_file,
            codex_vcf_converter_python_script = codex_vcf_converter_python_script,
            codex_docker = codex_docker
    }

    output {
        Array[Array[File]] coverage_files = CODEXSingleContig.coverage_files
        Array[File] qcmat_files = CODEXSingleContig.qcmat_file
        Array[File] choiceofK_files = CODEXSingleContig.choiceofK_file
        Array[File] segment_files = CODEXSingleContig.segment_file
        File output_vcf = ConvertSegmentFilesToVCF.output_vcf
    }
}

task ConvertSegmentFilesToVCF {
    Array[File] segment_files
    File sample_name_list_file
    String cohort_name
    File codex_vcf_header_file

    File codex_vcf_converter_python_script

    #Runtime parameters
    String codex_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    String output_vcf_filename = "${cohort_name}_CODEX_integer.vcf"

    command <<<
        set -e
        python3 ${codex_vcf_converter_python_script} \
            -segments ${sep=" -segments " segment_files} \
            -sample-name-list ${sample_name_list_file} \
            -codex-vcf-header ${codex_vcf_header_file} \
            -output-vcf ${output_vcf_filename}
    >>>

    runtime {
        docker: "${codex_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File output_vcf = output_vcf_filename
    }
}
