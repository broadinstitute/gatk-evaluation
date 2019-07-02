import "https://api.firecloud.org/ga4gh/v1/tools/slee:TH_Common_Tasks/versions/1/plain-WDL/descriptor" as CNVTasks

workflow CoverageMixingWorkflow {

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
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    Array[Int] mixing_tumor_percentages
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

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_fai, "GB"))
    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))
    Int normal_bam_size = ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB"))

    Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 20 + ceil(size(intervals, "GB")) + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0])

    Int preprocess_intervals_disk = ref_size + disk_pad
    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_preprocess_intervals,
            disk_space_gb = preprocess_intervals_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_counts_tumor_disk = tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsTumor {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            format = collect_counts_format,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_counts,
            disk_space_gb = collect_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_allelic_counts_tumor_disk = tumor_bam_size + ref_size + disk_pad
    call CNVTasks.CollectAllelicCounts as CollectAllelicCountsTumor {
        input:
            common_sites = common_sites,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            minimum_base_quality =  minimum_base_quality,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_allelic_counts,
            disk_space_gb = collect_allelic_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_counts_normal_disk = normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsNormal {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = normal_bam,
            bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            format = collect_counts_format,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_counts,
            disk_space_gb = collect_counts_normal_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_allelic_counts_normal_disk = normal_bam_size + ref_size + disk_pad
    call CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal {
        input:
            common_sites = common_sites,
            bam = normal_bam,
            bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            minimum_base_quality =  minimum_base_quality,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_allelic_counts,
            disk_space_gb = collect_allelic_counts_normal_disk,
            preemptible_attempts = preemptible_attempts
    }

    scatter (mixing_tumor_percentage in mixing_tumor_percentages) {
        call MixPair {
            input:
                normal_entity_id = CollectCountsNormal.entity_id,
                normal_counts_tsv = CollectCountsNormal.counts,
                normal_allelic_counts_tsv = CollectAllelicCountsNormal.allelic_counts,
                tumor_entity_id = CollectCountsTumor.entity_id,
                tumor_counts_tsv = CollectCountsTumor.counts,
                tumor_allelic_counts_tsv = CollectAllelicCountsTumor.allelic_counts,
                mixing_tumor_percentage = mixing_tumor_percentage,
                gatk_docker = gatk_docker
        }
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals

        File normal_read_counts_entity_id = CollectCountsNormal.entity_id
        File normal_read_counts = CollectCountsNormal.counts
        File normal_allelic_counts_entity_id = CollectAllelicCountsNormal.entity_id
        File normal_allelic_counts = CollectAllelicCountsNormal.allelic_counts

        File tumor_read_counts_entity_id = CollectCountsTumor.entity_id
        File tumor_read_counts = CollectCountsTumor.counts
        File tumor_allelic_counts_entity_id = CollectAllelicCountsTumor.entity_id
        File tumor_allelic_counts = CollectAllelicCountsTumor.allelic_counts

        Array[String] mixed_entity_ids = MixPair.entity_id
        Array[File] mixed_read_counts = MixPair.counts
        Array[File] mixed_allelic_counts = MixPair.allelic_counts
    }
}

task MixPair {
    String normal_entity_id
    File normal_counts_tsv
    File normal_allelic_counts_tsv
    String tumor_entity_id
    File tumor_counts_tsv
    File tumor_allelic_counts_tsv
    Int mixing_tumor_percentage
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    
    String entity_id_ = "N-${100 - mixing_tumor_percentage}-${normal_entity_id}-T-${mixing_tumor_percentage}-${tumor_entity_id}"

    command <<<
        set -e
        python - --normal_entity_id ${normal_entity_id} \
                 --normal_counts_tsv_path ${normal_counts_tsv} \
                 --normal_allelic_counts_tsv_path ${normal_allelic_counts_tsv} \
                 --tumor_entity_id ${tumor_entity_id} \
                 --tumor_counts_tsv_path ${tumor_counts_tsv} \
                 --tumor_allelic_counts_tsv_path ${tumor_allelic_counts_tsv} \
                 --mixing_tumor_percentage ${mixing_tumor_percentage} \
                 --output_path ${output_dir_} \
                 <<-'EOF'
        import argparse
        import os
        import pandas as pd
        import numpy as np

        def read_tsv(tsv_path):
            return pd.read_csv(tsv_path, sep='\t', comment='@', dtype={'CONTIG': str})

        def mix_counts(normal_counts_df, tumor_counts_df, mixing_tumor_percentage):
            assert np.all(normal_counts_df[['CONTIG', 'START', 'END']] == tumor_counts_df[['CONTIG', 'START', 'END']])
            mixed_counts_df = normal_counts_df.copy()
            mixed_counts_df['COUNT'] = np.random.binomial(normal_counts_df['COUNT'], 1 - 0.01 * mixing_tumor_percentage) + np.random.binomial(tumor_counts_df['COUNT'], 0.01 * mixing_tumor_percentage)
            return mixed_counts_df

        def mix_allelic_counts(normal_allelic_counts_df, tumor_allelic_counts_df, mixing_tumor_percentage):
            assert np.all(normal_allelic_counts_df[['CONTIG', 'POSITION']] == tumor_allelic_counts_df[['CONTIG', 'POSITION']])
            mixed_allelic_counts_df = normal_allelic_counts_df.copy() if mixing_tumor_percentage == 0 else tumor_allelic_counts_df.copy()
            mixed_allelic_counts_df['REF_COUNT'] = np.random.binomial(normal_allelic_counts_df['REF_COUNT'], 1 - 0.01 * mixing_tumor_percentage) + np.random.binomial(tumor_allelic_counts_df['REF_COUNT'], 0.01 * mixing_tumor_percentage)
            mixed_allelic_counts_df['ALT_COUNT'] = np.random.binomial(normal_allelic_counts_df['ALT_COUNT'], 1 - 0.01 * mixing_tumor_percentage) + np.random.binomial(tumor_allelic_counts_df['ALT_COUNT'], 0.01 * mixing_tumor_percentage)
            return mixed_allelic_counts_df

        class Pair:
            def __init__(self,
                         normal_entity_id, normal_counts_tsv_path, normal_allelic_counts_tsv_path,
                         tumor_entity_id, tumor_counts_tsv_path, tumor_allelic_counts_tsv_path):
                self.normal_entity_id = normal_entity_id
                self.normal_counts_df = read_tsv(normal_counts_tsv_path)
                self.normal_allelic_counts_df = read_tsv(normal_allelic_counts_tsv_path)

                self.tumor_entity_id = tumor_entity_id
                self.tumor_counts_df = read_tsv(tumor_counts_tsv_path)
                self.tumor_allelic_counts_df = read_tsv(tumor_allelic_counts_tsv_path)

                self.tumor_header = self._get_header(tumor_allelic_counts_tsv_path)

            def mix(self, mixing_tumor_percentage, output_path):
                assert isinstance(mixing_tumor_percentage, int) and 0 <= mixing_tumor_percentage <= 100
                mix_name = f'N-{100 - mixing_tumor_percentage}-{self.normal_entity_id}-T-{mixing_tumor_percentage}-{self.tumor_entity_id}'

                mixed_counts_df = mix_counts(self.normal_counts_df, self.tumor_counts_df, mixing_tumor_percentage)
                mixed_counts_tsv_path = os.path.join(output_path, mix_name + '.counts.tsv')
                mixed_counts_df.to_csv(mixed_counts_tsv_path, sep='\t', index=False)
                print('Mixed counts output to:', mixed_counts_tsv_path)

                mixed_allelic_counts_df = mix_allelic_counts(self.normal_allelic_counts_df, self.tumor_allelic_counts_df, mixing_tumor_percentage)
                mixed_allelic_counts_tsv_path = os.path.join(output_path, mix_name + '.allelicCounts.tsv')
                mixed_allelic_counts_df.to_csv(mixed_allelic_counts_tsv_path, sep='\t', index=False)
                print('Mixed allelic counts output to:', mixed_allelic_counts_tsv_path)

                #get, modify, and add headers
                mixed_header = self.tumor_header.replace(self.tumor_header.split(':')[-1], mix_name + '\n')
                self._prepend_header(mixed_counts_tsv_path, mixed_header)
                self._prepend_header(mixed_allelic_counts_tsv_path, mixed_header)

            def _get_header(self, path):
                header_lines = []
                with open(path, 'r') as f:
                    for line in f:
                        if line.startswith('@'):
                            header_lines.append(line)
                        else:
                            break
                return ''.join(header_lines)

            def _prepend_header(self, path, header):
                with open(path, 'r') as original:
                    data = original.read()
                with open(path, 'w') as modified:
                    modified.write(header + data)

        def main():
            parser = argparse.ArgumentParser(
                description='Mix counts and allelic counts from a matched pair.',
                formatter_class=argparse.RawTextHelpFormatter)

            parser.add_argument('--normal_entity_id',
                                type=str,
                                help='Normal entity ID')

            parser.add_argument('--normal_counts_tsv_path',
                                type=str,
                                help='Path to normal counts TSV.')

            parser.add_argument('--normal_allelic_counts_tsv_path',
                                type=str,
                                help='Path to normal allelic counts TSV.')

            parser.add_argument('--tumor_entity_id',
                                type=str,
                                help='Tumor entity ID')

            parser.add_argument('--tumor_counts_tsv_path',
                                type=str,
                                help='Path to tumor counts TSV.')

            parser.add_argument('--tumor_allelic_counts_tsv_path',
                                type=str,
                                help='Path to tumor allelic counts TSV.')

            parser.add_argument('--mixing_tumor_percentage',
                                type=int,
                                help='Mixing tumor percentage (an integer in [0, 100]).')

            parser.add_argument('--output_path',
                                type=str,
                                help='Path to output directory.')

            args = parser.parse_args()
            
            if not os.path.exists(args.output_path):
                os.mkdir(args.output_path)

            pair = Pair(args.normal_entity_id, args.normal_counts_tsv_path, args.normal_allelic_counts_tsv_path,
                        args.tumor_entity_id, args.tumor_counts_tsv_path, args.tumor_allelic_counts_tsv_path)

            pair.mix(args.mixing_tumor_percentage, args.output_path)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = entity_id_
        File counts = "${output_dir_}/${entity_id_}.counts.tsv"
        File allelic_counts = "${output_dir_}/${entity_id_}.allelicCounts.tsv"
    }
}