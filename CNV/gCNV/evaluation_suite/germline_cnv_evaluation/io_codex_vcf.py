import vcf
from typing import List, Optional, Dict, Tuple
from .core import GenericCNVCallSet, GenericCopyNumberVariant, get_overlapping_variants_set
from intervaltree_bio import GenomeIntervalTree
import logging

_logger = logging.getLogger(__name__)


def load_codex_vcf_file(codex_vcf_file: str,
                        codex_interval_list_bed_file: str,
                        max_records: Optional[int] = None,
                        log_frequency: int = 500,
                        min_overlap_fraction_for_variant_matching: float = 0.5) \
        -> Tuple[Dict[str, GenericCNVCallSet], GenomeIntervalTree]:

    allele_to_genotype_map = {'.': 'ref', '<DEL>': 'del', '<DUP>': 'dup'}

    # step 1. load VCF file
    codex_call_set_list: List[GenericCNVCallSet] = list()
    contig_set = set()
    with open(codex_vcf_file, 'r') as f:
        vcf_reader = vcf.Reader(f)
        sample_names = list(vcf_reader.samples)
        num_samples = len(sample_names)
        for sample_name in sample_names:
            sample_call_set = GenericCNVCallSet(sample_name, tags={"from CODEX segments VCF"})
            codex_call_set_list.append(sample_call_set)

        for record_num, record in enumerate(vcf_reader):

            if record_num % log_frequency == 0:
                _logger.info("Reading record number {0}...".format(record_num))

            if max_records is not None and record_num > max_records:
                break

            info = record.INFO
            contig = record.CHROM
            start = record.start + 1
            end = info['END']

            contig_set.add(contig)

            for si in range(num_samples):

                variant_frequency = float(info['AN']) / num_samples
                num_intervals = info['NT']

                original_genotype = record.samples[si]['GT']
                if original_genotype == '.':  # ref
                    genotype = 'ref'
                else:
                    genotype_index = int(original_genotype)
                    allele = str(record.alleles[genotype_index])
                    genotype = allele_to_genotype_map[allele]

                # Note: CODEX does not provide integer copy-number. We arbitrarily use 0 for deletion
                # and 3 for duplication
                if genotype == 'del':
                    var_copy_number = 0
                elif genotype == 'dup':
                    var_copy_number = 3
                else:
                    var_copy_number = 2

                original_quality = record.samples[si]['LRATIO']
                if original_quality is None:
                    quality = 0
                else:
                    quality = float(original_quality)

                # generate variant
                sample_var = GenericCopyNumberVariant(contig, start, end, var_copy_number, quality,
                                                      genotype=genotype,
                                                      num_intervals=num_intervals,
                                                      variant_frequency=variant_frequency,
                                                      variant_class=None)
                codex_call_set_list[si].add(sample_var)

    # step 2. included loci
    _logger.info("Obtaining included loci...")
    included_intervals = GenomeIntervalTree()
    with open(codex_interval_list_bed_file, 'r') as f:
        for line in f:
            if line[0] == '@':
                continue
            split_line = line.strip().split()
            contig = split_line[0]
            start = int(split_line[1])
            end = int(split_line[2])
            included_intervals[contig].addi(start, end)

    # step 3. non-variant "variants" are not needed anymore -- remove them from the interval tree
    _logger.info("Removing non-variant calls...")
    var_only_call_set_list: List[GenericCNVCallSet] = list()
    for sample_call_set in codex_call_set_list:
        var_only_call_set = GenericCNVCallSet(sample_call_set.sample_name, sample_call_set.tags)
        for contig in sample_call_set.contig_set:
            for variant in sample_call_set.iter_in_contig(contig):
                if variant.is_var:
                    var_only_call_set.add(variant)
        var_only_call_set.tags.add("Removed non-variant calls")
        var_only_call_set_list.append(var_only_call_set)

    # step 4. estimate variant frequency
    _logger.info("Calculating variant classes and frequencies...")
    num_samples = len(sample_names)
    for si in range(num_samples):
        for contig in contig_set:
            for variant in var_only_call_set_list[si].iter_in_contig(contig):
                all_variants = list()
                all_variants.append(variant)
                total_count = 1
                for osi in range(num_samples):
                    if osi == si:
                        continue
                    other_variants = get_overlapping_variants_set(
                        var_only_call_set_list[osi].genome_interval_tree,
                        variant, 'other', min_overlap_fraction_for_variant_matching)
                    total_count += len(other_variants) > 0
                    for other_variant, _ in other_variants:
                        all_variants.append(other_variant)
                variant.variant_frequency = float(total_count) / num_samples
                if all([_variant.is_dup for _variant in all_variants]) and len(all_variants) == num_samples:
                    variant.variant_class = 'dup'
                elif all([_variant.is_del for _variant in all_variants]) and len(all_variants) == num_samples:
                    variant.variant_class = 'del'
                else:
                    variant.variant_class = 'mixed'

    return {var_only_call_set_list[si].sample_name: var_only_call_set_list[si]
            for si in range(len(var_only_call_set_list))}, included_intervals
