import vcf
from typing import Optional, Dict, Tuple, List
from .core import GenericCNVCallSet, GenericCopyNumberVariant, get_overlapping_variants_set
from intervaltree_bio import GenomeIntervalTree
import logging

_logger = logging.getLogger(__name__)


def load_gcnv_segments_vcf_file(gcnv_segments_vcf_files: List[str],
                                interval_list_file: str,
                                max_records: Optional[int] = None,
                                log_frequency: int = 500,
                                quality_mode: str = "some",
                                min_overlap_fraction_for_variant_matching: float = 0.5) \
        -> Tuple[Dict[str, GenericCNVCallSet], GenomeIntervalTree]:
    """Loads gCNV-style segments VCF files and generates a dictionary of call sets.

    Args:
        gcnv_segments_vcf_files: a list of gCNV segments VCF files
        interval_list_file: interval list used to create the call-set
        max_records: (optional) maximum number of records to process
        log_frequency: log for every `log_frequency` processed records
        quality_mode: quality calculation mode. there are two acceptable options:
            "all": use "QA" from the VCF
            "some": use "QS" from the VCF
        min_overlap_fraction_for_variant_matching: used for finding overlapping variants (to set
            variant frequency and variant class)

    Returns:
        dict of call sets, included loci
    """

    _logger.warning("Assuming QS is normalized -- if you are using a later GATK, fix the evaluation code!")

    allele_to_genotype_map = {'N': 'ref', '<DEL>': 'del', '<DUP>': 'dup'}
    assert quality_mode in ["some", "all"]

    # step 1. load VCF file
    gcnv_call_set_list: List[GenericCNVCallSet] = list()
    sample_names = list()
    for gcnv_segments_vcf_file in gcnv_segments_vcf_files:
        with open(gcnv_segments_vcf_file, 'r') as f:
            vcf_reader = vcf.Reader(f)
            sample_name = vcf_reader.samples[0]
            _logger.info("Processing VCF for sample {0}...".format(sample_name))
            sample_call_set = GenericCNVCallSet(sample_name, tags={"from gCNV segments VCF"})
            gcnv_call_set_list.append(sample_call_set)
            sample_names.append(sample_name)

            for record_num, record in enumerate(vcf_reader):

                if record_num % log_frequency == 0:
                    _logger.info("Reading record number {0}...".format(record_num))

                if max_records is not None and record_num > max_records:
                    break

                info = record.INFO
                contig = record.CHROM
                start = record.start + 1
                end = info['END']
                var_copy_number = record.samples[0]['CN']
                genotype_index = int(record.samples[0]['GT'])
                allele = str(record.alleles[genotype_index])
                num_intervals = record.samples[0]['NP']
                genotype = allele_to_genotype_map[allele]

                if quality_mode == "some":
                    quality = num_intervals * record.samples[0]["QS"]
                elif quality_mode == "all":
                    quality = record.samples[0]["QA"]
                else:
                    raise Exception("Unknown quality calculation mode -- valid options are \"some\" and \"all\"")

                # generate variant
                sample_var = GenericCopyNumberVariant(contig, start, end, var_copy_number, quality,
                                                      genotype=genotype,
                                                      num_intervals=num_intervals,
                                                      variant_frequency=None,
                                                      variant_class=None)
                sample_call_set.add(sample_var)

            _logger.info("Done!")

    # step 2. included loci
    _logger.info("Obtaining included loci...")
    included_intervals = GenomeIntervalTree()
    with open(interval_list_file, 'r') as f:
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
    for sample_call_set in gcnv_call_set_list:
        var_only_call_set = GenericCNVCallSet(sample_call_set.sample_name, sample_call_set.tags)
        for contig in sample_call_set.contig_set:
            for variant in sample_call_set.iter_in_contig(contig):
                if variant.is_var:
                    var_only_call_set.add(variant)
        var_only_call_set.tags.add("Removed non-variant calls")
        var_only_call_set_list.append(var_only_call_set)

    # step 4. estimate variant frequency
    _logger.info("Calculating variant classes and frequencies...")
    contig_set = gcnv_call_set_list[0].contig_set
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
