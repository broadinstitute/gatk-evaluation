import vcf
from typing import List, Optional, Set, Dict, Tuple
from .core import GenericCNVCallSet, GenericCopyNumberVariant
from intervaltree_bio import GenomeIntervalTree, IntervalTree
import logging

_logger = logging.getLogger(__name__)


def _get_mode_int_list(array: List[int]):
    """Finds the mode of an integer list."""
    _dict = {v: 0 for v in range(max(array) + 1)}
    for entry in array:
        _dict[entry] += 1
    return sorted([(k, v) for k, v in _dict.items()], key=lambda x: x[1])[-1][0]


def load_genome_strip_vcf_file(genome_strip_vcf_file: str,
                               max_records: Optional[int] = None,
                               log_frequency: int = 500,
                               allosomal_contigs: Set[str] = {'X', 'Y'},
                               autosomal_ref_copy_number: int = 2) \
        -> Tuple[Dict[str, GenericCNVCallSet], GenomeIntervalTree]:
    """Loads Genome STRiP-style .VCF file and generates a dictionary of call sets, and infers
    an interval tree of genomic loci included in the analysis.

    Args:
        genome_strip_vcf_file: input Genome STRiP .VCF file
        max_records: (optional) maximum number of records to process
        log_frequency: log for every `log_frequency` processed records
        allosomal_contigs: set of allosomal contigs; used to determine ref copy number on such
            contigs
        autosomal_ref_copy_number: ref copy number for autosomal variants

    Note:
        The following replacements must be manually applied to any Genome STRiP VCF file in order
        to become compliant with the VCF specs:
            - replace '=NA' with '=.' everywhere
            - replace '"NA"' with '"."' everywhere
            - replace '= "NA"' with '= "."' everywhere
            - change GSCNCATEGORY type from Float => String

    Returns:
        dict of call sets, included loci
    """

    # step 1. load VCF file
    with open(genome_strip_vcf_file, 'r') as f:
        vcf_reader = vcf.Reader(f)
        sample_names = vcf_reader.samples
        num_samples = len(sample_names)
        gs_call_set_list = list()
        for sample_name in sample_names:
            gs_call_set_list.append(GenericCNVCallSet(sample_name, tags={"from Genome STRiP"}))

        for record_num, record in enumerate(vcf_reader):

            if record_num % log_frequency == 0:
                _logger.info("Reading record number {0}...".format(record_num))

            if max_records is not None and record_num > max_records:
                break

            info = record.INFO
            contig = record.CHROM
            start = record.start + 1
            end = info['END']
            num_variant_samples = info['GSNNONREF']
            variant_frequency = float(num_variant_samples) / num_samples

            for si, sample_record in enumerate(record.samples):
                # recalculate quality since the VCF caps it to 99
                # for some HOMDELs, CNP can be singleton; we set the quality to the segment quality
                copy_number_probs = sample_record.data.CNP
                quality = -10 * sorted(sample_record.data.CNP)[-2] \
                    if len(copy_number_probs) > 1 else info['GSCNQUAL']

                # generate variant
                sample_var = GenericCopyNumberVariant(contig, start, end, sample_record.data.CN, quality,
                                                      variant_frequency=variant_frequency,
                                                      variant_class=None,
                                                      genotype=None)
                gs_call_set_list[si].add(sample_var)

    # step 2. set the ref copy number
    #
    # Note:
    #   Since GS VCF contains variants discovered on any sample for all samples, the majority of
    #   variant contexts are REF for any given sample. We use this observation to determine the
    #   ploidy of sex chromosomes by identifying the most commonly called copy number state in
    #   sex chromosomes for each sample. For autosomal variants, we set `ref_copy_number` to
    #   `autosomal_ref_copy_number` (default=2).
    _logger.info("Determining reference copy numbers...")
    for gs_call_set in gs_call_set_list:
        _logger.info("Sample name: {0}".format(gs_call_set.sample_name))
        for contig in gs_call_set.contig_set:
            if contig in allosomal_contigs:
                inferred_ref_for_contig = _get_mode_int_list(
                    [iv.data.var_copy_number for iv in gs_call_set.get_contig_interval_tree(contig)])
                _logger.info("contig: {0}, inferred REF_CN: {1}".format(contig, inferred_ref_for_contig))
            else:
                inferred_ref_for_contig = autosomal_ref_copy_number
            for iv in gs_call_set.get_contig_interval_tree(contig):
                iv.data.set_genotype(inferred_ref_for_contig)

    # step 4. calculate variant classes
    _logger.info("Calculating variant classes...")
    contig_set = gs_call_set_list[0].contig_set
    for contig in contig_set:
        generators = [gs_call_set.iter_in_contig(contig) for gs_call_set in gs_call_set_list]
        while True:
            try:
                variants = [generator.__next__() for generator in generators]
                assert len(set(variants)) == 1  # assert that the variants indeed come from the same interval
                all_ref = all([not variant.is_var for variant in variants])
                all_dup = all([variant.is_dup for variant in variants])
                all_del = all([variant.is_del for variant in variants])
                if all_ref:
                    for variant in variants:
                        variant.variant_class = 'ref'
                elif all_dup:
                    for variant in variants:
                        variant.variant_class = 'dup'
                elif all_del:
                    for variant in variants:
                        variant.variant_class = 'del'
                else:
                    for variant in variants:
                        variant.variant_class = 'mixed'
            except StopIteration:
                break

    # step 5. infer GS included loci (from the first sample -- they are all the same)
    _logger.info("Calculate included loci...")
    included_intervals = GenomeIntervalTree()
    for contig in contig_set:
        contig_interval_tree: IntervalTree = included_intervals[contig]
        for iv in gs_call_set_list[0].get_contig_interval_tree(contig):
            contig_interval_tree.addi(iv.begin, iv.end)
        contig_interval_tree.merge_overlaps()

    # step 6. non-variant "variants" are not needed anymore -- remove them from the interval tree
    _logger.info("Removing non-variant calls...")
    var_only_gs_call_set_list = list()
    for gs_call_set in gs_call_set_list:
        var_only_gs_call_set = GenericCNVCallSet(gs_call_set.sample_name, gs_call_set.tags)
        for contig in gs_call_set.contig_set:
            for variant in gs_call_set.iter_in_contig(contig):
                if variant.is_var:
                    var_only_gs_call_set.add(variant)
        var_only_gs_call_set.tags.add("Removed non-variant calls")
        var_only_gs_call_set_list.append(var_only_gs_call_set)

    return {var_only_gs_call_set_list[si].sample_name: var_only_gs_call_set_list[si]
            for si in range(len(var_only_gs_call_set_list))}, included_intervals
