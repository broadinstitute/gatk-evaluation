import vcf
import numpy as np
from typing import Optional, Set, Dict, Tuple
from .core import GenericCNVCallSet, GenericCopyNumberVariant
from intervaltree_bio import GenomeIntervalTree
import gcnvkernel
import logging

_logger = logging.getLogger(__name__)


def load_sex_genotypes_tsv_file(sex_genotypes_tsv_file: str) -> Dict[str, str]:
    res = dict()
    with open(sex_genotypes_tsv_file, 'r') as f:
        for ln, line in enumerate(f):
            if ln == 0:
                continue
            parsed = line.strip().split()
            res[parsed[0]] = parsed[1]
    return res


def load_gcnv_spark_vcf_file(gcnv_spark_vcf_file: str,
                             sex_genotypes_tsv_file: str,
                             targets_tsv_file: str,
                             max_records: Optional[int] = None,
                             log_frequency: int = 500,
                             allosomal_contigs: Set[str] = {'X', 'Y'},
                             autosomal_ref_copy_number: int = 2) \
        -> Tuple[Dict[str, GenericCNVCallSet], GenomeIntervalTree]:
    """Loads gCNV-Spark .VCF file and generates a dictionary of call sets, and infers
    an interval tree of genomic loci included in the analysis.

    Args:
        gcnv_spark_vcf_file: input gCNV Spark .VCF file
        sex_genotypes_tsv_file: sex genotypes table
        targets_tsv_file: targets interval list
        max_records: (optional) maximum number of records to process
        log_frequency: log for every `log_frequency` processed records
        allosomal_contigs: set of allosomal contigs; used to determine ref copy number on such
            contigs
        autosomal_ref_copy_number: ref copy number for autosomal variants

    Returns:
        dict of call sets, included loci
    """

    # step 0. load sex genotypes .tsv file and targets .tsv file
    sex_genotypes_dict = load_sex_genotypes_tsv_file(sex_genotypes_tsv_file)
    targets_interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(targets_tsv_file)

    # map DEL/REF/DUP call to copy numbers
    call_to_cn_map = {0: 2, 1: 0, 2: 3}

    # map sex genotypes to contig ploidy
    sex_genotypes_ploidy_map = {
        'SEX_XX': {'X': 2, 'Y': 0},
        'SEX_XY': {'X': 1, 'Y': 1}
    }

    # step 1. load VCF file
    with open(gcnv_spark_vcf_file, 'r') as f:
        vcf_reader = vcf.Reader(f)
        sample_names = vcf_reader.samples
        gcnv_call_set_list = list()
        for sample_name in sample_names:
            gcnv_call_set_list.append(GenericCNVCallSet(sample_name, tags={"from gCNV-Spark"}))
    
        for record_num, record in enumerate(vcf_reader):
            if record_num % log_frequency == 0:
                _logger.info("Reading record number {0}...".format(record_num))
    
            if max_records is not None and record_num > max_records:
                break
    
            info = record.INFO
            contig = record.CHROM
            start = record.start
            end = info['END']
            num_intervals = info['NTARGETS']

            for si, sample_record in enumerate(record.samples):
                pl = sample_record.data.PL
                if pl is None:
                    var_copy_number = -1
                    quality = -1
                else:
                    var_copy_number = call_to_cn_map[int(np.argmin(pl))]
                    quality = sorted(pl)[1]

                # generate variant
                sample_var = GenericCopyNumberVariant(contig, start, end,
                                                      var_copy_number, quality,
                                                      genotype=None,
                                                      variant_frequency=None,
                                                      variant_class=None,
                                                      num_intervals=num_intervals)
                gcnv_call_set_list[si].add(sample_var)
    
    # step 2. set the ref copy number
    _logger.info("Determining reference copy numbers...")
    for gcnv_call_set in gcnv_call_set_list:
        _logger.info("Sample name: {0}".format(gcnv_call_set.sample_name))
        for contig in gcnv_call_set.contig_set:
            if contig in allosomal_contigs:
                inferred_ref_for_contig = sex_genotypes_ploidy_map[
                    sex_genotypes_dict[gcnv_call_set.sample_name]][contig]
                _logger.info("contig: {0}, inferred REF_CN: {1}".format(contig, inferred_ref_for_contig))
            else:
                inferred_ref_for_contig = autosomal_ref_copy_number
            for iv in gcnv_call_set.get_contig_interval_tree(contig):
                iv.data.set_genotype(inferred_ref_for_contig)
    
    # step 3. calculate variant frequencies
    _logger.info("Calculating variant frequencies...")
    contig_set = gcnv_call_set_list[0].contig_set
    for contig in contig_set:
        generators = [gcnv_call_set.iter_in_contig(contig) for gcnv_call_set in gcnv_call_set_list]
        while True:
            try:
                variants = [generator.__next__() for generator in generators]
                assert len(set(variants)) == 1  # assert that the variants indeed come from the same interval
                num_var_samples = sum([variant.is_var and variant.quality >= 0 for variant in variants])
                num_legit_vars = sum([variant.quality >= 0 for variant in variants])
                variant_frequency = float(num_var_samples) / num_legit_vars
                for variant in variants:
                    variant.variant_frequency = variant_frequency
            except StopIteration:
                break
    
    # step 4. calculate variant classes
    _logger.info("Calculating variant classes...")
    contig_set = gcnv_call_set_list[0].contig_set
    for contig in contig_set:
        generators = [gcnv_call_set.iter_in_contig(contig) for gcnv_call_set in gcnv_call_set_list]
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

    # step 5. included loci from the target list
    _logger.info("Calculate included loci...")
    included_intervals = GenomeIntervalTree()
    for interval in targets_interval_list:
        included_intervals[interval.contig].addi(interval.start, interval.end)

    # step 6. non-variant "variants" are not needed anymore -- remove them from the interval tree
    _logger.info("Removing non-variant calls...")
    var_only_gcnv_call_set_list = list()
    for gcnv_call_set in gcnv_call_set_list:
        var_only_gcnv_call_set = GenericCNVCallSet(gcnv_call_set.sample_name, gcnv_call_set.tags)
        for contig in gcnv_call_set.contig_set:
            for variant in gcnv_call_set.iter_in_contig(contig):
                if variant.is_var and variant.quality >= 0:
                    var_only_gcnv_call_set.add(variant)
        var_only_gcnv_call_set.tags.add("Removed non-variant calls")
        var_only_gcnv_call_set_list.append(var_only_gcnv_call_set)
    
    return {var_only_gcnv_call_set_list[si].sample_name: var_only_gcnv_call_set_list[si]
            for si in range(len(var_only_gcnv_call_set_list))}, included_intervals
