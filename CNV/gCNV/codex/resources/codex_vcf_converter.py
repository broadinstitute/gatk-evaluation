import re
from sys import argv
from collections import defaultdict, OrderedDict, Counter
from typing import List, Dict, Tuple

# This script aggregates multiple per-contig CODEX segment files into a single
# VCF file.  It should be invoked as follows:
#
#  python codex_vcf_converter.py \
#    -segments          segments-for-contig-1.seg \
#     ...
#    -segments          segments-for-contig-N.seg \
#    -sample-name-list  sample-names.txt \
#    -codex-vcf-header  codex-vcf-header.txt \
#    -output-vcf        output.vcf
#
# Each per-contig CODEX segment file should contain calls across all samples.
# We output all calls to a single VCF file.  Each variant line in the VCF
# file gives the calls for that variant across all samples.  The contig order
# follows that of the per-contig segment files and the sample-name order is
# taken from the sample-name list.

# Each call line in a segment file contains the tab-separated fields:
call_field_keys = ['sample_name', 'chr', 'cnv', 'st_bp', 'ed_bp', 'length_kb',
                   'st_exon', 'ed_exon', 'raw_cov', 'norm_cov', 'copy_no',
                   'lratio', 'mBIC']

# Each variant line in a VCF file contains the standard tab-separated fields:
variant_standard_field_keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                               'QUAL', 'FILTER', 'INFO', 'FORMAT']


def _convert_call_fields_to_variant(
        variant_id: str,
        call_fields_per_variant_sample: List[Dict[str, str]],
        sample_names: List[str]) -> Tuple[Dict[str, str], List[str]]:
    # parse call fields from first sample for standard variant fields
    contig = call_fields_per_variant_sample[0]['chr']
    start = call_fields_per_variant_sample[0]['st_bp']

    # construct variant standard fields
    num_samples = len(sample_names)
    info_field = _construct_info_field(call_fields_per_variant_sample,
                                       num_samples)
    variant_standard_fields = OrderedDict(zip(
        variant_standard_field_keys,
        [contig, start, variant_id, 'N', '<DEL>,<DUP>', '.', '.', info_field,
         'GT:GQ:PL']))

    # construct list of per-sample genotype fields
    # (we do not use a dict here to allow for duplicate sample names)
    # for all non-variant samples, the genotype field is '.'
    # for all variant samples, we construct the corresponding genotype field
    sample_to_genotype_field_map = defaultdict(
        lambda: '.',
        zip([call_fields['sample_name']
             for call_fields in call_fields_per_variant_sample],
            map(_construct_genotype_field, call_fields_per_variant_sample)))
    genotype_fields = [sample_to_genotype_field_map[sample_name]
                       for sample_name in sample_names]

    return variant_standard_fields, genotype_fields


def _construct_info_field(
        call_fields_per_variant_sample: List[Dict[str, str]],
        num_samples: int) -> str:
    # calculate ACs and AFs
    cnv_fields_per_variant_sample = \
        [call_fields['cnv']
         for call_fields in call_fields_per_variant_sample]
    cnv_counter = Counter(cnv_fields_per_variant_sample)
    ac_del = cnv_counter['del']
    ac_dup = cnv_counter['dup']
    an = len(cnv_fields_per_variant_sample)
    af_del = '{0:.3f}'.format(float(ac_del) / num_samples)
    af_dup = '{0:.3f}'.format(float(ac_dup) / num_samples)

    # parse call fields from first sample
    call_fields = call_fields_per_variant_sample[0]
    end = call_fields['ed_bp']
    start_exon = int(call_fields['st_exon'])
    end_exon = int(call_fields['ed_exon'])
    num_targets = end_exon - start_exon + 1

    # construct info field
    info_field = \
        'AC=' + str(ac_del) + ',' + str(ac_dup) + ';' \
              + 'AF=' + str(af_del) + ',' + str(af_dup) + ';' \
              + 'AN=' + str(an) + ';' \
              + 'END=' + str(end) + ';' \
              + 'NT=' + str(num_targets)

    return info_field


def _construct_genotype_field(
        call_fields: Dict[str, str]) -> str:
    # parse call fields
    cnv = call_fields['cnv']
    lratio = int(float(call_fields['lratio']))

    # construct genotype field
    if cnv == 'del':
        gt = '1'
        pl = str(lratio) + ',0,10000'
    elif cnv == 'dup':
        gt = '2'
        pl = str(lratio) + ',10000,0'
    else:
        raise Exception('CODEX \"cnv\" field must be \"del\" or \"dup\".')
    gq = min(99, lratio)
    genotype_field = \
        gt + ':' \
           + str(gq) + ':' \
           + pl

    return genotype_field


def write_vcf(segment_files: List[str],
              sample_name_list_file: str,
              codex_vcf_header_file: str,
              output_vcf_file: str):
    vcf = open(output_vcf_file, 'w')

    # write header to VCF file
    with open(codex_vcf_header_file, 'r') as f:
        vcf.write(f.read())

    # write column names (standard columns + sample names) to VCF file
    sample_names = [line.rstrip('\n')
                    for line in open(sample_name_list_file, 'r')]
    variant_column_names = list(variant_standard_field_keys)
    variant_column_names.extend(sample_names)
    vcf.write('#' + '\t'.join(variant_column_names) + '\n')

    for segment_file in segment_files:
        # make map of variant IDs -> call fields per variant sample
        variant_to_call_fields_map = defaultdict(list)
        with open(segment_file, 'r') as f:
            next(f)     # skip header line
            for call_line in f:
                call_fields = dict(zip(call_field_keys, call_line.split('\t')))
                contig = call_fields['chr']
                start = call_fields['st_bp']
                end = call_fields['ed_bp']
                variant_id = '_'.join(['CNV', contig, start, end])
                variant_to_call_fields_map[variant_id].append(call_fields)

        # construct list of variants
        # (this is a list of (variant_standard_fields, genotype_fields) tuples)
        variants = \
            [_convert_call_fields_to_variant(
                variant_id, call_fields_per_variant_sample, sample_names)
             for variant_id, call_fields_per_variant_sample
                in variant_to_call_fields_map.items()]

        # sort list of variant tuples by POS and then by END;
        # we use regex to recover END from the INFO field
        variants.sort(
            key=lambda variant:
                (int(variant[0]['POS']),
                 int(re.search(r'END=(.*?);', variant[0]['INFO'])
                       .group(1))))

        # write variants to VCF file
        for variant in variants:
            # variant[0] is an OrderedDict so that standard fields are in order
            vcf.write('\t'.join(list(variant[0].values()) + variant[1]) + '\n')

    vcf.close()


if __name__ == '__main__':
    # parse command-line arguments
    args = {}
    while argv:
        if argv[0][0] == '-':
            if argv[0] in args:
                args[argv[0]].append(argv[1])
            else:
                args[argv[0]] = [argv[1]]
        argv = argv[1:]
    segment_files = args['-segments']
    sample_name_list_file = args['-sample-name-list'][0]
    codex_vcf_header_file = args['-codex-vcf-header'][0]
    output_vcf_file = args['-output-vcf'][0]

    write_vcf(segment_files, sample_name_list_file, codex_vcf_header_file,
              output_vcf_file)
