##
# CLI for making plots of metrics and performance given results of a gCNV workflow(case or cohort) 
##
import argparse
import gzip
import shutil
import pybedtools

import plot_metrics
from callset import TruthCallset, GCNVCallset
from interval_collection import IntervalCollection
from evaluator import Evaluator
from reference_dictionary import ReferenceDictionary
import io_plt


def plot_gcnv_metrics(output_dir: str, gcnv_segment_vcfs: list, gcnv_model_dir: str, num_model_shards: int):
    plot_metrics.plot_number_of_events_distribution(output_dir, gcnv_segment_vcfs)
    plot_metrics.plot_quality_metric_distribution(output_dir, gcnv_segment_vcfs)
    plot_metrics.plot_event_size_distribution(output_dir, gcnv_segment_vcfs)
    plot_metrics.plot_ard_components(output_dir, gcnv_model_dir, num_model_shards)


def evaluate_performance_metrics_and_write_results(truth_calls: str, ref_dict_file: str, gcnv_segment_vcfs: list,
                                                   padded_interval_file: str, blacklisted_intervals_truth: str,
                                                   callset_filter_names: list, callset_filter_max_values: list,
                                                   callset_filter_num_bins: list, attribute_for_roc_creation: str,
                                                   output_dir: str, truth_allele_frequency_threshold: float):
    io_plt.log("Reading in callsets.")
    ref_dict = ReferenceDictionary.read_in(ref_dict_file)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_segment_vcfs, reference_dictionary=ref_dict)
    truth_callset = TruthCallset.read_in_callset(truth_file=truth_calls,
                                                 interval_file=padded_interval_file,
                                                 reference_dictionary=ref_dict,
                                                 allele_frequency_threshold=truth_allele_frequency_threshold)
    considered_intervals = IntervalCollection.read_interval_list(padded_interval_file)
    blacklisted_intervals_truth = IntervalCollection.read_interval_list(blacklisted_intervals_truth)
    io_plt.log("Evaluating the callset against the truth.")
    evaluator = Evaluator(evaluation_name="test_eval",
                          considered_intervals=considered_intervals,
                          blacklisted_intervals_truth=blacklisted_intervals_truth)
    result = evaluator.evaluate_callset(callset_truth=truth_callset,
                                        callset_to_evaluate=gcnv_callset,
                                        callset_filter_names=callset_filter_names,
                                        callset_filter_max_values=callset_filter_max_values,
                                        callset_filter_num_bins=callset_filter_num_bins)

    result.compute_f1_measures()
    result.write_performance_curves_to_file(output_dir, attribute_for_roc_creation)
    result.write_results(output_dir)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='output directory')

    parser.add_argument('--ref_dict', metavar='ReferenceDictionary', type=str,
                        help='Dict reference file')

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                        help='Segment VCFs output by gCNV')

    parser.add_argument('--sorted_truth_calls_bed', metavar='SortedTruthCallsBed', type=str,
                        help='Sorted bed file that contains truth calls on superset of samples that are being evaluated')

    parser.add_argument('--padded_intervals', metavar='PaddedIntervalsFile', type=str,
                        help='Not filtered, padded interval file')

    parser.add_argument('--blacklisted_intervals_truth', metavar='BlacklistTruth', type=str,
                        help='Blacklisted regions that were used to create truth callset')

    parser.add_argument('--callset_filter_names', metavar='CallsetFilterNames', type=str, nargs='+',
                        help='Genotype fields to filter on')

    parser.add_argument('--callset_filter_max_values', metavar='CallsetFilterMaxValues', type=float, nargs='+',
                        help='Max cutoff filter values for the corresponding genotype fields')

    parser.add_argument('--callset_filter_num_bins', metavar='CallsetFilterNumberOfBins', type=int, nargs='+',
                        help='Number of filter bins for the corresponding genotype field.')

    parser.add_argument('--attribute_for_roc_creation', metavar="AttributeROC", type=str,
                        help='Which genotype attribute for filtering to create ROC.')

    parser.add_argument('--gcnv_models_directory', metavar='gCNVModelDirectory', type=str,
                        help='Directory with subdirectories containing model shard results, where each '
                             'shard subdirectory is is named \shard-X where X is the index of the shard ')

    parser.add_argument('--num_model_shards', metavar='NumberOfModelShards', type=int,
                        help='Number of model shards')

    parser.add_argument('--truth_allele_frequency_threshold', metavar="TruthAlleleFrequencyThreshold", type=float,
                        help='Specify allele frequency threshold above which the calls will be filtered out in the'
                             'truth callset')

    ###################
    # Parse arguments #
    ###################
    args = parser.parse_args()

    # input arguments
    ref_dict_file = args.ref_dict
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    truth_calls = args.sorted_truth_calls_bed
    padded_intervals = args.padded_intervals
    blacklisted_intervals_truth = args.blacklisted_intervals_truth
    truth_allele_frequency_threshold = args.truth_allele_frequency_threshold

    # output arguments
    output_dir = args.output_dir

    # filtering arguments
    # TODO handle NoneType, i.e. if user does not provide these
    callset_filter_names = args.callset_filter_names
    callset_filter_max_values = args.callset_filter_max_values
    callset_filter_num_bins = args.callset_filter_num_bins
    attribute_for_roc_creation = args.attribute_for_roc_creation

    # arguments for model plots
    gcnv_model_dir = args.gcnv_models_directory
    num_model_shards = args.num_model_shards

    # Handle gzip-ed VCF gcnv files
    new_gcnv_vcfs_list = []
    for vcf_file in gcnv_segment_vcfs:
        if vcf_file.endswith(".gz"):
            new_vcf_file = vcf_file[:-3]
            with open(new_vcf_file, 'wt') as f_out, gzip.open(vcf_file, 'rt') as f_in:
                shutil.copyfileobj(f_in, f_out)
                vcf_file = new_vcf_file
        new_gcnv_vcfs_list.append(vcf_file)
    gcnv_segment_vcfs = new_gcnv_vcfs_list

    # Evaluate performance and produce metrics plots
    plot_gcnv_metrics(output_dir, gcnv_segment_vcfs, gcnv_model_dir, num_model_shards)
    evaluate_performance_metrics_and_write_results(truth_calls, ref_dict_file, gcnv_segment_vcfs, padded_intervals,
                                                   blacklisted_intervals_truth, callset_filter_names,
                                                   callset_filter_max_values, callset_filter_num_bins,
                                                   attribute_for_roc_creation, output_dir, truth_allele_frequency_threshold)
    io_plt.log("SUCCESS")


if __name__ == '__main__':
    main()
