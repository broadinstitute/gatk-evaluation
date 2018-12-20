##
# CLI for making plots of metrics and performance given results of a gCNV workflow(case or cohort) 
##
import argparse
import gzip
import shutil

from plot_metrics import plot_number_of_events_distribution
from callset import TruthCallset, GCNVCallset
from interval_collection import IntervalCollection
from evaluator import Evaluator
from filtering import BinningFilteringStrategy
import io_plt


def plot_gcnv_metrics(output_dir: str, gcnv_segment_vcfs: list):
    plot_number_of_events_distribution(output_dir + "event_number_distribution.png", gcnv_segment_vcfs)


def evaluate_performance_metrics_and_write_results(truth_calls: str, gcnv_segment_vcfs: list,
                                                   padded_interval_file: str, blacklisted_intervals_truth: str,
                                                   confusion_matrix_output_file: str, area_under_curve_output: str,
                                                   callset_filter_names: list, callset_filter_max_values: list,
                                                   callset_filter_num_bins: list, attribute_for_roc_creation: str):
    io_plt.log("Reading in callsets.")
    truth_callset = TruthCallset.read_in_callset(truth_file=truth_calls, interval_file=padded_interval_file)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_segment_vcfs)
    considered_intervals = IntervalCollection.read_interval_list(padded_interval_file)
    blacklisted_intervals_truth = IntervalCollection.read_interval_list(blacklisted_intervals_truth)
    binning_strategy = BinningFilteringStrategy(attributes=callset_filter_names,
                                                attributes_max_values=callset_filter_max_values,
                                                attributes_num_bins=callset_filter_num_bins)
    io_plt.log("Evaluating the callset against the truth.")
    evaluator = Evaluator(evaluation_name="test_eval",
                          considered_intervals=considered_intervals,
                          blacklisted_intervals_truth=blacklisted_intervals_truth,
                          samples_to_evaluate=gcnv_callset.sample_names,
                          binning_strategy=binning_strategy)
    result = evaluator.evaluate_callsets(callset_truth=truth_callset, callset_to_evaluate=gcnv_callset)
    result.compute_f1_measures()
    result.write_area_under_roc_to_file(area_under_curve_output,
                                        binning_strategy.get_single_attribute_filters(attribute_for_roc_creation))
    result.write_result(confusion_matrix_output_file)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='output directory')

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                        help='Segment VCFs output by gCNV')

    parser.add_argument('--sorted_truth_calls_bed', metavar='SortedTruthCallsBed', type=str,
                        help='Sorted bed file that contains truth calls on superset of samples that are being evaluated')

    parser.add_argument('--padded_intervals', metavar='PaddedIntervalsFile', type=str,
                        help='Not filtered, padded interval file')

    parser.add_argument('--confusion_matrix_output', metavar='ConfusionMatrixOutput', type=str,
                        help='A file to write the resulting confusion matrix')

    parser.add_argument('--area_under_curve_output', metavar='AUCOutput', type=str,
                        help='A file to write the area under ROC curve;'
                             ' first filter in list will be used to create ROC')

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

    # parser.add_argument('--gcnv_model_shards', metavar='gCNVModelShards', type=str, nargs='+',
    #                 help='Model shard directories')

    ###################
    # Parse arguments #
    ###################
    args = parser.parse_args()

    output_dir = args.output_dir
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    truth_calls = args.sorted_truth_calls_bed
    padded_intervals = args.padded_intervals
    confusion_matrix_output_file = args.confusion_matrix_output
    area_under_curve_output = args.area_under_curve_output
    blacklisted_intervals_truth = args.blacklisted_intervals_truth
    # filtering arguments
    # TODO handle NoneType, i.e. if user does not provide these
    callset_filter_names = args.callset_filter_names
    callset_filter_max_values = args.callset_filter_max_values
    callset_filter_num_bins = args.callset_filter_num_bins
    attribute_for_roc_creation = args.attribute_for_roc_creation

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
    plot_gcnv_metrics(output_dir, gcnv_segment_vcfs)
    evaluate_performance_metrics_and_write_results(truth_calls, gcnv_segment_vcfs, padded_intervals,
                                                   blacklisted_intervals_truth, confusion_matrix_output_file,
                                                   area_under_curve_output, callset_filter_names,
                                                   callset_filter_max_values, callset_filter_num_bins,
                                                   attribute_for_roc_creation)
    io_plt.log("SUCCESS")


if __name__ == '__main__':
    main()
