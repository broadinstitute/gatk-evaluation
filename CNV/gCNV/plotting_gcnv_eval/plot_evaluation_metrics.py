##
# CLI for making plots of metrics and performance given results of a gCNV workflow(case or cohort) 
##
import argparse
from plot_metrics import plot_number_of_events_distribution
from callset import TruthCallset, GCNVCallset
from interval_collection import IntervalCollection
from evaluator import Evaluator
from evaluation_results import EvaluationResult

def plot_gcnv_metrics(output_dir: str, gcnv_segment_vcfs: list):
    plot_number_of_events_distribution(output_dir + "event_number_distribution.png", gcnv_segment_vcfs)

def plot_performance_metrics_and_write_results(output_dir: str, truth_calls: str,
                             gcnv_segment_vcfs: list, padded_interval_file: str,
                             confusion_matrix_output_file: str, f_measure_output_file: str):
    print("Reading in callsets.")
    truth_callset = TruthCallset.read_in_callset(truth_file=truth_calls, interval_file=padded_interval_file)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_segment_vcfs)
    considered_intervals = IntervalCollection.read_interval_list(padded_interval_file)
    print("Evaluating the callset against the truth.")
    evaluator = Evaluator(evaluation_name="test_eval",
                          considered_intervals=considered_intervals,
                          samples_to_evaluate=gcnv_callset.sample_names)
    result = evaluator.evaluate_callsets(callset_truth=truth_callset, callset_to_evaluate=gcnv_callset)
    result.write_f1_measure(f_measure_output_file)
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
                     help='A file to write resulting confusion matrix to')

    parser.add_argument('--f_measure_output_file', metavar='FMeasureOutput', type=str,
                    help='A file to write F1 measure of the resulting confusion matric to')

    # parser.add_argument('--gcnv_interval_vcfs', metavar='gCNVIntervalVCF', type=str, nargs='+',
    #                 help='Interval VCFs output by gCNV')

    # parser.add_argument('--gcnv_model_shards', metavar='gCNVModelShards', type=str, nargs='+',
    #                 help='Model shard directories')

    args = parser.parse_args()
    output_dir = args.output_dir
    # truth_calls_file = args.truth_calls
    # padded_interval_file = args.padded_intervals
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    truth_calls = args.sorted_truth_calls_bed
    padded_intervals = args.padded_intervals
    confusion_matrix_output_file = args.confusion_matrix_output
    f_measure_output_file = args.f_measure_output_file

    plot_gcnv_metrics(output_dir, gcnv_segment_vcfs)
    plot_performance_metrics_and_write_results(output_dir, truth_calls, gcnv_segment_vcfs, padded_intervals, confusion_matrix_output_file, f_measure_output_file)


if __name__ == '__main__':
    main()
