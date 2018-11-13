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
    plot_number_of_events_distribution(output_dir, gcnv_segment_vcfs)

def plot_performance_metrics(output_dir: str, truth_calls: str, gcnv_segment_vcfs: list, padded_interval_file: str):
    print("Reading in callsets.")
    truth_callset = TruthCallset.read_in_callset(truth_file=truth_calls)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_segment_vcfs)
    considered_intervals = IntervalCollection.read_interval_list(padded_interval_file)
    print("Evaluating the callset against the truth.")
    evaluator = Evaluator(evaluation_name="test_eval",
                          considered_intervals=considered_intervals,
                          samples_to_evaluate=gcnv_callset.sample_names)
    result = evaluator.evaluate_callsets(callset_truth=truth_callset, callset_to_evaluate=gcnv_callset)
    print("Result is: \n")
    print(result)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                    help='output directory')

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                    help='Segment VCFs output by gCNV')

    parser.add_argument('--truth_calls', metavar='TruthCallsFile', type=str,
                    help='file that contains truth calls on superset of samples')

    parser.add_argument('--padded_intervals', metavar='PaddedIntervalsFile', type=str,
                     help='not filtered, padded interval file')

    # parser.add_argument('--gcnv_interval_vcfs', metavar='gCNVIntervalVCF', type=str, nargs='+',
    #                 help='Interval VCFs output by gCNV')

    # parser.add_argument('--gcnv_model_shards', metavar='gCNVModelShards', type=str, nargs='+',
    #                 help='Model shard directories')

    args = parser.parse_args()
    output_dir = args.output_dir
    # truth_calls_file = args.truth_calls
    # padded_interval_file = args.padded_intervals
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    truth_calls = args.truth_calls
    padded_intervals = args.padded_intervals

    plot_gcnv_metrics(output_dir, gcnv_segment_vcfs)
    plot_performance_metrics(output_dir, truth_calls, gcnv_segment_vcfs, padded_intervals)


if __name__ == '__main__':
    main()
