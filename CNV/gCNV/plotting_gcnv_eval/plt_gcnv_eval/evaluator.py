from callset import Callset
from call import EventType
from interval import Interval
from interval_collection import IntervalCollection
from evaluation_results import EvaluationResult
from typing import List

import io_plt


class Evaluator:
    """
    Computes a confusion matrix of a generated callset given a truth callset on an array of samples
    """

    def __init__(self, evaluation_name: str, considered_intervals: IntervalCollection,
                 blacklisted_intervals_truth: IntervalCollection):
        """

        Args:
            evaluation_name:
            considered_intervals:
            blacklisted_intervals_truth: intervals to exclude from evaluation
        """
        self.evaluation_name = evaluation_name
        self.considered_intervals = considered_intervals - blacklisted_intervals_truth

    def evaluate_callset(self, callset_truth: Callset, callset_to_evaluate: Callset, callset_filter_names: List[str],
                         callset_filter_max_values: List[float], callset_filter_num_bins: List[int]) -> EvaluationResult:
        """
        Note: we assume that everything marked NO_CALL in the truth callset is a reference allele

        Args:
            callset_truth: the truth callset that does not necessarily contains attributes
            callset_to_evaluate: callset to evaluate
            callset_filter_names: attribute names to filter on
            callset_filter_max_values: maximum values of corresponding attributes
            callset_filter_num_bins: number of bins to split each filtering interval by

        Returns:
            Evaluation result object that contains confusion matrices fo various filtering criteria
        """
        samples_to_evaluate = list(callset_truth.sample_names.intersection(callset_to_evaluate.sample_names))
        io_plt.log("Found %d out of %d provided samples in the truth file" % (len(samples_to_evaluate), len(callset_to_evaluate.sample_names)))

        evaluation_result = EvaluationResult(callset_filter_names, callset_filter_max_values, callset_filter_num_bins)

        for sample in samples_to_evaluate:
            for call_to_evaluate in callset_to_evaluate.sample_to_calls_map[sample]:
                for truth_interval, truth_call_event_type in callset_truth.find_intersection_with_interval(
                        call_to_evaluate.interval, sample):
                    for truncated_interval in self.considered_intervals.find_intersection_with_interval_and_truncate(
                            truth_interval):
                        Evaluator._evaluate_single_interval_and_update_results(interval=truncated_interval,
                                                                               truth_call=truth_call_event_type,
                                                                               eval_call=call_to_evaluate.event_type,
                                                                               call_attributes=call_to_evaluate.call_attributes,
                                                                               evaluation_result=evaluation_result)
        return evaluation_result

    @staticmethod
    def _evaluate_single_interval_and_update_results(interval: Interval, truth_call: EventType,
                                                     eval_call: EventType,
                                                     call_attributes: dict, evaluation_result: EvaluationResult):
        num_bases = interval.end - interval.start
        assert num_bases > 0
        if eval_call == EventType.NO_CALL:
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_no_call_match_bases(num_bases, call_attributes)
            else:
                evaluation_result.increase_fn(num_bases)
        elif eval_call == EventType.DUP or eval_call == EventType.DEL:
            if truth_call == eval_call:
                evaluation_result.increase_tp(num_bases, call_attributes)
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_fp(num_bases, call_attributes)
        elif eval_call == EventType.REF:
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_tn(num_bases)
            else:
                evaluation_result.increase_fn(num_bases)
