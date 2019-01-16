from callset import Callset
from callset import EventType
from interval import Interval
from interval_collection import IntervalCollection
from evaluation_results import EvaluationResult
from filtering import BinningFilteringStrategy


class Evaluator:
    """
    Computes a confusion matrix of a generated callset given a truth callset on an array of samples
    """

    def __init__(self, evaluation_name: str, considered_intervals: IntervalCollection,
                 blacklisted_intervals_truth: IntervalCollection, samples_to_evaluate: list,
                 binning_strategy: BinningFilteringStrategy):
        self.evaluation_name = evaluation_name
        self.considered_intervals = considered_intervals
        self.blacklisted_intervals_truth = blacklisted_intervals_truth
        self.samples = samples_to_evaluate
        self.binning_strategy = binning_strategy
        self.callset_filters = binning_strategy.filter_list

    def evaluate_callset(self, callset_truth: Callset, callset_to_evaluate: Callset) -> EvaluationResult:
        """
        Note: we assume that everything marked NO_CALL in the truth callset is a reference allele

        Args:
            callset_truth: the truth callset that does not necessarily contains attributes
            callset_to_evaluate: callset to evaluate

        Returns:
            Evaluation result object that contains confusion matrices for various filtering criteria
        """

        evaluation_result = EvaluationResult(self.callset_filters)

        for sample in self.samples:
            for call_to_evaluate in callset_to_evaluate.get_calls_with_filters_applied(sample):
                filter_index = self.binning_strategy.get_filter_index_by_values(call_to_evaluate.call_attributes)
                for truth_interval, truth_call_event_type in callset_truth.find_intersection_with_interval(
                        call_to_evaluate.interval, sample):
                    for truncated_interval in self.considered_intervals.find_intersection_with_interval_and_truncate(
                            truth_interval):
                        Evaluator._evaluate_single_interval_and_update_results(interval=truncated_interval,
                                                                               truth_call=truth_call_event_type,
                                                                               eval_call=call_to_evaluate.event_type,
                                                                               filter_index=filter_index,
                                                                               evaluation_result=evaluation_result)
        return evaluation_result

    @staticmethod
    def _evaluate_single_interval_and_update_results(interval: Interval, truth_call: EventType,
                                                     eval_call: EventType,
                                                     filter_index: int, evaluation_result: EvaluationResult):
        num_bases = interval.end - interval.start
        if eval_call == EventType.NO_CALL:
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_no_call_match_bases(num_bases, filter_index)
            else:
                evaluation_result.increase_fn(num_bases, filter_index)
        elif eval_call == EventType.DUP or eval_call == EventType.DEL:
            if truth_call == eval_call:
                evaluation_result.increase_tp(num_bases, filter_index)
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_fp(num_bases, filter_index)
        elif eval_call == EventType.REF:
            if truth_call == EventType.NO_CALL:
                evaluation_result.increase_tn(num_bases, filter_index)
            else:
                evaluation_result.increase_fn(num_bases, filter_index)
