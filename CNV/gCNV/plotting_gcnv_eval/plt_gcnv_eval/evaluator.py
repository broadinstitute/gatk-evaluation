from callset import Callset
from callset import EventType
from interval import Interval
from interval_collection import IntervalCollection
from evaluation_results import EvaluationResult

class Evaluator:
    """
    Computes a confusion matrix of a generated callset given a truth callset on an array of samples
    """

    def __init__(self, evaluation_name: str, considered_intervals: IntervalCollection, samples_to_evaluate: list):
        self.evaluation_name = evaluation_name
        self.considered_intervals = considered_intervals
        self.samples = samples_to_evaluate

    def evaluate_callsets(self, callset_truth: Callset, callset_to_evaluate: Callset):
        evaluation_result = EvaluationResult()
        # we assume that everything marked NO_CALL in the truth callset is a reference allele
        for sample in self.samples:
            for interval in self.considered_intervals.interval_list:

                truth_calls = callset_truth.find_intersection_with_interval(interval, sample) # this returns a list of pairs (interval, event_type)
                eval_calls = callset_to_evaluate.find_intersection_with_interval(interval, sample)
                assert (len(truth_calls) > 0 and len(eval_calls) > 0)
                assert (truth_calls[0][0].start == interval.start and eval_calls[0][0].start == interval.start)
                assert (truth_calls[-1][0].end == interval.end and eval_calls[-1][0].end == interval.end)
                current_position = interval.start
                truth_calls_current_index = 0
                eval_calls_current_index = 0

                while (current_position < interval.end):
                    truth_call = truth_calls[truth_calls_current_index][1]
                    eval_call = eval_calls[eval_calls_current_index][1]
                    truth_end = truth_calls[truth_calls_current_index][0].end
                    eval_end = eval_calls[eval_calls_current_index][0].end

                    if (truth_end == eval_end):
                        Evaluator.__evaluate_single_interval_and_update_results(interval=Interval(interval.chrom, current_position, truth_end),
                                                                                truth_call=truth_call,
                                                                                eval_call=eval_call,
                                                                                evaluation_result=evaluation_result)
                        current_position = truth_end
                        truth_calls_current_index += 1
                        eval_calls_current_index += 1
                    elif (truth_end < eval_end):
                        Evaluator.__evaluate_single_interval_and_update_results(interval=Interval(interval.chrom, current_position, truth_end),
                                                                                truth_call=truth_call,
                                                                                eval_call=eval_call,
                                                                                evaluation_result=evaluation_result)
                        current_position = truth_end
                        truth_calls_current_index += 1
                    elif (truth_end > eval_end):
                        Evaluator.__evaluate_single_interval_and_update_results(interval=Interval(interval.chrom, current_position, eval_end),
                                                                                truth_call=truth_call,
                                                                                eval_call=eval_call,
                                                                                evaluation_result=evaluation_result)
                        current_position = eval_end
                        eval_calls_current_index += 1
        return evaluation_result

    @staticmethod
    def __evaluate_single_interval_and_update_results(interval: Interval, truth_call: EventType,
                                                        eval_call: EventType, evaluation_result: EvaluationResult):
        num_bases = interval.end - interval.start
        if (eval_call == EventType.NO_CALL):
            if (truth_call == EventType.NO_CALL):
                evaluation_result.increase_no_call_match_bases(num_bases)
            else:
                evaluation_result.increase_fn(num_bases)
        elif (eval_call == EventType.DUP or eval_call == EventType.DEL):
            if (truth_call == eval_call):
                evaluation_result.increase_tp(num_bases)
            if (truth_call == EventType.NO_CALL):
                evaluation_result.increase_fp(num_bases)
        elif (eval_call == EventType.REF):
            if (truth_call == EventType.NO_CALL):
                evaluation_result.increase_tn(num_bases)
            else:
                evaluation_result.increase_fn(num_bases)