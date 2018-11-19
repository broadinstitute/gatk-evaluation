from enum import Enum
import pandas as pd
import numpy as np

class ConfusionType(Enum):
    TP_BASES = 0,
    TN_BASES = 1,
    FP_BASES = 2,
    FN_BASES = 3,
    NO_CALL_MATCH_BASES = 4

class EvaluationResult:
    """
    Class that represents confusion matrix of the evaluation
    """

    def __init__(self):
        num_columns = len(list(ConfusionType))
        column_names = [enumtype.name for enumtype in list(ConfusionType)]
        self.confusion_matrix_pd = pd.Series(index=column_names, data=np.zeros(num_columns), dtype='int64')

    def increase_tp(self, num_bases: int):
        self.confusion_matrix_pd[ConfusionType.TP_BASES.name] += num_bases

    def increase_tn(self, num_bases: int):
        self.confusion_matrix_pd[ConfusionType.TN_BASES.name] += num_bases

    def increase_fp(self, num_bases: int):
        self.confusion_matrix_pd[ConfusionType.FP_BASES.name] += num_bases

    def increase_fn(self, num_bases: int):
        self.confusion_matrix_pd[ConfusionType.FN_BASES.name] += num_bases

    def increase_no_call_match_bases(self, num_bases: int):
        self.confusion_matrix_pd[ConfusionType.NO_CALL_MATCH_BASES.name] += num_bases

    def get_recall(self):
        tp_bases = self.confusion_matrix_pd[ConfusionType.TP_BASES.name]
        fn_bases = self.confusion_matrix_pd[ConfusionType.FN_BASES.name]
        return tp_bases/(tp_bases + fn_bases)

    def get_precision(self):
        tp_bases = self.confusion_matrix_pd[ConfusionType.TP_BASES.name]
        fp_bases = self.confusion_matrix_pd[ConfusionType.FP_BASES.name]
        return tp_bases / (tp_bases + fp_bases)

    def write_f1_measure(self, file: str):
        recall = self.get_recall()
        precision = self.get_precision()

        f1_score = 2 * precision * recall/(precision + recall)

        with open(file, 'w') as output_file:
            output_file.write(str(f1_score))

    def write_result(self, file: str):
        self.confusion_matrix_pd.to_csv(path=file, index=True, sep='\t')

    def __str__(self):
        return ("TP: %d, TN: %d, FP: %d, FN: %d, NO_CALL_MATCH_BASES: %d" \
                 % (self.tp_bases, self.tn_bases, self.fp_bases, self.fn_bases, self.no_call_match_bases))