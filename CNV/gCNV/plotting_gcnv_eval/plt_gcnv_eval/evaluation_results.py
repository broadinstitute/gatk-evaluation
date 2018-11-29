from enum import Enum
import pandas as pd
import numpy as np

import constants

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

    def __init__(self, filters: list):
        num_cols = len(list(ConfusionType)) + 1 # extra column for the f1 score
        num_rows = len(filters)

        self.filter_names = [str(f) for f in filters]
        self.column_names = [enumtype.name for enumtype in list(ConfusionType)]
        self.f1_score_col = constants.F1_SCORE_COLUMN_NAME
        self.column_names.append(self.f1_score_col)
        self.confusion_matrix_pd = pd.DataFrame(index=self.filter_names, data=np.zeros((num_rows, num_cols)),
                                             columns=self.column_names)
        self.current_filter_name = None

    def increase_tp(self, num_bases: int):
        self.confusion_matrix_pd.loc[self.current_filter_name, ConfusionType.TP_BASES.name] += num_bases

    def increase_tn(self, num_bases: int):
        self.confusion_matrix_pd.loc[self.current_filter_name, ConfusionType.TN_BASES.name] += num_bases

    def increase_fp(self, num_bases: int):
        self.confusion_matrix_pd.loc[self.current_filter_name, ConfusionType.FP_BASES.name] += num_bases

    def increase_fn(self, num_bases: int):
        self.confusion_matrix_pd.loc[self.current_filter_name, ConfusionType.FN_BASES.name] += num_bases

    def increase_no_call_match_bases(self, num_bases: int):
        self.confusion_matrix_pd.loc[self.current_filter_name, ConfusionType.NO_CALL_MATCH_BASES.name] += num_bases

    def get_recall(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fn_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.FN_BASES.name]
        return tp_bases / (tp_bases + fn_bases)

    def get_precision(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.FP_BASES.name]
        return tp_bases / (tp_bases + fp_bases)

    def compute_f1_measures(self):
        for filter_name in self.filter_names:
            recall = self.get_recall(filter_name)
            precision = self.get_precision(filter_name)
            if (precision == 0 or recall == 0):
                f1_score = 0.0
            else:
                f1_score = 2 * precision * recall / (precision + recall)
            
            self.confusion_matrix_pd.loc[filter_name, self.f1_score_col] = f1_score

        #with open(file, 'w') as output_file:
        #    output_file.write(str(f1_score))

    def write_result(self, file: str):
        self.confusion_matrix_pd.to_csv(path_or_buf=file, index=True, sep='\t')

    def __str__(self):
        return ("TP: %d, TN: %d, FP: %d, FN: %d, NO_CALL_MATCH_BASES: %d" \
                 % (self.tp_bases, self.tn_bases, self.fp_bases, self.fn_bases, self.no_call_match_bases))