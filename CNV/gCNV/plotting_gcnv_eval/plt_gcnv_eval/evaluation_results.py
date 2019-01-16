from typing import List
from enum import Enum
import pandas as pd
import numpy as np

import constants
from filtering import CallsetFilter


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

    def __init__(self, filters: List[CallsetFilter]):
        # extra column for the f1 score
        num_cols = len(list(ConfusionType)) + 1
        num_rows = len(filters)
        self.filters = filters

        self.filter_names = [str(f) for f in filters]
        self.column_names = [t.name for t in list(ConfusionType)]
        self.f1_score_col = constants.F1_SCORE_COLUMN_NAME
        self.column_names.append(self.f1_score_col)
        self.confusion_matrix = np.zeros((num_rows, num_cols))
        self.confusion_matrix_pd = None

    def increase_tp(self, num_bases: int, filtering_bin_index: int):
        self.confusion_matrix[filtering_bin_index, ConfusionType.TP_BASES.value] += num_bases

    def increase_tn(self, num_bases: int, filtering_bin_index: int):
        self.confusion_matrix[filtering_bin_index, ConfusionType.TN_BASES.value] += num_bases

    def increase_fp(self, num_bases: int, filtering_bin_index: int):
        self.confusion_matrix[filtering_bin_index, ConfusionType.FP_BASES.value] += num_bases

    def increase_fn(self, num_bases: int, filtering_bin_index: int):
        self.confusion_matrix[filtering_bin_index, ConfusionType.FN_BASES.value] += num_bases

    def increase_no_call_match_bases(self, num_bases: int, filtering_bin_index: int):
        self.confusion_matrix[filtering_bin_index, ConfusionType.NO_CALL_MATCH_BASES.value] += num_bases

    def _get_recall(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fn_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.FN_BASES.name]
        overall = tp_bases + fn_bases
        return tp_bases / overall if overall else 0.0

    def _get_precision(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.FP_BASES.name]
        overall = tp_bases + fp_bases
        return tp_bases / overall if overall else 0.0

    def _get_false_positive_rate(self, filter_to_use: str):
        fp_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.FP_BASES.name]
        tn_bases = self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.TN_BASES.name] + \
            self.confusion_matrix_pd.loc[filter_to_use, ConfusionType.NO_CALL_MATCH_BASES.name]
        overall = fp_bases + tn_bases
        return fp_bases / overall if overall else 0.0

    def _initialize_pandas_dataframe(self):
        if self.confusion_matrix_pd is None:
            self.confusion_matrix_pd = pd.DataFrame(index=self.filter_names, data=self.confusion_matrix,
                                                    columns=self.column_names)

    def compute_f1_measures(self):
        self._initialize_pandas_dataframe()
        for filter_name in self.filter_names:
            recall = self._get_recall(filter_name)
            precision = self._get_precision(filter_name)
            if precision == 0 or recall == 0:
                f1_score = 0.0
            else:
                f1_score = 2 * precision * recall / (precision + recall)

            self.confusion_matrix_pd.loc[filter_name, self.f1_score_col] = f1_score

    def write_area_under_roc_to_file(self, file: str, filters: List[CallsetFilter]):
        self._initialize_pandas_dataframe()
        assert self.filters is not None and len(self.filters) > 0
        if filters is None:
            assert len(self.filters[0].attributes) == 1, "Ambiguous filtering strategy for calculating ROC chosen"
            filters = self.filters

        sensitivity_values = [self._get_recall(str(f)) for f in filters]
        false_positive_rate_values = [self._get_false_positive_rate(str(f)) for f in filters]
        # append the corner points of the ROC
        sensitivity_values.insert(0, 0.0)
        sensitivity_values.append(1.0)
        false_positive_rate_values.insert(0, 0.0)
        false_positive_rate_values.append(1.0)
        print(sensitivity_values)
        print(false_positive_rate_values)
        area_under_roc = np.trapz(y=sensitivity_values, x=false_positive_rate_values)
        with open(file, 'w') as output_file:
            output_file.write(str(area_under_roc))

    def write_result(self, file: str):
        self.confusion_matrix_pd.to_csv(path_or_buf=file, index=True, sep='\t')
