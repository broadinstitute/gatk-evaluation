from typing import List
import pandas as pd
from enum import Enum
import numpy as np

import constants
from filtering import FilteringBinCollection, LowerBoundedFilterBin


class ConfusionType(Enum):
    TP_BASES = 0
    TN_BASES = 1
    FP_BASES = 2
    FN_BASES = 3
    NO_CALL_MATCH_BASES = 4


class ClassificationResult:

    def __init__(self):
        self.tp_bases = 0
        self.tn_bases = 0
        self.fp_bases = 0
        self.fn_bases = 0
        self.no_call_match_bases = 0

    def get_ndarray(self):
        return np.array([self.tp_bases, self.tn_bases, self.fp_bases, self.fn_bases, self.no_call_match_bases])

    @staticmethod
    def get_name_header()->list:
        return [ct.name for ct in list(ConfusionType)]


class EvaluationResult:
    """
    Class that represents confusion matrix of the evaluation
    """

    def __init__(self, callset_filter_names, callset_filter_max_values, callset_filter_num_bins):
        self.filter_bin_collection = FilteringBinCollection(attributes=callset_filter_names,
                                                            attributes_max_values=callset_filter_max_values,
                                                            attributes_num_bins=callset_filter_num_bins)
        self.bounded_filters = self.filter_bin_collection.bounded_filter_list
        self.lower_bounded_filters = self.filter_bin_collection.lower_bounded_filter_list
        # extra column for the f1 score
        num_cols = len(list(ConfusionType)) + 1
        num_rows = len(self.bounded_filters)

        self.bounded_filter_names = [str(f) for f in self.bounded_filters]
        self.lower_bounded_filter_names = [str(f) for f in self.lower_bounded_filters]
        self.column_names = ClassificationResult.get_name_header()
        self.f1_score_col = constants.F1_SCORE_COLUMN_NAME
        self.column_names.append(self.f1_score_col)

        self.confusion_matrix_bounded_filters_pd = None
        self.confusion_matrix_lower_bounded_filters_pd = None

        self.confusion_matrix_bounded_filters = {key: ClassificationResult() for key in self.bounded_filters}

    def increase_tp(self, num_bases: int, attribute_to_value_map: dict):
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].tp_bases += num_bases

    def increase_tn(self, num_bases: int, attribute_to_value_map: dict):
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].tn_bases += num_bases

    def increase_fp(self, num_bases: int, attribute_to_value_map: dict):
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].fp_bases += num_bases

    def increase_fn(self, num_bases: int, attribute_to_value_map: dict):
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].fn_bases += num_bases

    def increase_no_call_match_bases(self, num_bases: int, attribute_to_value_map: dict):
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        assert filter_bin is not None, filter_bin
        self.confusion_matrix_bounded_filters[filter_bin].no_call_match_bases += num_bases

    def _get_recall(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fn_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.FN_BASES.name]
        overall = tp_bases + fn_bases
        return tp_bases / overall if overall else 0.0

    def _get_precision(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.TP_BASES.name]
        fp_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.FP_BASES.name]
        overall = tp_bases + fp_bases
        return tp_bases / overall if overall else 0.0

    def _get_false_positive_rate(self, filter_to_use: str):
        fp_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.FP_BASES.name]
        tn_bases = self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.TN_BASES.name] + \
            self.confusion_matrix_lower_bounded_filters_pd.loc[filter_to_use, ConfusionType.NO_CALL_MATCH_BASES.name]
        overall = fp_bases + tn_bases
        return fp_bases / overall if overall else 0.0

    def _initialize_pandas_dataframes(self):
        """

        Returns:

        """
        if self.confusion_matrix_bounded_filters_pd is None:
            confusion_matrix = np.zeros(shape=(len(self.bounded_filters), len(self.column_names)), dtype=np.int32)
            for idx, bounded_filter in enumerate(self.bounded_filters):
                confusion_entry_for_filter = self.confusion_matrix_bounded_filters[bounded_filter].get_ndarray()
                confusion_matrix[idx, :confusion_entry_for_filter.size] = confusion_entry_for_filter

            self.confusion_matrix_bounded_filters_pd = pd.DataFrame(index=self.bounded_filter_names, data=confusion_matrix,
                                                                    columns=self.column_names)
            self._compute_lower_bounded_confusion_matrix()

    def _compute_lower_bounded_confusion_matrix(self):
        """

        Returns:

        """
        if self.confusion_matrix_lower_bounded_filters_pd is None:
            confusion_matrix = np.zeros(shape=(len(self.lower_bounded_filters), len(self.column_names)), dtype=np.int32)
            for idx, lower_bounded_filter in enumerate(self.lower_bounded_filters):
                confusion_entry_sum = ClassificationResult().get_ndarray()
                for bounded_filter in self.filter_bin_collection.get_lower_bounded_filter_sum(lower_bounded_filter):
                    confusion_entry_for_filter = self.confusion_matrix_bounded_filters[bounded_filter].get_ndarray()
                    confusion_entry_sum += confusion_entry_for_filter
                confusion_matrix[idx, :confusion_entry_sum.size] = confusion_entry_sum
            self.confusion_matrix_lower_bounded_filters_pd = pd.DataFrame(index=self.lower_bounded_filter_names, data=confusion_matrix,
                                                                          columns=self.column_names)

    def compute_f1_measures(self):
        self._initialize_pandas_dataframes()
        for filter_name in self.lower_bounded_filter_names:
            recall = self._get_recall(filter_name)
            precision = self._get_precision(filter_name)
            if precision == 0 or recall == 0:
                f1_score = 0.0
            else:
                f1_score = 2 * precision * recall / (precision + recall)

            self.confusion_matrix_lower_bounded_filters_pd.loc[filter_name, self.f1_score_col] = f1_score

    def write_area_under_roc_to_file(self, file: str, attribute_for_filtering: str):
        self._initialize_pandas_dataframes()
        assert self.lower_bounded_filters is not None and len(self.lower_bounded_filters) > 0
        filters_list = self.filter_bin_collection.get_single_attribute_lower_bounded_filters(attribute_for_filtering)

        if filters_list is None:
            assert len(self.lower_bounded_filters[0].attributes) == 1, "Ambiguous filtering strategy for calculating ROC chosen"
            filters_list = self.lower_bounded_filters

        sensitivity_values = [self._get_recall(str(f)) for f in filters_list]
        false_positive_rate_values = [self._get_false_positive_rate(str(f)) for f in filters_list]
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
        self.confusion_matrix_bounded_filters_pd.to_csv(path_or_buf=file, index=True, sep='\t')
