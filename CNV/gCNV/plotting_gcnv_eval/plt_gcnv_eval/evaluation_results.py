import os
from typing import List
import pandas as pd
from enum import Enum
import numpy as np

import constants
from filtering import FilteringBinCollection
import plot_metrics


class ConfusionValueType(Enum):
    TP_BASES = 0
    TN_BASES = 1
    FP_BASES = 2
    FN_BASES = 3

    @staticmethod
    def get_name_header()->list:
        return [ct.name for ct in list(ConfusionValueType)]


class ClassificationResult:

    def __init__(self):
        """
        This class stores number of true and false positives calls (in bases) made in a specific filter bin
        """
        self.tp_bases = 0
        self.fp_bases = 0

    def __add__(self, other):
        assert type(other) == type(self)
        assert other.tp_bases >= 0 and other.fp_bases >= 0

        result = ClassificationResult()
        result.tp_bases = self.tp_bases + other.tp_bases
        result.fp_bases = self.fp_bases + other.fp_bases
        return result

    def __radd__(self, other):
        return self.__add__(other)

    def get_ndarray(self):
        return np.array([self.tp_bases, self.fp_bases])


class EvaluationResult:
    """
    Class that represents confusion matrix of the evaluation
    """

    def __init__(self, callset_filter_names: List[str], callset_filter_max_values: List[float], callset_filter_num_bins: List[int]):
        self.filter_bin_collection = FilteringBinCollection(attributes=callset_filter_names,
                                                            attributes_max_values=callset_filter_max_values,
                                                            attributes_num_bins=callset_filter_num_bins)
        self.bounded_filters = self.filter_bin_collection.bounded_filter_list
        self.lower_bounded_filters = self.filter_bin_collection.lower_bounded_filter_list

        self.bounded_filter_names = [str(f) for f in self.bounded_filters]
        self.lower_bounded_filter_names = [str(f) for f in self.lower_bounded_filters]
        self.f1_score_col = constants.F1_SCORE_COLUMN_NAME

        self.confusion_matrix_bounded_filters_df = None
        self.confusion_matrix_lower_bounded_filters_df = None

        self.confusion_matrix_bounded_filters = {key: ClassificationResult() for key in self.bounded_filters}

        self.overall_number_positive_bases = 0
        self.overall_number_negative_bases = 0

    def increase_tp(self, num_bases: int, attribute_to_value_map: dict):
        self.overall_number_positive_bases += num_bases
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].tp_bases += num_bases

    def increase_tn(self, num_bases: int):
        self.overall_number_negative_bases += num_bases

    def increase_fp(self, num_bases: int, attribute_to_value_map: dict):
        self.overall_number_negative_bases += num_bases
        filter_bin = self.filter_bin_collection.get_filter_bin_by_attribute_values(attribute_to_value_map)
        self.confusion_matrix_bounded_filters[filter_bin].fp_bases += num_bases

    def increase_fn(self, num_bases: int):
        self.overall_number_positive_bases += num_bases

    def compute_f1_measures(self):
        self._initialize_pandas_dataframes()
        for filter_name in self.lower_bounded_filter_names:
            recall = self._get_recall(filter_name)
            precision = self._get_precision(filter_name)
            if precision == 0 or recall == 0:
                f1_score = 0.0
            else:
                f1_score = 2 * precision * recall / (precision + recall)

            self.confusion_matrix_lower_bounded_filters_df.loc[filter_name, self.f1_score_col] = f1_score

    def write_area_under_roc_to_file(self, output_dir: str, attribute_for_filtering: str):
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
        area_under_roc = np.trapz(y=sensitivity_values, x=false_positive_rate_values)
        with open(os.path.join(output_dir, constants.AREA_UNDER_ROC_FILE_NAME), 'w') as output_file:
            output_file.write(str(area_under_roc))
        plot_metrics.plot_roc_curve(output_dir, false_positive_rate_values, sensitivity_values)

    def write_results(self, output_path: str):
        self.confusion_matrix_bounded_filters_df.to_csv(
            path_or_buf=open(os.path.join(output_path, constants.CONFUSION_MATRIX_WITH_BOUNDED_FILTERS), 'w'), index=True, sep='\t')
        self.confusion_matrix_lower_bounded_filters_df.to_csv(
            path_or_buf=open(os.path.join(output_path, constants.CONFUSION_MATRIX_WITH_LOWER_BOUNDED_FILTERS_NAME), 'w'), index=True, sep='\t')

    def _get_recall(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.TP_BASES.name]
        fn_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.FN_BASES.name]
        overall = tp_bases + fn_bases
        return tp_bases / overall if overall else 0.0

    def _get_precision(self, filter_to_use: str):
        tp_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.TP_BASES.name]
        fp_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.FP_BASES.name]
        overall = tp_bases + fp_bases
        return tp_bases / overall if overall else 0.0

    def _get_false_positive_rate(self, filter_to_use: str):
        fp_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.FP_BASES.name]
        tn_bases = self.confusion_matrix_lower_bounded_filters_df.loc[filter_to_use, ConfusionValueType.TN_BASES.name]
        overall = fp_bases + tn_bases
        return fp_bases / overall if overall else 0.0

    def _initialize_pandas_dataframes(self):
        """
        Initializes the full confusion matrices after processing and evaluating the callsets. In particular this method
        computes and fills in false negative and true negative entries of the confusion matrix for each specified
        filtering configuration.
        """
        self._compute_lower_bounded_confusion_matrix()
        column_names = ConfusionValueType.get_name_header()
        if self.confusion_matrix_bounded_filters_df is None:
            confusion_matrix = np.zeros(shape=(len(self.bounded_filters), len(column_names)), dtype=np.int64)
            self.confusion_matrix_bounded_filters_df = pd.DataFrame(index=self.bounded_filter_names,
                                                                    data=confusion_matrix,
                                                                    columns=column_names)
            for idx, bounded_filter in enumerate(self.bounded_filters):
                confusion_entry_for_filter = self.confusion_matrix_bounded_filters[bounded_filter]
                EvaluationResult._fill_in_dataframe_confusion_matrix_entry(
                    self.confusion_matrix_bounded_filters_df, bounded_filter, confusion_entry_for_filter,
                    self.overall_number_positive_bases, self.overall_number_negative_bases)

    def _compute_lower_bounded_confusion_matrix(self):
        """
        Helper method that initializes the confusion matrix with lower bounded filters only.
        """
        column_names = ConfusionValueType.get_name_header()
        column_names.append(constants.F1_SCORE_COLUMN_NAME)
        if self.confusion_matrix_lower_bounded_filters_df is None:
            confusion_matrix = np.zeros(shape=(len(self.lower_bounded_filters), len(column_names)), dtype=np.int64)
            self.confusion_matrix_lower_bounded_filters_df = pd.DataFrame(index=self.lower_bounded_filter_names,
                                                                          data=confusion_matrix,
                                                                          columns=column_names)
            for idx, lower_bounded_filter in enumerate(self.lower_bounded_filters):
                confusion_entry_sum = ClassificationResult()
                for bounded_filter in self.filter_bin_collection.get_lower_bounded_filter_sum(lower_bounded_filter):
                    confusion_entry_sum += self.confusion_matrix_bounded_filters[bounded_filter]

                EvaluationResult._fill_in_dataframe_confusion_matrix_entry(
                    self.confusion_matrix_lower_bounded_filters_df, lower_bounded_filter, confusion_entry_sum,
                    self.overall_number_positive_bases, self.overall_number_negative_bases)

    @staticmethod
    def _fill_in_dataframe_confusion_matrix_entry(confusion_df: pd.DataFrame, filter_, result: ClassificationResult,
                                                  overall_number_positive_bases, overall_number_negative_bases):
        confusion_df.loc[str(filter_), ConfusionValueType.TP_BASES.name] = result.tp_bases
        confusion_df.loc[str(filter_), ConfusionValueType.FP_BASES.name] = result.fp_bases
        confusion_df.loc[str(filter_), ConfusionValueType.FN_BASES.name] = overall_number_positive_bases - result.tp_bases
        confusion_df.loc[str(filter_), ConfusionValueType.TN_BASES.name] = overall_number_negative_bases - result.fp_bases
