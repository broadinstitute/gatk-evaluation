from typing import List
import itertools
from collections import OrderedDict
import numpy as np
import math


class CallsetFilter:

    def __init__(self, attribute_to_cutoff_map: map):
        """Constructor for a callset filter

        Args:
            attribute_to_cutoff_map: map from call attribute to its cutoff threshold
        """
        self.attribute_to_cutoff_map = attribute_to_cutoff_map
        self.attributes = list(attribute_to_cutoff_map.keys())
        self.filter_binary_lambda = self.__filter_binary_lambda
        self.representation = self.__name()

    def __name(self):
        representation = ""
        for attribute in self.attribute_to_cutoff_map.keys():
            representation += attribute
            representation += ">="
            representation += str(self.attribute_to_cutoff_map[attribute])
            representation += ";"
        return representation

    def __filter_binary_lambda(self, call):
        # TODO move Call class to a different file so that this ugly import is not needed here
        from callset import EventType
        if call.event_type == EventType.NO_CALL:
            return False
        for attribute in self.attribute_to_cutoff_map.keys():
            if call.call_attributes.get(attribute) < self.attribute_to_cutoff_map.get(attribute):
                return True
        return False

    def __str__(self):
        return self.representation


class BinningFilteringStrategy:

    def __init__(self, attributes: List[str],  attributes_max_values: List[float], attributes_num_bins: List[int]):
        """Constructor for a class that acts as a factory for creating callset filters
            Note that the minimum values for each filter is assumed to be 0

        Args:
            attributes: list of genotype fields to filter on
            attributes_max_values: upper boundaries of corresponding fields
            attributes_num_bins: number of bins for corresponding fields
        """
        assert len(attributes) == len(attributes_max_values) and len(attributes) == len(attributes_num_bins)
        self.attributes = attributes
        self.attributes_min_values = [0.0] * len(attributes)
        self.attributes_max_values = attributes_max_values
        self.attributes_num_bins = attributes_num_bins
        self.filter_list, self.filter_nd_array, self._single_attribute_filter_indices_map, self.size_steps =\
            self.__create_filters()

        self.current_idx = 0

    def __create_filters(self):
        cuttoffs_list = []
        size_steps = []
        for index, attribute in enumerate(self.attributes):
            step_size = (float(self.attributes_max_values[index]) - float(self.attributes_min_values[index])) \
                        / int(self.attributes_num_bins[index])
            size_steps.append(step_size)
            cutoffs = [i*step_size for i in range(self.attributes_num_bins[index])]
            cuttoffs_list.append(cutoffs)

        list_of_filter_bins = list(itertools.product(*cuttoffs_list))
        filter_nd_array = np.empty(shape=self.attributes_num_bins, dtype=CallsetFilter)
        filter_list = []
        single_attribute_filter_indices_map = {}
        for counter, nd_bin in enumerate(list_of_filter_bins):
            filter_map = OrderedDict()
            for index, attribute in enumerate(self.attributes):
                if nd_bin[index] == 0.0 or len(self.attributes) == 1:
                    single_attribute_filter_indices_map.setdefault(attribute, []).append(counter)
                filter_map[attribute] = nd_bin[index]

            nd_filter = CallsetFilter(filter_map)
            filter_list.append(nd_filter)
            filter_nd_array[np.unravel_index(counter, dims=self.attributes_num_bins)] = nd_filter

        return filter_list, filter_nd_array, single_attribute_filter_indices_map, size_steps

    def get_filter_index_by_values(self, attribute_to_value_map: dict):
        """

        Args:
            attribute_to_value_map:

        Returns:

        """
        if attribute_to_value_map is None:
            return 0
        indices = [0] * len(self.attributes_num_bins)
        for idx, attribute in enumerate(self.attributes):
            value = attribute_to_value_map.get(attribute)

            indices[idx] = int(math.ceil(min(value, self.attributes_max_values[idx]) / self.size_steps[idx])) - 1
        return np.ravel_multi_index(indices, dims=self.attributes_num_bins)

    def get_single_attribute_filters(self, attribute: str)->list:
        """

        Args:
            attribute: genotype field that was used to construct the filtering strategy

        Returns:
            a list of indices for filters that only filter based on that genotype field
            (e.g. [QS>0;QA>0, QS>10;QA>0, QS>20;QA>0] where QS is the input attribute)

        """
        assert attribute in self.attributes
        return [self.filter_list[i] for i in self._single_attribute_filter_indices_map.get(attribute)]

    def __iter__(self):
        self.current_idx = 0
        return self

    def __next__(self):
        if self.current_idx >= len(self.filter_list):
            raise StopIteration
        else:
            self.current_idx += 1
            return self.filter_list[self.current_idx - 1]
