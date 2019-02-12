from typing import List
import itertools
from collections import OrderedDict
import numpy as np


class LowerBoundedFilterBin:

    def __init__(self, attribute_to_cutoff_map: OrderedDict):
        """Constructor for a callset filter bin with no upper bound

        Args:
            attribute_to_cutoff_map: map from call attributes to their respective lower thresholds
        """
        self.attribute_to_cutoff_map = attribute_to_cutoff_map
        self.attributes = list(attribute_to_cutoff_map.keys())
        self.representation = self.__name()

    def __name(self):
        representation = ""
        for attribute in self.attribute_to_cutoff_map.keys():
            representation += attribute
            representation += ">="
            representation += str(self.attribute_to_cutoff_map[attribute])
            representation += ";"
        return representation

    def __str__(self):
        return self.representation


class BoundedFilterBin:

    def __init__(self, attribute_to_cutoff_map: OrderedDict,
                 attribute_to_step_size_map: OrderedDict,
                 attributes_max_values_map: OrderedDict):
        """Constructor for a callset filter bin with lower(inclusive) and upper(exclusive) bounds

        Args:
            attribute_to_cutoff_map: map from call attributes to their respective lower thresholds
        """
        self.attribute_to_cutoff_map = attribute_to_cutoff_map
        self.attributes = list(attribute_to_cutoff_map.keys())
        self.attribute_to_step_size_map = attribute_to_step_size_map
        self.attributes_max_values_map = attributes_max_values_map
        self.representation = self.__name()

    def __name(self):
        representation = ""
        for attribute in self.attribute_to_cutoff_map.keys():
            upper_bound = self.attribute_to_cutoff_map[attribute] + self.attribute_to_step_size_map[attribute]
            representation += str(upper_bound) if upper_bound < self.attributes_max_values_map[attribute] else "Infinity"
            representation += ">"
            representation += attribute
            representation += ">="
            representation += str(self.attribute_to_cutoff_map[attribute])
            representation += ";"
        return representation

    def __str__(self):
        return self.representation


class FilteringBinCollection:

    def __init__(self, attributes: List[str],  attributes_max_values: List[float], attributes_num_bins: List[int]):
        """Constructor for collection of callset filtering bins
            Note that the minimum values for each filter is assumed to be 0

        Args:
            attributes: list of genotype fields to filter on
            attributes_max_values: upper boundaries of corresponding fields
            attributes_num_bins: number of bins for corresponding fields
        """
        assert len(attributes) == len(attributes_max_values) and len(attributes) == len(attributes_num_bins)
        assert len(attributes) > 0
        for num_bins in attributes_num_bins:
            assert num_bins > 0
        self.attributes = attributes
        self.attributes_min_values = [0.0] * len(attributes)
        self.attributes_max_values_map = OrderedDict(zip(attributes, attributes_max_values))
        self.attributes_num_bins = attributes_num_bins
        self.attributes_num_bins_map = OrderedDict(zip(attributes, attributes_num_bins))
        self.bounded_filter_list, self.bounded_filter_nd_array, self.lower_bounded_filter_list, \
            self.lower_bounded_nd_array, self.attribute_to_step_size_map = self.__create_filters()

        self.current_idx = 0

    def __create_filters(self):
        cuttoffs_list = []
        attribute_to_step_size_map = OrderedDict()
        for index, attribute in enumerate(self.attributes):
            step_size = (float(self.attributes_max_values_map[attribute]) - float(self.attributes_min_values[index])) \
                        / int(self.attributes_num_bins[index])
            attribute_to_step_size_map[attribute] = step_size
            cutoffs = [i * step_size for i in range(self.attributes_num_bins[index])]
            cuttoffs_list.append(cutoffs)

        list_of_filter_bins = list(itertools.product(*cuttoffs_list))
        bounded_filter_nd_array = np.empty(shape=self.attributes_num_bins, dtype=BoundedFilterBin)
        lower_bounded_nd_array = np.empty(shape=self.attributes_num_bins, dtype=LowerBoundedFilterBin)
        bounded_filter_list = []
        lower_bounded_filter_list = []
        for counter, nd_bin in enumerate(list_of_filter_bins):
            attribute_to_cutoff_map = OrderedDict()
            for index, attribute in enumerate(self.attributes):
                attribute_to_cutoff_map[attribute] = nd_bin[index]

            bounded_filter = BoundedFilterBin(attribute_to_cutoff_map, attribute_to_step_size_map, self.attributes_max_values_map)
            bounded_filter_list.append(bounded_filter)
            bounded_filter_nd_array[np.unravel_index(counter, dims=self.attributes_num_bins)] = bounded_filter
            lower_bounded_filter = LowerBoundedFilterBin(attribute_to_cutoff_map)
            lower_bounded_filter_list.append(lower_bounded_filter)
            lower_bounded_nd_array[np.unravel_index(counter, dims=self.attributes_num_bins)] = lower_bounded_filter

        return bounded_filter_list, bounded_filter_nd_array, lower_bounded_filter_list, lower_bounded_nd_array, attribute_to_step_size_map

    def get_filter_bin_by_attribute_values(self, attribute_to_value_map: dict):
        """

        Args:
            attribute_to_value_map: map from attributes to the corresponding values

        Returns:
            A bounded filter bin whose boundaries contain the values given for corresponding attributes

        """
        if attribute_to_value_map is None:
            return 0
        indices = [0] * len(self.attributes_num_bins)
        for idx, attribute in enumerate(self.attributes):
            value = attribute_to_value_map.get(attribute)
            indices[idx] = min(int(value / self.attribute_to_step_size_map[attribute]), self.attributes_num_bins_map[attribute] - 1)

        return self.bounded_filter_nd_array[tuple(indices)]

    def get_single_attribute_lower_bounded_filters(self, attribute: str)->List[LowerBoundedFilterBin]:
        """

        Args:
            attribute: genotype field that was used to construct the filtering strategy

        Returns:
            a list of lower bounded filters that only filter based on that genotype field
            (e.g. [QS>0;QA>0, QS>10;QA>0, QS>20;QA>0] where QS is the input attribute)

        """
        assert attribute in self.attributes
        sl = ()
        for index, attr in enumerate(self.attributes):
            if attr == attribute:
                sl += (slice(None),)
            else:
                sl += (0,)
        return self.lower_bounded_nd_array[sl]

    def get_lower_bounded_filter_sum(self, lower_bounded_filter: LowerBoundedFilterBin)->List[BoundedFilterBin]:
        """
        Get list of bounded filters that sum up to the given lower bounded one

        Args:
            lower_bounded_filter: an instance of LowerBoundedFilterBin

        Returns:
            a list of BoundedFilterBin objects that sum up to the given bounded filter bin
        """
        index = np.where(self.lower_bounded_nd_array == lower_bounded_filter)
        assert len(index[0]) == 1
        sl = ()
        for idx in index:
            sl += (slice(idx[0], None, None),)
        return self.bounded_filter_nd_array[sl].flatten().tolist()
