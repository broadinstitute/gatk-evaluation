from typing import List
from itertools import product
from collections import OrderedDict


class CallsetFilter:

    def __init__(self, attribute_to_cutoff_map: map):
        """Constructor for a callset filter

        Args:
            attribute_to_cutoff_map: map from call attribute to its cutoff threshold
        """
        self.attribute_to_cutoff_map = attribute_to_cutoff_map
        self.attributes = attribute_to_cutoff_map.keys()
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
        self.attributes = attributes
        self.attributes_min_values = [0.0] * len(attributes)
        self.attributes_max_values = attributes_max_values
        self.attributes_num_bins = attributes_num_bins
        self.filter_list, self.single_attribute_filter_indices_map = self.__create_filters()

        self.current_idx = 0

    def __create_filters(self):
        cuttoffs_list = []
        for index, attribute in enumerate(self.attributes):
            step_size = (float(self.attributes_max_values[index]) - float(self.attributes_min_values[index])) / int(self.attributes_num_bins[index])
            cutoffs = [i*step_size for i in range(self.attributes_num_bins[index])]
            cuttoffs_list.append(cutoffs)

        list_of_filter_bins = list(product(*cuttoffs_list))
        
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

        return filter_list, single_attribute_filter_indices_map

    def get_single_attribute_filters(self, attribute: str):
        return [self.filter_list[i] for i in self.single_attribute_filter_indices_map.get(attribute)]

    def __iter__(self):
        self.current_idx = 0
        return self

    def __next__(self):
        if self.current_idx >= len(self.filter_list):
            raise StopIteration
        else:
            self.current_idx += 1
            return self.filter_list[self.current_idx - 1]
