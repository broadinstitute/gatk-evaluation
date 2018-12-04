from itertools import product
from collections import OrderedDict

class CallsetFilter():

    def __init__(self, attribute_to_cutoff_map: map):
        """Constructor for a callset filter

        Args:
            attribute_to_cutoff_map: map from call attribute to its cutoff threshold
        """
        self.attribute_to_cutoff_map = attribute_to_cutoff_map
        self.filter_binary_lambda = self.__filter_binary_lambda
        self.representation = self.__representation()

    def __representation(self):
        representation = ""
        for attribute in self.attribute_to_cutoff_map.keys():
            representation += attribute
            representation += ">="
            representation += str(self.attribute_to_cutoff_map[attribute])
            representation += ";"
        return representation

    def __filter_binary_lambda(self, call_attributes: list):
        for attribute in self.attribute_to_cutoff_map.keys():
            if call_attributes[attribute] < self.attribute_to_cutoff_map[attribute]:
                return True
        return False

    def __str__(self):
        return self.representation


class BinningFilteringStrategy():

    def __init__(self, attributes: list, attributes_min_values: list, attributes_max_values: list, attributes_num_bins: list):
        """Constructor for a class that acts as a factory for creating callse filters
            Note that the minimum values for each filter is assumed to be 0

        Args:
            attributes: list of genotype fields to filter on
            attributes_min_values: lower boundaries of corresponding fields
            attributes_max_values: upper boundaries of corresponding fields
            attributes_num_bins: number of bins for corresponding fields
        """
        self.attributes = attributes
        self.attributes_min_values = attributes_min_values
        self.attributes_max_values = attributes_max_values
        self.attributes_num_bins = attributes_num_bins

        self.filter_list = self.__create_filters()
        self.current_idx = 0

    def __create_filters(self):
        cuttoffs_list = []
        for index, attribute in enumerate(self.attributes):
            step_size = (float(self.attributes_max_values[index]) - float(self.attributes_min_values[index])) / int(self.attributes_num_bins[index])
            cutoffs = [i*step_size for i in range(self.attributes_num_bins[index])]
            cuttoffs_list.append(cutoffs)

        list_of_filter_bins = list(product(*cuttoffs_list))
        
        filter_list = []
        for nd_bin in list_of_filter_bins:
            filter_map = OrderedDict()
            for index, attribute in enumerate(self.attributes):
                filter_map[attribute] = nd_bin[index]
            nd_filter = CallsetFilter(filter_map)
            filter_list.append(nd_filter)

        return filter_list

    def __iter__(self):
        self.current_idx = 0
        return self

    def __next__(self):
        if self.current_idx >= len(self.filter_list):
            raise StopIteration
        else:
            self.current_idx += 1
            return self.filter_list[self.current_idx - 1]

