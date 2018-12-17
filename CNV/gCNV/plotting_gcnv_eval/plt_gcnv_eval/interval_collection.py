from interval import Interval
import io_plt
import pandas as pd
from collections import OrderedDict
from intervaltree import IntervalTree
from abc import ABC, abstractmethod


class LocatableCollection(ABC):
    """
    Abstract class for a collection of genomic intervals that allows search and corresponding data.
    """

    def find_intersection(self, interval: Interval) -> list:
        """

        Args:
            interval:

        Returns:

        """
        interval_tree_for_contig = self.get_interval_tree(interval.chrom)
        results = interval_tree_for_contig.search(interval.start, interval.end) if interval_tree_for_contig else []
        return [result.data for result in results]

    @abstractmethod
    def get_interval_tree(self, contig: str):
        """

        Args:
            contig:

        Returns:

        """
        pass

    @abstractmethod
    def get_intervals(self):
        """

        Returns: All intervals contained in this collection

        """
        pass


class IntervalCollection(LocatableCollection):
    """
    Concrete class for a collection of genomic intervals only
    """

    def __init__(self, interval_list: list, header: str = None):
        self.interval_list = interval_list
        self.header = header
        self.ordered_contigs = list(OrderedDict({t.chrom: None for t in self.interval_list}).keys())
        self.contig_to_intervals_map = {contig: IntervalTree() for contig in self.ordered_contigs}
        for interval in interval_list:
            self.contig_to_intervals_map[interval.chrom][interval.start:interval.end] = None

    @classmethod
    def read_interval_list(cls, interval_list_file):
        header = io_plt.read_comments_lines(interval_list_file)
        intervals_df = pd.read_csv(open(interval_list_file, 'r'), comment="@", delimiter="\t",
                                   usecols=[0, 1, 2], header=None, names=["CONTIG", "START", "END"],
                                   dtype={"CONTIG": str, "START": int, "END": int})
        intervals_series = intervals_df.apply(lambda x: Interval(str(x.CONTIG), int(x.START), int(x.END)), axis=1)
        interval_list = intervals_series.tolist()
        return cls(interval_list=interval_list, header=header)

    def get_interval_tree(self, contig: str):
        return self.contig_to_intervals_map.get(contig, [])

    def get_intervals(self):
        return self.interval_list


class FeatureCollection(LocatableCollection):
    """
    Concrete class for a collection of genomic intervals and corresponding genomic features
    that allows search.
    """

    def __init__(self, interval_to_features_dict: OrderedDict):
        """

        Args:
            interval_to_features_dict:
        """
        self.interval_list = interval_to_features_dict.keys()
        self.ordered_contigs = list(OrderedDict({t.chrom: None for t in self.interval_list}).keys())
        self.contig_to_intervals_map = {contig: IntervalTree() for contig in self.ordered_contigs}
        for interval in self.interval_list:
            self.contig_to_intervals_map[interval.chrom][interval.start:interval.end] = \
                interval_to_features_dict[interval]

    def get_interval_tree(self, contig: str):
        return self.contig_to_intervals_map.get(contig, [])

    def get_intervals(self):
        return self.interval_list

    def add_feature(self, interval: Interval, data):
        self.contig_to_intervals_map[interval.chrom][interval.start, interval.end] = data

    def find_feature_for_interval(self, interval: Interval):
        tree = self.contig_to_intervals_map[interval.chrom]
        tree_intervals = tree.search(interval.start, interval.end)
        assert len(tree_intervals) == 1 and tree_intervals[0].begin == interval.start \
            and tree_intervals[0].end == interval.end, "There should be exactly one entry in the collection"

        return tree_intervals[0].data

    def get_all_features_matching_criteria(self, criteria) -> list:
        features_matching_criteria = []
        for contig in self.ordered_contigs:
            for tree_interval in self.contig_to_intervals_map.get(contig):
                if criteria(tree_interval.data):
                    features_matching_criteria.append(tree_interval.data)
        return features_matching_criteria
