from abc import ABC, abstractmethod
import pandas as pd
from collections import OrderedDict
from intervaltree import IntervalTree

from interval import Interval
import io_plt


class LocatableCollection(ABC):
    """
    Abstract class for a collection of genomic intervals and their corresponding data.
    """

    def find_intersection(self, interval: Interval)->list:
        """
        Find all entries in the collection that overlap with a given interval

        Args:
            interval: a given genomic interval

        Returns: list of corresponding data attributes for entries overlapping with a given interval

        """
        interval_tree_for_contig = self._get_interval_tree(interval.chrom)
        results = interval_tree_for_contig.search(interval.start, interval.end) if interval_tree_for_contig else []
        return [result.data for result in results]

    @abstractmethod
    def _get_interval_tree(self, contig: str):
        """

        Args:
            contig: chromosome for which to retrieve a corresponding interval tree

        Returns:
            an interval tree containing all intervals for a particular chromosome
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
            self.contig_to_intervals_map[interval.chrom][interval.start:interval.end] = interval

    @classmethod
    def read_interval_list(cls, interval_list_file):
        header = io_plt.read_comments_lines(interval_list_file)
        intervals_df = pd.read_csv(open(interval_list_file, 'r'), comment="@", delimiter="\t",
                                   usecols=[0, 1, 2], header=None, names=["CONTIG", "START", "END"],
                                   dtype={"CONTIG": str, "START": int, "END": int})
        intervals_series = intervals_df.apply(lambda x: Interval(str(x.CONTIG), int(x.START), int(x.END)), axis=1)
        interval_list = intervals_series.tolist()
        return cls(interval_list=interval_list, header=header)

    def find_intersection_with_interval_and_truncate(self, interval: Interval)->list:
        """
        Args:
            interval: input interval

        Returns: list of intervals in the collection that intersect with a given intervals truncated by the end points
        of a given interval

        """
        intersecting_intervals = self.find_intersection(interval)

        result = sorted(intersecting_intervals)
        if not result:
            return []
        chrom = interval.chrom
        min_val = result[0].start
        max_val = result[-1].end
        result[0] = Interval(chrom, max(min_val, interval.start), result[0].end)
        result[-1] = Interval(chrom, result[-1].start, min(max_val, interval.end))
        return result

    def _get_interval_tree(self, contig: str):
        return self.contig_to_intervals_map.get(contig, [])

    def get_intervals(self):
        return self.interval_list


class FeatureCollection(LocatableCollection):
    """
    Concrete class for a collection of genomic intervals and corresponding genomic features/'p[;oliu
    that allows search.
    """

    def __init__(self, interval_to_features_dict: OrderedDict):
        self.interval_list = list(interval_to_features_dict.keys())
        self.ordered_contigs = list(OrderedDict({t.chrom: None for t in self.interval_list}).keys())
        self.contig_to_intervals_map = {contig: IntervalTree() for contig in self.ordered_contigs}
        for interval in self.interval_list:
            self.contig_to_intervals_map[interval.chrom][interval.start:interval.end] = \
                interval_to_features_dict[interval]

    def _get_interval_tree(self, contig: str):
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

    def get_interval_tree(self, contig: str)->IntervalTree:
        return self.contig_to_intervals_map.get(contig)

    def get_feature_list(self)->list:
        feature_list = []
        for contig in self.ordered_contigs:
            [feature_list.append(feature.data) for feature in sorted(self.contig_to_intervals_map.get(contig))]
        return feature_list

    def get_all_features_matching_criteria(self, criteria)->list:
        features_matching_criteria = []
        for contig in self.ordered_contigs:
            for tree_interval in self.contig_to_intervals_map.get(contig):
                if criteria(tree_interval.data):
                    features_matching_criteria.append(tree_interval.data)
        return features_matching_criteria
