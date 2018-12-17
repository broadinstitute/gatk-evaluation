from enum import Enum
import pandas as pd
from abc import ABC, abstractmethod
from collections import OrderedDict
import vcf
from interval import Interval
from intervaltree import IntervalTree
from intervaltree import Interval as TreeInterval
import pybedtools

from interval_collection import IntervalCollection, FeatureCollection
from filtering import CallsetFilter
import io_plt


class EventType(Enum):
    """Enumeration of possible alleles"""

    REF = 0
    DEL = 1
    DUP = 2
    NO_CALL = 3

    @classmethod
    def gcnv_genotype_to_event_type(cls, gcnv_call: int):
        return cls(gcnv_call)


class Call:
    """Stores an event type and call qualities for a single interval and single sample"""

    def __init__(self, interval: Interval, sample: str, event_type: EventType, call_attributes: map):
        self.interval = interval
        self.sample = sample
        self.event_type = event_type
        self.call_attributes = call_attributes

    def __eq__(self, other):
        return self.interval == other.interval and self.sample == other.sample \
               and self.event_type == other.event_type and self.call_attributes == other.call_attributes

    def __hash__(self) -> int:
        return super().__hash__()


class Callset(ABC):

    @abstractmethod
    def __init__(self, sample_to_calls_map: map):
        """Constructor for abstract callset representation

        Args:
            sample_to_calls_map: a map from samples to a FeatureCollection
        """

        assert len(sample_to_calls_map) > 0
        self.sample_names = sample_to_calls_map.keys()
        self.sample_to_calls_map = sample_to_calls_map
        self.filtered_calls = {}
        for sample in self.sample_names:
            self.filtered_calls[sample] = set()
        super().__init__()

    @classmethod
    @abstractmethod
    def read_in_callset(cls, **kwargs):
        pass

    def filter_callset(self, call_filter: CallsetFilter):
        """
        Filter callset given a binary lambda function that accepts call attributes map as argument

        Args:
            call_filter: filter to be applied to each of the calls
        """

        filter_lambda = call_filter.filter_binary_lambda
        for sample in self.sample_names:
            self.filtered_calls[sample] = \
                set(self.sample_to_calls_map.get(sample).get_all_features_matching_criteria(filter_lambda))

    def find_intersection_with_interval(self, interval: Interval, sample: str):
        """
        Given an interval find all overlapping calls in the callset and truncate them appropriately. In addition,
        fill in all missing gaps in callset with NO_CALL event types

        Args:
            interval: a given interval
            sample: sample from the callset

        Returns:
            A list of sorted, non-overlapping events that completely covers the given interval
        """

        calls = self.sample_to_calls_map.get(sample)
        intersecting_calls = calls.find_intersection(interval)
        filtered_intersecting_calls = [call for call in intersecting_calls if call
                                       not in self.filtered_calls.get(sample)]

        # check if there are any intersecting intervals
        if not filtered_intersecting_calls:
            return [(interval, EventType.NO_CALL)]
        else:
            result = IntervalTree([TreeInterval(call.interval.start, call.interval.end, call.event_type)
                                   for call in filtered_intersecting_calls])
            max_val = sorted(result)[-1].end
            min_val = sorted(result)[0].begin

            # fill potential gaps between calls that lie within interval with NO_CALL events
            result.addi(interval.start, interval.end, EventType.NO_CALL)
            result.split_overlaps()
            for call in filtered_intersecting_calls:
                result.remove_overlap(call.interval.center())
                result.addi(call.interval.start, call.interval.end, call.event_type)

            result.chop(interval.end, max(interval.end, max_val))
            result.chop(min(interval.start, min_val), interval.start)
            return [(Interval(interval.chrom, t.begin, t.end), t.data) for t in sorted(result)]


class TruthCallset(Callset):

    def __init__(self, sample_to_calls_map: map):
        super().__init__(sample_to_calls_map)

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "truth_file" in kwargs
        truth_file = kwargs["truth_file"]
        interval_file = kwargs["interval_file"]

        considered_interval_collection = IntervalCollection.read_interval_list(interval_file)

        truth_calls_pd = pd.read_csv(open(truth_file, 'r'), sep="\t", comment="#", header=None,
                                     names=["chrom", "start", "end", "name", "svtype", "samples"],
                                     dtype={"chrom": str, "start": int, "end": int, "name": str, "svtype": str,
                                            "samples": str})

        sample_to_calls_map = {}
        previous_interval_truth = None
        number_of_not_rescued_overlapping_events = 0
        number_of_overlapping_events_same_genotype = 0
        number_of_enveloped_events = 0
        overall_events = 0
        samples_set = set()
        for index, row in truth_calls_pd.iterrows():
            event_type = cls.__get_event_type_from_sv_type(row["svtype"])

            if event_type is None:
                continue

            interval = Interval(row["chrom"], int(row["start"]), int(row["end"]))
            if previous_interval_truth is not None and interval.chrom == previous_interval_truth.chrom \
                    and interval.start < previous_interval_truth.start:
                raise ValueError("Intervals Interval(%s) and Interval(%s) in truth callset are not in sorted order"
                                 % (previous_interval_truth, interval))
            if not considered_interval_collection.find_intersection(interval):
                continue

            sample_names = set(row["samples"].split(","))

            for sample_name in sample_names:
                samples_set.add(sample_name)
                call = Call(interval=interval, sample=sample_name, event_type=event_type, call_attributes=None)

                if sample_name in sample_to_calls_map:
                    overall_events += 1
                    if sample_to_calls_map[sample_name][0][-1].intersects_with(interval):
                        last_interval = sample_to_calls_map[sample_name][0][-1]
                        last_call = sample_to_calls_map[sample_name][1][-1]
                        if last_interval.end < interval.end and last_call.event_type == call.event_type:
                            # Merge overlapping events with the same call
                            sample_to_calls_map[sample_name][0][-1] = Interval(interval.chrom,
                                                                               last_interval.start,
                                                                               interval.end)
                            sample_to_calls_map[sample_name][1][-1].interval = \
                                sample_to_calls_map[sample_name][0][-1]
                            number_of_overlapping_events_same_genotype += 1
                        elif interval.end < last_interval.end:
                            # If one call is contained in another only keep the larger call
                            number_of_enveloped_events += 1
                            continue
                        else:
                            number_of_not_rescued_overlapping_events += 1
                            continue
                    sample_to_calls_map[sample_name][0].append(interval)
                    sample_to_calls_map[sample_name][1].append(call)
                else:
                    sample_to_calls_map[sample_name] = ([interval], [call])

            previous_interval_truth = interval

        for sample_name in samples_set:
            interval_to_call_map = OrderedDict()
            for index in range(len(sample_to_calls_map[sample_name][0])):
                interval_to_call_map[sample_to_calls_map[sample_name][0][index]] = \
                    sample_to_calls_map[sample_name][1][index]
            sample_to_calls_map[sample_name] = FeatureCollection(interval_to_call_map)

        io_plt.log("There are %d unique samples in truth set" % (len(samples_set)))
        io_plt.log("There are %d events for all samples in the truth call set" % overall_events)
        io_plt.log("There are %d intersecting events in truth set that were not rescued" %
                   number_of_not_rescued_overlapping_events)
        io_plt.log("There are %d overlapping events with the same genotype" %
                   number_of_overlapping_events_same_genotype)
        io_plt.log("There are %d enveloped events with different genotypes" % number_of_enveloped_events)
        return cls(sample_to_calls_map)

    @staticmethod
    def __get_event_type_from_sv_type(sv_type: str):
        """
        This method will return None if Structural Variation event type is not a Copy Number Variant
        """
        cnv_type_events = {"DUP", "DEL"}
        if sv_type not in cnv_type_events:
            return None
        else:
            return EventType[sv_type]


class GCNVCallset(Callset):

    def __init__(self, sample_to_calls_map: map):
        super().__init__(sample_to_calls_map)

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "gcnv_segment_vcfs" in kwargs
        gcnv_segment_vcfs = kwargs["gcnv_segment_vcfs"]
        sample_to_calls_map = {}
        for vcf_file in gcnv_segment_vcfs:

            vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            assert len(vcf_reader.samples) == 1
            sample_name = vcf_reader.samples[0]
            interval_to_call_map = OrderedDict()

            for record in vcf_reader:
                interval = Interval(record.CHROM, record.POS, record.INFO['END'])
                event_type = EventType.gcnv_genotype_to_event_type(int(record.genotype(sample_name)['GT']))
                attributes = {'QS': int(record.genotype(sample_name)['QS']),
                              'QA': int(record.genotype(sample_name)['QA']),
                              'NP': int(record.genotype(sample_name)['NP'])}
                call = Call(interval=interval, sample=sample_name, event_type=event_type, call_attributes=attributes)
                interval_to_call_map[interval] = call
            sample_to_calls_map[sample_name] = FeatureCollection(interval_to_call_map)
        return cls(sample_to_calls_map)
