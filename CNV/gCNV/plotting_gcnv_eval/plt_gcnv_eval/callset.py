from enum import Enum
import pandas as pd
from abc import ABC, abstractmethod
from collections import OrderedDict
import vcf 

from interval import Interval
from interval_collection import IntervalCollection
from filtering import CallsetFilter


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


class Callset(ABC):
    """Abstract callset class"""

    @abstractmethod
    def __init__(self, sample_to_calls_map: map):
        """Constructor for abstract callset representation

        Args:
            sample_names: list of sample names in the collection
            sample_to_calls_map: a map from samples to a pair (IntervalCollection, List[Call])
        """
        assert len(sample_to_calls_map) > 0
        self.sample_names = sample_to_calls_map.keys()
        self.sample_to_calls_map = sample_to_calls_map
        self.filtered_indices = {}
        for sample in self.sample_names:
            self.filtered_indices[sample] = set()
        super().__init__()

    @classmethod
    @abstractmethod
    def read_in_callset(cls, **kwargs):
        pass

    def filter_callset(call_filter: CallsetFilter):
        """Filter callset given a binary lambda function that accepts Call.call_info as argument"""
        for sample in sample_names:
            intervals, calls = self.sample_to_calls_map[sample]
            self.indices_filtered[sample] = set([i for i in range(len(intervals)) if call_filter(calls[i].call_attributes)])

    def find_intersection_with_interval(self, interval: Interval, sample: str):
        #return map from Interval to EventType
        intervals, calls = self.sample_to_calls_map[sample]
        unfiltered_indices = intervals.find_intersecting_interval_indices(interval)

        #Sanity check
        IntervalCollection.assert_interval_list_sorted([intervals.interval_list[index] for index in unfiltered_indices])

        filtered_indices = [index for index in unfiltered_indices if index not in self.filtered_indices[sample]]
        chromosome = interval.chrom

        #check if there are any intersecting intervals
        if (not filtered_indices):
            return [(interval, EventType.NO_CALL)]

        current_position = min(intervals.interval_list[filtered_indices[0]].start, interval.start)
        intersection = []
        # TODO refactor this
        for index in filtered_indices:
            call_interval_end = intervals.interval_list[index].end
            call_interval_start = intervals.interval_list[index].start
            if (call_interval_start <= interval.start):
                if (call_interval_end <= interval.end):
                    intersection.append((Interval(chromosome, interval.start, call_interval_end), calls[index].event_type))
                else:
                    intersection.append((Interval(chromosome, interval.start, interval.end), calls[index].event_type))
            else:
                #insert an event indicating that no call has been made on that part of the interval
                intersection.append((Interval(chromosome, current_position, call_interval_start), EventType.NO_CALL))
                if (call_interval_end <= interval.end):
                    intersection.append((Interval(chromosome, call_interval_start, call_interval_end), calls[index].event_type))
                else:
                    intersection.append((Interval(chromosome, call_interval_start, interval.end), calls[index].event_type))
            current_position = min(call_interval_end, interval.end)

        if (current_position < interval.end):
            intersection.append((Interval(chromosome, current_position, interval.end), EventType.NO_CALL))

        return intersection


class TruthCallset(Callset):
    """Callset containing truth calls from SV pipeline"""

    def __init__(self, sample_to_calls_map: map):
        super().__init__(sample_to_calls_map)

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "truth_file" in kwargs
        truth_file = kwargs["truth_file"]
        interval_file = kwargs["interval_file"]

        considered_interval_collection = IntervalCollection.read_interval_list(interval_file)

        truth_calls_pd = pd.read_csv(open(truth_file, 'r'), sep="\t", comment="#", header=None,
                                     names=["chrom", "start", "end", "name" , "svtype", "samples"],
                                     dtype={"chrom": str, "start": int, "end": int, "name": str, "svtype": str, "samples": str})

        sample_to_calls_map = {}
        previous_interval_truth = None
        number_of_not_rescued_overlapping_events = 0
        number_of_overlapping_events_same_genotype = 0
        number_of_enveloped_events = 0
        overall_events = 0 
        samples_set = set()
        for index, row in truth_calls_pd.iterrows():
            event_type = cls.__get_event_type_from_svtype(row["svtype"])

            if (event_type == None):
                continue

            interval = Interval(row["chrom"], int(row["start"]), int(row["end"]))
            if (previous_interval_truth is not None and interval.chrom == previous_interval_truth.chrom and interval.start < previous_interval_truth.start):
                raise ValueError("Intervals Interval(%s) and Interval(%s) in truth callset are not in sorted order" % (previous_interval_truth, interval))

            sample_names = set(row["samples"].split(","))

            for sample_name in sample_names:
                call = Call(interval=interval, sample=sample_name, event_type=event_type, call_attributes=None)
                if sample_name in sample_to_calls_map:
                    samples_set.add(sample_name)
                    intersect_indices = considered_interval_collection.find_intersecting_interval_indices(interval)
                    if (len(intersect_indices) > 0):
                        overall_events += 1
                        if (len(sample_to_calls_map[sample_name][0].find_intersecting_interval_indices(interval)) > 0):
                            last_interval = sample_to_calls_map[sample_name][0].interval_list[-1]
                            last_call = sample_to_calls_map[sample_name][1][-1]
                            if (last_interval.end < interval.end and last_call.event_type == call.event_type):
                                # Merge overlapping events with the same call
                                sample_to_calls_map[sample_name][0].interval_list[-1] = Interval(interval.chrom, last_interval.start, interval.end)
                                sample_to_calls_map[sample_name][1][-1].interval = sample_to_calls_map[sample_name][0].interval_list[-1]
                                number_of_overlapping_events_same_genotype += 1
                            elif (interval.end < last_interval.end):
                                #If one call is contained in another only keep the larger call
                                number_of_enveloped_events += 1
                                continue
                            else:
                                number_of_not_rescued_overlapping_events += 1
                                continue

                    sample_to_calls_map[sample_name][0].interval_list.append(interval)
                    sample_to_calls_map[sample_name][1].append(call)
                else:
                    sample_to_calls_map[sample_name] = (IntervalCollection([interval], None), [call])

            previous_interval_truth = interval

        print("There are %d unique samples in truth set" % (len(samples_set)))
        print("There are %d events for all samples in the truth call set" % (overall_events))
        print("There are %d intersecting events in truth set that were not rescued" % (number_of_not_rescued_overlapping_events))
        print("There are %d overlapping events with the same genotype" % (number_of_overlapping_events_same_genotype))
        print("There are %d enveloped events with different genotypes" % (number_of_enveloped_events))
        return cls(sample_to_calls_map)

    @staticmethod
    def __get_event_type_from_svtype(sv_type: str):
        """
        This method will return None if Structural Variation event type is not a Copy Number Variant
        """
        CNV_TYPE_EVENTS={"DUP", "DEL"}
        if(not sv_type in CNV_TYPE_EVENTS):
            return None
        else:
            return EventType[sv_type]


class GCNVCallset(Callset):
    """Callset containing segment gCNV calls"""
    
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
            sample_to_calls_map[sample_name] = (IntervalCollection([], None), [])

            for record in vcf_reader:
                interval = Interval(record.CHROM, record.POS, record.INFO['END'])
                event_type = EventType.gcnv_genotype_to_event_type(int(record.genotype(sample_name)['GT']))
                attributes = {'QS': int(record.genotype(sample_name)['QS']),
                              'QA': int(record.genotype(sample_name)['QA']),
                              'NP': int(record.genotype(sample_name)['NP'])}
                call = Call(interval=interval, sample=sample_name, event_type=event_type, call_attributes=attributes)
                sample_to_calls_map[sample_name][0].interval_list.append(interval)
                sample_to_calls_map[sample_name][1].append(call)
        return cls(sample_to_calls_map)

