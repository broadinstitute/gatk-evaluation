from abc import ABC, abstractmethod
from collections import OrderedDict

import pandas as pd
import vcf
from intervaltree import Interval as TreeInterval
from intervaltree import IntervalTree

import io_plt
from call import Call, EventType
from interval import Interval
from interval_collection import IntervalCollection, FeatureCollection
from reference_dictionary import ReferenceDictionary


class Callset(ABC):

    @abstractmethod
    def __init__(self, sample_to_calls_map: dict, ref_dict: ReferenceDictionary):
        """Constructor for abstract callset representation

        Args:
            sample_to_calls_map: a map from samples to a FeatureCollection
            ref_dict: reference dictionary
        """

        assert len(sample_to_calls_map) > 0
        self.ref_dict = ref_dict
        self.sample_names = set(sample_to_calls_map.keys())
        self.sample_to_calls_map = sample_to_calls_map
        self.__preprocess()

        # TODO add a check to make sure the callset is not malformed, i.e. the calls don't intersect and
        # TODO the intervals in the featurecollections equal to the intervals stored in their corresponding calls
        # TODO Also make sure that code doesn't break if one of the contigs is not in the callset
        super().__init__()

    @classmethod
    @abstractmethod
    def read_in_callset(cls, **kwargs):
        pass

    def __preprocess(self):
        """
        Preprocess the callset by filling the regions with no calls with EventType.NO_CALL events, thereby assigning
        an event to every single base.

        """

        for sample in self.sample_names:
            interval_to_call_map = OrderedDict()
            for contig in self.ref_dict.contigs:
                contig_interval = self.ref_dict.get_contig_interval_for_chrom_name(contig)
                events_on_contig = self.sample_to_calls_map.get(sample).get_interval_tree(contig)
                if events_on_contig is None:
                    continue

                result = events_on_contig.copy()
                # TODO make code aware of 1-based representation
                # TODO i.e. right now some events overlap by a single base
                # This hacky code fills potential gaps between calls that lie within interval with NO_CALL events
                result.addi(contig_interval.start, contig_interval.end, Call(interval=contig_interval,
                                                                             sample=sample,
                                                                             event_type=EventType.NO_CALL,
                                                                             call_attributes={"QS": 0, "QA": 0}))
                result.split_overlaps()
                for interval in events_on_contig.items():
                    result.remove_overlap(interval.begin, interval.end)
                for interval in events_on_contig.items():
                    result.addi(interval.begin, interval.end, Call.deep_copy(interval.data))

                for t in sorted(result):
                    if t.end - t.begin == 1 and t.data.event_type == EventType.NO_CALL:
                        # intervaltree.split_overlaps will create single base regions which we want to discard
                        continue
                    call = Call.deep_copy(t.data)
                    if t.data.event_type == EventType.NO_CALL:
                        call.interval = Interval(contig, t.begin, t.end)
                    interval_to_call_map[Interval(contig, t.begin, t.end)] = call
            self.sample_to_calls_map[sample] = FeatureCollection(interval_to_call_map)

    def find_intersection_with_interval(self, interval: Interval, sample: str):
        """
        Given an interval find all overlapping calls in the callset and truncate them appropriately.
        Note: we assume that the calls in the callset do not overlap for a single sample.

        Args:
            interval: a given interval
            sample: sample from the callset

        Returns:
            A list of sorted, non-overlapping events that completely cover a given interval
        """

        assert sample in self.sample_names, "Sample %s is not in the callset" % sample

        calls = self.sample_to_calls_map.get(sample)
        intersecting_calls = calls.find_intersection(interval)

        if not intersecting_calls:
            return [(interval, EventType.NO_CALL)]
        else:
            result = IntervalTree([TreeInterval(call.interval.start, call.interval.end, call.event_type)
                                   for call in intersecting_calls])
            max_val = sorted(result)[-1].end
            min_val = sorted(result)[0].begin
            result.chop(interval.end, max(interval.end, max_val))
            result.chop(min(interval.start, min_val), interval.start)
            return [(Interval(interval.chrom, t.begin, t.end), t.data) for t in sorted(result)]

    def to_string_sample(self, sample):
        print("#sample=%s" % sample)
        callset_feature_collection = self.sample_to_calls_map.get(sample)
        for contig in callset_feature_collection.ordered_contigs:
            for tree_interval in sorted(callset_feature_collection.get_interval_tree(contig)):
                print(str(Interval(contig, tree_interval.begin, tree_interval.end)) + '\t' + str(tree_interval.data))
                print(str(Interval(contig, tree_interval.begin, tree_interval.end)) + '\t' + str(tree_interval.data))


class TruthCallset(Callset):

    def __init__(self, sample_to_calls_map: map, ref_dict: ReferenceDictionary):
        super().__init__(sample_to_calls_map, ref_dict)

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "truth_file" in kwargs
        truth_file = kwargs["truth_file"]
        interval_file = kwargs["interval_file"]
        ref_dict = kwargs["reference_dictionary"]

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
                    if sample_to_calls_map.get(sample_name)[-1].interval.intersects_with(interval):
                        last_interval = sample_to_calls_map.get(sample_name)[-1].interval
                        last_call = sample_to_calls_map.get(sample_name)[-1]
                        if last_interval.end <= interval.end and last_call.event_type == call.event_type:
                            # Merge overlapping events with the same call
                            new_interval = Interval(interval.chrom, last_interval.start, interval.end)
                            sample_to_calls_map.get(sample_name)[-1].interval = new_interval
                            number_of_overlapping_events_same_genotype += 1
                        elif interval.end < last_interval.end:
                            # If one call is contained in another only keep the larger call
                            number_of_enveloped_events += 1
                        else:
                            number_of_not_rescued_overlapping_events += 1
                        continue
                    sample_to_calls_map.get(sample_name).append(call)
                else:
                    sample_to_calls_map[sample_name] = [call]

            previous_interval_truth = interval

        for sample_name in samples_set:
            interval_to_call_map = OrderedDict()
            for index in range(len(sample_to_calls_map.get(sample_name))):
                interval_to_call_map[sample_to_calls_map.get(sample_name)[index].interval] = \
                    sample_to_calls_map.get(sample_name)[index]
            sample_to_calls_map[sample_name] = FeatureCollection(interval_to_call_map)

        io_plt.log("There are %d unique samples in truth set" % (len(samples_set)))
        io_plt.log("There are %d events for all samples in the truth call set" % overall_events)
        io_plt.log("There are %d intersecting events in truth set that were not rescued" %
                   number_of_not_rescued_overlapping_events)
        io_plt.log("There are %d overlapping events with the same genotype" %
                   number_of_overlapping_events_same_genotype)
        io_plt.log("There are %d enveloped events with different genotypes" % number_of_enveloped_events)
        return cls(sample_to_calls_map, ref_dict)

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

    def __init__(self, sample_to_calls_map: map, ref_dict: ReferenceDictionary):
        super().__init__(sample_to_calls_map, ref_dict)

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "gcnv_segment_vcfs" in kwargs
        gcnv_segment_vcfs = kwargs["gcnv_segment_vcfs"]
        ref_dict = kwargs["reference_dictionary"]
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
        return cls(sample_to_calls_map, ref_dict)
