import numpy as np
import pickle
from gcnvkernel import Interval
from typing import Optional, Set, Generator, Dict, Tuple, Callable, Union
from intervaltree_bio import GenomeIntervalTree, IntervalTree
from collections import namedtuple
import logging

_logger = logging.getLogger(__name__)


class GenericCopyNumberVariant(Interval):
    """This class represent a generic CNV locus.

    Note:
        Equality testing and hashing is done similarly to 'Interval', i.e. calls and other
        annotations are ignored.
    """
    def __init__(self,
                 contig: str, start: int, end: int,
                 var_copy_number: int,
                 quality: float,
                 genotype: Optional[str] = None,
                 variant_frequency: Optional[float] = None,
                 num_intervals: Optional[int] = 1,
                 variant_class: Optional[str] = None):
        """Initializer.

        Args:
            contig: contig
            start: start
            end: end
            var_copy_number: var copy number
            quality: quality score
            genotype: (optional) genotype: 'ref', 'dup', 'del'
            variant_frequency: (optional) variant frequency in the cohort
            num_intervals: (optional) number of intervals that are merged to create the variant
            variant_class: (optional) variant class in a population: 'ref', 'dup', 'del', 'mixed'
        """
        super().__init__(contig, start, end)

#        assert var_copy_number >= 0
        assert quality >= 0.
        assert genotype is None or genotype in ['ref', 'dup', 'del', 'cnv']
        assert variant_frequency is None or 0. <= variant_frequency <= 1.
        assert variant_class is None or variant_class in ['ref', 'dup', 'del', 'mixed']
        assert num_intervals is None or num_intervals >= 1

        self.var_copy_number = var_copy_number
        self.genotype = genotype
        self.quality = quality
        self.variant_frequency = variant_frequency
        self.num_intervals = num_intervals
        self.variant_class = variant_class

    def set_genotype(self, ref_copy_number):
        if self.var_copy_number == ref_copy_number:
            self.genotype = 'ref'
        elif self.var_copy_number < ref_copy_number:
            self.genotype = 'del'
        else:
            self.genotype = 'dup'

    @property
    def is_dup(self):
        assert self.genotype is not None, "Genotype for {0} is not available. " \
                                          "Cannot determine whether the variant is DUP.".format(self)
        return self.genotype == 'dup'

    @property
    def is_del(self):
        assert self.genotype is not None, "Genotype for {0} is not available. " \
                                          "Cannot determine whether the variant is DEL.".format(self)
        return self.genotype == 'del'

    @property
    def is_var(self):
        assert self.genotype is not None, "Genotype for {0} is not available. " \
                                          "Cannot determine whether the variant is VAR (non-REF).".format(self)
        return self.genotype != 'ref'

    def variant_frequency_lte_value(self, max_variant_frequency: float):
        assert self.variant_frequency is not None, "Variant {0} has no variant frequency annotation".format(self)
        return self.variant_frequency <= max_variant_frequency

    def variant_frequency_gte_value(self, min_variant_frequency: float):
        assert self.variant_frequency is not None, "Variant {0} has no variant frequency annotation".format(self)
        return self.variant_frequency >= min_variant_frequency

    def quality_lte_value(self, max_quality: float):
        return self.quality <= max_quality

    def quality_gte_value(self, min_quality: float):
        return self.quality >= min_quality

    def length_lte_value(self, max_length: int):
        return (self.end - self.start) <= max_length

    def length_gte_value(self, min_length: int):
        return (self.end - self.start) >= min_length

    def get_padded(self, padding: int, keep_annotations=False):
        return GenericCopyNumberVariant(self.contig, self.start - padding, self.end + padding,
                                        self.var_copy_number, self.quality,
                                        genotype=self.genotype,
                                        variant_frequency=self.variant_frequency,
                                        variant_class=self.variant_class,
                                        num_intervals=self.num_intervals)

    @staticmethod
    def get_variant_filter(min_quality: Optional[float] = None,
                           max_quality: Optional[float] = None,
                           min_variant_frequency: Optional[float] = None,
                           max_variant_frequency: Optional[float] = None,
                           min_length: Optional[int] = None,
                           max_length: Optional[int] = None,
                           included_variant_classes: Optional[Set[str]] = None,
                           included_contigs: Optional[Set[str]] = None) -> Callable[['GenericCopyNumberVariant'], bool]:

        def min_quality_pass(variant: GenericCopyNumberVariant):
            return min_quality is None or variant.quality_gte_value(min_quality)

        def max_quality_pass(variant: GenericCopyNumberVariant):
            return max_quality is None or variant.quality_lte_value(max_quality)

        def min_variant_frequency_pass(variant: GenericCopyNumberVariant):
            return min_variant_frequency is None or variant.variant_frequency_gte_value(min_variant_frequency)

        def max_variant_frequency_pass(variant: GenericCopyNumberVariant):
            return max_variant_frequency is None or variant.variant_frequency_lte_value(max_variant_frequency)

        def min_length_pass(variant: GenericCopyNumberVariant):
            return min_length is None or variant.length_gte_value(min_length)

        def max_length_pass(variant: GenericCopyNumberVariant):
            return max_length is None or variant.length_lte_value(max_length)

        def included_variant_classes_pass(variant: GenericCopyNumberVariant):
            return included_variant_classes is None or variant.variant_class in included_variant_classes

        def included_contigs_pass(variant: GenericCopyNumberVariant):
            return included_contigs is None or variant.contig in included_contigs

        def all_pass(variant: GenericCopyNumberVariant):
            return (min_quality_pass(variant) and max_quality_pass(variant) and
                    min_variant_frequency_pass(variant) and max_variant_frequency_pass(variant) and
                    min_length_pass(variant) and max_length_pass(variant) and
                    included_variant_classes_pass(variant) and
                    included_contigs_pass(variant))

        return all_pass

    @staticmethod
    def get_overlap_fraction(first: Interval, second: Interval, ref: str) -> float:
        """Calculates the overlap fraction between self and other.

        Args:
            first: first `Interval`
            second: second `Interval`
            ref: overlap fraction calculation reference:
                'self': first is the reference
                'other': second is the reference
                'symmetric': symmetric
                'larger': between first and second, choose the one that has the larger overlap fraction
        Returns:
            a float in range [0, 1]
        """
        assert ref in ['self', 'other', 'larger', 'symmetric']
        if first.contig != first.contig:
            return 0.0
        else:
            if ref == 'self':
                interval_length = first.end - first.start
            elif ref == 'other':
                interval_length = second.end - second.start
            elif ref == 'larger':
                interval_length = min(first.end - first.start, second.end - second.start)
            else:
                interval_length = max(first.end, second.end) - min(first.start, second.start)

            overlap_length = min(first.end, second.end) - max(first.start, second.start)
            if interval_length == 0:
                if overlap_length == 0:
                    return 1.0
                else:
                    return 0.0
            overlap_length = max(0, overlap_length)
            return float(overlap_length) / interval_length

    def __repr__(self):
        return "({0}, {1}, {2}), VAR_CN: {3}, GT: {4}, GQ: {5:.3f}, VF: {6}, NI: {7}, class: {8}".format(
            self.contig, self.start, self.end,
            self.var_copy_number,
            "N/A" if self.genotype is None else self.genotype,
            self.quality,
            "N/A" if self.variant_frequency is None else "{0:.3f}".format(self.variant_frequency),
            "N/A" if self.num_intervals is None else self.num_intervals,
            "N/A" if self.variant_class is None else self.variant_class)

    __str__ = __repr__


def get_overlapping_variants_set(genome_interval_tree: GenomeIntervalTree,
                                 interval: Interval,
                                 ref: Optional[str],
                                 min_overlap_fraction: Optional[float]) \
        -> Set[Union[GenericCopyNumberVariant, Tuple[GenericCopyNumberVariant, float]]]:
    assert (ref is not None and min_overlap_fraction is not None) or (ref is None and min_overlap_fraction is None)
    tree_query = genome_interval_tree[interval.contig].search(interval.start, interval.end)
    found_variants = {found_interval.data for found_interval in tree_query}

    if min_overlap_fraction is None:
        return found_variants
    else:
        variant_and_fraction_set = set()
        for variant in found_variants:
            variant_and_fraction_set.add(
                (variant, GenericCopyNumberVariant.get_overlap_fraction(variant, interval, ref)))
        return {pair for pair in variant_and_fraction_set if pair[1] >= min_overlap_fraction}


def overlaps(genome_interval_tree: GenomeIntervalTree,
             interval: Interval,
             ref: Optional[str],
             min_overlap_fraction: Optional[float]) -> bool:
    if min_overlap_fraction is None:
        return genome_interval_tree[interval.contig].overlaps(interval.start, interval.end)
    else:
        assert ref is not None
        tree_query = genome_interval_tree[interval.contig].search(interval.start, interval.end)
        return any(
            GenericCopyNumberVariant.get_overlap_fraction(
                Interval(interval.contig, found_interval.begin, found_interval.end),
                interval, ref) >= min_overlap_fraction
            for found_interval in tree_query)


GenericCNVCallSetPickleBundle = namedtuple(
    'GenericCNVCallSetPickleBundle', 'sample_name, tags, genome_interval_tree')


class GenericCNVCallSet:
    """This class represents a generic CNV call set and is used for standardizing the output of different
    tools for fair comparison."""
    def __init__(self, sample_name: str, tags: Set[str] = set()):
        self.sample_name = sample_name
        self.tags = tags
        self.genome_interval_tree = GenomeIntervalTree()

    @staticmethod
    def from_genome_interval_tree(sample_name: str,
                                  genome_interval_tree: GenomeIntervalTree,
                                  tags: Set[str] = set()):
        call_set = GenericCNVCallSet(sample_name, tags=tags)
        call_set.genome_interval_tree = genome_interval_tree
        return call_set

    @staticmethod
    def from_pickle(pickle_file: str) -> Tuple[Dict[str, 'GenericCNVCallSet'], GenomeIntervalTree]:
        with open(pickle_file, 'rb') as f:
            unpickler = pickle.Unpickler(f)
            included_loci = unpickler.load()
            call_set_dict = dict()
            while True:
                try:
                    pickle_bundle = unpickler.load()
                    call_set_dict[pickle_bundle.sample_name] = GenericCNVCallSet.from_genome_interval_tree(
                        pickle_bundle.sample_name, pickle_bundle.genome_interval_tree, pickle_bundle.tags)
                except EOFError:
                    break
            return call_set_dict, included_loci

    @staticmethod
    def to_pickle(pickle_file: str,
                  call_set_dict: Dict[str, 'GenericCNVCallSet'],
                  included_loci: GenomeIntervalTree):
        with open(pickle_file, 'wb') as f:
            pickle.dump(included_loci, f)
            for call_set in call_set_dict.values():
                pickle_bundle = GenericCNVCallSetPickleBundle(
                    call_set.sample_name,
                    call_set.tags,
                    call_set.genome_interval_tree)
                pickle.dump(pickle_bundle, f)

    def add(self, variant: GenericCopyNumberVariant):
        self.genome_interval_tree.addi(variant.contig, variant.start, variant.end, data=variant)

    def get_overlapping_variants_set(self,
                                     interval: Interval,
                                     ref: Optional[str],
                                     min_overlap_fraction: Optional[float]) \
            -> Set[Union[GenericCopyNumberVariant, Tuple[GenericCopyNumberVariant, float]]]:
        return get_overlapping_variants_set(self.genome_interval_tree, interval, ref, min_overlap_fraction)

    def overlaps(self, interval: Interval) -> bool:
        return self.genome_interval_tree[interval.contig].overlaps(interval.start, interval.end)

    def iter_in_contig(self, contig: str) -> Generator[GenericCopyNumberVariant, None, None]:
        for entry in self.genome_interval_tree[contig].iter():
            yield entry.data

    def get_contig_interval_tree(self, contig: str) -> IntervalTree:
        return self.genome_interval_tree[contig]

    @property
    def contig_set(self) -> Set[str]:
        return set(self.genome_interval_tree.keys())

    def merge_overlapping_variants(self, padding: int,
                                   interval_consensus_strategy: str = 'highest_quality',
                                   call_consensus_strategy: str = 'highest_quality') -> 'GenericCNVCallSet':
        """Returns a new call set in which nearby variants are merged together.

        Args:
            padding: (positive integer) the amount by which the intervals are to be padded
                symmetrically before detecting overlapping sets
            interval_consensus_strategy: strategy for obtaining the consensus interval.
                the choices are are follows:
                `highest_quality`: select the highest quality variant in the set
                `envelope`: yields the smallest interval that envelopes all variants in the
                overlapping set
            call_consensus_strategy: strategy for obtaining the consensus call (i.e.
                var copy number, quality, etc.). the choices are as follows:
                'highest_quality`: calls from the highest quality variant will be copied
                `average`: use average of the calls in the overlapping set
        """

        assert interval_consensus_strategy in {'highest_quality', 'envelope'}, \
            "Unrecognized interval consensus strategy"
        assert call_consensus_strategy in {'highest_quality', 'average'}, \
            "Unrecognized call consensus strategy"

        def generate_consensus_variant(_overlapping_variants_set: Set[GenericCopyNumberVariant]) \
                -> GenericCopyNumberVariant:
            """From a set of overlapping variants, chooses the variant with highest quality as
            a representative, and returns a new variant that envelopes the overlapping set with
            the same attributes as the highest quality variant.

            Args:
                _overlapping_variants_set: a set of overlapping `GenericCopyNumberVariant`, each
                    padded symmetrically by `padding`
            """
            assert len(_overlapping_variants_set) > 0

            if len(_overlapping_variants_set) == 1:
                return _overlapping_variants_set.__iter__().__next__()

            _contig = next(iter(_overlapping_variants_set)).contig
            num_intervals = sum([_variant.num_intervals for _variant in _overlapping_variants_set])

            highest_quality_variant = sorted(_overlapping_variants_set, key=lambda _variant: _variant.quality)[-1]

            if interval_consensus_strategy == 'highest_quality':
                start = highest_quality_variant.start + padding
                end = highest_quality_variant.end - padding
            elif interval_consensus_strategy == 'envelope':
                start = min([_variant.start for _variant in _overlapping_variants_set]) + padding
                end = max([_variant.end for _variant in _overlapping_variants_set]) - padding
            else:
                raise Exception("Should not reach here")

            if call_consensus_strategy == 'highest_quality':
                genotype = highest_quality_variant.genotype
                var_copy_number = highest_quality_variant.var_copy_number
                quality = highest_quality_variant.quality
                variant_frequency = highest_quality_variant.variant_frequency
                variant_class = highest_quality_variant.variant_class
            elif call_consensus_strategy == 'average':
                raise NotImplementedError
            else:
                raise Exception("Should not reach here")

            return GenericCopyNumberVariant(
                _contig, start, end,
                var_copy_number, quality,
                genotype=genotype,
                variant_frequency=variant_frequency,
                num_intervals=num_intervals,
                variant_class=variant_class)

        # create an empty genome interval tree
        merged_genome_interval_tree = GenomeIntervalTree()

        # merge variants in each contig if necessary
        num_merge_operations = 0
        for contig in self.contig_set:
            sorted_var_list = [iv.data for iv in sorted(list(self.get_contig_interval_tree(contig)))]

            merged_interval_tree = IntervalTree()
            current_overlapping_set = set()
            current_enveloping_interval = None
            for variant in sorted_var_list:
                padded_variant = variant.get_padded(padding)
                if len(current_overlapping_set) == 0:
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig, padded_variant.start, padded_variant.end)
                elif current_enveloping_interval.overlaps_with(padded_variant):
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig,
                        min(padded_variant.start, current_enveloping_interval.start),
                        max(padded_variant.end, current_enveloping_interval.end))
                else:
                    merged_variant = generate_consensus_variant(current_overlapping_set)
                    if len(current_overlapping_set) > 1:
                        num_merge_operations += 1
                    merged_interval_tree.addi(merged_variant.start, merged_variant.end, data=merged_variant)
                    current_overlapping_set.clear()
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig, padded_variant.start, padded_variant.end)
            if len(current_overlapping_set) > 0:
                merged_variant = generate_consensus_variant(current_overlapping_set)
                if len(current_overlapping_set) > 1:
                    num_merge_operations += 1
                merged_interval_tree.addi(merged_variant.start, merged_variant.end,
                                          data=merged_variant)

            merged_genome_interval_tree[contig] = merged_interval_tree

        tags = self.tags.copy()
        tags.add("Merged overlapping intervals with padding={0}, "
                 "interval_consensus_strategy={1}, "
                 "call_consensus_strategy={2}, "
                 "num_merge_operations={3}".format(
                    padding, interval_consensus_strategy, call_consensus_strategy, num_merge_operations))
        return GenericCNVCallSet.from_genome_interval_tree(
            self.sample_name, merged_genome_interval_tree, tags)


class TruthAndTrialVariants:
    """Stores a par of truth and trial variants."""
    def __init__(self,
                 truth_variant: GenericCopyNumberVariant,
                 trial_variant: GenericCopyNumberVariant,
                 overlap_fraction: float,
                 num_overlaps: int,
                 ref: str):
        self.truth_variant = truth_variant
        self.trial_variant = trial_variant
        self.overlap_fraction = overlap_fraction
        self.num_overlaps = num_overlaps
        self.ref = ref

    def __repr__(self):
        return "[truth: {0}\ntrial: {1}\noverlap fraction: {2:f}\n" \
               "number of overlaps: {3}\nreference: {4}]".format(self.truth_variant, self.trial_variant,
                                                                 self.overlap_fraction, self.num_overlaps, self.ref)


class CNVCallSetAnalysisSummary:
    """Stores a summary of CNV call set analysis."""

    def __init__(self):
        # overlaps with a truth variant, has the same copy number call
        self.exact_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are DUP
        self.qualitative_dup_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are DEL
        self.qualitative_del_matches: Set[TruthAndTrialVariants] = set()

        # truth is dup(del) but call is del(dup)
        self.opposite_calls: Set[TruthAndTrialVariants] = set()

        # dup truth variants that are not called (the set contains such truth variants)
        self.missed_dup_calls: Set[GenericCopyNumberVariant] = set()

        # del truth variants that are not called (the set contains such truth variants)
        self.missed_del_calls: Set[GenericCopyNumberVariant] = set()

        # truth is ref but call is dup (the set contains such trial variants)
        self.false_dup_calls: Set[GenericCopyNumberVariant] = set()

        # truth is ref but call is del (the set contains such trial variants)
        self.false_del_calls: Set[GenericCopyNumberVariant] = set()

        # trial variant excluded in truth analysis (the set contains such truth variants)
        self.trial_excluded_calls: Set[GenericCopyNumberVariant] = set()

        # truth variant excluded in trial analysis (the set contains such trial variants)
        self.truth_excluded_calls: Set[GenericCopyNumberVariant] = set()

    def copy(self):
        summary_copy = CNVCallSetAnalysisSummary()
        summary_copy.exact_matches = self.exact_matches.copy()
        summary_copy.qualitative_dup_matches = self.qualitative_dup_matches.copy()
        summary_copy.qualitative_del_matches = self.qualitative_del_matches.copy()
        summary_copy.opposite_calls = self.opposite_calls.copy()
        summary_copy.missed_dup_calls = self.missed_dup_calls.copy()
        summary_copy.missed_del_calls = self.missed_del_calls.copy()
        summary_copy.false_dup_calls = self.false_dup_calls.copy()
        summary_copy.false_del_calls = self.false_del_calls.copy()
        summary_copy.trial_excluded_calls = self.trial_excluded_calls.copy()
        summary_copy.truth_excluded_calls = self.truth_excluded_calls.copy()
        return summary_copy

    def get_filtered_summary(self,
                             truth_filter: Callable[[GenericCopyNumberVariant], bool],
                             trial_filter: Callable[[GenericCopyNumberVariant], bool]):
        filtered_summary = CNVCallSetAnalysisSummary()

        # variants that overlap between truth and trial
        for original_match_set, filtered_match_set in zip([
                self.exact_matches,
                self.qualitative_dup_matches,
                self.qualitative_del_matches,
                self.opposite_calls], [
                filtered_summary.exact_matches,
                filtered_summary.qualitative_dup_matches,
                filtered_summary.qualitative_del_matches,
                filtered_summary.opposite_calls]):

            for joint_vars in original_match_set:
                truth_pass = truth_filter(joint_vars.truth_variant)
                trial_pass = trial_filter(joint_vars.trial_variant)

                if truth_pass:
                    if trial_pass:  # count as a match
                        filtered_match_set.add(joint_vars)
                    else:  # count as a missed variant
                        if joint_vars.truth_variant.is_del:
                            filtered_summary.missed_del_calls.add(joint_vars.truth_variant)
                        elif joint_vars.truth_variant.is_dup:
                            filtered_summary.missed_dup_calls.add(joint_vars.truth_variant)
                        else:
                            raise Exception("Should not reach here!")
                else:  # if a truth variant (that overlaps a trial) does not pass, exclude both
                    if trial_pass:  # count as false positive
                        if joint_vars.trial_variant.is_del:
                            filtered_summary.false_del_calls.add(joint_vars.trial_variant)
                        elif joint_vars.trial_variant.is_dup:
                            filtered_summary.false_dup_calls.add(joint_vars.trial_variant)
                        else:
                            raise Exception("Should not reach here!")
                    else:  # dismiss
                        filtered_summary.truth_excluded_calls.add(joint_vars.truth_variant)
                        filtered_summary.trial_excluded_calls.add(joint_vars.trial_variant)

        # missed truth variants
        for original_missed_set, filtered_missed_set in zip(
                [self.missed_dup_calls, self.missed_del_calls],
                [filtered_summary.missed_dup_calls, filtered_summary.missed_del_calls]):

            for var in original_missed_set:
                truth_pass = truth_filter(var)
                if truth_pass:  # if truth passes, it is still a false negative
                    filtered_missed_set.add(var)
                else:  # otherwise, exclude this truth variant
                    filtered_summary.truth_excluded_calls.add(var)

        # false trial variants
        for original_false_set, filtered_false_set in zip(
                [self.false_dup_calls, self.false_del_calls],
                [filtered_summary.false_dup_calls, filtered_summary.false_del_calls]):

            for var in original_false_set:
                trial_pass = trial_filter(var)
                if trial_pass:  # if trial passes, it is still a false positive
                    filtered_false_set.add(var)
                else:  # otherwise, exclude this trial variant
                    filtered_summary.trial_excluded_calls.add(var)

        # excluded
        for variant in self.truth_excluded_calls:
            filtered_summary.truth_excluded_calls.add(variant)

        for variant in self.trial_excluded_calls:
            filtered_summary.trial_excluded_calls.add(variant)

        return filtered_summary

    @property
    def all_matches(self):
        return self.exact_matches.union(self.qualitative_del_matches).union(self.qualitative_dup_matches)

    @property
    def all_misses(self):
        return self.missed_del_calls.union(self.missed_dup_calls)

    @property
    def all_falsities(self):
        return self.false_del_calls.union(self.false_dup_calls)

    @property
    def true_positive_count(self):
        exact_matches_count = len(self.exact_matches)
        qualitative_dup_matches_count = len(self.qualitative_dup_matches)
        qualitative_del_matches_count = len(self.qualitative_del_matches)
        return exact_matches_count + qualitative_dup_matches_count + qualitative_del_matches_count

    @property
    def false_positive_count(self):
        false_dup_calls_count = len(self.false_dup_calls)
        false_del_calls_count = len(self.false_del_calls)
        return false_dup_calls_count + false_del_calls_count

    @property
    def false_negative_count(self):
        missed_dup_calls_count = len(self.missed_dup_calls)
        missed_del_calls_count = len(self.missed_del_calls)
        return missed_dup_calls_count + missed_del_calls_count

    @property
    def opposite_calls_count(self):
        return len(self.opposite_calls)

    @property
    def trial_excluded_calls_count(self):
        return len(self.trial_excluded_calls)

    @property
    def truth_excluded_calls_count(self):
        return len(self.truth_excluded_calls)

    @property
    def sensitivity(self):
        total_truth_calls = self.true_positive_count + self.false_negative_count
        if total_truth_calls > 0:
            return float(self.true_positive_count) / total_truth_calls
        else:
            return np.nan

    @property
    def precision(self):
        total_trial_calls = self.true_positive_count + self.false_positive_count + self.opposite_calls_count
        if total_trial_calls > 0:
            return float(self.true_positive_count) / total_trial_calls
        else:
            return np.nan

    def __str__(self):
        output = ""
        output += "- Number of exact matches: {0}\n".format(len(self.exact_matches))
        output += "- Number of qualitative DUP matches: {0}\n".format(len(self.qualitative_dup_matches))
        output += "- Number of qualitative DEL matches: {0}\n".format(len(self.qualitative_del_matches))
        output += "- Number of missed DUP variants: {0}\n".format(len(self.missed_dup_calls))
        output += "- Number of missed DEL variants: {0}\n".format(len(self.missed_del_calls))
        output += "- Number of false DUP variants: {0}\n".format(len(self.false_dup_calls))
        output += "- Number of false DEL variants: {0}\n".format(len(self.false_del_calls))
        output += "- Number of opposite calls: {0}\n".format(len(self.opposite_calls))
        output += "- Number of excluded truth variants: {0}\n".format(len(self.truth_excluded_calls))
        output += "- Number of excluded trial variants: {0}\n".format(len(self.trial_excluded_calls))
        output += "* Sensitivity: {0:.3f}\n".format(self.sensitivity)
        output += "* Precision: {0:.3f}\n".format(self.precision)
        return output

    __repr__ = __str__


class CNVTrialCallSetEvaluator:
    def __init__(self,
                 truth_call_set: GenericCNVCallSet,
                 trial_call_set: GenericCNVCallSet,
                 trial_included_loci: GenomeIntervalTree,
                 min_overlap_fraction_for_truth_filtration: float,
                 min_overlap_fraction_for_variant_matching: float,
                 truth_overlapping_set_selection_strategy: str = "largest_overlap"):
        """Evaluates a trial CNV call set against a truth call set in genomic regions analyzed by both
        the truth and the trial callers.

        Args:
            truth_call_set: truth call set
            trial_call_set: trial call set
            trial_included_loci: genomic loci included by the trial caller
            min_overlap_fraction_for_variant_matching: consider two variants overlapping only if
                their overlap fraction is higher than this value. If not provided, any overlap
                is considered.
            truth_overlapping_set_selection_strategy: "largest_overlap" or "highest_quality"
        """
        if truth_call_set.sample_name != trial_call_set.sample_name:
            _logger.warning("CallSets have different samples names; truth: {0}, trial: {1}".format(
                truth_call_set.sample_name, trial_call_set.sample_name))

        self.sample_name = truth_call_set.sample_name
        self.truth_call_set = truth_call_set
        self.trial_call_set = trial_call_set
        self.trial_included_loci = trial_included_loci
        self.min_overlap_fraction_for_variant_matching = min_overlap_fraction_for_variant_matching
        self.truth_overlapping_set_selection_strategy = truth_overlapping_set_selection_strategy
        self.all_contigs = truth_call_set.contig_set.union(trial_call_set.contig_set)

        # filter truth
        self.relevant_truth_call_set = GenericCNVCallSet(truth_call_set.sample_name, truth_call_set.tags)
        for contig in self.all_contigs:
            for truth_variant in truth_call_set.iter_in_contig(contig):
                if overlaps(trial_included_loci, truth_variant, 'self', min_overlap_fraction_for_truth_filtration):
                    self.relevant_truth_call_set.add(truth_variant)

    def __call__(self) -> CNVCallSetAnalysisSummary:
        """Generates the analysis summary."""

        def perform_analysis_for_contig(_contig: str, _summary: CNVCallSetAnalysisSummary) -> None:
            """Performs the evaluation for a given contig and updates the analysis summary."""

            # keeps track of truth variants that do not overlap with trial variants
            remaining_truth_variants = self.relevant_truth_call_set.get_contig_interval_tree(_contig).copy()

            # first, iterate over trial variants
            for trial_variant in self.trial_call_set.iter_in_contig(_contig):

                # get a set of overlapping truth variants
                truth_overlapping_variants_and_overlap_fractions = get_overlapping_variants_set(
                    self.relevant_truth_call_set.genome_interval_tree,
                    trial_variant,
                    'self',
                    self.min_overlap_fraction_for_variant_matching)

                if len(truth_overlapping_variants_and_overlap_fractions) == 0:  # not in truth

                    if trial_variant.is_dup:
                        _summary.false_dup_calls.add(trial_variant)
                    else:
                        _summary.false_del_calls.add(trial_variant)

                else:

                    # choose a truth variant from the overlapping set
                    num_overlaps = len(truth_overlapping_variants_and_overlap_fractions)
                    if self.truth_overlapping_set_selection_strategy == "largest_overlap":
                        best = max(
                            truth_overlapping_variants_and_overlap_fractions,
                            key=lambda variant_and_fraction_pair: variant_and_fraction_pair[1])
                        best_truth_variant = best[0]
                        overlap_fraction = best[1]
                    elif self.truth_overlapping_set_selection_strategy == "highest_quality":
                        best = max(
                            truth_overlapping_variants_and_overlap_fractions,
                            key=lambda variant_and_fraction_pair: variant_and_fraction_pair[0].quality)
                        best_truth_variant = best[0]
                        overlap_fraction = best[1]
                    else:
                        raise Exception("Should not reach here!")

                    # remove matched truth variants from the remaining truth variants
                    try:
                        remaining_truth_variants.removei(
                            best_truth_variant.start, best_truth_variant.end, data=best_truth_variant)
                    except ValueError:
                        # a variant could be already removed and it is fine
                        # (note that we have queried truth_call_set)
                        pass

                    joint_variant = TruthAndTrialVariants(best_truth_variant, trial_variant,
                                                          overlap_fraction,
                                                          num_overlaps,
                                                          'self')

                    # exact match
                    if best_truth_variant.var_copy_number == trial_variant.var_copy_number:
                        _summary.exact_matches.add(joint_variant)
                        continue

                    # qualitative dup match
                    if best_truth_variant.is_dup and trial_variant.is_dup:
                        _summary.qualitative_dup_matches.add(joint_variant)
                        continue

                    # qualitative del match
                    if best_truth_variant.is_del and trial_variant.is_del:
                        _summary.qualitative_del_matches.add(joint_variant)
                        continue

                    # opposite
                    _summary.opposite_calls.add(joint_variant)

            # next, iterate over remaining truth variants
            for truth_entry in remaining_truth_variants.iter():
                truth_variant: GenericCopyNumberVariant = truth_entry.data

                if truth_variant.is_dup:
                    _summary.missed_dup_calls.add(truth_variant)
                else:
                    _summary.missed_del_calls.add(truth_variant)

        summary = CNVCallSetAnalysisSummary()

        for contig in self.all_contigs:
            perform_analysis_for_contig(contig, summary)
        return summary
