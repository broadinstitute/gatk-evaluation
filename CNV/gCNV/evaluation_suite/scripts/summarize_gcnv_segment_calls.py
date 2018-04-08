import os
import json
import argparse
import numpy as np
import matplotlib.pylab as plt
import pickle

import germline_cnv_evaluation as cnv_eval
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')

_logger = logging.getLogger("summarize_gcnv_segment_calls")

parser = argparse.ArgumentParser(description="gCNV vs. Genome STRiP analysis")

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--gcnv_segments_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Path to gCNV segments")

group.add_argument("--gcnv_workflow_json_file",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="JSON file corresponding to the gCNV run")

group.add_argument("--gcnv_intervals_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="A path containing the interval_list file mentioned in the JSON file")

group.add_argument("--gs_pkl_file",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Genome STRiP pickled call-set")

group.add_argument("--output_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path")

args = parser.parse_args()

run_prefix = args.gcnv_workflow_json_file.split('/')[-1].replace(".json", "")

gcnv_segments_vcf_files = [os.path.join(args.gcnv_segments_path, f)
                           for f in os.walk(args.gcnv_segments_path).__next__()[2]
                           if "genotyped-segments" in f]

with open(args.gcnv_workflow_json_file, 'r') as f:
    job_config = json.load(f)
job_interval_list_filename = job_config['CNVGermlineCohortWorkflow.intervals'].split('/')[-1]

trial_call_set_dict, trial_included_loci = cnv_eval.load_gcnv_segments_vcf_file(
    gcnv_segments_vcf_files,
    os.path.join(args.gcnv_intervals_path, job_interval_list_filename),
    quality_mode="some",
    min_overlap_fraction_for_variant_matching=0.5)

cnv_eval.GenericCNVCallSet.to_pickle(
    os.path.join(args.output_path, run_prefix + ".pkl"), trial_call_set_dict, trial_included_loci)





# load pickled results
truth_call_set_dict, truth_included_loci = cnv_eval.GenericCNVCallSet.from_pickle(args.gs_pkl_file)

# for GS calls: no need to pad segments before merging (by design, the best segment for each sample
# is already there)
truth_merge_padding = 0

# for GS calls: choose the highest quality segment from the overlapping set (it is reasonable
# to assume that the highest quality segment for each sample is the right call for that sample)
truth_merge_interval_consensus_strategy = 'highest_quality'

# for GS calls: the intervals end-points are the same as the highest quality segment
truth_merge_call_consensus_strategy = 'highest_quality'

# merge overlapping segments
truth_merged_call_set_dict = dict()
for sample_name in truth_call_set_dict.keys():
    truth_merged_call_set_dict[sample_name] = truth_call_set_dict[sample_name].merge_overlapping_variants(
        truth_merge_padding,
        interval_consensus_strategy=truth_merge_interval_consensus_strategy,
        call_consensus_strategy=truth_merge_call_consensus_strategy)

truth_sample_names = set(truth_call_set_dict.keys())
mutual_samples = truth_sample_names.intersection(trial_call_set_dict.keys())
mutual_samples = list(mutual_samples)

summary_dict = dict()
for sample_name in mutual_samples:
    summary_dict[sample_name] = cnv_eval.CNVTrialCallSetEvaluatorTargetResolved(
        truth_merged_call_set_dict[sample_name],
        trial_call_set_dict[sample_name],
        trial_included_loci)()

# pickle summary
summary_pkl_file = os.path.join(args.output_path, run_prefix + "-gs-concordance-summary.pkl")
with open(summary_pkl_file, 'w') as f:
    pickle.dump(summary_dict, f)


def get_concordance_from_summary_dict(truth_filter, trial_filter):
    TP = 0.
    FP = 0.
    FN = 0.
    TN = 0.
    for _sample_name in mutual_samples:
        c = summary_dict[_sample_name].get_filtered_summary(truth_filter, trial_filter).target_level_concordance
        TP += c.TP
        FP += c.FP
        FN += c.FN
        TN += c.TN
    return cnv_eval.core.Concordance(TP=TP, FP=FP, FN=FN, TN=TN)


def get_roc_curve(truth_min_variant_frequency,
                  truth_max_variant_frequency,
                  truth_min_length,
                  truth_min_quality,
                  trial_max_quality):

    truth_filter = cnv_eval.GenericCopyNumberVariant.get_variant_filter(
        min_quality=truth_min_quality,
        min_variant_frequency=truth_min_variant_frequency,
        max_variant_frequency=truth_max_variant_frequency,
        included_variant_classes={'mixed'},
        min_length=truth_min_length)

    _fdr_list = []
    _fpr_list = []
    _tpr_list = []

    for trial_min_quality in np.linspace(0, trial_max_quality, 20):

        trial_filter = cnv_eval.GenericCopyNumberVariant.get_variant_filter(
            min_quality=trial_min_quality)

        concordance = get_concordance_from_summary_dict(truth_filter, trial_filter)

        _fdr_list.append(1. - concordance.PPV)
        _tpr_list.append(concordance.TPR)
        _fpr_list.append(concordance.FPR)

    return _fdr_list, _tpr_list, _fpr_list


fdr_list, tpr_list, fpr_list = get_roc_curve(0., 1., 0, truth_min_quality=0, trial_max_quality=99)

fig = plt.figure()
ax = plt.gca()
ax.plot(fpr_list, tpr_list, label=run_prefix)
ax.set_ylim((0, 1))
ax.set_xlabel('false positive rate')
ax.set_ylabel('true positive rate')
ax.legend()
plt.savefig(os.path.join(args.output_path, run_prefix + "-ROC.pdf"))

