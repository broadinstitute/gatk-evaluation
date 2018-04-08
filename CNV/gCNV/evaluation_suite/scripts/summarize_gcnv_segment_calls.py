import os
import json
import argparse

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

call_set, included_intervals = cnv_eval.load_gcnv_segments_vcf_file(
    gcnv_segments_vcf_files,
    os.path.join(args.gcnv_intervals_path, job_interval_list_filename),
    quality_mode="some",
    min_overlap_fraction_for_variant_matching=0.5)


cnv_eval.GenericCNVCallSet.to_pickle(
    os.path.join(args.output_path, run_prefix + ".pkl", call_set, included_intervals))

