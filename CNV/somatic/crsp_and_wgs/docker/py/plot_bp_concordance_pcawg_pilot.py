import math

import matplotlib
# Force matplotlib to not use any Xwindows backend.
from pandas.core.frame import DataFrame
from pandas.core.series import Series

CONTIG_COL = "CONTIG"
START_COL = "START"
END_COL = "END"
LOG2_COPY_RATIO_POSTERIOR_90_COL = "LOG2_COPY_RATIO_POSTERIOR_90"
LOG2_COPY_RATIO_POSTERIOR_10_COL = "LOG2_COPY_RATIO_POSTERIOR_10"
LARGE_EVENT_TYPE_COL = "type"
NUM_BASES_COL = "num_bases"
POSSIBLE_GERMLINE_COL = "POSSIBLE_GERMLINE"

STAR_COL = "star"

PURITY_PLOIDY_INPUT_PURITY_CONF_MAD_COL = "purity_conf_mad"
PURITY_PLOIDY_INPUT_PLOIDY_COL = "ploidy"
PURITY_PLOIDY_INPUT_PURITY_COL = "purity"
PURITY_PLOIDY_INPUT_SAMPLENAME_COL = "samplename"

FINAL_RESULTS_STARS_COL = "stars"
FINAL_RESULTS_MAE_NO_GERMLINE_COL = "mae_no_germline"
FINAL_RESULTS_MAE_COL = "mae"
FINAL_RESULTS_BASES_EVALUATED_NO_GERMLINE_COL = "num_bases_evaluated_no_germline"
FINAL_RESULTS_BASES_EVALUATED_COL = "num_bases_evaluated"
FINAL_RESULTS_PROP_SUCCESS_NO_GERMLINE_COL = "prop_success_no_germline"
FINAL_RESULTS_PROP_SUCCESS_COL = "prop_success"
FINAL_RESULTS_DISTANCE_SUCCESS_COL = "distance_success"
FINAL_RESULTS_PLOIDY_COL = "ploidy"
FINAL_RESULTS_PURITY_COL = "purity"
FINAL_RESULTS_SAMPLE_COL = "samplename"


matplotlib.use('Agg')

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
import matplotlib.colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy
import pandas
from sklearn.metrics import mean_absolute_error as mae

CR_GUESS_COLUMN = "LOG2_COPY_RATIO_POSTERIOR_50"
PURITY_PLOIDY_FILENAME = "consensus.20170217.purity.ploidy.txt"
UUID_BARCODE_FILE = "UUID-BARCODE.tsv"
GT_CN_COLUMN_NAME = "final_total_cn"
GT_CR_COLUMN_NAME = GT_CN_COLUMN_NAME.replace("cn", "cr")
MIN_SEGMENT_LENGTH_BP = 0
MIN_STARS = 2.5
PREVALENCE_MIN_FOR_BLACKLIST = 0.1
FRACMATCH_MIN = 0.95
FINAL_RESULTS_COLUMN_LIST = [FINAL_RESULTS_SAMPLE_COL, FINAL_RESULTS_PURITY_COL, FINAL_RESULTS_PLOIDY_COL,
                             FINAL_RESULTS_DISTANCE_SUCCESS_COL,
                             FINAL_RESULTS_PROP_SUCCESS_COL, FINAL_RESULTS_PROP_SUCCESS_NO_GERMLINE_COL,
                             FINAL_RESULTS_BASES_EVALUATED_COL,
                             FINAL_RESULTS_BASES_EVALUATED_NO_GERMLINE_COL, FINAL_RESULTS_MAE_COL, FINAL_RESULTS_MAE_NO_GERMLINE_COL,
                             FINAL_RESULTS_STARS_COL]


def calculate_mae_with_weights(comparison_df):
    # type: (DataFrame) -> float
    comparison_df_mad = comparison_df[~(comparison_df[CR_GUESS_COLUMN].isnull()) & ~(comparison_df[GT_CR_COLUMN_NAME].isnull())]
    if len(comparison_df_mad) == 0:
        return 0.0
    return mae(comparison_df_mad[CR_GUESS_COLUMN], comparison_df_mad[GT_CR_COLUMN_NAME],
               sample_weight=comparison_df_mad[NUM_BASES_COL])


def calculate_mae_no_weights(comparison_df):
    # type: (DataFrame) -> float
    comparison_df_mad = comparison_df[~(comparison_df[CR_GUESS_COLUMN].isnull()) & ~(comparison_df[GT_CR_COLUMN_NAME].isnull())]
    if len(comparison_df_mad) == 0:
        return 0.0
    return mae(comparison_df_mad[CR_GUESS_COLUMN], comparison_df_mad[GT_CR_COLUMN_NAME])


def plot_bp_seg_concordance(output_dir, df_to_plot, short_name, other_lines_for_title=None, title_plot=None):
    # type: (str, DataFrame, str, list[str], str) -> None
    """
    Save a png file of the concordance plot between WGS ground truth and the test case.

    :param output_dir: Directory to save the output plots.
    :param df_to_plot: DataFrame to plot.  Includes columns for both the ground truth and the test case.
    :param short_name: a string that can be used as an identifier in the filename of the plot.
    :param other_lines_for_title: Arbitrary strings to be included in the title.
    :param title_plot:
    :return:
    """
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    if title_plot is None:
        title_plot = "Concordance"
    title_plot += "\nSegment stars: " + str(df_to_plot["star"].unique().tolist())

    if all([math.isnan(x) for x in list(df_to_plot[CR_GUESS_COLUMN])]):
        print("All values are nan.  Cannot plot...")
        return

    if other_lines_for_title is None:
        str_other_lines_for_title = ""
    else:
        str_other_lines_for_title = "\n" + "\n".join(other_lines_for_title)

    h = plt.figure()
    mae_bp = calculate_mae_with_weights(df_to_plot)
    hi_cn = 3.5
    step_cn = 0.1
    plt.hist2d(df_to_plot[CR_GUESS_COLUMN], df_to_plot[GT_CR_COLUMN_NAME],
               bins=[numpy.arange(0, hi_cn, step_cn),numpy.arange(0, hi_cn, step_cn)], cmin=0.0001,
               norm=matplotlib.colors.LogNorm(), cmap='viridis_r',
               weights=df_to_plot[NUM_BASES_COL])
    plt.xlabel("GATK CNV 1.5 Model Segments Copy Ratio")
    plt.ylabel("Calculated Ground Truth Copy Ratio")
    plt.title(title_plot + "  bp count.  MAD: " + ("%2.3f" % mae_bp) + str_other_lines_for_title)
    plt.colorbar()
    h.gca().add_line(plt.Line2D([0, hi_cn], [0, hi_cn], ls=":", lw=1, color='g'))
    h.gca().add_patch(patches.Rectangle((0.9, 0.9), 0.2, 0.2, facecolor="#aaaaaa", edgecolor="#000000"))
    h.set_size_inches(9.5, 5.5)
    h.savefig(output_dir + short_name + "_weighted_by_bp.png", dpi=90)
    plt.close(h)

    h = plt.figure()
    mae_seg = calculate_mae_no_weights(df_to_plot)
    plt.hist2d(df_to_plot[CR_GUESS_COLUMN], df_to_plot[GT_CR_COLUMN_NAME],
               bins=[numpy.arange(0, hi_cn, step_cn),numpy.arange(0, hi_cn, step_cn)], cmin=0.0001,
               norm=matplotlib.colors.LogNorm(), cmap='viridis_r')
    plt.xlabel("GATK CNV 1.5 Model Segments Copy Ratio")
    plt.ylabel("Calculated Ground Truth Copy Ratio")
    plt.title(title_plot + "  seg count.  MAD: " + ("%2.3f" % mae_seg) + str_other_lines_for_title)
    plt.colorbar()
    h.gca().add_line(plt.Line2D([0, hi_cn], [0, hi_cn], ls=":", lw=1, color='g'))
    h.gca().add_patch(patches.Rectangle((0.9, 0.9), 0.2, 0.2, facecolor="#aaaaaa", edgecolor="#000000"))
    h.set_size_inches(9.5, 5.5)
    h.savefig(output_dir + short_name + "_weighted_by_seg.png", dpi=90)
    plt.close(h)


def retrieve_purity_ploidy(gt_fn):
    # type: (str) -> (float,float,float)
    """
    Retrieve the purity and ploidy from the provided consensus file.  There is a hack in that this depends on the
     filename.

     Purity/ploidy files uses the same sample names as PCAWG, not TCGA, so the samplename column will match the ground truth segment files.

    :param gt_fn: Filename of PCAWG consensus ground truth
    :return:
    """
    purity_ploidy_df = pandas.read_csv(PURITY_PLOIDY_FILENAME, sep="\t")  # type: DataFrame

    base_filename = os.path.basename(gt_fn)
    final_search_string = base_filename.split(".")[0]

    purity = purity_ploidy_df[purity_ploidy_df[PURITY_PLOIDY_INPUT_SAMPLENAME_COL] == final_search_string][
        PURITY_PLOIDY_INPUT_PURITY_COL].tolist()[0]
    ploidy = purity_ploidy_df[purity_ploidy_df[PURITY_PLOIDY_INPUT_SAMPLENAME_COL] == final_search_string][
        PURITY_PLOIDY_INPUT_PLOIDY_COL].tolist()[0]
    purity_conf_mad = purity_ploidy_df[purity_ploidy_df[PURITY_PLOIDY_INPUT_SAMPLENAME_COL] == final_search_string][
        PURITY_PLOIDY_INPUT_PURITY_CONF_MAD_COL].tolist()[0]
    return purity, ploidy, purity_conf_mad


def create_segments_df_for_comparison(input_tsv, purity, ploidy):
    # type: (str, float, float) -> DataFrame
    segs_df = pandas.read_csv(input_tsv, sep="\t", comment="@")

    cr_gt = 1 + (purity * ((segs_df[GT_CN_COLUMN_NAME] / ploidy) - 1))
    cr_gt.rename(GT_CR_COLUMN_NAME, inplace=True)
    cr = 2 ** segs_df[CR_GUESS_COLUMN]
    cr10 = 2 ** segs_df[LOG2_COPY_RATIO_POSTERIOR_10_COL]
    cr90 = 2 ** segs_df[LOG2_COPY_RATIO_POSTERIOR_90_COL]
    weight = segs_df[END_COL] - segs_df[START_COL]
    weight.rename(NUM_BASES_COL, inplace=True)

    # TODO: Get rid of this next statement, since it will cut additional columns that are added in later veresions
    comparison = pandas.concat([segs_df[CONTIG_COL], segs_df[START_COL], segs_df[END_COL], cr_gt, cr, weight, segs_df[STAR_COL],
                                cr10, cr90, segs_df[POSSIBLE_GERMLINE_COL], segs_df[LARGE_EVENT_TYPE_COL]], axis=1)
    comparison_pruned = comparison[~(((comparison[CR_GUESS_COLUMN] < 1.1) & (comparison[CR_GUESS_COLUMN] > 0.9)) &
                               ((comparison[GT_CR_COLUMN_NAME] < 1.1) & (comparison[GT_CR_COLUMN_NAME] > 0.9)))]

    comparison_pruned = comparison_pruned[comparison_pruned[STAR_COL] > MIN_STARS]
    print("Stars being considered: " + str(comparison_pruned[STAR_COL].unique().astype('int32').tolist()))

    comparison_pruned = comparison_pruned[comparison_pruned[NUM_BASES_COL] > MIN_SEGMENT_LENGTH_BP]
    print("Removing segments that are less than minimum  of " + str(MIN_SEGMENT_LENGTH_BP) + " bp")

    return comparison_pruned


def create_segments_df_for_comparison_germline_removed(segments_df):
    # type: (DataFrame) -> (DataFrame)
    """
    Prune the segments tagged as possible germline.  I.e. the POSSIBLE_GERMLINE column is non-zero.

    :param segments_df: DataFrame with the POSSIBLE_GERMLINE column
    :return: DataFrame
    """
    comparison_edit_germline_removed = segments_df

    # Remove centromeres
    comparison_edit_germline_removed = comparison_edit_germline_removed[
        (comparison_edit_germline_removed["type"].isnull())]

    # Remove possible germline events
    comparison_edit_germline_removed[POSSIBLE_GERMLINE_COL] = comparison_edit_germline_removed[POSSIBLE_GERMLINE_COL].astype('str')
    comparison_edit_germline_removed = comparison_edit_germline_removed[
        (comparison_edit_germline_removed[POSSIBLE_GERMLINE_COL] == "0" or comparison_edit_germline_removed[POSSIBLE_GERMLINE_COL] == "0.0" or comparison_edit_germline_removed[POSSIBLE_GERMLINE_COL] == "0.")]

    return comparison_edit_germline_removed


def parse_options():
    epilog = """ This tool will create plots of the concordance with PCAWG pilot ground truth files. This tool is very fussy about the incoming data being in its expected format.
    """
    desc = "Create PCAWG performance to ground truth for GATK CNV comparison"
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('-I', '--input_files', action='append', help='Input seg files that must be merged with gt.', required=True)
    parser.add_argument('-G', '--gt_files', action='append', help='Ground truth seg files.  Must be in corresponding order to input_files', required=True)
    # python arg.py -I 1234 -I 2345 -I 3456 -I 4567
    parser.add_argument("-O", "--output_dir", type=str, required=True, help="Output dir.  Will be created if it does not exist.")
    # Process arguments
    args_parsed = parser.parse_args()

    return args_parsed

if __name__ == '__main__':
    args = parse_options()
    output_dir = args.output_dir
    input_tsvs = args.input_files
    gt_tsvs = args.gt_files

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    #######################

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    final_results_df = DataFrame(columns=FINAL_RESULTS_COLUMN_LIST)
    for i, input_tsv in enumerate(input_tsvs):
        print("\n============" + input_tsv)
        purity_ploidy_tuple = retrieve_purity_ploidy(gt_tsvs[i])
        purity = purity_ploidy_tuple[0]
        ploidy = purity_ploidy_tuple[1]

        print("purity " + str(purity) + "   ploidy " + str(ploidy))

        comparison_edit = create_segments_df_for_comparison(input_tsv, purity, ploidy)
        if len(comparison_edit) == 0:
            print(input_tsv + " has no segments that remain after initial pruning by star count.")
            continue
        print("Copy neutral region not considered.")
        print("Total number of bases (in Mbp): " + str((comparison_edit[NUM_BASES_COL].astype('float') / 1e6).sum()))

        comparison_edit_germline_removed = create_segments_df_for_comparison_germline_removed(comparison_edit)

        print("Removing common CNV segments, centromeres, and segDupeSegments.")
        print("   Proportion of bases remaining: " + str(float(comparison_edit_germline_removed[NUM_BASES_COL].sum()) / float(comparison_edit[
                                                                                                                          NUM_BASES_COL].sum())))
        print("Total number of bases remaining (in Mbp): " + str((comparison_edit_germline_removed[NUM_BASES_COL].astype('float') / 1e6).sum()))
        print("Median length of segment (breakpoints merged with ground truth): " + str(comparison_edit[
                                                                                            NUM_BASES_COL].median()))
        print("Threshold (min) length of segment (breakpoints merged with ground truth): " + str(MIN_SEGMENT_LENGTH_BP))

        sample_file = os.path.splitext(os.path.basename(input_tsv))[0]

        comparison_edit.to_csv(output_dir + os.path.basename(input_tsv) + "_comparison_edit.txt", sep="\t")
        comparison_edit_germline_removed.to_csv(output_dir + os.path.basename(input_tsv) + "_comparison_edit_germline_removed.txt", sep="\t")

        distance_success = 0.1
        distance_df = abs(comparison_edit_germline_removed[CR_GUESS_COLUMN] - comparison_edit_germline_removed[GT_CR_COLUMN_NAME])
        comparison_edit_germline_removed_success = comparison_edit_germline_removed[distance_df < distance_success]

        distance_all_df = abs(comparison_edit[CR_GUESS_COLUMN] - comparison_edit[GT_CR_COLUMN_NAME])
        comparison_edit_success = comparison_edit_germline_removed[distance_all_df < distance_success]

        success_bp = comparison_edit_germline_removed_success[NUM_BASES_COL].sum()
        all_bp = comparison_edit_germline_removed[NUM_BASES_COL].sum()
        print(str(success_bp))
        print(str(all_bp))
        print("Proportion of bp within " + str(distance_success) + ": " + str(float(success_bp)/float(all_bp)))

        success_bp_all = comparison_edit_success[NUM_BASES_COL].sum()
        all_bp_all = comparison_edit[NUM_BASES_COL].sum()
        print("With germline events unfiltered... Proportion of bp within " + str(distance_success) + ": " + str(float(success_bp_all)/float(all_bp_all)))

        line_dict = {FINAL_RESULTS_SAMPLE_COL: sample_file, FINAL_RESULTS_PURITY_COL: purity,
                     FINAL_RESULTS_PLOIDY_COL: ploidy, FINAL_RESULTS_DISTANCE_SUCCESS_COL: distance_success,
                     FINAL_RESULTS_PROP_SUCCESS_COL: float(success_bp_all) / float(all_bp_all),
                     FINAL_RESULTS_PROP_SUCCESS_NO_GERMLINE_COL: float(success_bp) / float(all_bp),
                     FINAL_RESULTS_BASES_EVALUATED_COL: all_bp_all,
                     FINAL_RESULTS_BASES_EVALUATED_NO_GERMLINE_COL: all_bp,
                     FINAL_RESULTS_MAE_COL: round(calculate_mae_with_weights(comparison_edit), 3),
                     FINAL_RESULTS_MAE_NO_GERMLINE_COL: round(calculate_mae_with_weights(comparison_edit_germline_removed), 3),
                     FINAL_RESULTS_STARS_COL: ",".join(comparison_edit[STAR_COL].unique().astype('str').tolist())}
        final_results_df = final_results_df.append(Series(name=sample_file, data=line_dict))

        other_lines_unfiltered = ["Proportion of bp w/in CR of " + str(distance_success) + ": " + str(float(success_bp_all)/float(all_bp_all)),
                                  "Total bp: " + str(all_bp_all)]
        other_lines_filtered = ["Proportion of bp w/in CR of " + str(distance_success) + ": " + str(float(success_bp)/float(all_bp)),
                                "Total bp: " + str(all_bp)]

        plot_bp_seg_concordance(output_dir, comparison_edit, short_name=sample_file + "Germline_Unfiltered_bp_concordance",
                                title_plot="Concordance Germline Unfiltered", other_lines_for_title=other_lines_unfiltered)
        plot_bp_seg_concordance(output_dir, comparison_edit_germline_removed, short_name=sample_file + "Germline_Filtered_bp_concordance",
                                title_plot="Concordance Germline Filtered", other_lines_for_title=other_lines_filtered)

    final_results_filename = output_dir + "final_results.txt"
    print("Writing final results: " + final_results_filename)
    final_results_df.to_csv(final_results_filename, sep="\t")

