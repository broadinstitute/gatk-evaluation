
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from os.path import basename

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas
from pandas import DataFrame
from pandas import Series

SMALL_SEG_THRESHOLD = 7

from clopper_pearson import clopper_pearson

GUESS_CALL_COLUMN = "CALL"
GT_CALL_COLUMN = "Segment_Call"
GT_SEG_MEAN_COLUMN = "Segment_Mean"

CALL_TO_ENGLISH = {"+": "Amplification", "-": "Deletion"}


def parse_options():
    epilog = """ This tool will sensitivity plots of the ground truth. This tool assumes that CombineSegmentBreakpoints has already been run between the candidate file and its ground truth.
    """
    desc = ""
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('-I', '--input_files', action='append', help='Input seg files that have already been merged with gt.', required=True)
    parser.add_argument("-O", "--output_dir", type=str, required=True, help="Output dir.  Will be created if it does not exist.")

    # Process arguments
    args_parsed = parser.parse_args()

    return args_parsed


def create_sensitivity_df(list_of_lists_tp_fn, sample_names):
    # type: (list[list[int]], list[str]) -> DataFrame
    """

    :param list_of_lists_tp_fn:
    :param sample_names:
    :return:
    """
    if len(sample_names) != len(list_of_lists_tp_fn):
        raise Exception("list of lists must match number of sample names")

    final_df = DataFrame(index=sample_names, columns=["sensitivity", "sens_lo", "sens_hi", "sens_N"])

    for i, p in enumerate(list_of_lists_tp_fn):
        n = p[0] + p[1]
        if n == 0:
            return
        sensitivity = float(p[0]) / n
        sens_binom_interval = clopper_pearson(p[0], n)
        final_df.ix[i] = Series(
            {"sensitivity": sensitivity, "sens_lo": sens_binom_interval[0], "sens_hi": sens_binom_interval[1],
             "sens_N": n})
    return final_df


def plot_clinical_sensitivity(sensitivity_df, output_dir, plot_title, min_sens_val=0.8,
                              file_prefix='clinical_'):
    # type: (DataFrame, str, str, float, str) -> None
    """
    Create a clinical sensitivity png file.  The filename will be:
    output_dir + '/' + file_prefix + plot_title + '.png'

    :param sensitivity_df: DataFrame to plot.  Must have columns: sensitivity, sens_lo, sens_hi, sens_N
    :param output_dir: directory to drop the plots.
    :param plot_title: Title string
    :param min_sens_val: Where to plot a cutoff minimum sensitivity
    :param file_prefix: prefix for the output file.
    :return:
    """
    h = plt.figure()
    h.hold(True)
    xvals = range(0, len(sensitivity_df.index))

    if len(sensitivity_df) == 0:
        plt.text(0.25, 0.25, "There are no events fitting this criteria")
    else:

        es = plt.errorbar(xvals, list(sensitivity_df.sensitivity),
                          yerr=[list(sensitivity_df.sensitivity - sensitivity_df.sens_lo),
                                list(sensitivity_df.sens_hi - sensitivity_df.sensitivity)],
                          barsabove=True, marker='x', mew=2, lw=4, capsize=10, color='b', ls='None')
        min_sens = plt.Line2D([-0.25, 100], [min_sens_val, min_sens_val], ls=":", lw=2, color='b')
        h.gca().add_line(min_sens)
        plt.legend([es, min_sens],
                   ['Sensitivity', 'Minimum Sensitivity'],
                   loc='lower right')

    xlabels = [x + " (N=" + str(int(sensitivity_df.loc[x].sens_N)) + ")" for x
               in sensitivity_df.index]

    plt.xticks(xvals, xlabels, rotation=45, ha='right', size='x-small')
    plt.ylim(0, 1.0)

    plt.title(plot_title)
    plt.savefig(output_dir + '/' + file_prefix + plot_title + '.png', dpi=200, bbox_inches='tight')

if __name__ == '__main__':
    args = parse_options()
    output_dir = args.output_dir
    input_files = args.input_files

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    #######################

    dfs = [pandas.read_csv(f, sep="\t", comment="@") for f in input_files]
    sample_names = [basename(f) for f in input_files]
    for gt_call_of_interest in ["+", "-"]:

        target_final_df = DataFrame(columns=["sensitivity", "sens_lo", "sens_hi", "sens_N"])
        seg_final_df = DataFrame(columns=["sensitivity", "sens_lo", "sens_hi", "sens_N"])
        small_seg_final_df = DataFrame(columns=["sensitivity", "sens_lo", "sens_hi", "sens_N"])

        for j, df in enumerate(dfs):

            # Get sample name
            sample_name = sample_names[j]

            # relevant is where we have a ground truth call and we have at least one target.
            #  Please note that we generally do not count ground truth with only one overlapping target.
            #  But that pruning would be done later.
            relevant_df = df[(~df["Segment_Call"].isnull()) & (~df["NAME"].isnull())]

            # keys must be dependent only on GT columns
            relevant_df["keys"] = relevant_df["Segment_Call"] + "__" + relevant_df["Segment_Mean"].astype(str)

            gt_seg_value_counts = relevant_df["keys"].value_counts()

            # [TP, FN]
            small_totals_by_target = [0, 0]
            totals_by_target = [0, 0]
            small_totals_by_seg = [0, 0]
            totals_by_seg = [0, 0]

            for i, k in enumerate(gt_seg_value_counts):

                # Only consider segments with at least 2 targets in the gt
                if k < 2:
                    continue

                overlapping_segs = relevant_df[(relevant_df["keys"] == gt_seg_value_counts.index[i]) &
                                               (relevant_df[GT_CALL_COLUMN] == gt_call_of_interest)]
                is_small = (k <= SMALL_SEG_THRESHOLD)

                tp_by_target = (overlapping_segs[GT_CALL_COLUMN] == overlapping_segs[GUESS_CALL_COLUMN]).sum()
                fn_by_target = (overlapping_segs[GT_CALL_COLUMN] != overlapping_segs[GUESS_CALL_COLUMN]).sum()

                totals_by_target[0] = totals_by_target[0] + tp_by_target
                totals_by_target[1] = totals_by_target[1] + fn_by_target

                tp_by_seg = 1 if (tp_by_target > 0) else 0
                fn_by_seg = 1 if (fn_by_target > 0) else 0

                totals_by_seg[0] = totals_by_seg[0] + tp_by_seg
                totals_by_seg[1] = totals_by_seg[1] + fn_by_seg

                if is_small:
                    small_totals_by_target[0] = small_totals_by_target[0] + tp_by_target
                    small_totals_by_target[1] = small_totals_by_target[1] + fn_by_target
                    small_totals_by_seg[0] = small_totals_by_seg[0] + tp_by_seg
                    small_totals_by_seg[1] = small_totals_by_seg[1] + fn_by_seg

            target_final_df = pandas.concat([target_final_df, create_sensitivity_df([totals_by_target], sample_names=[sample_name])])
            seg_final_df = pandas.concat([seg_final_df, create_sensitivity_df([totals_by_seg], sample_names=[sample_name])])
            small_seg_final_df = pandas.concat([small_seg_final_df, create_sensitivity_df([small_totals_by_seg], sample_names=[sample_name])])
        call_of_interest_english = CALL_TO_ENGLISH[gt_call_of_interest]
        plot_clinical_sensitivity(target_final_df, output_dir, call_of_interest_english + " sensitivity by target")
        plot_clinical_sensitivity(seg_final_df, output_dir, call_of_interest_english + " sensitivity by segment")
        plot_clinical_sensitivity(small_seg_final_df, output_dir, call_of_interest_english + " sensitivity by small segments ( max " + str(SMALL_SEG_THRESHOLD) + " targets)")