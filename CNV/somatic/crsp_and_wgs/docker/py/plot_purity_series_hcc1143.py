import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import os
from argparse import RawDescriptionHelpFormatter, ArgumentParser
from CNV.somatic.crsp_and_wgs.docker.py.clopper_pearson import clopper_pearson
import matplotlib.pyplot as plt

import pandas
from pandas import DataFrame
from pandas import Series

ploidy = 3.7
"""sample_id	purity
SM-74NEG	0
SM-74P2T	10
SM-74P35	30
SM-74P3J	40
SM-74P3K	70
SM-74P3M	60
SM-74P4M	100
SM-74P51	80
SM-74P56	90"""

# purities = [0.0, 0.1, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]
# input_tsvs = ["SM-74NEG.bam.calls.tsv", "SM-74P2T.bam.calls.tsv", "SM-74P35.bam.calls.tsv", "SM-74P3J.bam.calls.tsv",
#               "SM-74P3M.bam.calls.tsv", "SM-74P3K.bam.calls.tsv", "SM-74P51.bam.calls.tsv", "SM-74P56.bam.calls.tsv",
#               "SM-74P4M.bam.calls.tsv"]

gt_cn_column_name = "cn"
gt_cr_column_name = "gt_cr"
guess_cr_column_name = "guess_cr"
# input_guess_cr_column_name = "LOG2_COPY_RATIO_POSTERIOR_50"
input_guess_cr_column_name = "MEAN_LOG2_COPY_RATIO"
is_log2_guess_cr = True


def moreThanOneValue(v):
    if not isinstance(v, str):
        return True
    return v.find("__") == -1

purity_dict = {"SM-74NEG": 0.0,
               "SM-74P2T": 0.1,
               "SM-74P35": 0.3,
               "SM-74P3J": 0.4,
               "SM-74P3M": 0.6,
               "SM-74P3K": 0.7,
               "SM-74P51": 0.8,
               "SM-74P4M": 1.0,
               "SM-74NF5": 1.0,
               "SM-74P56": 0.9}


def find_purity_from_filename(fn):
    for k in purity_dict.keys():
        if fn.find(k) != -1:
            return purity_dict[k]
    return None

def find_sample_from_filename(fn):
    for k in purity_dict.keys():
        if fn.find(k) != -1:
            return k
    return None


def plot_purity_series(output_dir, df_to_plot, plot_title, min_sensitivity, min_precision, min_supported_purity, is_show_line=True):
    h = plt.figure()
    h.hold(True)

    # Useful link: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.errorbar
    if is_show_line:
        linestyle = 'solid'
    else:
        linestyle = 'None'

    es = plt.errorbar(list(df_to_plot.purity), list(df_to_plot.sensitivity),
                      yerr=[list(df_to_plot.sensitivity - df_to_plot.sens_lo),
                            list(df_to_plot.sens_hi - df_to_plot.sensitivity)],
                      barsabove=True, marker='x', mew=2, lw=2, capsize=10, color='b', ls=linestyle)
    ep = plt.errorbar(list(df_to_plot.purity), list(df_to_plot.precision),
                      yerr=[list(df_to_plot.precision - df_to_plot.prec_lo),
                            list(df_to_plot.prec_hi - df_to_plot.precision)],
                      barsabove=True, marker='x', mew=2, lw=2, capsize=10, color='g', ls=linestyle)
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)
    plt.title(plot_title)
    plt.xlabel("purity")
    min_sens = plt.Line2D([0, 100], [min_sensitivity, min_sensitivity], ls=":", lw=2, color='b')
    min_prec = plt.Line2D([0, 100], [min_precision, min_precision], ls=":", lw=2, color='g')
    min_purity = plt.Line2D([min_supported_purity, min_supported_purity], [0, 1.0], ls=":", lw=2, color='r')
    plt.legend([es, ep, min_sens, min_prec, min_purity],
               ['Sensitivity', 'Precision', 'Minimum Sensitivity', 'Minimum Precision', 'Min Supported Purity'],
               loc='lower right')
    h.gca().add_line(min_purity)
    h.gca().add_line(min_sens)
    h.gca().add_line(min_prec)
    plt.savefig(output_dir + 'purity_series_' + plot_title + '.png', dpi=200, bbox_inches='tight')


def is_passing(hi_sens, hi_prec, purity, min_sensitivity, min_precision, min_supported_purity):
    if purity <= min_supported_purity:
        return True
    if hi_sens < min_sensitivity:
        return False
    if hi_prec < min_precision:
        return False
    return True

def run_purity_plotting(input_tsvs, output_dir):
    result_cols = ["sensitivity", "sens_lo", "sens_hi", "sens_N", "precision", "prec_lo", "prec_hi", "prec_N", "purity", "pass"]
    amp_results_df = DataFrame(columns=result_cols)
    del_results_df = DataFrame(columns=result_cols)

    min_sensitivity = 0.85
    min_precision = 0.8
    min_supported_purity = 0.39

    for i, input_tsv in enumerate(input_tsvs):
        purity = find_purity_from_filename(input_tsv)
        print(input_tsv + "  purity: " + str(purity))

        if purity is None:
            print("The file " + input_tsv + " is unrecognized as being a HCC1143T purity file, so it is being skipped.  Please see the src code here if you believe this is an error.")
            continue

        sample = find_sample_from_filename(input_tsv)

        segs_df_tmp = pandas.read_csv(input_tsv, sep="\t", comment="@")

        # Clean up by removing all locations where there was more than one ground truth value for copy number/ratio
        segs_df = segs_df_tmp[segs_df_tmp[gt_cn_column_name].apply(moreThanOneValue)]
        tmp = segs_df[gt_cn_column_name]
        tmp = pandas.to_numeric(tmp, errors='coerce', downcast='integer')
        segs_df[gt_cn_column_name] = tmp

        cr_gt = 1 + (purity * ((segs_df[gt_cn_column_name] / ploidy) - 1))
        cr_gt.rename(gt_cr_column_name, inplace=True)

        if is_log2_guess_cr:
            cr_guess = 2 ** segs_df[input_guess_cr_column_name]
        else:
            cr_guess = segs_df[input_guess_cr_column_name]

        cr_guess.rename(guess_cr_column_name, inplace=True)

        segs_df[gt_cr_column_name] = cr_gt
        segs_df[guess_cr_column_name] = cr_guess
        segs_gt_to_consider = segs_df[~segs_df["CALL"].isnull() & (segs_df["CONTIG"] != "2")]

        ## Amps
        tp = segs_gt_to_consider[(segs_gt_to_consider["CALL"] == "+") & (segs_gt_to_consider[gt_cn_column_name] >= 5)]
        all_gt_amp = segs_gt_to_consider[segs_gt_to_consider[gt_cn_column_name] >= 5]
        sens_amps = float(len(tp)) / float(len(all_gt_amp))
        sens_amps_ci = clopper_pearson(len(tp), len(all_gt_amp))
        sens_amps_N = len(all_gt_amp)

        fp = segs_gt_to_consider[(segs_gt_to_consider["CALL"] == "+") & (segs_gt_to_consider[gt_cn_column_name] <= 4)]
        prec_amps = float(len(tp)) / float(len(tp) + len(fp))
        prec_amps_ci = clopper_pearson(len(tp), (len(tp) + len(fp)))
        prec_amps_N = len(tp) + len(fp)

        amp_result = Series(name=sample, data={result_cols[0]: sens_amps, result_cols[1]: sens_amps_ci[0],
                                               result_cols[2]: sens_amps_ci[1], result_cols[3]: sens_amps_N,
                                               result_cols[4]: prec_amps, result_cols[5]: prec_amps_ci[0],
                                               result_cols[6]: prec_amps_ci[1], result_cols[7]: prec_amps_N,
                                               result_cols[8]: purity,
                                               result_cols[9]: is_passing(sens_amps_ci[1], prec_amps_ci[1], purity, min_sensitivity, min_precision,
                                                          min_supported_purity)})

        amp_results_df = amp_results_df.append(amp_result)
        amp_results_df.sort_values(result_cols[8], inplace=True)

        print("Amp sensitivity: " + str(sens_amps) + "  " + str(sens_amps_ci))
        print("Amp precision: " + str(prec_amps) + "  " + str(prec_amps_ci))

        ## Dels
        tp_del = segs_gt_to_consider[
            (segs_gt_to_consider["CALL"] == "-") & (segs_gt_to_consider[gt_cn_column_name] <= 2)]
        all_gt_del = segs_gt_to_consider[segs_gt_to_consider[gt_cn_column_name] <= 2]
        sens_dels = float(len(tp_del)) / float(len(all_gt_del))
        sens_dels_ci = clopper_pearson(len(tp_del), len(all_gt_del))
        sens_dels_N = len(all_gt_del)

        fp_del = segs_gt_to_consider[
            (segs_gt_to_consider["CALL"] == "-") & (segs_gt_to_consider[gt_cn_column_name] > 2)]
        prec_dels = float(len(tp_del)) / float(len(tp_del) + len(fp_del))
        prec_dels_ci = clopper_pearson(len(tp_del), (len(tp_del) + len(fp_del)))
        prec_dels_N = len(tp_del) + len(fp_del)

        del_result = Series(name=sample, data={result_cols[0]: sens_dels, result_cols[1]: sens_dels_ci[0],
                                               result_cols[2]: sens_dels_ci[1], result_cols[3]: sens_dels_N,
                                               result_cols[4]: prec_dels, result_cols[5]: prec_dels_ci[0],
                                               result_cols[6]: prec_dels_ci[1], result_cols[7]: prec_dels_N,
                                               result_cols[8]: purity,
                                               result_cols[9]: is_passing(sens_dels_ci[1], prec_dels_ci[1], purity,
                                                                          min_sensitivity, min_precision,
                                                                          min_supported_purity)})
        del_results_df = del_results_df.append(del_result)
        del_results_df.sort_values(result_cols[8], inplace=True)

        print("Del sensitivity: " + str(sens_dels) + "  " + str(sens_dels_ci))
        print("Del precision: " + str(prec_dels) + "  " + str(prec_dels_ci))

    if len(amp_results_df) > 0 and len(del_results_df) > 0:
        plot_purity_series(output_dir, amp_results_df, "Amplifications", min_sensitivity, min_precision, min_supported_purity)
        plot_purity_series(output_dir, del_results_df, "Deletions", min_sensitivity, min_precision, min_supported_purity)
    amp_results_df.to_csv(output_dir + "Amplifications_table.tsv", sep="\t")
    del_results_df.to_csv(output_dir + "Deletions_table.tsv", sep="\t")


def parse_options():
    epilog = """ This tool will create plots of the sensitivity and precision.  This script can only be used w/ HCC1143T purity series results
        produced by the sequencing platform (ICE).  Ground truth must be from JaBbA. 
    """
    desc = "Create purity performance for ground truth to GATK CNV comparison"
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("input_files", nargs="+", type=str, help="Files that are the called segments merged with the ground truth.")
    parser.add_argument("-O", "--output_dir", type=str, required=True, help="Output dir.  Will be created if it does not exist.")
    # Process arguments
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_options()
    input_tsvs = args.input_files
    output_dir = args.output_dir

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    try:
        os.makedirs(output_dir)
    except OSError:
        pass
    run_purity_plotting(input_tsvs, output_dir)
