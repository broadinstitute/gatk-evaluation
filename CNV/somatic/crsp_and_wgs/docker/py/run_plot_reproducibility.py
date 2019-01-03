import os
# Force matplotlib to not use any Xwindows backend.
from math import ceil, sqrt

import matplotlib
from pandas.core.frame import DataFrame
from pandas.core.series import Series

matplotlib.use('Agg')
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import numpy
import pandas

# When the parameter below is True, the validation assumes the output of CallModeledSegments. If it is False, then it assumes the output of CallCopyRatioSegments.
ASSUME_CALL_MODELED_SEGMENTS_OUTPUT = True

# shard 5 and shard 10
def plot_reproducibility(dfj, output_dir, sample1, sample2, threePThresh, title, step_in_plot, hi_lim=40.0,
                         file_prefix='reproducibility_', label="Copy Number", col1="guess_cn_1", col2="guess_cn_2",
                         is_plot_3p=True):
    # type: (DataFrame, str, str, str, float, str, float, float, str, str, str, str, bool) -> None
    """

    Create a log 2D histogram from a DataFrame (dfj) and two columns (col1, col2), which will be x, y, respectively.
    Writes a png to the output directory:  output_dir + '/' + file_prefix + "Reproducibility.png"

    :param dfj: DataFrame with the data to be counted in the 2D histogram.
    :param output_dir: Name of the output directory.
    :param sample1: Name associated with the first column.  This is usually a sample name.  This will be prefix of the x-axis label
    :param sample2: Name associated with the second column.  This is usually a sample name.  This will be prefix of the y-axis label
    :param threePThresh: Value for the 3*ploidy, since this plots in copy number space.
    :param title: Title in the plot.
    :param step_in_plot: Steps for x and y axis.  I.e. bin sizes.
    :param hi_lim: Maximum copy number to plot.
    :param file_prefix: Starting string for the png filename.
    :param label: axis label that will be appended to sample1 and sample2 strings.
    :param col1: column name in dfj to use for the x-axis
    :param col2: column name in dfj to use for the y-axis
    :param is_plot_3p: Whether to plot the three * ploidy line.
    :return:
    """

    h = plt.figure()

    plt.hist2d(dfj[col1], dfj[col2],
               bins=[numpy.arange(0, hi_lim, step_in_plot), numpy.arange(0, hi_lim, step_in_plot)], cmin=0.0001,
               cmap='viridis_r', norm=matplotlib.colors.LogNorm())
    plt.xlabel(sample1 + " " + label)
    plt.ylabel(sample2 + " " + label)
    plt.colorbar()
    plt.title(title)
    h.gca().add_line(plt.Line2D([0, hi_lim], [0, hi_lim], ls=":", lw=1, color='g'))

    if is_plot_3p:
        h.gca().add_line(plt.Line2D([threePThresh, threePThresh], [0, threePThresh], ls=":", lw=2, color='m'))
        h.gca().add_line(plt.Line2D([0, threePThresh], [threePThresh, threePThresh], ls=":", lw=2, color='m'))

    plt.savefig(output_dir + '/' + file_prefix + "Reproducibility.png", dpi=200, bbox_inches='tight')


def parse_options():
    epilog = """ This tool will create plots for HCC1143T (100% purity) reproducibility.  This can be used for other 
    samples, but you will probably want to specify the `--cr` parameter. 
    """
    desc = "Create reproducibility performance for GATK CNV comparison"
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("input_file", type=str, help="File of first sample.  From called seg file combined with other sample called seg file -- no ground truth involved.")
    parser.add_argument("sample1_name", type=str, help="Name of first sample (x-axis in plots).")
    parser.add_argument("sample2_name", type=str, help="Name of second sample (y-axis in plots).")
    parser.add_argument("output_dir", type=str, help="Output dir.  Will be created if it does not exist.")
    parser.add_argument("ploidy", type=str, help="ploidy.  Ignored if --cr is specified.", default=3.7)
    parser.add_argument("-l", "--log_filename", type=str, default="run_plot_reproducibility.log", help="log filename.  Will be overwritten if it exists.")
    parser.add_argument("--cr", action='store_true', required=False, help="Plot as if input was copy ratio.  A dummy ploidy will be used in some places, overriding what is specified.")
    # Process arguments
    args = parser.parse_args()

    return args


def main():

    args = parse_options()
    log_filename = args.log_filename
    input_file = args.input_file
    ploidy = float(args.ploidy)
    output_dir = args.output_dir
    sample1_name = args.sample1_name
    sample2_name = args.sample2_name

    is_cr = args.cr
    if (is_cr):
        ploidy = 1.0

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    sample_df = pandas.read_csv(input_file, sep="\t", comment="@")

    # Back out a copy ratio (not log'd) from the output of the CNV tool.
    if "LOG2_COPY_RATIO_POSTERIOR_50_1" in sample_df.columns:
        sample_df["guess_cr_1"] = 2 ** sample_df["LOG2_COPY_RATIO_POSTERIOR_50_1"]
        sample_df["guess_cr_2"] = 2 ** sample_df["LOG2_COPY_RATIO_POSTERIOR_50_2"]
    else:
        sample_df["guess_cr_1"] = 2 ** sample_df["MEAN_LOG2_COPY_RATIO_1"]
        sample_df["guess_cr_2"] = 2 ** sample_df["MEAN_LOG2_COPY_RATIO_2"]

    # Back out a copy number.
    #  Important note:  This assumes a purity of 100%
    sample_df["guess_cn_1"] = sample_df["guess_cr_1"] * ploidy
    sample_df["guess_cn_2"] = sample_df["guess_cr_2"] * ploidy

    sample_df_targets_only = sample_df[(~sample_df["LOG2_COPY_RATIO"].isnull()) &
                                       (~sample_df["guess_cn_1"].isnull()) &
                                       (~sample_df["guess_cn_2"].isnull())]

    # Calculate all information needed to do the reproducibility and plot it.
    rmse = pow(sample_df_targets_only.guess_cn_1 - sample_df_targets_only.guess_cn_2, 2)
    sample_df_targets_only["rmse"] = rmse

    # Uncomment when minor allele fraction is in the call file
    # rmse_af = pow(sample_df_targets_only.MINOR_ALLELE_FRACTION_POSTERIOR_50_1 - sample_df_targets_only.MINOR_ALLELE_FRACTION_POSTERIOR_50_2, 2)
    # sample_df_targets_only["rmse_af"] = rmse_af
    # rmse_af = sqrt(sum(rmse_af) / len(rmse_af))

    rmse = sqrt(sum(rmse) / len(rmse))
    low_indices = (sample_df_targets_only.guess_cn_1 < (3.0 * ploidy)) & (sample_df_targets_only.guess_cn_2 < (3.0 * ploidy))
    rmse_low = pow(sample_df_targets_only.guess_cn_1[low_indices] - sample_df_targets_only.guess_cn_2[low_indices], 2)
    rmse_low = sqrt(sum(rmse_low) / len(rmse_low))

    # Write the target level information to a table on disk.
    sample_df_targets_only.to_csv(output_dir + "sample_df_targets_only.tsv", sep="\t")

    # Plot the reproducibility
    base_step_in_plot = 0.5
    if ploidy <= 2.0:
        base_step_in_plot = 0.2

    if not is_cr:
        plot_reproducibility(sample_df_targets_only, output_dir, sample1_name, sample2_name, 3.0 * ploidy, step_in_plot=base_step_in_plot,
                             title="RMSE by Target Seg Mean: (3*ploidy/all) " + ("%0.4f" % rmse_low) + "/" + ("%0.4f" % rmse))
        plot_reproducibility(sample_df_targets_only, output_dir, sample1_name, sample2_name, 3.0 * ploidy, hi_lim=ceil(3.0 * ploidy),
                             step_in_plot=base_step_in_plot/2.0, file_prefix='reproducibility_zoom_',
                             title="RMSE by Target Seg Mean: (3*ploidy/all) " + ("%0.4f" % rmse_low) + "/" + ("%0.4f" % rmse))
    else:
        plot_reproducibility(sample_df_targets_only, output_dir, sample1_name, sample2_name, 3.0 * ploidy,
                             hi_lim=ceil(3.0 * ploidy),
                             step_in_plot=base_step_in_plot / 2.0, file_prefix='reproducibility_zoom_',
                             title="RMSE by Target Seg Mean: " + (
                             "%0.4f" % rmse), is_plot_3p=False, label="Copy Ratio")

    # Uncomment when minor allele fraction is in the call file
    # plot_reproducibility(sample_df_targets_only, rmse_af, rmse_af, output_dir, sample1_name, sample2_name, 3.0 * ploidy, hi_lim=0.5,
    #                      step_in_plot=0.025, file_prefix='reproducibility_af_', label="Minor Allelic Fraction",
    #                      title="RMSE AF by Target Seg Mean: " + ("%0.4f" % rmse_af), is_plot_3p=False,
    #                      col1="MINOR_ALLELE_FRACTION_POSTERIOR_50_1",
    #                      col2 = "MINOR_ALLELE_FRACTION_POSTERIOR_50_2")

    # Create and write a confusion matrix for the caller.
    conf_mat = DataFrame(index=['amp','none','del'], columns=['amp', 'none', 'del'])
    conf_mat.ix['amp']['amp'] = sum((sample_df_targets_only.CALL_1 == sample_df_targets_only.CALL_2) & (sample_df_targets_only.CALL_1 == "+"))
    conf_mat.ix['none']['none'] = sum((sample_df_targets_only.CALL_1 == sample_df_targets_only.CALL_2) & (sample_df_targets_only.CALL_1 == "0"))
    conf_mat.ix['del']['del'] = sum((sample_df_targets_only.CALL_1 == sample_df_targets_only.CALL_2) & (sample_df_targets_only.CALL_1 == "-"))
    conf_mat.ix['amp']['del'] = sum((sample_df_targets_only.CALL_1 == "+") & (sample_df_targets_only.CALL_2 == "-"))
    conf_mat.ix['none']['del'] = sum((sample_df_targets_only.CALL_1 == "0") & (sample_df_targets_only.CALL_2 == "-"))
    conf_mat.ix['none']['amp'] = sum((sample_df_targets_only.CALL_1 == "0") & (sample_df_targets_only.CALL_2 == "+"))
    conf_mat.ix['del']['amp'] = sum((sample_df_targets_only.CALL_1 == "-") & (sample_df_targets_only.CALL_2 == "+"))
    conf_mat.ix['del']['none'] = sum((sample_df_targets_only.CALL_1 == "-") & (sample_df_targets_only.CALL_2 == "0"))
    conf_mat.ix['amp']['none'] = sum((sample_df_targets_only.CALL_1 == "+") & (sample_df_targets_only.CALL_2 == "0"))

    conf_mat.to_csv(output_dir + "call_conf_mat.tsv", sep="\t")
    call_concordance = (conf_mat.ix['amp']['amp'] + conf_mat.ix['none']['none'] + conf_mat.ix['del']['del']) / conf_mat.sum().sum()

    # Create and write a performance summary.
    #sample1 	sample2 	TargetsMatch? 	CN_RMSE 	Call_Concordance 	CN<11 RMSE 	isPass 	RMSE_CN<3p_Thold 	Call_Concordance_Thold
    #SM-74NF5-100 	SM-74P4M-100 	1 	0.1598 	0.9762 	0.1414 	YES 0.1500 	0.9700
    call_concordance_df = DataFrame(columns=["sample1", "sample2", "CN_RMSE", "Call_Concordance", "CN<3p RMSE", "isPass", "RMSE_CN<3p_Thold", "Call_Concordance_Thold"])
    info = {c: "" for c in call_concordance_df.columns}

    info["sample1"] = sample1_name
    info["sample2"] = sample2_name
    info["CN_RMSE"] = rmse
    info["Call_Concordance"] = call_concordance
    info["CN<3p RMSE"] = rmse_low

    info["RMSE_CN<3p_Thold"] = 0.15
    info["Call_Concordance_Thold"] = 0.97
    info["isPass"] = str(call_concordance > info["Call_Concordance_Thold"] and rmse_low < info["RMSE_CN<3p_Thold"])

    call_concordance_df.ix["Results"] = Series(info)
    call_concordance_df.to_csv(output_dir + "call_concordance.tsv", sep="\t")

if __name__ == "__main__":
    main()
