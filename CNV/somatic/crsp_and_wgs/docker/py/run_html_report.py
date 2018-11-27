import logging
import math
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from glob import glob

import django
import pandas
from django.conf import settings
from django.template import Context, Template
from pandas.core.frame import DataFrame
from pandas.core.series import Series

import CliUtils


def parse_options():
    epilog = """ This tool will create a directory and some html (with images) for a final report of performance.  Many image locations are hardcoded.
    
    Formatting of the tables is generated dynamically in this script and then directly imported by the html template.
    """
    desc = "Create little website for evaluation results"
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("--purity-dir", type=str, help="Directory of the purity plots.  Should be a relative path to the output_dir specified below.")
    parser.add_argument("--clinical-dir", type=str, help="Directory of the clinical plots.  Should be a relative path to the output_dir specified below.")
    parser.add_argument("--reproducibility-dir", type=str, help="Directory of the reproducibility plots.  Should be a relative path to the output_dir specified below.")
    parser.add_argument("--bp-concordance-dir", type=str, help="Directory of the WGS bp concordance plots.  Should be a relative path to the output_dir specified below.")
    parser.add_argument("--gatk-docker", type=str, help="GATK docker image used.  This is purely for display in the result.  There is no way to enforce whether results were generated with this image.")
    parser.add_argument("--eval-docker", type=str, help="GATK eval image used.  This is purely for display in the result.  There is no way to enforce whether results were generated with this image.")
    parser.add_argument("--mad-threshold", required=False, type=float, default=0.1, help="Threshold for passing WGS concordance.")
    parser.add_argument("--html_template", type=str, default="html/aggregate_template.html", help="Template to use.")
    parser.add_argument("output_dir", type=str, help="Output dir.  Will be created if it does not exist.")
    parser.add_argument("-l", "--log_filename", type=str, default="html_aggregate.log", help="log filename.  Will be overwritten if it exists.")
    # Process arguments
    args = parser.parse_args()

    return args


def float_to_pretty_str(f):
    # type float -> str
    """
    Standardize how we print floating point numbers.
    :param f: The float to be printed.
    :return: string representation of the float.
    """
    return "%3.3f" % f


def color_row_value_below(df_in, index_name, val, color='red'):
    # type: (DataFrame, str, float, str) -> DataFrame
    """
    Method used to color the cells of a single row of a table based on whether a specific column meets a threshold.
    Returns a format string for each row that can be used in the html.

    :param df_in: DataFrame of values to drive table color determination.
    :param index_name: Column name that is compared to the threshold val.
    :param val: min threshold for a background color change.
    :param color: color to use for the background.
    :return: lists of strings corresponding to each entry in Series s.  Empty string if the value did not meet min
        threshold.

    Inspired by https://stackoverflow.com/questions/47469478/how-to-color-whole-rows-in-python-pandas-based-on-value
    """
    #copy df to new - original data are not changed
    df = df_in.copy()
    #set by condition
    mask = df[index_name] < val
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


def color_row_str_equals(x, index_name, val, color='red'):
    # type: (DataFrame, str, str, str) -> DataFrame
    """
    See :func color_row_value_below:, but val is a string that must match the column value for a
    background color change.

    :param x:
    :param index_name:
    :param val:
    :param color:
    :return:
    """
    #copy df to new - original data are not changed
    df = x.copy()
    #set by condition
    mask = df[index_name].astype(str) == str(val)
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


def color_row_value_above(x, index_name, val, color='red'):
    # type: (DataFrame, str, float, str) -> DataFrame
    """
    See :func color_row_value_below:, but val is the max threshold for a background color change.

    :param x:
    :param index_name:
    :param val:
    :param color:
    :return:
    """
    #copy df to new - original data are not changed
    df = x.copy()
    #set by condition
    mask = df[index_name] > val
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


def render_purity_summary_html(purity_summary_df):
    # type: (DataFrame) -> str
    """
    Create a html block that formats the table for the purity summary.  Colors any row where the "pass" field is "False"

    :param purity_summary_df: DataFrame of purity results.
    :return: string of html.
    """
    return purity_summary_df.style.set_table_styles([
        {'selector': 'td', 'props': [
            ('padding', '4px'),
            ('text-align', 'center'), ("font-size", "11px")]},
        {'selector': 'th', 'props': [
            ('padding', '4px'),
            ('text-align', 'left'), ("font-size", "11px")]},
        {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#ddd')]}
        ]).apply(color_row_str_equals, axis=None, index_name='pass', val="False", color='#f99').render()


def render_reproducibility_summary_html(reproducibility_summary):
    # type: (DataFrame) -> str
    """
    Create a html block that formats the table for the reproducibility summary.  Colors any row where the "isPass" field is "False"

    :param reproducibility_summary: DataFrame of purity results.
    :return: string of html.
    """
    return reproducibility_summary.style.set_table_styles([
        {'selector': 'td', 'props': [
            ('padding', '4px'),
            ('text-align', 'center'), ("font-size", "12px")]},
        {'selector': 'th', 'props': [
            ('padding', '4px'),
            ('text-align', 'left'), ("font-size", "12px")]},
        {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#ddd')]}
        ]).apply(color_row_str_equals, axis=None, index_name='isPass', val="False", color='#f99').render()


def main():

    args = parse_options()
    log_filename = args.log_filename
    gatk_docker = args.gatk_docker
    eval_docker = args.eval_docker
    purity_plots_dir = args.purity_dir
    reproducibility_plots_dir = args.reproducibility_dir
    clinical_plots_dir = args.clinical_dir
    bp_concordance_dir = args.bp_concordance_dir
    html_template = args.html_template
    output_dir = args.output_dir
    CliUtils.setup_logging(log_filename, is_verbose_stdout=True)
    mae_threshold = args.mad_threshold
    logging.info(str(args))

    # Setup the django template engine.
    settings.configure(DEBUG=True, TEMPLATES=[{'BACKEND': 'django.template.backends.django.DjangoTemplates'}])
    django.setup()

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    clinical_amp_plots = glob(clinical_plots_dir + "/clinical_Amp*target*.png")
    clinical_del_plots = glob(clinical_plots_dir + "/clinical_Del*target*.png")
    reproducibility_plots = glob(reproducibility_plots_dir + "/reproducibility_*.png")
    bp_concordance_filtered_plots = sorted(glob(bp_concordance_dir + "/*bp_concordance*bp.png"))
    bp_concordance_maf_plots = sorted(glob(bp_concordance_dir + "/*bp_maf_concordance*bp.png"))
    bp_concordance_input_files = [".".join(os.path.basename(p).split(".")[0:2]) for p in bp_concordance_filtered_plots]
    bp_concordance_raw_data_files = glob(bp_concordance_dir + "/*.seg_comparison_edit.txt")

    # https://pandas.pydata.org/pandas-docs/stable/style.html

    ### WGS bp Concordance table
    bp_concordance_summary = pandas.read_csv(bp_concordance_dir + "/final_results.txt", sep="\t", index_col=0)
    bp_concordance_summary = bp_concordance_summary.drop(['samplename'], axis=1)
    bp_concordance_summary.sort_index(inplace=True)

    ### Reproducibility
    reproducibility_call_conf_mat = pandas.read_csv(reproducibility_plots_dir + "/call_conf_mat.tsv",
                                                                   sep="\t", index_col=0)
    reproducibility_call_concordance_summary = pandas.read_csv(reproducibility_plots_dir + "/call_concordance.tsv",
                                                                  sep="\t", index_col=0)

    ### Purity Table
    purity_amp_summary = pandas.read_csv(purity_plots_dir + "/Amplifications_table.tsv", sep="\t", index_col=0)
    purity_del_summary = pandas.read_csv(purity_plots_dir + "/Deletions_table.tsv", sep="\t", index_col=0)

    ## Construct the context dictionary
    context_dict = {"gatk_docker": gatk_docker, "eval_docker": eval_docker,
                    "clinical_amp_plots": clinical_amp_plots,
                    "clinical_del_plots": clinical_del_plots,
                    "zip_clinical_plots": zip(clinical_amp_plots, clinical_del_plots),
                    "bp_concordance_plots": bp_concordance_filtered_plots,
                    "bp_concordance_mad_germline_filtered_threshold": str(mae_threshold),
                    "bp_concordance_summary": bp_concordance_summary.style.set_table_styles([
                        {'selector': 'td', 'props': [
                                     ('padding', '4px'),
                                     ('text-align', 'center'), ("font-size", "11px")]},
                        {'selector': 'th', 'props': [
                            ('padding', '4px'),
                            ('text-align', 'left'),("font-size", "11px")]},
                        {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#ddd')]}
                        ])
                        .apply(color_row_value_above, axis=None, index_name='mae', val=mae_threshold, color='#f99')
                        .apply(color_row_value_below, axis=None, index_name='mae', val=1e-16, color='red')
                        .render(),
                    "bp_concordance_raw_data_files": bp_concordance_raw_data_files,
                    "reproducibility_plots": reproducibility_plots,
                    "reproducibility_call_conf_mat": reproducibility_call_conf_mat.to_html(),
                    "reproducibility_call_concordance_summary": render_reproducibility_summary_html(reproducibility_call_concordance_summary),
                    "reproducibility_link": reproducibility_plots_dir + "/sample_df_targets_only.tsv",
                    "zip_zip_reproducibility_names": zip(bp_concordance_input_files, bp_concordance_filtered_plots, bp_concordance_raw_data_files, bp_concordance_maf_plots),
                    "purity_amp_plot": purity_plots_dir + "/purity_series_Amplifications.png",
                    "purity_del_plot": purity_plots_dir + "/purity_series_Deletions.png",
                    "purity_amp_summary": render_purity_summary_html(purity_amp_summary),
                    "purity_del_summary": render_purity_summary_html(purity_del_summary)
                    }

    ## Render and write the actual html
    html_template_str = file(html_template).read()
    out_fp = file(output_dir + "report.html", 'w')
    out_fp.write(Template(html_template_str).render(Context(context_dict)))

    out_fp.close()

if __name__ == "__main__":
    main()