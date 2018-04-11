import logging
import math
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from glob import glob

import django
import pandas
from django.conf import settings
from django.template import Context, Template

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
    parser.add_argument("--gatk-docker", type=str, help="GATK docker image used.  This is purely for display in the result.  Tehre is no way to enforce whether results were generated with this image.")
    parser.add_argument("--eval-docker", type=str, help="GATK eval image used.  This is purely for display in the result.  Tehre is no way to enforce whether results were generated with this image.")
    parser.add_argument("--mad-threshold", required=False, type=float, default=0.1, help="Threshold for passing WGS concordance.")
    parser.add_argument("--html_template", type=str, default="aggregate_template.html", help="Template to use.")
    parser.add_argument("output_dir", type=str, help="Output dir.  Will be created if it does not exist.")
    parser.add_argument("-l", "--log_filename", type=str, default="html_aggregate.log", help="log filename.  Will be overwritten if it exists.")
    # Process arguments
    args = parser.parse_args()

    return args


def float_to_pretty_str(f):
    return "%3.3f" % f

def color_values_above(s, val, color='red'):
    return ['background-color: ' + color if (si > val or math.isnan(si)) else '' for si in s]


def color_values_below(s, val, color='red'):
    return ['background-color: ' + color if (si < val or math.isnan(si)) else '' for si in s]


def color_row_value_below(x, index_name, val, color='red'):
    #copy df to new - original data are not changed
    df = x.copy()
    #set by condition
    mask = df[index_name] < val
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


def color_row_str_equals(x, index_name, val, color='red'):
    #copy df to new - original data are not changed
    df = x.copy()
    #set by condition
    mask = df[index_name].astype(str) == str(val)
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


#TODO: Get citation for this method, so that you are not rude.
def color_row_value_above(x, index_name, val, color='red'):
    #copy df to new - original data are not changed
    df = x.copy()
    #set by condition
    mask = df[index_name] > val
    df.loc[mask, :] = 'background-color: ' + color
    df.loc[~mask,:] = 'background-color: ""'
    return df


def render_purity_summary_html(purity_summary_df):
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
    mad_threshold = args.mad_threshold
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
    bp_concordance_filtered_plots = sorted(glob(bp_concordance_dir + "/*seg*_Filtered_*bp.png"))
    bp_concordance_unfiltered_plots = sorted(glob(bp_concordance_dir + "/*seg*_Unfiltered_*bp.png"))
    bp_concordance_input_files = [".".join(os.path.basename(p).split(".")[0:2]) for p in bp_concordance_filtered_plots]
    bp_concordance_raw_data_files_germline_removed = glob(bp_concordance_dir + "/*.seg_comparison_edit_germline_removed.txt")
    bp_concordance_raw_data_files = glob(bp_concordance_dir + "/*.seg_comparison_edit.txt")
    zip_bp_concordance_plots = zip(bp_concordance_filtered_plots, bp_concordance_unfiltered_plots)

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
                    "zip_bp_concordance_plots": zip_bp_concordance_plots,
                    "bp_concordance_mad_germline_filtered_threshold": str(mad_threshold),
                    "bp_concordance_summary": bp_concordance_summary.style.set_table_styles([
                        {'selector': 'td', 'props': [
                                     ('padding', '4px'),
                                     ('text-align', 'center'), ("font-size", "11px")]},
                        {'selector': 'th', 'props': [
                            ('padding', '4px'),
                            ('text-align', 'left'),("font-size", "11px")]},
                        {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#ddd')]}
                        ])
                        .apply(color_row_value_above, axis=None, index_name='mad_no_germline', val=mad_threshold, color='#f99')
                        .apply(color_row_value_below, axis=None, index_name='mad_no_germline', val=1e-16, color='red')
                        .render(),
                    "bp_concordance_raw_data_files_germline_removed": bp_concordance_raw_data_files_germline_removed,
                    "bp_concordance_raw_data_files": bp_concordance_raw_data_files,
                    "reproducibility_plots": reproducibility_plots,
                    "reproducibility_call_conf_mat": reproducibility_call_conf_mat.to_html(),
                    "reproducibility_call_concordance_summary": render_reproducibility_summary_html(reproducibility_call_concordance_summary),
                    "reproducibility_link": reproducibility_plots_dir + "/sample_df_targets_only.tsv",
                    "zip_zip_reproducibility_names": zip(bp_concordance_input_files, zip_bp_concordance_plots, bp_concordance_raw_data_files_germline_removed, bp_concordance_raw_data_files),
                    "purity_amp_plot": purity_plots_dir + "/purity_series_Amplifications.png",
                    "purity_del_plot": purity_plots_dir + "/purity_series_Deletions.png",
                    "purity_amp_summary": render_purity_summary_html(purity_amp_summary),
                    "purity_del_summary": render_purity_summary_html(purity_del_summary)
                    }

    html_template_str = file(html_template).read()
    out_fp = file(output_dir + "report.html", 'w')
    out_fp.write(Template(html_template_str).render(Context(context_dict)))

    out_fp.close()

if __name__ == "__main__":
    main()