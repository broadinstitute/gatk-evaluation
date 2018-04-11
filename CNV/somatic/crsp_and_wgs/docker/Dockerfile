FROM ubuntu:16.04

RUN apt-get update && \
    apt-get install -y python-pip curl libxft-dev libfreetype6 libfreetype6-dev libblas-dev liblapack-dev liblapacke-dev gfortran && \
    pip install --upgrade pip
RUN pip install numpy==1.11.1
RUN pip install pandas==0.20.1
RUN pip install scipy==0.19.1
RUN pip install matplotlib==2.0.2
RUN pip install scikit-learn==0.19.1

# Get OpenJDK8 on here.
RUN apt-get update && \
apt-get install -y openjdk-8-jdk

# For HTML rendering, install django
RUN pip install Jinja2
RUN pip install Django==1.11.11

COPY consensus.20170217.purity.ploidy.txt /root/consensus.20170217.purity.ploidy.txt
COPY UUID-BARCODE.tsv /root/UUID-BARCODE.tsv
COPY plot_purity_series_hcc1143.py /root/plot_purity_series_hcc1143.py
COPY plot_bp_concordance_pcawg_pilot.py /root/plot_bp_concordance_pcawg_pilot.py
COPY run_plot_reproducibility.py /root/run_plot_reproducibility.py
COPY final_centromere_hg19.seg /root/final_centromere_hg19.seg
COPY clopper_pearson.py /root/clopper_pearson.py
COPY plot_clinical_sensitivity.py /root/plot_clinical_sensitivity.py
COPY run_html_report.py /root/run_html_report.py
COPY style.css /root/style.css
COPY CliUtils.py /root/CliUtils.py
COPY aggregate_template.html /root/aggregate_template.html
