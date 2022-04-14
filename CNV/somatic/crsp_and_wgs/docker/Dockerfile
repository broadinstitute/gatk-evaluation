FROM ubuntu:20.04

RUN apt-get update
RUN apt-get install -y python3-pip
#curl libxft-dev libfreetype6 libfreetype6-dev libblas-dev liblapack-dev liblapacke-dev gfortran
RUN pip install --upgrade pip
RUN pip install numpy==1.20.3
RUN pip install pandas==1.3.2
RUN pip install scipy==1.6.2
RUN pip install matplotlib==3.4.2
RUN pip install scikit-learn==0.24.2

COPY py/* /root/

# Get OpenJDK8 on here.
#RUN apt-get update && \
#apt-get install -y openjdk-8-jdk

# For HTML rendering, install django (and Jinja)
#RUN pip install Jinja2
#RUN pip install Django==1.11.11

#COPY data/consensus.20170217.purity.ploidy.txt /root/consensus.20170217.purity.ploidy.txt
#COPY data/UUID-BARCODE.tsv /root/UUID-BARCODE.tsv
#COPY data/final_centromere_hg19.seg /root/final_centromere_hg19.seg

#COPY py/plot_purity_series_hcc1143.py /root/plot_purity_series_hcc1143.py
#COPY py/plot_bp_concordance_pcawg_pilot.py /root/plot_bp_concordance_pcawg_pilot.py
#COPY py/run_plot_reproducibility.py /root/run_plot_reproducibility.py
#COPY py/clopper_pearson.py /root/clopper_pearson.py
#COPY py/plot_clinical_sensitivity.py /root/plot_clinical_sensitivity.py
#COPY py/run_html_report.py /root/run_html_report.py
#COPY py/CliUtils.py /root/CliUtils.py

#COPY html/style.css /root/style.css
#COPY html/aggregate_template.html /root/aggregate_template.html
