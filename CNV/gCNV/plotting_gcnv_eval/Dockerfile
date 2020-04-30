FROM continuumio/miniconda3 

RUN apt-get -y update
RUN apt-get -y install --reinstall build-essential 
RUN apt-get install libz-dev

ADD environment.yml .
RUN conda env create -f environment.yml 

ADD ./plot_evaluation_metrics.py /root/
ADD ./plt_gcnv_eval/ /root/
RUN export PYTHONPATH="${PYTHONPATH}:/root/plt_gcnv_eval/"

RUN echo "source activate gcnv-eval" > ~/.bashrc
ENV PATH /opt/conda/envs/gcnv-eval/bin:$PATH

EXPOSE 80 
