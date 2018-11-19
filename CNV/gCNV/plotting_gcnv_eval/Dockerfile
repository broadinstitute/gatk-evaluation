FROM continuumio/miniconda3 

ADD environment.yml .
RUN conda env create -f environment.yml 

ADD ./plot_evaluation_metrics.py /root/
ADD ./plt_gcnv_eval/ /root/
RUN export PYTHONPATH="${PYTHONPATH}:/root/plt_gcnv_eval/"

EXPOSE 80

CMD [ "source activate gcnv-eval" ] 
