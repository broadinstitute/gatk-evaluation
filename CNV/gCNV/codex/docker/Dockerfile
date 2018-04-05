FROM bioconductor/release_base2
WORKDIR /root

RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("CODEX", suppressUpdates=TRUE)'

EXPOSE 80
