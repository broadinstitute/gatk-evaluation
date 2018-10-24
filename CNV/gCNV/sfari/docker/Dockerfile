FROM bioconductor/release_base2
WORKDIR /root

RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("rtracklayer", "VariantAnnotation", "stringr"), suppressUpdates=TRUE)'

EXPOSE 80
