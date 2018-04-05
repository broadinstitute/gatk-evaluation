#!/usr/bin/env Rscript

#This script collects coverage at targets on a specified contig from a BAM.
#It requires arguments in the following order:
#
#   args[1] - BAM file
#   args[2] - BED file containing targets
#   args[3] - sample name (can differ from SM tags in BAM)
#   args[4] - contig (e.g., 22)
#   args[5] - minimum mapping quality

suppressMessages(library(CODEX))

#parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Exactly 5 arguments must be supplied")
}
bam_file <- args[1]
target_bed_file <- args[2]
sample_name <- as.matrix(args[3])
contig <- args[4]
min_mapping_quality <- as.numeric(args[5])
output_dir <- args[6]

#collect coverage
bambed_obj <- getbambed(bamdir=bam_file, bedFile=target_bed_file, sampname=sample_name, projectname="coverage", contig)
coverage_obj <- getcoverage(bambed_obj, mapqthres=min_mapping_quality)
targets <- bambed_obj$ref
coverage_df <- data.frame(sample_name, bambed_obj$chr, start(targets), end(targets), coverage_obj$Y)
colnames(coverage_df) <- c("SAMPLE", "CONTIG", "START", "END", "Y")
write.table(coverage_df, file=paste(sample_name, "_", contig, "_coverage.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
