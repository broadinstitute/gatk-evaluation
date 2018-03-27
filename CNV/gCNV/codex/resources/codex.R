#!/usr/bin/env Rscript

#This script requires arguments in the following order:
#
#   args[1] - file containing list of coverage TSVs 
#             (produced by codex-coverage.R, all targets must be identical)
#   args[2] - contig (e.g., 22)
#   args[3] - cohort name (e.g., TCGA_samples)
#   args[4] - output directory
#   args[5] - per-target QC minimum median coverage 
#   args[6] - per-target QC maximum median coverage
#   args[7] - per-target QC minimum length
#   args[8] - per-target QC maximum length
#   args[9] - per-target QC minimum GC
#   args[10] - per-target QC maximum GC
#   args[11] - per-target QC minimum mappability
#   args[12] - maximum number of latent factors

suppressMessages(library(CODEX))

#parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 12) {
  stop("Exactly 12 arguments must be supplied")
}
coverage_tsv_files <- scan(args[1], what="", quiet=TRUE)
contig <- args[2]
cohort_name <- args[3]
output_dir <- args[4]
min_median_coverage <- as.numeric(args[5])
max_median_coverage <- as.numeric(args[6])
min_length <- as.numeric(args[7])
max_length <- as.numeric(args[8])
min_GC <- as.numeric(args[9])
max_GC <- as.numeric(args[10])
min_mappability <- as.numeric(args[11])
max_num_latent_factors <- as.numeric(args[12])

#aggregate coverage TSVs and get targets
message("Aggregating coverage files...")
num_samples = length(coverage_tsv_files)
sample_names <- vector("list", num_samples)
if (num_samples < max_num_latent_factors) {
    stop("Maximum number of latent factors must be less than the number of samples.")
}
for (i in seq_along(coverage_tsv_files)) {
    coverage_df <- read.csv(coverage_tsv_files[i], sep="\t", stringsAsFactors=FALSE)
    if (i == 1) {
        targets <- coverage_df[c("CONTIG", "START", "END")]
        ref <- IRanges(start=targets[["START"]], end=targets[["END"]])
        Y <- matrix(nrow=nrow(targets), ncol=num_samples)
    } else if (!all.equal(coverage_df[c("CONTIG", "START", "END")], targets)) {
        stop("Targets must be identical across all coverage files.")
    }
    sample_names[i] <- coverage_df[[1, "SAMPLE"]]
    Y[, i] <- coverage_df[["Y"]]
}

#the following is nearly verbatim from the CODEX example, with some tidying up and additional logging

message("Computing GC content and mappability...")
gc <- getgc(contig, ref)
mapp <- getmapp(contig, ref)

message("Performing per-sample and per-target QC and writing QC matrix...")
qcObj <- qc(Y, sample_names, contig, ref, mapp, gc,
            cov_thresh=c(min_median_coverage, max_median_coverage),
            length_thresh=c(min_length, max_length),
            mapp_thresh=min_mappability,
            gc_thresh=c(min_GC, max_GC))
Y_qc <- qcObj$Y_qc
sampname_qc <- qcObj$sampname_qc
gc_qc <- qcObj$gc_qc
ref_qc <- qcObj$ref_qc
qcmat <- qcObj$qcmat
write.table(qcmat, file=paste(output_dir, "/", cohort_name, "_", contig, "_qcmat", ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

message("Fitting Poisson latent-factor model...")
normObj <- normalize(Y_qc, gc_qc, K=1:max_num_latent_factors)
Yhat <- normObj$Yhat

message("Determining number of latent factors and writing latent-factors plot...")
AIC <- normObj$AIC
BIC <- normObj$BIC
RSS <- normObj$RSS
K <- normObj$K
optK = K[which.max(BIC)]
choiceofK(AIC, BIC, RSS, K, filename=paste(output_dir, "/", cohort_name, "_", contig, "_choiceofK", ".pdf", sep=""))

message("Performing segmentation and writing segments...")
finalcall <- segment(Y_qc, Yhat, optK=optK, K=K, sampname_qc, ref_qc, contig, lmax=200, mode="integer")     #use integer mode for germline CNV calling
write.table(finalcall, file=paste(output_dir, "/", cohort_name, "_", contig, "_", optK, "_CODEX_integer.seg", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
