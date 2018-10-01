#!/usr/bin/env Rscript

################################################################################################################################################
### This script is to convert all calls (truth and gCNV) into bin-based coordinates for calculating statistics.
### [Achieved by mapping genomic coordinates into IRanges (of bin numbers)]
###
### It requires arguments in the following order:
###   args[1] - Path to the gCNV calls
###   args[2] - Path to the filtered, annotated intervals TSV file
###   args[3] - Path to the truth bed file
###   args[4] - Path to the sample-ID map TSV with columns (SFARI.ID, Sample.ID)
###   args[5] - QS threshold
###   args[6] - Number of cores to use
###   args[7] - Output directory
###
#############################################################################################################################################

### Parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
  stop("Exactly 7 arguments must be supplied")
}

gcnv_calls_path <- args[1]
filtered_annotated_intervals_file <- args[2]
truth_bed_file <- args[3]
sample_id_map_file <- args[4]
threshold_QS <- as.numeric(args[5])
num_cores <- as.numeric(args[6])
out_dir <- args[7]

################################################################################################################################################
### Functions
################################################################################################################################################
### Function for converting calls into bin-space
toBinSpace <- function(calls, bins, gcnv=FALSE){
    ol <- findOverlaps(calls, bins)
    ir <- GRanges()
    if(length(ol)>0){
        ir <- by(subjectHits(ol), queryHits(ol), function(x) GRanges(0, IRanges(min(x), max(x))))
            names(ir) <- NULL
            ir <- do.call(c, ir)
        if(gcnv==TRUE){
            ir$QA <- calls$QA[unique(queryHits(ol))]
            ir$QS <- calls$QS[unique(queryHits(ol))]
            ir$QSS <- calls$QSS[unique(queryHits(ol))]
            ir$QSE <- calls$QSE[unique(queryHits(ol))]
        }
    }
    return(ir)
}
### Function for creating a GRanges object for each subject [samp] that has that subject's DEL and DUP truth calls.
extractCallsTruth <- function(samp, truth_calls){
    ind <- which(str_detect(truth_calls[,6], samp))
    out <- GRanges(truth_calls[,1], IRanges(truth_calls[,2], truth_calls[,3]), call=truth_calls[,5], sample=samp, af=str_count(truth_calls[,6], ",")/2076)[ind]    
    strand(out) <- "-"; strand(out)[out$call=="DUP"] <- "+"
    return(out)
}
### Function for creating a GRanges object for each subject [samp] that has that subject's DEL and DUP gCNV calls.
extractCallsgCNV <- function(file){
    vcf <- readVcf(file)    
    out <- rowRanges(vcf)
        end(out) <- info(vcf)$END
    out$QA <- geno(vcf)$QA
    out$QS <- geno(vcf)$QS
    out$QSS <- geno(vcf)$QSS
    out$QSE <- geno(vcf)$QSE
    out$CN <- as.numeric(unlist(geno(vcf)$CN))
    out <- out[out$CN!=2]
    out$call <- "DEL"
        out$call[out$CN>2] <- "DUP"
    strand(out) <- "-"
        strand(out)[out$CN>2] <- "+"
    names(out) <- NULL
    
    out <- out[seqnames(out) %in% 1:22]
    return(out)
}
### Compile recipricol overlap information on bin-space
getROL <- function(i, bs_gcnv, bs_truth){
    ol <- findOverlaps(bs_gcnv[[i]], bs_truth[[i]])
    sample <- names(bs_gcnv)[i]
    if(length(ol)>0){
        gcnv <- bs_gcnv[[i]][queryHits(ol)]
        truth <- bs_truth[[i]][subjectHits(ol)]
        int <- intersect(gcnv, truth)
        rol <- round(width(int)/apply(cbind(width(gcnv), width(truth)), 1, max), 3)

        info <- data.frame(gcnv_start=start(gcnv), gcnv_end=end(gcnv), gcnv_wid=width(gcnv),
            truth_start=start(truth), truth_end=end(truth), truth_wid=width(truth),
            rol=rol, QA=gcnv$QA, QS=gcnv$QS)

        ## Fill in gcnv and truth sets that do not have correpsonding in the other
        gcnv2 <- bs_gcnv[[i]][-queryHits(ol),,drop=FALSE]
        truth2 <- bs_truth[[i]][-subjectHits(ol),,drop=FALSE]
        info2 <- NULL
        info3 <- NULL

        if(length(gcnv2)>0)
            info2 <- data.frame(gcnv_start=start(gcnv2), gcnv_end=end(gcnv2), gcnv_wid=width(gcnv2),
                truth_start=0, truth_end=0, truth_wid=0, rol=0, QA=gcnv2$QA, QS=gcnv2$QS)
        if(length(truth2)>0)
            info3 <- data.frame(gcnv_start=0, gcnv_end=0, gcnv_wid=0,
                truth_start=start(truth2), truth_end=end(truth2), truth_wid=width(truth2), rol=0, QA=0, QS=0)

        out <- rbind(info, info2, info3)
    }else{
        info2 <- NULL
        info3 <- NULL
        gcnv <- bs_gcnv[[i]]
        truth <- bs_truth[[i]]

        if(length(gcnv)>0)
            info2 <- data.frame(gcnv_start=start(gcnv), gcnv_end=end(gcnv), gcnv_wid=width(gcnv),
                truth_start=0, truth_end=0, truth_wid=0, rol=0, QA=gcnv$QA, QS=gcnv$QS)
        if(length(truth)>0)
            info3 <- data.frame(gcnv_start=0, gcnv_end=0, gcnv_wid=0,
                truth_start=start(truth), truth_end=end(truth), truth_wid=width(truth), rol=0, QA=0, QS=0)
        out <- rbind(info2, info3)
    }
    out$sample <- sample
    return(out)
}
### Compile performance
calculatePerformance <- function(rol_set, threshold_rol, threshold_bin){
    results <- NULL
    for(i in threshold_bin){
        results <- rbind(results, c(
        ### Sensitivity
        mean(rol_set$rol[rol_set$truth_wid>=i] > threshold_rol),
        ### n_var_truth
        sum(rol_set$truth_wid==i),
        ### PPV
        mean(rol_set$rol[rol_set$gcnv_wid>=i] > threshold_rol),
        ### n_var_gcnv
        sum(rol_set$gcnv_wid==i)))

    }
    results <- rbind(results, c(
    ### Sensitivity
    mean(rol_set$rol[rol_set$truth_wid>i] > threshold_rol),
    ### n_var_truth
    sum(rol_set$truth_wid>i),
    ### PPV
    mean(rol_set$rol[rol_set$gcnv_wid>i] > threshold_rol),
    ### n_var_gcnv
    sum(rol_set$gcnv_wid>i)))

    results <- round(results, 2)
    rownames(results) <- c(threshold_bin, paste0(">", i))
    colnames(results) <- c("sens", "n_var_truth", "ppv", "n_var_gcnv")
    return(results)
}
### Plot for evaluating sens, spec, and var sizes
performancePlot <- function(matrix_overall, matrix_del, matrix_dup, wids, label, path){
    pdf(paste0(path, "/", label, ".pdf"), width=12, height=4)
    par(mfrow=c(1,3))
    par(mar=rep(4,4))
    xs <- 1:dim(matrix_overall)[1]
    plot(x=c(0, dim(matrix_overall)[1]), y=c(0, 1), ty='n', pch='', xaxt='n', ylab = 'sensitivity', xlab = 'number of bins')
        points(matrix_overall[,1]~xs, pch=19, cex=0.5)
        points(matrix_del[,1]~xs, pch=19, col=2, cex=0.5)
        points(matrix_dup[,1]~xs, pch=19, col=4, cex=0.5)
        legend(20, y=0.2, pch=19, legend=c("Overall", "Dels", "Dups"), col=c(1, 2, 4))
        axis(side=3, at=xs, labels=matrix_overall[,2], las=2)
        axis(side=1, at=xs, labels=rownames(matrix_overall), las=2)
    plot(x=c(0, dim(matrix_overall)[1]), y=c(0, 1), ty='n', pch='', xaxt='n', ylab = 'ppv', xlab = 'number of bins', main="")
        points(matrix_overall[,3]~xs, pch=19, cex=0.5)
        points(matrix_del[,3]~xs, pch=19, col=2, cex=0.5)
        points(matrix_dup[,3]~xs, pch=19, col=4, cex=0.5)
        axis(side=3, at=xs, labels=matrix_overall[,4], las=2)
        axis(side=1, at=xs, labels=rownames(matrix_overall), las=2)
    wids_use <- wids
    wids_use[wids[,2]>(max(xs)-1),2] <- 26
    boxplot(log(wids_use[,1],10)~wids_use[,2], xlab="number of bins", ylab="log 10 variant size", main=label, xaxt="n")
        axis(side=1, at=xs, labels=rownames(matrix_overall), las=2)
        abline(h=4, col=2)
    dev.off()
}

################################################################################################################################################
### Setting up necessary resources
##########################################################################################
library(rtracklayer); library(VariantAnnotation); library(stringr)

bins <- read.table(filtered_annotated_intervals_file, comment.char="@", header=TRUE)
bins <- GRanges(bins[,1], IRanges(bins[,2], bins[,3]), gc=bins$GC_CONTENT, segdup=bins$SEGMENTAL_DUPLICATION_CONTENT)
bins <- bins[seqnames(bins) %in% 1:22]
bins <- bins[bins$segdup==0]
sample_id_map <- read.table(sample_id_map_file, header=TRUE)

##########################################################################################
### Determining which families have full quartests
##########################################################################################
files <- list.files(gcnv_calls_path, pattern=".vcf.gz", full.names=TRUE)
samples <- str_replace(str_replace(basename(files), "genotyped-segments-", ""), ".vcf.gz", "")
ids <- as.character(sample_id_map$SFARI.ID[match(samples, sample_id_map$Sample.ID)])

################################################################################################################################################
### Parsing the truth callset and converting to bin-space
### There are many more samples in the table than the smaller exome set
##########################################################################################
truth_full <- read.table(truth_bed_file)
truth_full <- truth_full[truth_full[,1] %in% c(1:22),]
truth_del <- truth_full[truth_full[,5]=="DEL",]
truth_dup <- truth_full[truth_full[,5]=="DUP",]

calls_truth_del <- mclapply(ids, extractCallsTruth, truth_del, mc.cores=num_cores)
names(calls_truth_del) <- ids
calls_truth_dup <- mclapply(ids, extractCallsTruth, truth_dup, mc.cores=num_cores)
names(calls_truth_dup) <- ids

bs_truth_del <- mclapply(calls_truth_del, toBinSpace, bins, mc.cores=num_cores)
bs_truth_dup <- mclapply(calls_truth_dup, toBinSpace, bins, mc.cores=num_cores)

##########################################################################################
### Parsing the gCNV callset and converting to bin-space.
##########################################################################################
calls_gcnv <- mclapply(files, extractCallsgCNV, mc.cores=num_cores)
names(calls_gcnv) <- ids
calls_gcnv_bins <- mclapply(calls_gcnv, function(x){x$numBins=countOverlaps(x, bins); return(x)}, mc.cores=num_cores)

calls_gcnv_del <- mclapply(calls_gcnv_bins, function(x) x[x$call=="DEL"], mc.cores=num_cores)
calls_gcnv_dup <- mclapply(calls_gcnv_bins, function(x) x[x$call=="DUP"], mc.cores=num_cores)

bs_gcnv_del <- mclapply(calls_gcnv_del, toBinSpace, bins, gcnv=TRUE, mc.cores=num_cores)
bs_gcnv_dup <- mclapply(calls_gcnv_dup, toBinSpace, bins, gcnv=TRUE, mc.cores=num_cores)

##########################################################################################
### Comparing callsets by calculating bin-space recipricol overlap
##########################################################################################
rols_del <- mclapply(1:length(bs_gcnv_del), getROL, bs_gcnv_del, bs_truth_del, mc.cores=num_cores)
rols_dup <- mclapply(1:length(bs_gcnv_dup), getROL, bs_gcnv_dup, bs_truth_dup, mc.cores=num_cores)

rols_del <- do.call(rbind, rols_del)
rols_del$var <- "-"
rols_dup <- do.call(rbind, rols_dup)
rols_dup$var <- "+"
rols_all <- rbind(rols_del, rols_dup)
rownames(rols_all) <- NULL
wids_all <- data.frame(var_wid=end(bins)[rols_all[,2]]-start(bins)[rols_all[,1]], bins=rols_all$gcnv_wid[rols_all$gcnv_wid!=0])

threshold_bin <- 1:25
performance <- calculatePerformance(rols_all, threshold_rol=0, threshold_bin=1:25)
performance_del <- calculatePerformance(rols_all[rols_all$var=="-",], threshold_rol=0, threshold_bin=1:25)
performance_dup <- calculatePerformance(rols_all[rols_all$var=="+",], threshold_rol=0, threshold_bin=1:25)

##########################################################################################
### Evaluating performance (Rare)
### Removing variants more than 20% covered by bins called in 4 or more parents
### in either the truth set or the gcnv set
##########################################################################################
threshold_ac <- 4
threshold_ol <- 0.2
ir_tots <- IRanges(rols_all[,4], rols_all[,5])
ir_parents <- ir_tots[stringr::str_detect(rols_all$sample, "(mo)|(fa)")]
cluster_parents <- as.numeric(coverage(ir_tots))
filter_parents <- which(cluster_parents>=threshold_ac)
filter_wgs <- IRanges(filter_parents, filter_parents)

ir_gcnv <- IRanges(rols_all[,1], rols_all[,2])
ir_parents <- ir_gcnv[stringr::str_detect(rols_all$sample, "(mo)|(fa)")]
cluster_parents <- as.numeric(coverage(ir_gcnv))
filter_parents <- which(cluster_parents>=threshold_ac)
filter_gcnv <- IRanges(filter_parents, filter_parents)

col_wgs <- countOverlaps(ir_tots, filter_wgs)/width(ir_tots)
col_gcnv <- countOverlaps(ir_gcnv, filter_gcnv)/width(ir_gcnv)

rols_rare <- rols_all[col_wgs<=threshold_ol & col_gcnv<=threshold_ol,]
wids_rare <- data.frame(var_wid=end(bins)[rols_rare[,2]]-start(bins)[rols_rare[,1]], bins=rols_rare$gcnv_wid[rols_rare$gcnv_wid!=0])

performance_rare <- calculatePerformance(rols_rare, threshold_rol=0, threshold_bin=1:25)
performance_rare_del <- calculatePerformance(rols_rare[rols_rare$var=="-",], threshold_rol=0, threshold_bin=1:25)
performance_rare_dup <- calculatePerformance(rols_rare[rols_rare$var=="+",], threshold_rol=0, threshold_bin=1:25)

##########################################################################################
### Evaluating performance (Rare, good batch, QS)
##########################################################################################
rols_rare_batch_QS <- rols_rare
rols_rare_batch_QS[which(rols_rare_batch_QS$gcnv_wid>0 & rols_rare_batch_QS$QS<threshold_QS), c(1:3, 7)] <- 0
rols_rare_batch_QS <- rols_rare_batch_QS[-which(apply(rols_rare_batch_QS[,1:6], 1, sum)==0),]
wids_rare_batch_QS <- data.frame(var_wid=end(bins)[rols_rare_batch_QS[,2]]-start(bins)[rols_rare_batch_QS[,1]], bins=rols_rare_batch_QS$gcnv_wid[rols_rare_batch_QS$gcnv_wid!=0])

performance_rare_batch_QS <- calculatePerformance(rols_rare_batch_QS, threshold_rol=0, threshold_bin=1:25) 
performance_rare_batch_QS_del <- calculatePerformance(rols_rare_batch_QS[rols_rare_batch_QS$var=="-",], threshold_rol=0, threshold_bin=1:25)
performance_rare_batch_QS_dup <- calculatePerformance(rols_rare_batch_QS[rols_rare_batch_QS$var=="+",], threshold_rol=0, threshold_bin=1:25)

##########################################################################################
### Visualize Sensitivity, PPV, and number of variants
##########################################################################################
dir.create(out_dir)
performancePlot(performance, performance_del, performance_dup, wids_all, "performance_overall", out_dir)
performancePlot(performance_rare, performance_rare_del, performance_rare_dup, wids_rare, "performance_rare", out_dir)
performancePlot(performance_rare_batch_QS, performance_rare_batch_QS_del, performance_rare_batch_QS_dup, wids_rare_batch_QS, "performance_QS_filtered", out_dir)
