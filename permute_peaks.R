knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(GenomicRanges)
library(rslurm)
library(regioneR)
library(effectsize)
source("util/intersect_functions.R")
source("util/_setup.R")

# IMPORT peaks and features to test overlaps during the randomizations(permutations) ####
peak_file_list <- list.files("/scratch/Shares/rinnclass/data/consensus_peaks/", full.names = T)
consensus_peaks <- lapply(peak_file_list, rtracklayer::import)
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/data/consensus_peaks//|.bed", "", peak_file_list)
# Filter the consensus peaks with too few peaks
peak_list <- consensus_peaks[sapply(consensus_peaks, length) > 250]

base_path <- "/scratch/Shares/rinnclass/JR/CLASS_2021/analysis/00_consensus_peaks/results/gene_annotations/"
lncrna_promoters <- rtracklayer::import(paste0(base_path, "lncrna_promoters.gtf"))
mrna_promoters <- rtracklayer::import(paste0(base_path, "mrna_promoters.gtf"))
lncrna_mrna_promoters <- rtracklayer::import(paste0(base_path, "lncrna_mrna_promoters.gtf"))

lncrna_matrix <- count_peaks_per_feature(lncrna_promoters, peak_list, type = "occurrence")
mrna_matrix <- count_peaks_per_feature(mrna_promoters, peak_list, type = "occurrence")
lncrna_mrna_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, peak_list, type = "occurrence")

stopifnot(file.exists("/Shares/rinn_class/data/k562_chip/util/rmsk.txt"))
rmsk <- import_repeatmasker("/Shares/rinn_class/data/k562_chip/util/rmsk.txt")
rmsk_family <- subset_rmsk(rmsk, rep_level = "family")
names(rmsk_family) <- paste(names(rmsk_family), "family", sep = "_")
rmsk_class <- subset_rmsk(rmsk, rep_level = "class")
names(rmsk_class) <- paste(names(rmsk_class), "class", sep = "_")

hg38 <- getGenome("hg38")

region_list <- c("lncrna_promoters" = list(lncrna_promoters), 
                 "mrna_promoters" = list(mrna_promoters), 
                 "lncrna_mrna_promoters" = list(lncrna_mrna_promoters),
                 rmsk_family, rmsk_class)


canonical_chr <- as.character(unique(seqnames(region_list[[1]])))
# Sanitize RMSK to only those repeats on canonical chromosomes
for(i in 1:length(region_list)) {
  region_list[[i]] <- region_list[[i]][which(seqnames(region_list[[i]]) %in% canonical_chr)]
}

#### DEFINE Permutation TESTS ####
pars <- expand.grid(1:length(region_list), 1:length(peak_list)) %>% 
  as.data.frame()
names(pars) <- c("region_index", "peak_index")



#### TEST FUNCTION ####
perm_test <- function(region_index, peak_index, npermutations = 1000) {
  
  set.seed(12044593)
  region <- names(region_list)[[region_index]]
  tf <- names(peak_list)[[peak_index]]
  
  cat("Running overlap test for ", region, "  & ", tf, "\n\n")
  a_regions <- region_list[[region_index]]
  b_regions <- peak_list[[peak_index]]
  
  suppressWarnings(pt <- overlapPermTest(A = a_regions, 
                                         B = b_regions, 
                                         ntimes = npermutations, 
                                         non.overlapping = FALSE, 
                                         verbose = FALSE,
                                         genome = hg38,
                                         alternative =  "auto", 
                                         mc.cores = 1))
  
  ptdf <- data.frame("region" = region,
                     "tf" = tf,
                     "pval" = pt$numOverlaps$pval,
                     "zscore" = pt$numOverlaps$zscore,
                     "nperm" = pt$numOverlaps$ntimes,
                     "alternative" = pt$numOverlaps$alternative,
                     "observed" = pt$numOverlaps$observed,
                     "permuted" = paste(pt$numOverlaps$permuted, collapse = ";"))
  return(ptdf)
}

# Note this took about 400 hours to run
# TODO: re-write to be faster
res_files <- list.files("_rslurm_perm_overlaps/", full.names = T, pattern = "results")
if(length(res_files) == 0) {
  sjob <- slurm_apply(perm_test, pars, jobname = 'perm_overlaps',
                      add_objects = c("region_list", "peak_list", "hg38", "overlapPermTest"),
                      nodes = 15, cpus_per_node = 30, 
                      slurm_options = list(time = '400:00:00'),
                      submit = FALSE)
}