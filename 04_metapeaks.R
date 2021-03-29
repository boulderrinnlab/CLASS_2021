options(stringsAsFactors = FALSE)
knitr::opts_chunk$set(echo = TRUE)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


data_path <- "/Shares/rinn_class/data/CLASS_2021/CLASS_2021"

gencode_gr <- rtracklayer::import(file.path(data_path,"data/gencode.v32.annotation.gtf"))
consensus_peaks <- import_peaks(file.path(data_path, "results/filtered_peaks"))

midbump_promoters <- read_csv("analysis/02_humps/results/bump2_promoters_geneids.csv")

genes <- gencode_gr[gencode_gr$type == "gene"]
# This function is a convenience function built into GenomicRanges
all_promoters_gr <- promoters(genes, upstream = 3e3, downstream = 3e3)
table(width(all_promoters_gr))



profile_tss <- function(peaks, 
                        promoters_gr,
                        upstream = 3e3,
                        downstream = 3e3) {
  
  
  peak_coverage <- coverage(peaks)
  
  coverage_length <- elementNROWS(peak_coverage)
  coverage_gr <- GRanges(seqnames = names(coverage_length),
                         IRanges(start = rep(1, length(coverage_length)), 
                                 end = coverage_length))
  
  promoters_gr <- subsetByOverlaps(promoters_gr, 
                                   coverage_gr, 
                                   type="within", 
                                   ignore.strand=TRUE)
  chromosomes <- intersect(names(peak_coverage), 
                           unique(as.character(seqnames(promoters_gr))))
  peak_coverage <- peak_coverage[chromosomes]
  
  promoters_ir <- as(promoters_gr, "IntegerRangesList")[chromosomes]
  
  promoter_peak_view <- Views(peak_coverage, promoters_ir)
  
  promoter_peak_view <- lapply(promoter_peak_view, function(x) t(viewApply(x, as.vector)))
  promoter_peak_matrix <- do.call("rbind", promoter_peak_view)
  
  minus_idx <- which(as.character(strand(promoters_gr)) == "-")
  promoter_peak_matrix[minus_idx,] <- promoter_peak_matrix[minus_idx,
                                                           ncol(promoter_peak_matrix):1]
  
  promoter_peak_matrix <- promoter_peak_matrix[rowSums(promoter_peak_matrix) > 1,]
  
  peak_sums <- colSums(promoter_peak_matrix)
  peak_dens <- peak_sums/sum(peak_sums)
  
  metaplot_df <- data.frame(x = -upstream:(downstream-1),
                            dens = peak_dens)
  
  return(metaplot_df)
}

ZNF34_peaks <- consensus_peaks[["ZNF34"]]
ZNF343_peaks <- consensus_peaks[["ZNF343"]]
ATF5_peaks <- consensus_peaks[["ATF5"]]

ZNF34_metaplot_df <- profile_tss(ZNF34_peaks,all_promoters_gr)
ZNF343_metaplot_df <- profile_tss(ZNF343_peaks,all_promoters_gr)
ATF5_metaplot_df <- profile_tss(ATF5_peaks,all_promoters_gr)

ZNF34_metaplot_df["dbp"] <- "ZNF34"
ZNF343_metaplot_df["dbp"] <- "ZNF343"
ATF5_metaplot_df["dbp"] <- "ATF5"

all_df <- bind_rows(ZNF343_metaplot_df,ZNF34_metaplot_df,ATF5_metaplot_df)

ggplot(all_df, aes(x = x, y = dens, color = dbp)) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = .50) + 
  ggtitle("ZNF343 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency")
