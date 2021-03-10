```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ggpubr)
library(Gviz)
source("util/_setup.R")
source("util/intersect_functions.R")
source("util/plotting_functions.R")
```

```{r}
Group1consensus_peaks <- create_consensus_peaks('/scratch/Shares/rinnclass/Chelsea/second_chipseq/results_data_replicates')

#Peaks from replicates are merged

names(Group1consensus_peaks) <- sapply(Group1consensus_peaks, function(x){
  unlist(strsplit(x$name, "_"))[[1]]
})

#Peaks are renamed to something readable via sapply string splitting

num_peaks_df <- data.frame("dbp" = names(Group1consensus_peaks),
                           "num_peaks" = sapply(Group1consensus_peaks, length))
#We are collecting the number of peaks for each of our 4 proteins

g <- ggplot(num_peaks_df, aes(x = num_peaks)) + 
  geom_histogram(bins = 70)

show(g)

for(i in 1:length(Group1consensus_peaks)) {
  rtracklayer::export(Group1consensus_peaks[[i]], paste0("results/filtered_peaks/", names(Group1consensus_peaks)[i], 
                                                            "_consensus_peaks_filter.bed"))
}

num_peaks_df$total_peak_length <- sapply(Group1consensus_peaks, function(x) sum(width(x)))


gencode_gr <- rtracklayer::import("/Shares/rinn_class/data/genomes/human/gencode/v32/gencode.v32.annotation.gtf")

lncrna_mrna_promoters <- get_promoter_regions(gencode_gr, biotype = c("lncRNA", "protein_coding"))
rtracklayer::export(lncrna_mrna_promoters, "results/lncrna_mrna_promoters.gtf")

lncrna_promoters <- get_promoter_regions(gencode_gr, biotype = "lncRNA")
rtracklayer::export(lncrna_promoters, "results/lncrna_promoters.gtf")

mrna_promoters <- get_promoter_regions(gencode_gr, biotype = "protein_coding")
rtracklayer::export(mrna_promoters, "results/mrna_promoters.gtf")

promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, Group1consensus_peaks, type = "occurrence")

num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_promoters$gene_id])

num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_promoters$gene_id])
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/Olivia/CLASS_2021/results/num_peaks_df.csv")

```