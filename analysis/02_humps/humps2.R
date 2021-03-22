knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(tidyverse)
source("./util/plotting_functions.R")
source("./util/_setup.R")
source("./util/intersect_functions.R")



#so, lets find the proteins that occur in the 2nd peak! 

bump2_promoters <- read_csv('results/bump2_promoters_geneids.csv')
num_peaks_df <- read_csv('./analysis/00_consensus_peaks/results/num_peaks_df.csv')

#change column names
names(bump2_promoters)[names(bump2_promoters) == "bump2_promoters$gene_id"] <- "gene_id"

#read in the lncrna mrna peak occurrence
lncrna_mrna_promoters_peaks <- 
  read_table('/Users/maca9088/CLASS/Maria/CLASS_2021/analysis/00_consensus_peaks/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv')

#opa we'll do that later, let's look at zinc fingers

lncrna_mrna_promoters <- rtracklayer::import("./analysis/00_consensus_peaks/promoters/lncrna_mrna_promoters.gtf")

zinc_fingers <- filter(num_peaks_df, num_peaks_df$dbd == 'C2H2 ZF')

zf_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, zinc_fingers, 
                                                   type = "occurrence")
genebody_peak_occurence <- count_peaks_per_feature(lncrna_mrna_genebody, filtered_consensus_peaks, 
                                                   type = "occurrence")



peak_occurence_df <- read_csv('./analysis/00_consensus_peaks/results/peak_occurence_dataframe.csv')
g <- ggplot(peak_occurence_df, aes(x = number_of_dbp,))
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes")