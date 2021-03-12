knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(tidyverse)
source("util/plotting_functions.R")
source("util/_setup.R")


#load in our consensus peaks and counted peak occurrence
num_peaks_df <- read_csv('./analysis/00_consensus_peaks/results/num_peaks_df.csv')
peak_occurence_df <- read_csv('./analysis/00_consensus_peaks/results/peak_occurence_dataframe.csv')

#load genebodies of lncrna and mrna
lncrna_mrna_genebody <- rtracklayer::import("./analysis/00_consensus_peaks/genebodies/lncrna_mrna_genebody.gtf")
mrna_genebody <- rtracklayer::import("./analysis/00_consensus_peaks/genebodies/mrna_genebody.gtf")
lncrna_mrna_genebody <- rtracklayer::import("./analysis/00_consensus_peaks/genebodies/lncrna_mrna_genebody.gtf")

#load promoters of lncrna and mrna
lncrna_mrna_promoters <- rtracklayer::import("./analysis/00_consensus_peaks/promoters/lncrna_mrna_promoters.gtf")
lncrna_promoters <- rtracklayer::import("./analysis/00_consensus_peaks/promoters/lncrna_promoters.gtf")
mrna_promoters <- rtracklayer::import("./analysis/00_consensus_peaks/promoters/mrna_promoters.gtf")

#g ranges to data frames
mrna_promoters_df <- as.data.frame(mrna_promoters)
lncrna_promoters_df <- as.data.frame(lncrna_promoters)

lncrna_genebodies_df <- as.data.frame(lncrna_genebody)
lncrna_genebodies_df <- as.data.frame(mrna_genebody)
lncrna_mrna_genebodies_df <-as.data.frame(lncrna_mrna_genebody)


#filter out the lncrnas, etc
peak_occurence_lncrna <- filter(peak_occurence, gene_type == "lncRNA")
peak_occurence_mrna <- filter(peak_occurence, gene_type == "protein_coding")


#plotting lncrna peak occurence
g <- ggplot(peak_occurence_df, aes(x = number_of_dbp, fill = gene_type, color = gene_type))
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes")


#plotting lncrna peak occurence
g <- ggplot(peak_occurence_lncrna, aes(x = number_of_dbp), fill=gene_type)
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes")


#plotting mrna peak occurence
g <- ggplot(peak_occurence_mrna, aes(x = number_of_dbp), fill=gene_type)
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes")



genebody_peak_occurence <- count_peaks_per_feature(lncrna_mrna_genebody, filtered_consensus_peaks, 
                                                   type = "occurrence")
# Output to promoter_peak_occurecne_matrix
write.table(genebody_peak_occurence, "./analysis/00_consensus_peaks/results/genebodies_peak_occurence_matrix.tsv")
# Now we want to make into a data frame using the promoter annotations as rows and attributes as columns.
# We will use lncrna_mrna_promoters to index "all promoters"
# First make sure promoter_peak_occurence and lncrna_mrna_promoters are in the same order
    stopifnot(all(colnames(genebody_peak_occurence) == lncrna_mrna_promoters$gene_id))



# We are going to use the promoter peak occurence matrix above to essentially
# recreate a working version of num_peaks_df. However we will now organize it more
# this is a good example of how to set up a bunch of columns in dataframe. using the
# data.frame() fucntion.
# essentially we will index values from the objects we created and make a .CSV 
# to keep adding onto in future classes.

peak_occurence_genebody_df <- data.frame("gene_id" = colnames(genebody_peak_occurence),
                                "gene_name" = lncrna_mrna_genebody$gene_name,
                                "gene_type" = lncrna_mrna_genebody$gene_type,
                                "chr" = lncrna_mrna_genebody@seqnames,   
                                "3kb_up_tss_start" = lncrna_mrna_genebody@ranges@start,
                                "strand" = lncrna_mrna_genebody@strand,
                                "number_of_dbp" = colSums(genebody_peak_occurence))
# This is the CSV file we will start building upon adding columns of properties as we analyze them
# The output file name will change based on what is added later, but the "peak_occurence_df" will be used throughout.
write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")

#plotting mrna peak occurence
g <- ggplot(peak_occurence_genebody_df, aes(x = number_of_dbp, fill = gene_type, color = gene_type))
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Genebody binding events",
          subtitle = "mRNA and lncRNA gene bodies")

