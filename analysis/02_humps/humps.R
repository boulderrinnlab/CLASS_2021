knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(tidyverse)
source("util/plotting_functions.R")
source("util/_setup.R")


#load in our consensus peaks and counted peak occurrence

filtered_consensus_peaks_files <- list.files("/Users/maca9088/CLASS/Maria/CLASS_2021/analysis/00_consensus_peaks/results/filtered_consensus_peaks",
                                             pattern = "*.bed",
                                             full.names = TRUE)
filtered_consensus_peaks <- lapply(filtered_consensus_peaks_files, rtracklayer::import)
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


#count the peaks per lncrna and mrna genebody
genebody_peak_occurence <- count_peaks_per_feature(lncrna_mrna_genebody, filtered_consensus_peaks, 
                                                   type = "occurrence")
# Output to promoter_peak_occurecne_matrix

#write.table(genebody_peak_occurence, "./analysis/02_humps/results/genebodies_peak_occurence_matrix.tsv")
# Now we want to make into a data frame using the promoter annotations as rows and attributes as columns.




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
#write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")

#plotting mrna peak occurence
g <- ggplot(peak_occurence_genebody_df, aes(x = number_of_dbp, fill = gene_type, color = gene_type))
g + geom_density(alpha = 0.2, color = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Genebody binding events",
          subtitle = "mRNA and lncRNA gene bodies")


#filter promoter peaks into two groups: above 200 dbp and below
peak_occurence_promoters_df <- peak_occurence_df

#filter genes above 290 DBPs and less than 310 DBPs
bump2_promoters <- filter(peak_occurence_df, peak_occurence_df$number_of_dbp > 200)

bump2_promoters_genenames <- as.data.frame(bump2_promoters$gene_name)
bump2_promoters_geneids <- as.data.frame(bump2_promoters$gene_id)

mid_bump2_promoters <- filter(peak_occurence_df, peak_occurence_df$number_of_dbp > 290)
mid_bump2_promoters <- filter(mid_bump2_promoters, mid_bump2_promoters$number_of_dbp < 310)

mid_bump2_promoters_geneids <-as.data.frame(mid_bump2_promoters$gene_id)

bump1_promoters <- filter(peak_occurence_df, peak_occurence_df$number_of_dbp < 100)
bump1_promoters_genenames <- as.data.frame(bump1_promoters$gene_name)
bump1_promoters_geneids <- as.data.frame(bump1_promoters$gene_id)


#write the bump1 and bump2 gene names to csv file

write_csv(bump1_promoters_genenames, "results/bump1_promoters_genenames.csv")
write_csv(bump2_promoters_genenames, "results/bump2_promoters_genenames.csv")
write_csv(bump1_promoters_geneids, "results/bump1_promoters_geneids.csv")
write_csv(bump2_promoters_geneids, "results/bump2_promoters_geneids.csv")

write_csv(mid_bump2_promoters_geneids, "results/mid_bump2_promoters_geneids.csv")

#remove the version
#sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' file_with_version > file_without_version 

#percentage of mRNA in bump2 

n_bump2_mrna_promoters = filter(bump2_promoters, bump2_promoters$gene_type == 'protein_coding')
n_bump2_lncrna_promoters = filter(bump2_promoters, bump2_promoters$gene_type == 'lncRNA')

n_bump2_mrna_promoters <- nrow(n_bump2_mrna_promoters)
n_bump2_lncrna_promoters <- nrow(n_bump2_lncrna_promoters)

percent_mrna_bump2 <- n_bump2_mrna_promoters/(n_bump2_mrna_promoters+n_bump2_lncrna_promoters)

#percent mrna in bump1
n_bump1_mrna_promoters = filter(bump1_promoters, bump1_promoters$gene_type == 'protein_coding')
n_bump1_lncrna_promoters = filter(bump1_promoters, bump1_promoters$gene_type == 'lncRNA')

n_bump1_mrna_promoters <- nrow(n_bump1_mrna_promoters)
n_bump1_lncrna_promoters <- nrow(n_bump1_lncrna_promoters)

percent_mrna_bump1 <- n_bump1_mrna_promoters/(n_bump1_mrna_promoters+n_bump1_lncrna_promoters)


#interesting!!! The percent mrna in bump1 is 39.53 % and in bump2 is 76.06%
