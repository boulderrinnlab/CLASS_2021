knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(ggdendro)
library(GenomicRanges)
source("util/intersect_functions.R")
source("util/plotting_functions.R")
source("util/_setup.R")

# Let's grab both the lncRNA and mRNA promoters
promoters <- rtracklayer::import("results/lncrna_mrna_promoters.gtf")

peak_occurence_matrix <- read.table("results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Let's filter again to make sure we have at least 250 peaks:

peak_occurence_matrix <- peak_occurence_matrix[rowSums(peak_occurence_matrix) >250, ]

# ok now we have all the things we need to start clustering !
# First we will use hclust

# ?hclust
# bin_hier
# R is so savy with statistics that this is actually a base R function!
# There are lots of methods to cluster by let's start with binary.

# Hierarchical clustering with binary distance measure
bin_hier <- hclust(dist(peak_occurence_matrix, method = "binary"))


bin_hier$labels[bin_hier$order]
# First we are going to order the samples by the correlation matrix we just made.

# Now we can plot
ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3, 
                       theme_dendro = TRUE) +
  coord_flip() +
  # scale_y_continuous() +
  # scale_x_continuous(position = "top") +
  scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
                     labels = bin_hier$labels[bin_hier$order], position = "top",
                     expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, hjust  = 1)) + 
  theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
  scale_y_reverse(expand = c(0.01, 0)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

pdf("./figures/dbp_hclust_dendro.pdf", height = 12, width = 70)
plot(bin_hier)
dev.off()


#tree cutting:
clust <- cutree(bin_hier,h=.52)
table(clust)
clust_df <- data.frame(dbp = names(clust), cluster = clust)

#choos dbps in cluster11
cluster11 <- filter(clust_df,cluster == 11)

#peak occurrences of cluster 11 dbps
po_clust11 <- peak_occurence_matrix[row.names(peak_occurence_matrix) %in% row.names(cluster11),]

#promoters that ALL dpbs in cluster 11 bind (there are 10 dpbs)
pro_cluster11 <- po_clust11[,colSums(po_clust11) > 9]

#pull out the promoter names that all dbps in cluster11 bound
pro_clus11 <- data.frame(promoters = colnames(pro_cluster11))




write_csv(pro_clus11,"./promoterID_cluster11_2")
