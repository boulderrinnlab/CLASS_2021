#playing around
knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(tximport)
library(DESeq2)

source("util/plotting_functions.R")
source("util/_setup.R")

# We want to include the gene annotation info in the DEseq2 object,
# so we'll need to read that in.
gencode_gtf <- rtracklayer::import("/scratch/Shares/rinnclass/data/genomes/gencode.v32.annotation.gtf")
genes <- gencode_gtf[gencode_gtf$type == "gene"]

# Salmon's primary level of quantifying is at the transcript-level.
# Since we're doing differential expression at the gene-level
# Salmon needs to be able to map transcripts to genes.
tx2gene <- gencode_gtf[gencode_gtf$type == "transcript"] %>% 
  as.data.frame() %>%
  dplyr::select(transcript_id, gene_id)

# We'll also create an object that maps gene IDs to gene symbols
# This will come in handy in looking at our results
g2s <- genes %>% as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

samplesheet <- read_csv("rnaseq/samplesheet.csv")

# We of course need to swap the labels for total and insoluble cytoplasmic
# From looking in the browser at genes that should be nuclear: NEAT1 for ex.
# We have concluded that 
# hepg2_R2 -- is whole cell / total
# hepg2_insoluble_cytoplasmic_fraction_R2 -- is whole cell / total
# hepg2_R1 -- is insoluble_cytoplasmic
# hepg2_insoluble_cytoplasmic_fraction_R1 -- is insoluble_cytoplasmic

# The cleanest thing to do would be to fix it on the design file and re-run the nf-core pipeline,
# but for now we can swap it here.
samplesheet[which(samplesheet$sample_name == "hepg2_R1"), "condition"] <- "insoluble_cytoplasmic_fraction"
samplesheet[which(samplesheet$sample_name == "hepg2_insoluble_cytoplasmic_fraction_R2"), "condition"] <- "total"
```


```{r}
# We now won't filter the sample sheet and we'll read in all the data
files <- file.path("rnaseq/results/star_salmon",
                   samplesheet$sample_name, 
                   "quant.sf")
names(files) <- samplesheet$sample_name

txi <- tximport(files, type = c("salmon"), tx2gene = tx2gene)
```


```{r}
# We'll do the same housekeeping as before to make sure the genes and the samples are in the same order as
# the count data.
names(genes) <- genes$gene_id
genes <- genes[rownames(txi$counts)]
stopifnot(all(genes$gene_id == rownames(txi$counts)))

# Create rownames for the sample sheet and then use that to re-order
samplesheet <- samplesheet %>%
  mutate(row_names = sample_name) %>%
  column_to_rownames("row_names")
samplesheet <- samplesheet[colnames(txi$counts),]
stopifnot(all(rownames(samplesheet) == colnames(txi$counts)))
```

```{r}
# In order to make the total (whole cell) condition be the one that everything is compared back to
# we will want to set the factor levels with "total" first
samplesheet$condition <- factor(samplesheet$condition, levels = c("total", "membrane_fraction", 
                                                                  "insoluble_cytoplasmic_fraction", 
                                                                  "cytosolic_fraction", "nuclear_fraction"))

dds <- DESeqDataSetFromTximport(txi,
                                design = ~ condition,
                                colData = samplesheet,
                                rowData = genes)

# Run the DESeq stats
dds <- DESeq(dds)
```

```{r}
# We now have a bunch more results
resultsNames(dds)

# Let's just look at one of the results


res <- results(dds, name = "condition_membrane_fraction_vs_total")

res_shrunken <- lfcShrink(dds, coef = "condition_membrane_fraction_vs_total",  res = res)




#doing our own stuff

res_membrane <- results(dds, name = "condition_membrane_fraction_vs_total")
res_s_membrane <- lfcShrink(dds, coef = "condition_membrane_fraction_vs_total",  res = res)

res_df_mem <- res_membrane %>% as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  merge(g2s) %>%
  mutate(result_name = "condition_membrane_fraction_vs_total")

res_ins <- results(dds, name = "condition_insoluble_cytoplasmic_fraction_vs_total")
res_s_ins <- lfcShrink(dds, coef = "condition_insoluble_cytoplasmic_fraction_vs_total",  res = res)

res_df_ins <- res_ins %>% as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  merge(g2s) %>%
  mutate(result_name = "condition_insoluble_cytoplasmic_fraction_vs_total")

res_cyt <- results(dds, name = "condition_cytosolic_fraction_vs_total")
res_s_cyt <- lfcShrink(dds, coef = "condition_cytosolic_fraction_vs_total",  res = res)

res_df_cyt <- res_cyt %>% as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  merge(g2s) %>%
  mutate(result_name = "condition_cytosolic_fraction_vs_total")

res_nuc <- results(dds, name = "condition_nuclear_fraction_vs_total")
res_s_nuc <- lfcShrink(dds, coef = "condition_nuclear_fraction_vs_total",  res = res)

res_df_nuc <- res_nuc %>% as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  merge(g2s) %>%
  mutate(result_name = "condition_nuclear_fraction_vs_total")



res_df <- rbind(res_df_mem,res_df_ins,res_df_cyt,res_df_nuc)


lfc_matrix <- res_df %>% dplyr::select(gene_id, log2FoldChange, result_name) %>% pivot_wider(names_from = "result_name", values_from = "log2FoldChange") 
lfc_matrix <- na.omit(lfc_matrix)
lfc_matrix <- lfc_matrix %>% column_to_rownames("gene_id") %>% as.matrix()
#Heatmap plot


Heatmap(lfc_matrix,show_row_names = FALSE,show_column_names = TRUE)
