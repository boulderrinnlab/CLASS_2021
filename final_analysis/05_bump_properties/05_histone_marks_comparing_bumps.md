``` r
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/data/gencode.v32.annotation.gtf")

fl <- list.files("/scratch/Shares/rinnclass/data/consensus_peaks", full.names = TRUE)
consensus_peaks <- lapply(fl, rtracklayer::import)
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/data/consensus_peaks/|.bed", "", fl)
```

defining first and second bumps from plot of \# of genes vs \# of dbps bound in 02
==================================================================================

Read in the peak occurrence matrix and filter out the promoters in bump 1 (defined as 1-46 DBPs) and bump 2 (defined as 230-352 DBPs)

``` r
peak_occurrence_df <- read_csv('../01_global_peak_properties/results/peak_occurence_dataframe.csv')
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   gene_id = col_character(),
    ##   gene_name = col_character(),
    ##   gene_type = col_character(),
    ##   chr = col_character(),
    ##   X3kb_up_tss_start = col_double(),
    ##   strand = col_character(),
    ##   number_of_dbp = col_double()
    ## )

``` r
bump1_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp < 46)
bump1_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp > 1)

bump2_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp > 230)
bump2_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp < 352)
```

make promoter windows
=====================

``` r
#pull out gene annotations
genes <- gencode_gr[gencode_gr$type == "gene"]
# pull out TSS locations and define promoter window as +-3kb (function built into GenomicRanges)
all_promoters <- promoters(genes, upstream = 3e3, downstream = 3e3)

#mRNA genes
mrna_genes <- genes[genes$gene_type == "protein_coding"]
mrna_promoters <- promoters(mrna_genes, upstream = 3e3, downstream = 3e3)

# lncRNA promoter profiles
lncrna_genes <- genes[genes$gene_type == "lncRNA"]
lncrna_promoters <- promoters(lncrna_genes, upstream = 3e3, downstream = 3e3)

#first bump promoter profiles
#turn bump1_promoters back into GRanges
list_bump1_genes <- bump1_promoters$gene_id
bump1_genes_gr <- genes[genes$gene_id %in% list_bump1_genes]
#check that the genes in GRanges are the same as in the bump1 list
identical(bump1_genes_gr$gene_id, list_bump1_genes)
```

    ## [1] TRUE

``` r
#make windows
bump1_promoters_gr <- promoters(bump1_genes_gr, upstream = 3e3, downstream = 3e3)

#second bump promoter profiles
list_bump2_genes <- bump2_promoters$gene_id
bump2_genes_gr <- genes[genes$gene_id %in% list_bump2_genes]
#check that the genes in GRanges are the same as in the bump2 list
identical(bump2_genes_gr$gene_id, list_bump2_genes)
```

    ## [1] TRUE

``` r
#make windows
bump2_promoters_gr <- promoters(bump2_genes_gr, upstream = 3e3, downstream = 3e3)
```

define profile\_TSS function
============================

``` r
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
```

use profile\_tss for each histone mark
======================================

``` r
#H4K20me1 is associated with transcriptional activation. The most highly transcribed group of genes tend to have H4K20me1 present in addition to the core group of modifications at active promoters (Wang et al., 2008).

H4K20me1_allGenes_tss_profile <- profile_tss(consensus_peaks[["H4K20me1"]], all_promoters)

#plot all genes
ggplot(H4K20me1_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H4K20me1 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H4K20me1%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H4K20me1_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H4K20me1_mRNA_tss_profile <- profile_tss(consensus_peaks[["H4K20me1"]], mrna_promoters)
H4K20me1_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H4K20me1"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H4K20me1_mRNA_tss_profile$gene_type <- "mRNA"
H4K20me1_lncRNA_tss_profile$gene_type <- "lncRNA"
H4K20me1_combined_metaplot_profile <- bind_rows(H4K20me1_mRNA_tss_profile, H4K20me1_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H4K20me1_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H4K20me1 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H4K20me1%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H4K20me1_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H4K20me1_bump1_tss_profile <- profile_tss(consensus_peaks[["H4K20me1"]], bump1_promoters_gr)
H4K20me1_bump2_tss_profile <- profile_tss(consensus_peaks[["H4K20me1"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H4K20me1_bump1_tss_profile$gene_type <- "bump_1"
H4K20me1_bump2_tss_profile$gene_type <- "bump_2"
H4K20me1_bumps_metaplot_profile <- bind_rows(H4K20me1_bump1_tss_profile, H4K20me1_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H4K20me1_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H4K20me1 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H4K20me1%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H4K20me1_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repression. indicates the tri-methylation at the 9th lysine residue of the histone H3 protein and is often associated with heterochromatin.

H3K9me3_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K9me3"]], all_promoters)

#plot all genes
ggplot(H3K9me3_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9me3 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9me3%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K9me3_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K9me3_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K9me3"]], mrna_promoters)
H3K9me3_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K9me3"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K9me3_mRNA_tss_profile$gene_type <- "mRNA"
H3K9me3_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K9me3_combined_metaplot_profile <- bind_rows(H3K9me3_mRNA_tss_profile, H3K9me3_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K9me3_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9me3 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9me3%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K9me3_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K9me3_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K9me3"]], bump1_promoters_gr)
H3K9me3_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K9me3"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K9me3_bump1_tss_profile$gene_type <- "bump_1"
H3K9me3_bump2_tss_profile$gene_type <- "bump_2"
H3K9me3_bumps_metaplot_profile <- bind_rows(H3K9me3_bump1_tss_profile, H3K9me3_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K9me3_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9me3 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9me3%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K9me3_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#Activated if this mark is acetylated (and silences them if methylated, above)

H3K9ac_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K9ac"]], all_promoters)

#plot all genes
ggplot(H3K9ac_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9ac Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9ac%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K9ac_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K9ac_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K9ac"]], mrna_promoters)
H3K9ac_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K9ac"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K9ac_mRNA_tss_profile$gene_type <- "mRNA"
H3K9ac_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K9ac_combined_metaplot_profile <- bind_rows(H3K9ac_mRNA_tss_profile, H3K9ac_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K9ac_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9ac Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9ac%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K9ac_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K9ac_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K9ac"]], bump1_promoters_gr)
H3K9ac_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K9ac"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K9ac_bump1_tss_profile$gene_type <- "bump_1"
H3K9ac_bump2_tss_profile$gene_type <- "bump_2"
H3K9ac_bumps_metaplot_profile <- bind_rows(H3K9ac_bump1_tss_profile, H3K9ac_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K9ac_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K9ac Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K9ac%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K9ac_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#marks promoters, and distinguishes from enhancers. If it is H3K4me1, the region is an enhancer, and if it is H3K4me3, the region is a promoter

H3K4me3_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K4me3"]], all_promoters)

#plot all genes
ggplot(H3K4me3_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me3 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me3%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K4me3_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K4me3_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me3"]], mrna_promoters)
H3K4me3_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me3"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me3_mRNA_tss_profile$gene_type <- "mRNA"
H3K4me3_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K4me3_combined_metaplot_profile <- bind_rows(H3K4me3_mRNA_tss_profile, H3K4me3_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me3_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me3 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me3%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K4me3_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K4me3_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K4me3"]], bump1_promoters_gr)
H3K4me3_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K4me3"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me3_bump1_tss_profile$gene_type <- "bump_1"
H3K4me3_bump2_tss_profile$gene_type <- "bump_2"
H3K4me3_bumps_metaplot_profile <- bind_rows(H3K4me3_bump1_tss_profile, H3K4me3_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me3_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me3 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me3%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K4me3_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#at promoters, of transcriptionally active genes as well as genes primed for future expression during cell development in higher eukaryotes

H3K4me2_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K4me2"]], all_promoters)

#plot all genes
ggplot(H3K4me2_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me2 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me2%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K4me2_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K4me2_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me2"]], mrna_promoters)
H3K4me2_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me2"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me2_mRNA_tss_profile$gene_type <- "mRNA"
H3K4me2_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K4me2_combined_metaplot_profile <- bind_rows(H3K4me2_mRNA_tss_profile, H3K4me2_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me2_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me2 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me2%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K4me2_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K4me2_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K4me2"]], bump1_promoters_gr)
H3K4me2_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K4me2"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me2_bump1_tss_profile$gene_type <- "bump_1"
H3K4me2_bump2_tss_profile$gene_type <- "bump_2"
H3K4me2_bumps_metaplot_profile <- bind_rows(H3K4me2_bump1_tss_profile, H3K4me2_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me2_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me2 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me2%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K4me2_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#marks enhancers and distingushes from promoters. If it is H3K4me1, the region is an enhancer, and if it is H3K4me3, the region is a promoter

H3K4me1_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K4me1"]], all_promoters)

#plot all genes
ggplot(H3K4me1_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me1 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me1%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K4me1_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K4me1_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me1"]], mrna_promoters)
H3K4me1_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K4me1"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me1_mRNA_tss_profile$gene_type <- "mRNA"
H3K4me1_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K4me1_combined_metaplot_profile <- bind_rows(H3K4me1_mRNA_tss_profile, H3K4me1_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me1_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me1 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me1%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K4me1_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K4me1_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K4me1"]], bump1_promoters_gr)
H3K4me1_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K4me1"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K4me1_bump1_tss_profile$gene_type <- "bump_1"
H3K4me1_bump2_tss_profile$gene_type <- "bump_2"
H3K4me1_bumps_metaplot_profile <- bind_rows(H3K4me1_bump1_tss_profile, H3K4me1_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K4me1_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K4me1 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K4me1%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K4me1_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#in gene bodies, deposited in the wake of Pol II in yeast, could play roles in defining exons and splicing, also linked to DDR
#H3K36 is like a fine wine: complex, intriguing, and an active source of interest among researchers. https://epigenie.com/key-epigenetic-players/histone-proteins-and-modifications/histone-h3k36/
#H3K36me3 serves as the postman to send chromatin information to DNA damage repair processors. (Z Sun et al. 2020)

H3K36me3_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K36me3"]], all_promoters)

#plot all genes
ggplot(H3K36me3_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K36me3 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K36me3%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K36me3_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K36me3_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K36me3"]], mrna_promoters)
H3K36me3_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K36me3"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K36me3_mRNA_tss_profile$gene_type <- "mRNA"
H3K36me3_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K36me3_combined_metaplot_profile <- bind_rows(H3K36me3_mRNA_tss_profile, H3K36me3_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K36me3_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K36me3 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K36me3%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K36me3_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K36me3_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K36me3"]], bump1_promoters_gr)
H3K36me3_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K36me3"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K36me3_bump1_tss_profile$gene_type <- "bump_1"
H3K36me3_bump2_tss_profile$gene_type <- "bump_2"
H3K36me3_bumps_metaplot_profile <- bind_rows(H3K36me3_bump1_tss_profile, H3K36me3_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K36me3_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K36me3 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K36me3%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K36me3_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#marks silencers and cell-type specific genes (generally repressed unless in that cell type). (Cai et al. 2021)

H3K27me3_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K27me3"]], all_promoters)

#plot all genes
ggplot(H3K27me3_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27me3 Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27me3%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K27me3_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K27me3_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K27me3"]], mrna_promoters)
H3K27me3_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K27me3"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K27me3_mRNA_tss_profile$gene_type <- "mRNA"
H3K27me3_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K27me3_combined_metaplot_profile <- bind_rows(H3K27me3_mRNA_tss_profile, H3K27me3_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K27me3_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27me3 Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27me3%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K27me3_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K27me3_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K27me3"]], bump1_promoters_gr)
H3K27me3_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K27me3"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K27me3_bump1_tss_profile$gene_type <- "bump_1"
H3K27me3_bump2_tss_profile$gene_type <- "bump_2"
H3K27me3_bumps_metaplot_profile <- bind_rows(H3K27me3_bump1_tss_profile, H3K27me3_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K27me3_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27me3 Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27me3%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K27me3_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#marks active enhancers. H3K27ac Distinguishes Active from Inactive Enhancers. H3K4me1 alone marks inactive enhancers (Creyghton et al. 2010)

H3K27ac_allGenes_tss_profile <- profile_tss(consensus_peaks[["H3K27ac"]], all_promoters)

#plot all genes
ggplot(H3K27ac_allGenes_tss_profile, 
       aes(x = x, y = dens)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27ac Promoter Metaplot - All Genes") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27ac%20TSS%20profiles-1.png)

``` r
ggsave("./figures/H3K27ac_allGenes_tss_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for lncRNA and mRNA
H3K27ac_mRNA_tss_profile <- profile_tss(consensus_peaks[["H3K27ac"]], mrna_promoters)
H3K27ac_lncRNA_tss_profile <- profile_tss(consensus_peaks[["H3K27ac"]], lncrna_promoters)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K27ac_mRNA_tss_profile$gene_type <- "mRNA"
H3K27ac_lncRNA_tss_profile$gene_type <- "lncRNA"
H3K27ac_combined_metaplot_profile <- bind_rows(H3K27ac_mRNA_tss_profile, H3K27ac_lncRNA_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K27ac_combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27ac Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27ac%20TSS%20profiles-2.png)

``` r
ggsave("./figures/H3K27ac_combined_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

``` r
#repeat for first and second bump
H3K27ac_bump1_tss_profile <- profile_tss(consensus_peaks[["H3K27ac"]], bump1_promoters_gr)
H3K27ac_bump2_tss_profile <- profile_tss(consensus_peaks[["H3K27ac"]], bump2_promoters_gr)
# We can row bind these dataframes so that we can plot mRNA vs lncRNA on the same plot
H3K27ac_bump1_tss_profile$gene_type <- "bump_1"
H3K27ac_bump2_tss_profile$gene_type <- "bump_2"
H3K27ac_bumps_metaplot_profile <- bind_rows(H3K27ac_bump1_tss_profile, H3K27ac_bump2_tss_profile)

#plot mRNA and lncRNA
ggplot(H3K27ac_bumps_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type)) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(size = 1.5) + 
  ggtitle("H3K27ac Promoter Metaplot - bump 1 vs bump 2") + 
  scale_x_continuous(breaks = c(-3000, 0, 3000),
                     labels = c("-3kb", "TSS", "+3kb"),
                     name = "") + 
  ylab("Peak frequency") + 
  scale_color_manual(values = c("#424242","#a8404c"))
```

![](05_histone_marks_comparing_bumps_files/figure-markdown_github/H3K27ac%20TSS%20profiles-3.png)

``` r
ggsave("./figures/H3K27ac_bumps_metaplot_profile.pdf")
```

    ## Saving 7 x 5 in image

make summary table
==================

``` r
#make a df with all the histone marks
summary_histone_marks <- data.frame(Histone_marks = c('H4K20me1', 'H3K9me3', 'H3K9ac', 'H3K4me3', 'H3K4me2', 'H3K4me1', 'H3K36me3', 'H3K27me3', 'H3K27ac'), 
                                    mark_function = c('promoter activation', 'promoter repression', 'promoter activation', 'identifies promoters', 'promoter activation', 'identifies enhancers', 'gene bodies', 'silencers', 'enhancer activation'), 
                                    difference_in_RNA_type = c('Yes', 'maybe', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'maybe'), 
                                    difference_in_bumps = c('No', 'Yes', 'No', 'No', 'No', 'No', 'maybe', 'maybe', 'No')
                                    )

summary_histone_marks$bumps_notes <- c('', 'present in bump2 promoters only (matching all genes) but absent at bump1 promoters', '', '', '', '', 'slightly less depleted at bump2 promoters (potentially due to smaller % of lncRNA - flat line across lncRNA TSS)', 'depleted directly at TSS for both but different patters maybe +-1kb out from TSS', '')

summary_histone_marks
```

    ##   Histone_marks        mark_function difference_in_RNA_type difference_in_bumps
    ## 1      H4K20me1  promoter activation                    Yes                  No
    ## 2       H3K9me3  promoter repression                  maybe                 Yes
    ## 3        H3K9ac  promoter activation                    Yes                  No
    ## 4       H3K4me3 identifies promoters                    Yes                  No
    ## 5       H3K4me2  promoter activation                    Yes                  No
    ## 6       H3K4me1 identifies enhancers                    Yes                  No
    ## 7      H3K36me3          gene bodies                    Yes               maybe
    ## 8      H3K27me3            silencers                    Yes               maybe
    ## 9       H3K27ac  enhancer activation                  maybe                  No
    ##                                                                                                        bumps_notes
    ## 1                                                                                                                 
    ## 2                               present in bump2 promoters only (matching all genes) but absent at bump1 promoters
    ## 3                                                                                                                 
    ## 4                                                                                                                 
    ## 5                                                                                                                 
    ## 6                                                                                                                 
    ## 7 slightly less depleted at bump2 promoters (potentially due to smaller % of lncRNA - flat line across lncRNA TSS)
    ## 8                                 depleted directly at TSS for both but different patters maybe +-1kb out from TSS
    ## 9
