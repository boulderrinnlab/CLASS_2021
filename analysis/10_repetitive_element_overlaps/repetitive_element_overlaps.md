``` r
# Read in the results of the permutation tests.
# fl <- list.files("../../_rslurm_perm_overlaps/", pattern = "results_", full.names = T)
# res_list <- lapply(fl, read_rds)
# 
# res_df <- bind_rows(res_list[[1]])
# for(i in 2:length(res_list)) {
#   tmp_res_df <- bind_rows(res_list[[i]])
#   res_df <- bind_rows(res_df, tmp_res_df)
# }
# 
# write_csv(res_df, "results/repetitive_elements_permutation_results.csv")
# write_csv(res_df, "~/_rslurm_perm_overlaps/repetitive_elements_permutation_results.csv")
res_df <- read_csv("results/repetitive_elements_permutation_results.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   region = col_character(),
    ##   tf = col_character(),
    ##   pval = col_double(),
    ##   zscore = col_double(),
    ##   nperm = col_double(),
    ##   alternative = col_character(),
    ##   observed = col_double(),
    ##   permuted = col_character()
    ## )

``` r
hist(res_df$pval)
```

![](repetitive_element_overlaps_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# sig_res_df <- res_df %>% filter(pval <= min(pval))
# min(res_df$pval)
# nrow(sig_res_df) / nrow(res_df)
# Get the ERVH
# repeat_names <- as.character(unique(sig_res_df$region))
# repeat_names[grep("LTR", repeat_names)]
```

Class level repeat overlaps with DBPs
=====================================

``` r
class_ov_matrix <- res_df %>%
  # filter(!(grepl("promoters", region))) %>%
  filter(grepl("_class", region)) %>%
  dplyr::select(region, tf, zscore) %>%
  pivot_wider(names_from = tf, values_from = zscore) %>%
  column_to_rownames("region") %>%
  as.matrix()
class_ov_matrix <- class_ov_matrix[complete.cases(class_ov_matrix),]
pdf("figures/class_ov_heatmap.pdf", height = 10, width = 40)
pheatmap::pheatmap(class_ov_matrix, fontsize_col = 7)
dev.off()
```

    ## pdf 
    ##   3

``` r
png("figures/class_ov_heatmap.png", height = 10, width = 40, units = "in", res = 300)
pheatmap::pheatmap(class_ov_matrix, fontsize_col = 7)
dev.off()
```

    ## pdf 
    ##   3

![class\_heatmap](figures/class_ov_heatmap.png)

Family level repeat overlaps with DBPs
======================================

``` r
family_ov_matrix <- res_df %>%
  # filter(!(grepl("promoters", region))) %>%
  filter(grepl("_family", region)) %>%
  dplyr::select(region, tf, zscore) %>%
  pivot_wider(names_from = tf, values_from = zscore) %>%
  column_to_rownames("region") %>%
  as.matrix()
family_ov_matrix <- family_ov_matrix[complete.cases(family_ov_matrix),]
pdf("figures/family_ov_heatmap.pdf", height = 15, width = 40)
pheatmap::pheatmap(family_ov_matrix, fontsize_col = 7)
dev.off()
```

    ## pdf 
    ##   3

``` r
png("figures/family_ov_heatmap.png", height = 15, width = 40, units = "in", res = 300)
pheatmap::pheatmap(family_ov_matrix, fontsize_col = 7)
dev.off()
```

    ## pdf 
    ##   3

![family\_heatmap](figures/family_ov_heatmap.png)

``` r
# sva_enriched <- sig_res_df %>% filter(grepl("SVA", region), zscore > 0) 
# 
# 
# class_res <- sig_res_df %>% filter(grepl("class", region))
# family_res <- sig_res_df %>% filter(grepl("family", region))
# class_region_clust <- hclust(dist(class_ov_matrix))
# class_tf_clust <- hclust(dist(t(class_ov_matrix)))
# plot(class_region_clust)
```
