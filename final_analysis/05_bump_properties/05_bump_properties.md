Read in the peak occurrence matrix.

``` r
num_peaks_df <- read_csv('../01_global_peak_properties/results/num_peaks_df.csv')
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   dbp = col_character(),
    ##   num_peaks = col_double(),
    ##   total_peak_length = col_double(),
    ##   peaks_overlapping_promoters = col_double(),
    ##   peaks_overlapping_lncrna_promoters = col_double(),
    ##   peaks_overlapping_mrna_promoters = col_double(),
    ##   peaks_overlapping_genebody = col_double(),
    ##   peaks_overlapping_lncrna_genebody = col_double(),
    ##   peaks_overlapping_mrna_genebody = col_double(),
    ##   ensembl_id = col_character(),
    ##   dbd = col_character(),
    ##   tf = col_character()
    ## )

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
peak_occurrence_lncrna <- filter(peak_occurrence_df, gene_type == "lncRNA")
peak_occurrence_mrna <- filter(peak_occurrence_df, gene_type == "protein_coding")
```

Now filter out the promoters in bump 1 (defined as 1-50 DBPs) and bump 2 (defined as 250-350 DBPs)

``` r
bump1_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp < 46)
bump1_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp > 1)

bump2_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp > 230)
bump2_promoters <- filter(peak_occurrence_df, peak_occurrence_df$number_of_dbp < 352)
```

Determine the percent of mRNA and lncRNA in bump 1 and bump2.

``` r
#percent of mRNA and lncRNA bump 1
bump1_mrna_promoters = filter(bump1_promoters, bump1_promoters$gene_type == 'protein_coding')
bump1_lncrna_promoters = filter(bump1_promoters, bump1_promoters$gene_type == 'lncRNA')

n_bump1_mrna_promoters <- nrow(bump1_mrna_promoters)
n_bump1_lncrna_promoters <- nrow(bump1_lncrna_promoters)

percent_mrna_bump1 <- n_bump1_mrna_promoters/(n_bump1_mrna_promoters+n_bump1_lncrna_promoters)
percent_lncrna_bump1 <- n_bump1_lncrna_promoters/(n_bump1_mrna_promoters+n_bump1_lncrna_promoters)

print(paste('The percent mRNA in bump1 (low binders) is ',percent_mrna_bump1*100))
```

    ## [1] "The percent mRNA in bump1 (low binders) is  59.5035460992908"

``` r
print(paste('The percent lncRNA in bump1 (low binders) is ',percent_lncrna_bump1*100))
```

    ## [1] "The percent lncRNA in bump1 (low binders) is  40.4964539007092"

``` r
#percent of mRNA bump 1
bump2_mrna_promoters = filter(bump2_promoters, bump2_promoters$gene_type == 'protein_coding')
bump2_lncrna_promoters = filter(bump2_promoters, bump2_promoters$gene_type == 'lncRNA')

n_bump2_mrna_promoters <- nrow(bump2_mrna_promoters)
n_bump2_lncrna_promoters <- nrow(bump2_lncrna_promoters)

percent_mrna_bump2 <- n_bump2_mrna_promoters/(n_bump2_mrna_promoters+n_bump2_lncrna_promoters)
percent_lncrna_bump2 <- n_bump2_lncrna_promoters/(n_bump2_mrna_promoters+n_bump2_lncrna_promoters)


print(paste('The percent mRNA in bump2 (high binders) is ',percent_mrna_bump2*100))
```

    ## [1] "The percent mRNA in bump2 (high binders) is  53.440276127426"

``` r
print(paste('The percent lncRNA in bump2 (high binders) is ',percent_lncrna_bump2*100))
```

    ## [1] "The percent lncRNA in bump2 (high binders) is  46.559723872574"

Let's explore gene length!

``` r
#read out the promoter gene ids for each bump
bump1_gene_ids <- bump1_promoters$gene_id %>% as.data.frame()
bump2_gene_ids <- bump2_promoters$gene_id %>% as.data.frame()

bump1_gene_ids$gene_id <- bump1_promoters$gene_id 
bump2_gene_ids$gene_id <- bump2_promoters$gene_id 


#read in all the gencode data and make it a data frame! 
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/data/gencode.v32.annotation.gtf")
gene_df <- as.data.frame(gencode_gr)

#extract some useful information (we want gene width!)
gencode_length <- data.frame("gene_id" = gencode_gr@elementMetadata$gene_id,
                           "gene_type" = gencode_gr@elementMetadata$gene_type,
                           "type" = gencode_gr@elementMetadata$type,
                           "start" = gencode_gr@ranges@start,
                           "width" = gencode_gr@ranges@width)

#let's just look at genes (no exons or anything else!)                     
gene_length <- filter(gencode_length, type == "gene")

#filter based on bump1 or bump2
bump1_gene_lengths <- semi_join(gene_length, bump1_gene_ids, by = "gene_id")
bump2_gene_lengths <- semi_join(gene_length, bump2_gene_ids, by = "gene_id")

bump1_gene_lengths$bump <- 'blue'
bump2_gene_lengths$bump <- 'red'

both_bumps_df <- rbind(bump1_gene_lengths,bump2_gene_lengths)
```

Now let's look at the average length of genes for promoters in the two bumps.

``` r
#calculate average gene length for the promoters in bump 1 and 2
av_width_b1 <- mean(bump1_gene_lengths[,"width"])
av_width_b2 <- mean(bump2_gene_lengths[,"width"])

print(paste('The average length of promoter in bump1 (low binders) is ',av_width_b1))
```

    ## [1] "The average length of promoter in bump1 (low binders) is  54181.2397163121"

``` r
print(paste('The average length of promoter in bump2 (high binders) is ',av_width_b2))
```

    ## [1] "The average length of promoter in bump2 (high binders) is  52119.2027669326"

Some cool histogram of gene length in the two bumps.

``` r
#calculate average gene length for the promoters in bump 1 and 2

widths_b1 <- bump1_gene_lengths$width
widths_b2 <- bump2_gene_lengths$width

hist(widths_b1,breaks = 200, main='Gene Lengths Bump 1')
```

![](05_bump_properties_files/figure-markdown_github/gene%20length%20histograms-1.png)

``` r
hist(widths_b2,breaks = 200,main='Gene Lengths Bump 2')
```

![](05_bump_properties_files/figure-markdown_github/gene%20length%20more%20histograms-1.png)
