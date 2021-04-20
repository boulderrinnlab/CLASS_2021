set working directory to 01\_global\_peak\_properties

import filtered consensus peaks and annotations:

``` r
filtered_consensus_peaks_files <- list.files("../00_consensus_peaks/results/chipseq/filtered_consensus", 
                                             pattern = "*.bed",
                                             full.names = TRUE)
filtered_consensus_peaks <- lapply(filtered_consensus_peaks_files, rtracklayer::import)

names(filtered_consensus_peaks) <- gsub("../00_consensus_peaks/results/chipseq/filtered_consensus/|_filtered_consensus_peaks.bed", "", filtered_consensus_peaks_files)

# Import annotations
lncrna_mrna_promoters <- rtracklayer::import("../00_consensus_peaks/results/lncrna_mrna_promoters.gtf")
lncrna_mrna_genebody <- rtracklayer::import("../00_consensus_peaks/results/lncrna_mrna_genebody.gtf")

lncrna_gene_ids <- lncrna_mrna_genebody$gene_id[lncrna_mrna_genebody$gene_type == "lncRNA"]
mrna_gene_ids <- lncrna_mrna_genebody$gene_id[lncrna_mrna_genebody$gene_type == "protein_coding"]
```

make num\_peaks\_df to hold our data:

``` r
num_peaks_df <- data.frame("dbp" = names(filtered_consensus_peaks),
                           "num_peaks" = sapply(filtered_consensus_peaks, length))



num_peaks_df$total_peak_length <- sapply(filtered_consensus_peaks, function(x) sum(width(x)))
```

add DBP info to num\_peaks\_df this counts the number of peaks in each promoter (if there are 2 peaks for one DBP in one promoter that will count as +2)

``` r
#make peak count matrix for each promoter (+- 3kb from TSS)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_peaks, type = "counts")

#add column to num_peaks_df that has the count of peaks in each promoter
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

#add column that has the number of peaks in lncRNA promoters only
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

#add column that has the number of peaks in mRNA promoters only
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])
```

same as above but for gene bodies instead of promoters

``` r
#make peak count matrix for if a DBP binds in each gene body
genebody_peak_counts <- count_peaks_per_feature(lncrna_mrna_genebody, 
                                                filtered_consensus_peaks, 
                                                type = "counts")
#add info to num_peaks_df like above
num_peaks_df$peaks_overlapping_genebody <- rowSums(genebody_peak_counts)
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])
```

adding if the DBPs were classified as TFs in another study and adding that to num\_peaks\_df

``` r
# The human TFs
# https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx
human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/data/mmc2.xlsx",
                                sheet = 2, skip = 1)
```

    ## New names:
    ## * `` -> ...4

``` r
names(human_tfs)[4] <- "is_tf"

length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 438

``` r
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)
```

write out num\_peaks\_df

``` r
write_csv(num_peaks_df, "./results/num_peaks_df.csv")
```

make peak occurence matrix and add to num\_peaks\_df like above but using the occurance rather than the counts this is a binary output for if a DBP binds a promoter (if there are two peaks in a promoter for one DBP that will only report a 1)

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_peaks, 
                                               type = "occurrence")
# Output to promoter_peak_occurecne_matrix
write.table(promoter_peak_occurence, "./results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
# Now we want to make into a data frame using the promoter annotations as rows and attributes as columns.
# We will use lncrna_mrna_promoters to index "all promoters"
# First make sure promoter_peak_occurence and lncrna_mrna_promoters are in the same order
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))


#make a df to put our info in
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "3kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))
# This is the CSV file we will start building upon adding columns of properties as we analyze them
# The output file name will change based on what is added later, but the "peak_occurence_df" will be used throughout.
write_csv(peak_occurence_df, "./results/peak_occurence_dataframe.csv")
```

looking into super promoters (over 350 DBPs bound) and zinc fingers

``` r
#reload occurence matrix and num_peaks_df if needed and didn't run it in previous chunks
#occurence_matrix <- read.table("./results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
#num_peaks_df <- read.csv("./results/num_peaks_df.csv")

#making matrix of just super promoters
super_promoter <- occurence_matrix[, colSums(occurence_matrix) > 350]

#make df of all DBPs that bind at super promoters
super_promoter_dbps <- data.frame(dbp = rownames(super_promoter), num_super_promoter_overlaps = rowSums(super_promoter))

#add num_peaks_df info for DBPs that bind at super promoters, and add ZF data
super_prom_features <- merge(super_promoter_dbps, num_peaks_df)

super_prom_features <- dplyr::select(super_prom_features, dbp, num_super_promoter_overlaps, everything() )

super_prom_features <- super_prom_features %>%
  mutate(percent_super_prom_ov = num_super_promoter_overlaps/peaks_overlapping_promoters * 100, 
         percent_of_total_peaks_at_super_prom = num_super_promoter_overlaps/num_peaks * 100,
         zincfinger = ifelse(grepl(" ZF",dbd), "ZF", "Other"))

#check how many of the DBPs in our dataset are 
table(super_prom_features$zincfinger)
```

    ## 
    ## Other    ZF 
    ##   261   199

``` r
#plot likelyhood of ZFs to be in super pomoters vs other DBPs
ggplot(super_prom_features, aes(x = percent_super_prom_ov, color = zincfinger)) +
  geom_density()
```

![](01_global_peak_properties_files/figure-markdown_github/super%20promoters%20and%20ZFs-1.png)

``` r
#ploting percent of total peaks at super-prom vs percent super prom ov
ggplot(super_prom_features, aes(x = percent_of_total_peaks_at_super_prom, y = percent_super_prom_ov)) + 
  geom_point()
```

![](01_global_peak_properties_files/figure-markdown_github/super%20promoters%20and%20ZFs-2.png)

``` r
write.csv(super_prom_features, "./results/super_promoter_features.csv")
```