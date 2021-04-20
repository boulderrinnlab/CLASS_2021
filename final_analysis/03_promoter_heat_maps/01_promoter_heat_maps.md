read.c--- title: "03\_promoter\_heat\_maps" author: "James Pratt" date: "4/14/2021" output: html\_document --- Determine the binding patterns of DNA binding proteins that binds to each gene on a gene by gene basis. Helpful to determine where DBPS are binding when you have a high number of binders (&gt;200).

Read in gencode annotations as Granges object, read in consensus peaks

``` r
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/data/gencode.v32.annotation.gtf")

fl <- list.files("/scratch/Shares/rinnclass/data/consensus_peaks", full.names = TRUE)
consensus_peaks <- lapply(fl, rtracklayer::import)
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/data/consensus_peaks/|.bed", "", fl)

num_peaks_threshold <- 250

filtered_consensus_peaks <- consensus_peaks[sapply(consensus_peaks, length) > num_peaks_threshold]
```

Create promoters for all genes - here defined as 3kb up and downstream of called peak.

``` r
genes <- gencode_gr[gencode_gr$type == "gene"]
all_promoters_gr <- promoters(genes, upstream = 3e3, downstream = 3e3)
```

To make the actual promoter heat map (as in, show peaks for all DNA binding proteins that have peaks called for a genes annotated promoter).

``` r
# Example gene promoter heat map, can search for any gene name that is in the all_promoters_gr - use the num_peaks_df to determine genes of interest. 
ZNF470_promoter <- all_promoters_gr[all_promoters_gr$gene_name == "ZNF470"]
ZNF470_promoter_matrix <- make_promoter_binding_matrix(filtered_consensus_peaks, ZNF470_promoter)
plot_promoter_peak_matrix(ZNF470_promoter_matrix, gene_name = "ZNF470", save_pdf = TRUE)
```

    ## RStudioGD 
    ##         2