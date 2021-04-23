The goal for today is to "cluster" the 438 DBPs to find which ones have similar binding profiles across the genome -- and those that are different.

Let's start by refreshing our memory of what all went down in 01\_global\_peak\_properties.

Briefly we did 2 really good things that we will need:

1.  We made a list of (i) all promoters (ii) lncRNA promoters only (iii) mRNA promoters
2.  We made a peak\_occurance\_dataframe (.tsv) that is comprised of columns for each dbp. Rows are the DBP and filled 1 or 0 if that DBP binds to a given promoter. Or:

Promoter 1, 2, 3, .... DBP 0, 0, 1,....

Thus, each row has ~36,000 1/0 for a given DBP. This is binary data (1/0) and not contious data. We need to consider that when clustering as these similarity measurements will be influenced by continious or binary data.

We exported this file to "results" folder and same for the GTFs. So we are all set to start clustering.

There are numerous ways to do clustering, but the basic principle is that the algorithm will start from the bottom and move up. This means it will find to samples with highest similarity and then find the next and the next.. in agglormative manner.

Really all this is a correlation (e.g. Pearson) between two vectors of ~60,000 values of 1/0.

Let's get started by clustering all promoters first, then later we can see how this changes when we seperate lncRNAs and mRNAs.

``` r
# Let's grab both the lncRNA and mRNA promoters
# promoters <- rtracklayer::import("results/lncrna_mrna_promoters.gtf")
peak_occurence_matrix <- read.table("../01_global_peak_properties/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Let's filter again to make sure we have at least 250 peaks
# This is a different filtering because we're now filtering on at least 250 peaks on 
# promoter regions rather than 250 peaks across the genome.
peak_occurence_matrix <- peak_occurence_matrix[rowSums(peak_occurence_matrix) > 250, ]

# ok now we have all the things we need to start clustering !
# First we will use hclust

# ?hclust
# bin_hier
# R is so savy with statistics that this is actually a base R function!
# There are lots of methods to cluster by let's start with binary.

# Clustering requires a distance matrix which represents each row in the matrix.
# In this case rows are DBPs. We could also cluster on the columns (promoters) but
# for now we'll start with the DBPs.
# If we just want to calculate the euclidean distance between the first two rows,
# we can use the dist function.
peak_occurence_matrix[1:2,1:10]
```

    ##      ENSG00000243485.5 ENSG00000237613.2 ENSG00000186092.6 ENSG00000238009.6
    ## ADNP                 0                 0                 0                 0
    ## AFF4                 0                 0                 0                 0
    ##      ENSG00000239945.1 ENSG00000239906.1 ENSG00000241860.7 ENSG00000241599.1
    ## ADNP                 0                 0                 0                 0
    ## AFF4                 0                 0                 0                 0
    ##      ENSG00000286448.1 ENSG00000236601.2
    ## ADNP                 0                 0
    ## AFF4                 0                 0

``` r
dist(peak_occurence_matrix[1:2,], method = "euclidean")
```

    ##          ADNP
    ## AFF4 98.74715

``` r
# Since this is a binary vector, let's use the binary distance measure.
# The binary measure only considers those elements which are not zero in both vectors.
# It is then the proportion of elements which are not shared (i.e. 0 in vector 1, and 1 in vector 2).
dist(peak_occurence_matrix[1:4,], method = "binary")
```

    ##            ADNP      AFF4     AHDC1
    ## AFF4  0.8088089                    
    ## AHDC1 0.8358010 0.6518106          
    ## AHR   0.8238244 0.2774043 0.7021701

``` r
# We can create a distance matrix for each pairwise vector comparison for the whole matrix now.
# This will result in a 460x460 matrix that is symmetric on either side of the diagonal.
peak_occurence_dist <- dist(peak_occurence_matrix, method = "binary")

# Now that we have that, we can do hierarchical clustering. In brief, this will cluster more closely
# vectors that have small distances, and will iterate that process until everything is in the same cluster.
# Hierarchical clustering with binary distance measure
bin_hier <- hclust(peak_occurence_dist, method = "complete")

plot(bin_hier)
```

![](00_global_clustering_files/figure-markdown_github/Clustering%20all%20Promoters-1.png)

``` r
# To make this bigger, we can plot it as a pdf
pdf("figures/dbp_hclust_dendro.pdf", height = 12, width = 70)
plot(bin_hier)
dev.off()
```

    ## png 
    ##   2

``` r
bin_hier$labels[bin_hier$order]
```

    ##   [1] "BRCA1"           "ZBTB33"          "NRF1"            "SUZ12"          
    ##   [5] "SP2"             "ZNF460"          "ADNP"            "ZMYM3"          
    ##   [9] "RFX3"            "POLR2A"          "POLR2AphosphoS2" "FUBP1"          
    ##  [13] "ASH2L"           "GMEB2"           "YY1"             "HCFC1"          
    ##  [17] "SIN3B"           "TBL1XR1"         "ELF1"            "SIN3A"          
    ##  [21] "LRRFIP1"         "ZNF526"          "ZC3H8"           "TCF3"           
    ##  [25] "ZNF678"          "GABPA"           "MXI1"            "RFX5"           
    ##  [29] "CHD2"            "RCOR1"           "ZNF629"          "ZNF644"         
    ##  [33] "H3K4me1"         "POGZ"            "ZNF598"          "ZNF33B"         
    ##  [37] "ZNF784"          "ZNF786"          "LIN54"           "HDAC1"          
    ##  [41] "TBX2"            "MIER3"           "HMG20B"          "ZCCHC11"        
    ##  [45] "NFKBIZ"          "KLF11"           "ZNF511"          "KDM3A"          
    ##  [49] "ZNF639"          "DMAP1"           "KAT7"            "AHR"            
    ##  [53] "MEF2D"           "TFE3"            "ZSCAN9"          "SMAD3"          
    ##  [57] "ZBTB25"          "ATF2"            "IKZF5"           "DRAP1"          
    ##  [61] "HMGXB3"          "ZNF335"          "MIER2"           "PRDM10"         
    ##  [65] "RBPJ"            "KDM4B"           "MXD4"            "ZNF205"         
    ##  [69] "ZNF580"          "THRA"            "KLF6"            "KMT2B"          
    ##  [73] "ZNF792"          "ZBED4"           "ZNF124"          "FOXP4"          
    ##  [77] "ZNF217"          "PAXIP1"          "SOX6"            "LCOR"           
    ##  [81] "LCORL"           "CEBPD"           "CBX5"            "PITX1"          
    ##  [85] "ZNF883"          "ZNF274"          "BCL3"            "ZNF547"         
    ##  [89] "ZBTB7B"          "KMT2A"           "PHF8"            "ETV5"           
    ##  [93] "FOXO1"           "TFAP4"           "ZNF331"          "ARID4B"         
    ##  [97] "SAP130"          "HNF1B"           "FOXA3"           "ARID5B"         
    ## [101] "GATAD2A"         "SOX13"           "BCL6"            "SMAD4"          
    ## [105] "ELF3"            "ERF"             "IRF2"            "ZGPAT"          
    ## [109] "KAT8"            "THAP11"          "GATAD1"          "KLF16"          
    ## [113] "RARA"            "RXRB"            "RREB1"           "NFIA"           
    ## [117] "NR2F6"           "TEAD1"           "TEAD3"           "ZNF614"         
    ## [121] "NFIL3"           "PPARG"           "HOMEZ"           "MIXL1"          
    ## [125] "CEBPA"           "CEBPG"           "HMG20A"          "RCOR2"          
    ## [129] "STAG1"           "NRL"             "ZNF280B"         "ZNF891"         
    ## [133] "HIVEP1"          "ZBTB38"          "ZNF407"          "ZNF687"         
    ## [137] "ZNF572"          "H3K4me2"         "HOXA3"           "ZNF574"         
    ## [141] "MNX1"            "UBTF"            "GABPB1"          "MXD3"           
    ## [145] "ZNF501"          "EGR1"            "PATZ1"           "KDM2A"          
    ## [149] "ZBTB26"          "ZNF318"          "ZNF768"          "THAP9"          
    ## [153] "ZNF225"          "ZNF350"          "DNMT1"           "ZBED5"          
    ## [157] "AKAP8"           "TIGD6"           "ZNF276"          "ZNF619"         
    ## [161] "TRAFD1"          "PAF1"            "ZNF326"          "ZNF652"         
    ## [165] "PHF21A"          "RXRA"            "RBAK"            "ZKSCAN5"        
    ## [169] "ZBTB49"          "ZHX3"            "ZNF608"          "GPN1"           
    ## [173] "EEA1"            "WIZ"             "ZNF710"          "ZNF776"         
    ## [177] "EED"             "ZNF788"          "NFKB2"           "ZNF709"         
    ## [181] "ZBTB44"          "HOXD1"           "KLF12"           "IRF5"           
    ## [185] "ZNF142"          "ZSCAN31"         "ZIK1"            "ZNF616"         
    ## [189] "DZIP1"           "ZNF138"          "ZNF720"          "STAT5B"         
    ## [193] "ZFP36L2"         "RERE"            "SIX1"            "MEF2A"          
    ## [197] "SP140L"          "ZNF333"          "ZNF143"          "RFXAP"          
    ## [201] "ZNF264"          "DLX6"            "ZBTB43"          "TBP"            
    ## [205] "MAX"             "MAZ"             "ZC3H4"           "HDAC2"          
    ## [209] "NFYA"            "NFYB"            "NFYC"            "NFATC3"         
    ## [213] "SOX18"           "KLF9"            "ZNF772"          "TOE1"           
    ## [217] "FOXK1"           "SIX4"            "MLX"             "ZBTB21"         
    ## [221] "PHF20"           "H3K4me3"         "H3K9ac"          "TAF1"           
    ## [225] "H3K27ac"         "POLR2AphosphoS5" "NR2C2"           "ZNF230"         
    ## [229] "ATF7"            "ESRRA"           "ZNF761"          "HNF1A"          
    ## [233] "ZNF483"          "ARNTL"           "ZNF543"          "TFDP1"          
    ## [237] "SP5"             "ZFP64"           "SMAD9"           "ZNF697"         
    ## [241] "YEATS4"          "HOXA5"           "ISL2"            "YEATS2"         
    ## [245] "TSC22D2"         "GLI4"            "GMEB1"           "HMGXB4"         
    ## [249] "MYRF"            "MYNN"            "RELA"            "ZNF12"          
    ## [253] "ZSCAN29"         "ZNF317"          "SFPQ"            "ZNF414"         
    ## [257] "GZF1"            "ZNF280D"         "ZKSCAN8"         "AFF4"           
    ## [261] "ATF6"            "RFXANK"          "ZNF296"          "ZNF800"         
    ## [265] "ZFP41"           "LBX2"            "NR0B2"           "KDM5B"          
    ## [269] "MXD1"            "ZNF777"          "NFAT5"           "ZNF607"         
    ## [273] "SAFB2"           "E2F4"            "E2F5"            "ZNF263"         
    ## [277] "HINFP"           "ZNF556"          "ZNF691"          "ZNF782"         
    ## [281] "ARID4A"          "ZNF781"          "ZNF451"          "SPEN"           
    ## [285] "ZFY"             "IKZF4"           "CBFB"            "ETV6"           
    ## [289] "TFDP2"           "ATAD3A"          "KLF13"           "ZNF550"         
    ## [293] "ZNF281"          "MTA1"            "ZNF609"          "FUBP3"          
    ## [297] "ZNF548"          "ZNF510"          "ZZZ3"            "ATF5"           
    ## [301] "ZNF34"           "ZNF343"          "ZBTB7A"          "MATR3"          
    ## [305] "ZBTB14"          "ZNF740"          "GTF3A"           "ZNF660"         
    ## [309] "DMTF1"           "IRF1"            "ZBTB24"          "ZSCAN22"        
    ## [313] "TIGD3"           "ZNF224"          "ZNF512"          "ZNF512B"        
    ## [317] "NFE2L1"          "ZNF83"           "ZNF839"          "SMAD1"          
    ## [321] "ZNF713"          "ZNF737"          "ZNF576"          "SNAI1"          
    ## [325] "SNAPC4"          "PRDM15"          "ZNF766"          "CREB3"          
    ## [329] "SATB2"           "ZBTB42"          "ZNF557"          "ZNF232"         
    ## [333] "ZNF25"           "ZNF256"          "MAF1"            "E2F1"           
    ## [337] "ZSCAN20"         "NKX3-1"          "MED1"            "MED13"          
    ## [341] "ZNF503"          "ELK1"            "ZFP37"           "ZNF485"         
    ## [345] "ZNF44"           "ZNF446"          "ZNF362"          "ZNF384"         
    ## [349] "NR2F1"           "THRB"            "NR5A1"           "KDM6A"          
    ## [353] "DPF2"            "HNF4A"           "ONECUT1"         "ONECUT2"        
    ## [357] "AHDC1"           "GATA2"           "MAFF"            "MAFK"           
    ## [361] "SMC3"            "CTCF"            "RAD21"           "MAFG"           
    ## [365] "NFE2L2"          "BHLHE40"         "ATF3"            "USF1"           
    ## [369] "USF2"            "BRD4"            "GTF2F1"          "ETS1"           
    ## [373] "ZKSCAN1"         "ZBTB39"          "FOXJ3"           "TP53"           
    ## [377] "DDIT3"           "STAT6"           "ZNF337"          "ELF4"           
    ## [381] "NAIF1"           "H4K20me1"        "CREM"            "E2F2"           
    ## [385] "ZNF513"          "REST"            "ZBTB4"           "ZBTB46"         
    ## [389] "PAWR"            "ZNF383"          "RARG"            "ZNF775"         
    ## [393] "ZFP1"            "ZFP14"           "THAP8"           "ZNF234"         
    ## [397] "CERS6"           "DBP"             "PREB"            "ZFP82"          
    ## [401] "ZXDC"            "ZNF30"           "ZNF773"          "FOXM1"          
    ## [405] "ISX"             "TEF"             "CAMTA2"          "ZNF674"         
    ## [409] "ZBTB1"           "ZBTB10"          "ZNF180"          "ZNF549"         
    ## [413] "KIAA2018"        "ZNF569"          "ATRX"            "MTA3"           
    ## [417] "NR3C1"           "HSF2"            "ZNF367"          "XBP1"           
    ## [421] "AKNA"            "MTF1"            "ZNF615"          "ZNF490"         
    ## [425] "ZNF570"          "ZBTB8A"          "ZNF101"          "ZNF20"          
    ## [429] "ZNF778"          "ZNF703"          "PLSCR1"          "ZNF530"         
    ## [433] "NFIC"            "SP1"             "SRF"             "JUN"            
    ## [437] "KDM1A"           "CEBPB"           "JUND"            "ARID3A"         
    ## [441] "FOSL2"           "FOXP1"           "HNF4G"           "FOXA1"          
    ## [445] "FOXA2"           "H3K9me3"         "EZH2"            "MTF2"

Nice so we now have object "bin\_hier" that will contain all the cross correlation values across all the samples.

Now we get to PLOT !! We will use ggdenro package that will plot the branch lengths that indicate how similar two samples are -- the longer the branch the more different they are -- closer more similar.

Here is an alternative way to plot it:

``` r
# First we are going to order the samples by the correlation matrix we just made.

# You can also use ggdendro to plot with ggplot syntax -- which looks nice but for this many 
# DBPs, it's a little bit easier to visualize with base R.

# Now we can plot
ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3, 
                        theme_dendro = TRUE) +
   coord_flip() +
    scale_y_continuous() +
    scale_x_continuous(position = "top") +
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
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.
    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

![](00_global_clustering_files/figure-markdown_github/Plot%20CLuster%20using%20ggdendro-1.png)

``` r
# Go over the plotting aspects as we discussed in 11_intro_plotting_R.Rmd
# Especially thje scale_x_continious part -- we saw that in the 11_ class 

# So this is very difficult to see since we have SO much data :)
# one thing you can see is "phospho S5" and "phospho S2".
# These are modifications on Pol II that indicate two different states:
# S2 = distal or elongating Pol II so less on promoter more in genebody
# S5 = Indicative of Pol II at promoter proximilar regions.
# Note how these don't cluster anywhere near each other -- that is good news!


# Let's save this as a figure -- we want to keep organized so we made a ChIPseq
# folder and figures inside that folder. Let's put the figure there:

ggsave("figures/all_hclust_binary_dist.pdf", height = 44, width = 5.0)

# Nice, this will make it much easier to read :)
# use file transfer to bring to desktop and take a look!
# This is looking biologically relevant! Notice Pol II phospho 2 is
# similar to all other compents of Pol II -- that is exactly what we expect.
# also Pol II S5 is associated with H3K27ac which is a mark of active promoters
# again what we would expect of promoter proximal Pol II.
# MED DBPs all near each other too!
# Sweet -- with a simple spot check we can see all of this data is making sense!
```

Here's some code to cluster mRNA promoters and lncRNA promoters separately. They end up looking pretty similar.

``` r
# # Let's start with lncRNA promoters -- why not?
# 
# # We first need to load the .GTF we made for lncRNAs and mRNAs.
# 
lncrna_promoters <- rtracklayer::import("../00_consensus_peaks/results/lncrna_promoters.gtf")
mrna_promoters <- rtracklayer::import("../00_consensus_peaks/results/mrna_promoters.gtf")
 
# Cool now we can do the same as above:
 
lncrna_peak_occurence <- peak_occurence_matrix[,lncrna_promoters$gene_id]
 
bin_hier_lncrna <- hclust(dist(lncrna_peak_occurence, method = "binary"))
 
ggdendro::ggdendrogram(bin_hier_lncrna, rotate = TRUE,  size = 3)
```

![](00_global_clustering_files/figure-markdown_github/Cluster%20mRNA%20promoters%20and%20lncRNA%20promoters%20separately-1.png)

``` r
ggsave("figures/lncrna_hclust_binary_dist.pdf", height = 44, width = 6)
 
#Cool we can see the Pol II S2 and S5 diff right away again!
 
 
#Now for mRNA
 
mrna_peak_occurence <- peak_occurence_matrix[,mrna_promoters$gene_id]
 
bin_hier_mrna <- hclust(dist(lncrna_peak_occurence, method = "binary"))
 
ggdendro::ggdendrogram(bin_hier, rotate = TRUE,  size = 3)
```

![](00_global_clustering_files/figure-markdown_github/Cluster%20mRNA%20promoters%20and%20lncRNA%20promoters%20separately-2.png)

``` r
ggsave("figures/mrna_hclust_binary_dist.pdf", height = 44, width = 6)


# Do a file transfer and open up lncRNA and mRNA clusters side by side:
# They are almost identical !! That is another good sign.
```
