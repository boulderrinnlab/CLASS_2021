# This file will eventually knit all the RMarkdowns
library(rmarkdown)
rmarkdown::render(file.path("01_global_peak_properties",
                            "01_global_peak_properties.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("02_content_of_peaks_per_promoter/",
                            "02_content_of_peaks_per_promoter_OL.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("03_global_clustering/",
                            "03_gloabal_clustering.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("04_metaplot/",
                            "04_metaplots.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("05_RNA-seq_expression/",
                            "05_01_Retrieving_RNA-seq_data.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("05_RNA-seq_expression/",
                            "05_02_RNA-seq_differential_expression.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("06_HEPG2_binding_versus_expression/",
                            "06_HEPG2_binding_vs_expression_OL.Rmd"), 
                  md_document(variant = "markdown_github"))



rmarkdown::render(file.path("08_genebody_vs_promoter/",
                            "08_01_global_peak_properties_genebody_NR.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("08_genebody_vs_promoter/",
                            "08_03_high_binders_promoter_vs_genebody_NR.Rmd"), 
                  md_document(variant = "markdown_github"))
