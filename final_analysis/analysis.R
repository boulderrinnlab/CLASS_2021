# analysis.R

library(rmarkdown)
rmarkdown::render(file.path("00_consensus_peaks",
                            "00_consensus_peaks.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("01_global_peak_properties",
                            "01_global_peak_properties.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("02_exploring_bumps",
                            "02_exploring_bumps.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("03_promoter_heat_maps",
                            "01_promoter_heat_maps.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("03_promoter_heat_maps",
                            "00_global_clustering.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("04_resevoir_analysis",
                            "04_01_RNAseq_Design.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("04_resevoir_analysis",
                            "04_resevoir_analysis.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("05_bump_properties",
                            "05_histone_marks_comparing_bumps.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("05_bump_properties",
                            "05_01_gene_length_by_bump.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("05_bump_properties",
                            "05_bump_properties.Rmd"), 
                  md_document(variant = "markdown_github"))