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