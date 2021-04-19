# analysis.R

library(rmarkdown)
rmarkdown::render(file.path("00_consensus_peaks",
                            "00_consensus_peaks.Rmd"), 
                  md_document(variant = "markdown_github"))

rmarkdown::render(file.path("01_global_peak_properties",
                            "01_global_peak_properties.Rmd"), 
                  md_document(variant = "markdown_github"))
