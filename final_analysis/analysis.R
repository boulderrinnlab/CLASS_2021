# analysis.R

library(rmarkdown)
rmarkdown::render(file.path("00_consensus_peaks",
                            "00_consensus_peaks.Rmd"), 
                  md_document(variant = "markdown_github"))
