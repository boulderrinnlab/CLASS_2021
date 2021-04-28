library(rmarkdown)

rmarkdown::render("09_high_binding_properties/high_binding_properties.Rmd",
                  md_document(variant = "markdown_github"))

rmarkdown::render("06_promoter_conservation/promoter_conservation.Rmd",
                  md_document(variant = "markdown_github"))

rmarkdown::render("10_repetitive_element_overlaps/repetitive_element_overlaps.Rmd",
                  md_document(variant = "markdown_github"))


rmarkdown::render("07_conservative_reservoirs/conservative_reservoirs.Rmd",
                  md_document(variant = "markdown_github"))
                  