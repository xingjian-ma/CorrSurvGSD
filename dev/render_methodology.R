options(
  rmarkdown.html_dependency.header_attr = FALSE,
  rmarkdown.html_vignette.check_title = FALSE
)

rmarkdown::render(
  "vignettes/methodology.Rmd",
  output_dir = "inst/shiny/www/methodology",
  quiet = TRUE
)
