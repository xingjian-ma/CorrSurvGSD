# run_shiny.R — launch the development Shiny application.

# Run this script with CorrSurvGSD/ as the working directory.
devtools::load_all(".")

shiny::runApp("inst/shiny")
