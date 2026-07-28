# run_shiny.R — launch the development Shiny application.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../R/utils.R")
source("../R/per_subject_moments.R")
source("../R/calendar_cutoff.R")
source("../R/joint_cor_matrix.R")
source("../R/closed_gsd_os_and_pfs.R")
source("../R/simulation.R")
source("../R/pipeline.R")

shiny::runApp("../inst/shiny")
