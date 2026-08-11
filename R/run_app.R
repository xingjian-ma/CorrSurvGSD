#' Run the CorrSurvGSD Shiny application.
#'
#' @param ... Arguments passed to [shiny::runApp()].
#'
#' @return The result returned by [shiny::runApp()].
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The Shiny application requires the 'shiny' package.")
  }

  app_dir <- system.file("shiny", package = "CorrSurvGSD")
  if (!nzchar(app_dir)) {
    stop("The CorrSurvGSD Shiny application files were not found.")
  }

  shiny::runApp(appDir = app_dir, ...)
}
