#' @param ... passed to shiny::runApp()
#' @export
run_shiny_app <- function(...) {
  app_dir <- system.file("shiny", package = "InteracDiagnosis")
  if (app_dir == "") stop("Shiny app directory not found. Did you install the package?")
  shiny::runApp(app_dir, ...)
}











