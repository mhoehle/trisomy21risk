#' Shiny app demoing the pkg.
#' 
#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "t21app", package = "trisomy21risk")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `trisomy21risk`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
