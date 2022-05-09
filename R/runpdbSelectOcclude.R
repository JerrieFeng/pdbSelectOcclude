#' Launch Shiny App for pdbSelectOcclude
#'
#' A function that launches the Shiny app for pdbSelectOcclude.
#' The purpose of this app is only to illustrate how a Shiny
#' app works. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' pdbSelectOcclude::runpdbSelectOcclude(pdb, "1bm8")
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runpdbSelectOcclude <- function(pdbFile, name) {

  #infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  appDir <- system.file("shiny-scripts",
                        package = "pdbSelectOcclude")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
# [END]
