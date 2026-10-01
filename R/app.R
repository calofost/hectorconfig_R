# Shiny checker ---------------------------------------------------------------

#' Create the HECTOR checker Shiny application
#'
#' @param constants Optional named list of constants.
#' @param constants_path Optional YAML constants path.
#' @return A Shiny application object.
#' @export
hector_checker_app <- function(constants = NULL, constants_path = NULL) {
  engine <- .hector_new_engine(
    constants = constants,
    constants_path = constants_path,
    visualise = FALSE
  )
  app_file <- system.file("shiny", "hector_checker_app.R",
                          package = "hectorconfig")
  if (!nzchar(app_file)) {
    stop("The HECTOR checker app could not be located.", call. = FALSE)
  }
  sys.source(app_file, envir = engine)
  engine$app
}

#' Run the HECTOR checker Shiny application
#'
#' @inheritParams hector_checker_app
#' @param ... Additional arguments passed to `shiny::runApp()`.
#' @export
hector_run_app <- function(constants = NULL, constants_path = NULL, ...) {
  shiny::runApp(
    hector_checker_app(constants = constants, constants_path = constants_path),
    ...
  )
}
