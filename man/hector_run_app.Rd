% Generated manually for the initial package scaffold.
\name{hector_run_app}
\alias{hector_run_app}
\title{Run the HECTOR checker Shiny application}
\description{Creates and launches the HECTOR checker application.}
\usage{
hector_run_app(constants = NULL, constants_path = NULL, ...)
}
\arguments{
  \item{constants}{Optional named list of constants.}
  \item{constants_path}{Optional YAML constants path.}
  \item{...}{Additional arguments passed to \code{shiny::runApp}.}
}
\value{Invisibly returns the result of the Shiny application.}
\seealso{\code{hector_checker_app}}
