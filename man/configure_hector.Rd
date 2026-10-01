% Generated manually for the initial package scaffold.
\name{configure_hector}
\alias{configure_hector}
\title{Configure a HECTOR field}
\description{Reads distortion-corrected HECTOR tile and guide tables, configures galaxy, standard-star, and guide probes, and optionally writes configured tables and a diagnostic plot.}
\usage{
configure_hector(tile_file, guide_file, hexafile_out = NULL,
  guidefile_out = NULL, plot_file = NULL, constants = NULL,
  constants_path = NULL, visualise = FALSE, plot_pause = 0.1,
  plot_device = c("current", "quartz"))
}
\arguments{
  \item{tile_file}{Path to the distortion-corrected hexabundle CSV.}
  \item{guide_file}{Path to the distortion-corrected guide-star CSV.}
  \item{hexafile_out}{Optional output path for the configured hexabundle CSV.}
  \item{guidefile_out}{Optional output path for the configured guide CSV.}
  \item{plot_file}{Optional output path for a PDF diagnostic plot.}
  \item{constants}{Optional named list of constants.}
  \item{constants_path}{Optional YAML constants path.}
  \item{visualise}{Whether to draw intermediate plots during configuration.}
  \item{plot_pause}{Seconds to pause after each visual redraw.}
  \item{plot_device}{Graphics device for intermediate plots: \code{"current"}
  uses the active device; \code{"quartz"} opens a separate macOS Quartz
  window.}
}
\value{A list containing the configuration, output tables, prepared inputs, and private engine environment.}
\seealso{\code{hector_read_tables}, \code{hector_checker_app}}
