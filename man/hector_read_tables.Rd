% Generated manually for the initial package scaffold.
\name{hector_read_tables}
\alias{hector_read_tables}
\title{Read HECTOR tile and guide tables}
\description{Reads and validates the two CSV tables required by the configuration workflow.}
\usage{
hector_read_tables(tile_file, guide_file)
}
\arguments{
  \item{tile_file}{Path to the hexabundle CSV.}
  \item{guide_file}{Path to the guide-star CSV.}
}
\value{A list containing \code{tile} and \code{guide} data frames.}
