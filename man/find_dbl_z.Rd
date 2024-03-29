% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deisotope.R
\name{find_dbl_z}
\alias{find_dbl_z}
\title{Searches for a possible doubled charge state.}
\usage{
find_dbl_z(
  ch = 2L,
  p = 2L,
  sta,
  mass,
  moverzs,
  max_charge = 4L,
  ppm = 5L,
  f = `+`
)
}
\arguments{
\item{ch}{The initial charge state.}

\item{p}{The position in \code{moverzs}.}

\item{sta}{The position of start.}

\item{mass}{The current m-over-z.}

\item{moverzs}{A vector of m-over-z values.}

\item{max_charge}{The maximum charge state for consideration.}

\item{ppm}{Mass error tolerance.}

\item{f}{A function of \code{+} or \code{-} for forward or backward
searching.}
}
\description{
Not used.
}
