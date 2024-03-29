% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deisotope.R
\name{find_charge_state}
\alias{find_charge_state}
\title{Finds charge state.}
\usage{
find_charge_state(
  mass,
  imax,
  mint,
  moverzs,
  msxints,
  n_ms1s,
  lenm,
  max_charge = 4L,
  n_fwd = 20L,
  ms_lev = 1L,
  is_dda = TRUE,
  ppm = 8
)
}
\arguments{
\item{mass}{The current mass.}

\item{imax}{The index of the most intense peak.}

\item{mint}{The maximum intensity (corresponding to the peak at \code{imax}).}

\item{moverzs}{MS1 or MS2 Mass-to-charge ratios. The inputs are typically at
weighted-mean statistics.}

\item{msxints}{MS1 or MS2 peak intensities. The inputs are typically at mean
statistics.}

\item{n_ms1s}{The underlying counts that have contributed to \code{moverzs}
and \code{maxints}.}

\item{lenm}{The length of \code{moverzs}.}

\item{max_charge}{The maximum charge state.}

\item{n_fwd}{Forward looking up to \code{n_fwd} mass entries. The default is
20 for MS1 and 10 for MS2.}

\item{ms_lev}{MS level.}

\item{is_dda}{Logical; is DDA or not.}

\item{ppm}{Allowance in mass error when deisotoping.}
}
\description{
Finds charge state.
}
