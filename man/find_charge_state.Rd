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
  max_charge = 4L,
  n_fwd = 15L,
  n_bwd = 20L,
  ms_lev = 1L,
  is_dda = TRUE,
  tol = 8e-06
)
}
\arguments{
\item{mass}{The mass of a peak for questing its charge state.}

\item{imax}{The index of the peak, which is often the most intense.}

\item{mint}{Not yet used. The intensity of the peak at \code{imax}.}

\item{moverzs}{MS1 or MS2 Mass-to-charge ratios. The inputs are typically at
weighted-mean statistics.}

\item{msxints}{MS1 or MS2 peak intensities. The inputs are typically at mean
statistics.}

\item{n_ms1s}{The underlying counts that have contributed to \code{moverzs}
and \code{maxints}.}

\item{max_charge}{The maximum charge state.}

\item{n_fwd}{Forward looking up to \code{n_fwd} mass entries when determining
the charge state.}

\item{n_bwd}{Backward looking up to \code{n_bwd} mass entries when
determining the charge state.}

\item{ms_lev}{MS level.}

\item{is_dda}{Logical; is DDA or not.}

\item{tol}{The tolerance in mass errors.}
}
\description{
Finds charge state.
}
