% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mzml.R
\name{traceDIA}
\alias{traceDIA}
\title{Traces DIA-MS2 against MS1.}
\usage{
traceDIA(
  df = NULL,
  icenter = 1L,
  ws = NULL,
  n_dia_scans = 4L,
  n_dia_ms2bins = 1L,
  min_mass = 200L,
  min_ms2mass = 115L,
  step1 = 1e-05,
  step2 = 1e-05,
  temp_dir = NULL,
  join_ms = FALSE,
  spread_ohw = FALSE,
  delayed_diams2_tracing = FALSE
)
}
\arguments{
\item{df}{A data frame.}

\item{icenter}{The index of an isolation center.}

\item{ws}{Weights.}

\item{n_dia_scans}{The number of adjacent MS scans for constructing a peak
profile and thus for determining the apex scan number of an moverz value
along LC.}

\item{n_dia_ms2bins}{Allowance in adjacent DIA-MS2 bins for retention time
alignment with a DIA-MS1.}

\item{min_mass}{A positive integer; the minimum precursor mass for
interrogation. The default is an arbitrarily low value (the primary guard
against low molecular-weight precursors is \code{min_len}).}

\item{min_ms2mass}{A positive integer; the minimum MS2 mass for
interrogation. The default is 110.}

\item{step1}{A step size for MS1.}

\item{step2}{A step size for MS2.}

\item{temp_dir}{A temporary directory.}

\item{join_ms}{Logical; if TRUE, combine adjacent entries.}

\item{spread_ohw}{logical; if TRUE, spread one-hit wonders to the nearest
neighbors.}

\item{delayed_diams2_tracing}{Logical; if TRUE, delays the tracing of MS1-MS2
after ion matches.}
}
\description{
For a single MS1 isolation range, e.g., 592.5 to 604.5.
}
