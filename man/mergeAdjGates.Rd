% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lfq.R
\name{mergeAdjGates}
\alias{mergeAdjGates}
\title{Merge adjacent gates between columns ka and kb}
\usage{
mergeAdjGates(
  xs_ka,
  ys_ka,
  xs_kb,
  ys_kb,
  rngs_ka,
  rngs_kb,
  rts_ka,
  rts_kb,
  scans_ka,
  scans_kb
)
}
\arguments{
\item{xs_ka}{A vector of X values under column ka.}

\item{ys_ka}{A vector of Y values under column ka.}

\item{xs_kb}{A vector of X values under column kb.}

\item{ys_kb}{A vector of Y values under column kb.}

\item{rngs_ka}{A vector of scan range values under column ka.}

\item{rngs_kb}{A vector of scan range values under column kb.}

\item{rts_ka}{A vector of retention time values under column ka.}

\item{rts_kb}{A vector of retention time values under column kb.}

\item{scans_ka}{A vector of scan number values under column ka.}

\item{scans_kb}{A vector of scan number values under column kb.}
}
\description{
Exact tangent gates for merging.
}
