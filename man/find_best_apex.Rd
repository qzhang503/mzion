% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lfq.R
\name{find_best_apex}
\alias{find_best_apex}
\title{Find the apex under a mass column}
\usage{
find_best_apex(
  xvs,
  yvs,
  xbars,
  yints,
  aps,
  rts,
  rngs,
  xref,
  scan_ms1,
  ss,
  step_tr = 6e-06,
  min_y = 0,
  min_n1 = 10L,
  min_n2 = 20L,
  min_n3 = 15L
)
}
\arguments{
\item{xvs}{X values.}

\item{yvs}{Y values.}

\item{xbars}{Averaged m-over-zs under each gate.}

\item{yints}{Peak areas under under each gate.}

\item{aps}{All apexes scan numbers under a mass column.}

\item{rts}{All apexes retention times under a mass column.}

\item{rngs}{The scan ranges of apexes under column k.}

\item{xref}{The experimental X value for tracing.}

\item{scan_ms1}{The preceding MS1 scan number of an MS2 event at \code{xref}.}

\item{ss}{All of the scan numbers.}

\item{step_tr}{A step size.}

\item{min_y}{The cut-off of intensity values.}
}
\description{
Find the apex under a mass column
}
