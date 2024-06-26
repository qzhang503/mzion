% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{hadd_primatches}
\alias{hadd_primatches}
\title{Helper of \link{add_primatches}}
\usage{
hadd_primatches(
  out_path = NULL,
  is_notched = FALSE,
  add_ms2theos = FALSE,
  add_ms2theos2 = FALSE,
  add_ms2moverzs = FALSE,
  add_ms2ints = FALSE
)
}
\arguments{
\item{out_path}{A file path of outputs.}

\item{is_notched}{Logical; is a search with MS1 notches or not.}

\item{add_ms2theos}{Logical. If true, adds the sequence of primary
  theoretical MS2 m/z values (\code{pep_ms2_theos}). The sequence order at a
  given \code{type_ms2ions} is:

  \tabular{ll}{ \strong{Type}   \tab \strong{Sequence}\cr by \tab \eqn{b1,
  b2..., y1, y2...} \cr ax \tab \eqn{a1, a2..., x1, x2...} \cr cz \tab
  \eqn{c1, c2..., z1, z2...} \cr }}

\item{add_ms2theos2}{Logical. If true, adds the sequence of secondary
  theoretical MS2 m/z values (\code{pep_ms2_theos2}). The sequence order at a
  given \code{type_ms2ions} is:

  \tabular{ll}{ \strong{Type}   \tab \strong{Order of sequences}\cr by \tab
  \eqn{b2, b*, b*2, b0, b02, y2, y*, y*2, y0, y02} \cr ax \tab \eqn{a2, a*,
  a*2, a0, a02, x2} \cr cz \tab \eqn{c2, z2} \cr }}

\item{add_ms2moverzs}{Logical; if TRUE, adds the sequence of experimental
\eqn{m/z} values (\code{pep_ms2_moverzs}).}

\item{add_ms2ints}{Logical; if TRUE, adds the sequence of experimental MS2
intensity values (\code{pep_ms2_ints}).}
}
\description{
Helper of \link{add_primatches}
}
