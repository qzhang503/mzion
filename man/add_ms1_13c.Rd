% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{add_ms1_13c}
\alias{add_ms1_13c}
\title{Adds Carbon-13 masses.}
\usage{
add_ms1_13c(peps, n_13c = 0L, max_mass = 4500L)
}
\arguments{
\item{peps}{A named vector of peptide sequences. Sequences in names and
masses in values.}

\item{n_13c}{Number(s) of 13C off-sets in precursor masses, for example, over
the range of \code{-1:2}. The default is 0.}

\item{max_mass}{A positive integer; the maximum precursor mass for
interrogation.}
}
\description{
Not currently used.
}
