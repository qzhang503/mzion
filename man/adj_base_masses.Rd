% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{adj_base_masses}
\alias{adj_base_masses}
\title{Adjusts the masses for the all-fixed mode}
\usage{
adj_base_masses(
  fwd_peps_1,
  aa_masses_0,
  aa_masses_1,
  min_mass = 200L,
  max_mass = 4500L,
  digits = 4L
)
}
\arguments{
\item{fwd_peps_1}{The first list of forward peptides with masses.}

\item{aa_masses_0}{The amino-acid masses look-up table with the coercion of
fixed to variable modification.}

\item{aa_masses_1}{The amino-acid masses look-up table without the coercion
of fixed to variable modification.}

\item{min_mass}{A positive integer; the minimum precursor mass for
interrogation. The default is an arbitrarily low value (the primary guard
against low molecular-weight precursors is \code{min_len}).}

\item{max_mass}{A positive integer; the maximum precursor mass for
interrogation.}

\item{digits}{A non-negative integer; the number of decimal places to be
used. The default is 4.}
}
\description{
The \code{fwd_peps_1} was at first based on the coerced aa_masses_0;
}
