% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vmod_ms1_labels.R
\name{make_ms1_vmodsets}
\alias{make_ms1_vmodsets}
\title{Makes the sets of labels of variable modifications.}
\usage{
make_ms1_vmodsets(
  aa_masses_all = NULL,
  maxn_vmods_per_pep = 5L,
  maxn_sites_per_vmod = 3L
)
}
\arguments{
\item{aa_masses_all}{All of the amino acid lookup tables.}

\item{maxn_vmods_per_pep}{A non-negative integer; the maximum number of
\code{Anywhere} (non-terminal) variable modifications per peptide. The
default is 5.}

\item{maxn_sites_per_vmod}{A non-negative integer; the maximum number of
  combinatorial \code{Anywhere} (non-terminal) variable modifications per
  site in a peptide sequence. The default is 3.

  For instance, variable modifications of \code{Carbamyl (M)} and
  \code{Oxidation (M)} both have site \code{M}. In order to have a
  combination of two \code{Carbamyl (M)} and two \code{Oxidation (M)} being
  considered, the value of \code{maxn_sites_per_vmod} needs to be four or
  greater.}
}
\description{
No position permutation (for MS1 masses).
}
