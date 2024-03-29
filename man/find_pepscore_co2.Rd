% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{find_pepscore_co2}
\alias{find_pepscore_co2}
\title{Finds the cut-off in peptide scores for a given \code{target_fdr}.}
\usage{
find_pepscore_co2(td, target_fdr = 0.01)
}
\arguments{
\item{td}{A data frame of target-decoy results at a given peptide length.}

\item{target_fdr}{A numeric; the targeted false-discovery rate (FDR) at the
levels of PSM, peptide or protein. The default is 0.01. See also argument
\code{fdr_type}.}
}
\description{
Assume log-normal for decoy scores.
}
