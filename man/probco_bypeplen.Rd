% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{probco_bypeplen}
\alias{probco_bypeplen}
\title{Helper of \link{calc_pepfdr}.}
\usage{
probco_bypeplen(
  len,
  td,
  fdr_type = "protein",
  target_fdr = 0.01,
  min_pepscores_co = 30,
  fct_score = 5,
  is_long_len = FALSE,
  sco_long_len = 30,
  ref_probco = NULL,
  out_path = NULL
)
}
\arguments{
\item{len}{Numeric; the length of peptides.}

\item{td}{A target-decoy pair.}

\item{fdr_type}{A character string; the type of FDR control. The value is in
  one of c("protein", "peptide", "psm"). The default is \code{protein}.

  Note that \code{fdr_type = protein} is comparable to \code{fdr_type =
  peptide} with the additional filtration of data at \code{prot_tier == 1}.}

\item{target_fdr}{A numeric; the targeted false-discovery rate (FDR) at the
levels of PSM, peptide or protein. The default is 0.01. See also argument
\code{fdr_type}.}

\item{min_pepscores_co}{A non-negative numeric; the lower limit in the
cut-offs of peptide scores for discriminating significant and insignificant
identities.}

\item{fct_score}{A factor to convert p-values to scores.}

\item{sco_long_len}{Empirical score threshold at one-percent FDR for
sequences
\eqn{\ge 35} residues.}

\item{ref_probco}{Reference cut-offs of probabilities.}

\item{out_path}{A file path of outputs.}

\item{is_long_len;}{Logical; is a long peptide length or not.}
}
\description{
Calculates the probability cut-off for target-decoy pairs at a given peptide
length.
}
