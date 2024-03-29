% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{post_pepfdr}
\alias{post_pepfdr}
\title{Post processing after peptide FDR.}
\usage{
post_pepfdr(prob_cos = NULL, n_13c = 0L, out_path = NULL, fct_score = 5L)
}
\arguments{
\item{prob_cos}{Probability cut-offs (in data frame).}

\item{n_13c}{Number(s) of 13C off-sets in precursor masses, for example, over
the range of \code{-1:2}. The default is 0.}

\item{out_path}{An output path.}

\item{fct_score}{A factor to convert p-values to scores.}
}
\description{
Post processing after peptide FDR.
}
