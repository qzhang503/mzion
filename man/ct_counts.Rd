% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{ct_counts}
\alias{ct_counts}
\title{Counts the number of trailing residues from C-term for the replacment of "-".}
\usage{
ct_counts(max_miss = 2L)
}
\arguments{
\item{max_miss}{The maximum number of cleavages.}
}
\description{
For full-enzymes: n(i+1) = n(i) + (i+1). Not applicable for semi-enzymes.
}
