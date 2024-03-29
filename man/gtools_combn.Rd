% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vmod_ms1_labels.R
\name{gtools_combn}
\alias{gtools_combn}
\title{Finds combination.}
\usage{
gtools_combn(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
}
\arguments{
\item{n}{The number of balls.}

\item{r}{The number of slots.}

\item{v}{The labels of the n balls.}

\item{set}{Logical; duplicates should be removed from the source vector v or
not.}

\item{repeats.allowed}{Logical; vector v may include duplicated values or
not.}
}
\description{
From gtools::combinations
}
