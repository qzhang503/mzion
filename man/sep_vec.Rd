% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_engine.R
\name{sep_vec}
\alias{sep_vec}
\title{Separates a vector into chunks.}
\usage{
sep_vec(vec, fold = 5L)
}
\arguments{
\item{vec}{A vector.}

\item{fold}{The number of folds.}
}
\value{
A vector of fold indexes.
}
\description{
Separates a vector into chunks.
}
\examples{
vec <- sort(rep(LETTERS[1:5], 1:5))
mzion:::sep_vec(vec, 6L)
}
