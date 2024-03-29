% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{find_probco_valley}
\alias{find_probco_valley}
\title{Find the pep_len that yields the lowest probability.}
\usage{
find_probco_valley(prob_cos, guess = 12L)
}
\arguments{
\item{prob_cos}{A vector of probability cut-offs.}

\item{guess}{An integer of guessed valley.}
}
\description{
Closest to 13L; favors smaller "left" if with a "right" tie.
}
