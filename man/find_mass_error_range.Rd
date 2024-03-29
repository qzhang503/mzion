% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_engine.R
\name{find_mass_error_range}
\alias{find_mass_error_range}
\title{Finds the error range of a number.}
\usage{
find_mass_error_range(x = 500L, ppm = 20L)
}
\arguments{
\item{x}{A numeric value.}

\item{ppm}{Numeric; the ppm allowed from \code{x}.}
}
\value{
The lower and the upper bound to \eqn{x} by \eqn{ppm}.
}
\description{
Assumes \eqn{x} is positive without checking.
}
