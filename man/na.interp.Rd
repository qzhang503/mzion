% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{na.interp}
\alias{na.interp}
\title{Interpolates missing values in a time series.}
\usage{
na.interp(
  x,
  lambda = NULL,
  linear = (frequency(x) <= 1 | sum(!is.na(x)) <= 2 * frequency(x))
)
}
\arguments{
\item{x}{Time series}

\item{lambda}{Box-Cox transformation parameter. If lambda="auto", then a
transformation is automatically selected using BoxCox.lambda. The
transformation is ignored if NULL. Otherwise, data transformed before model
is estimated.}

\item{linear}{Linear.}
}
\description{
From \code{forecast}.
}
