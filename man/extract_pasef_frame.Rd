% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pasefreader.R
\name{extract_pasef_frame}
\alias{extract_pasef_frame}
\title{Extracts MS1 or MS2 data from a PASEF frame.}
\usage{
extract_pasef_frame(data, ms_lev = 1L, ymin = 50, ymax = 1e+07)
}
\arguments{
\item{data}{A PASEF frame (with multiple slices).}

\item{ms_lev}{The level of MS data.}

\item{ymin}{The cut-off of intensity.}

\item{ymax}{The maximum intensity.}
}
\description{
Extracts MS1 or MS2 data from a PASEF frame.
}
