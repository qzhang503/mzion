% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pasefreader.R
\name{extract_pasef_ms2}
\alias{extract_pasef_ms2}
\title{Extracts PASEF MS1 data.}
\usage{
extract_pasef_ms2(ms2data, lens2, iso_info, keys, ymin = 10, step = 1.6e-05)
}
\arguments{
\item{ms2data}{A list of \emph{MS2} PASEF frames. Each list entry
corresponding to one MS2 frame.}

\item{lens2}{The number of lines in each MS2 frame.}

\item{iso_info}{The information of MS2 isolation windows etc.}

\item{keys}{The key names of output lists.}

\item{ymin}{The minimum intensity values for considerations.}

\item{step}{A step size for mass binning.}
}
\description{
Extracts PASEF MS1 data.
}
