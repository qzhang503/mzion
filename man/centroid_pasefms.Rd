% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pasefreader.R
\name{centroid_pasefms}
\alias{centroid_pasefms}
\title{Sum PASEF MS1 peak area}
\usage{
centroid_pasefms(xs, ys, reso = 60000, maxn = 2000L, ymin = 100L, tol = 0.1)
}
\arguments{
\item{xs}{A vector of ascending m-over-z values.}

\item{ys}{A vector of intensity values.}

\item{reso}{The resolution of a peak.}

\item{maxn}{The maximum number of peaks.}

\item{ymin}{The minimum Y values for considering in peak centroiding.}

\item{tol}{The tolerance of Y for defining a peak profile.}
}
\description{
Sum PASEF MS1 peak area
}
\examples{
# example code
mzion:::centroid_pasefms(c(500), c(10))
mzion:::centroid_pasefms(c(500, 500.01), c(1, 10))
mzion:::centroid_pasefms(c(500, 500.01, 500.2), c(1, 10, 20))
mzion:::centroid_pasefms(c(500, 500.01, 500.2), c(10, 5, 2))

# ignore the trailing half peak
mzion:::centroid_pasefms(500 + .01 * 0:3, c(1, 10, 2, 5))

mzion:::centroid_pasefms(500 + .01 * 0:4, c(1, 10, 2, 5, 3))

# trailing max at the last position
mzion:::centroid_pasefms(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,33,23,22,41))

mzion:::centroid_pasefms(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,15,23,30,41))
}
