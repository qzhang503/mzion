% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_masses.R
\name{binTheoSeqs}
\alias{binTheoSeqs}
\title{Helper of \link{bin_theoseqs}.}
\usage{
binTheoSeqs(
  idxes = NULL,
  res = NULL,
  min_mass = 200L,
  max_mass = 4500L,
  ppm_ms1 = 10L,
  enzyme = "trypsin_p",
  out_path = NULL
)
}
\arguments{
\item{idxes}{A set of indexes, e.g. "1" for \code{pepmasses_1.rds} and
"rev_1" for \code{pepmasses_rev_1.rds}.}

\item{res}{Lists of data containing theoretical peptides and masses from
\link{readRDS}.}

\item{min_mass}{Numeric; the minimum MS1 mass.}

\item{max_mass}{Numeric; the maximum MS1 mass.}

\item{ppm_ms1}{Numeric; (half of) the error tolerance of MS1 mass in ppm.}

\item{enzyme}{The assume enzyme activity.}

\item{out_path}{The output path.}
}
\value{
Lists of theoretical peptides binned by MS1 masses. The lists
  correspond to the lists of \code{res}.
}
\description{
Helper of \link{bin_theoseqs}.
}
\examples{
\donttest{
library(mzion)

# res <- readRDS("~/mzion/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
# theopeps <- mzion:::binTheoSeqs(res)
}
}
