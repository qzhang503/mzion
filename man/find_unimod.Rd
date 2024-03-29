% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unimods.R
\name{find_unimod}
\alias{find_unimod}
\title{Finds a Unimod.}
\usage{
find_unimod(
  unimod = "Carbamidomethyl (C)",
  xml_files = c("master.xml", "custom.xml")
)
}
\arguments{
\item{unimod}{The name of a \href{https://www.unimod.org/}{Unimod} modification.}

\item{xml_files}{Name(s) of Unimod ".xml" files. The file path is a system
setting of \code{system.file("extdata", xml_file, package = "mzion")}.}
}
\description{
Finds the mono-isotopic mass, position, site and neutral losses of a
modification.
}
\details{
In the field of \code{position_site}, \code{position} is the name and
\code{site} is the value.
}
\examples{
\donttest{
library(mzion)

x1 <- find_unimod("Carbamidomethyl (C)")
x2 <- find_unimod("Carbamidomethyl (M)")
x3 <- find_unimod("Acetyl (Protein N-term)")
x4 <- find_unimod("Gln->pyro-Glu (N-term = Q)")
x5 <- find_unimod("Hex(5)HexNAc(2) (N)")
}

\dontrun{
# Prefer an "=" sign between 'N-term' and 'Q'
x <- find_unimod("Gln->pyro-Glu (N-term Q)")
}
}
\seealso{
\link{table_unimods}, \link{parse_unimod}.
}
