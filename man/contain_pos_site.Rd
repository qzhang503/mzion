% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispatch.R
\name{contain_pos_site}
\alias{contain_pos_site}
\title{Function factories for checking the existence of an amino-acid residue at a
position.}
\usage{
contain_pos_site(pos)
}
\arguments{
\item{pos}{The Unimod position of a variable modification, which must be one
of \code{c("Protein N-term", "Protein C-term", "Any N-term", "Any C-term",
"Anywhere")}. No need to check its value for the way the function is called
(by the developer).}
}
\value{
A function for checking the existence of a residue at the position
  specified by the argument \code{pos}. For each function, it takes a list of
  variable modifications specified by argument \code{vmods}, and the
  corresponding postitions and counts as inputs.
}
\description{
No need to check \code{is.null(p)} and \code{is.null(posns)}.
}
