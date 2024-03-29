% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_os.R
\name{reloc_col_after}
\alias{reloc_col_after}
\title{Relocate column "to_move" immediately after column "col_before".}
\usage{
reloc_col_after(df, to_move = "after_anchor", col_before = "anchor")
}
\arguments{
\item{df}{The original data frame.}

\item{to_move}{The column to be moved.}

\item{col_before}{The anchor column to which the \code{to_move} will be moved
after.}
}
\description{
Relocate column "to_move" immediately after column "col_before".
}
