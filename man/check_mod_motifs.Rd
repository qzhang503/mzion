% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{check_mod_motifs}
\alias{check_mod_motifs}
\title{Checks mod_motifs}
\usage{
check_mod_motifs(mod_motifs, mods)
}
\arguments{
\item{mod_motifs}{The motifs to restrict \code{Anywhere} variable
modification. For example, provided the \code{Anywhere} variable
modifications containing \code{c("Oxidation (M)", "Deamidated (N)")} and
\code{mod_motifs = list(`Deamidated (N)` = c("NG", "NM"), `Oxidation (M)` =
c("NM", "MP"))},
variable modifications will only be considered at sites that satisfy the
motifs.}

\item{mods}{A concatenated list of fixed and variable modifications.}
}
\description{
Checks mod_motifs
}
