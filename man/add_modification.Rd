% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unimods.R
\name{add_modification}
\alias{add_modification}
\title{Adds a new node of modification.}
\usage{
add_modification(
  node = NULL,
  title = "Acetyl",
  full_name = "",
  site = "N-term",
  position = "Protein N-term",
  mod_mono_mass = "42.010565",
  mod_avge_mass = "42.0367",
  mod_composition = "H(2) C(2) O",
  neuloss_mono_mass = "0",
  neuloss_avge_mass = "0",
  neuloss_composition = "0",
  record_id = "0"
)
}
\arguments{
\item{node}{A node of top-level \code{modification} (one of the top-4).}

\item{title}{The title of a modification.}

\item{full_name}{The full-name description of a modification.}

\item{site}{The site of a modification.}

\item{position}{The position of a modification.}

\item{mod_mono_mass}{The mono-isotopic mass delta of a modification.}

\item{mod_avge_mass}{The average mass delta of a modification.}

\item{mod_composition}{The chemical composition of a modification.}

\item{neuloss_mono_mass}{The mono-isotopic mass delta of neutral loss in a
positive value.}

\item{neuloss_avge_mass}{The average mass delta of neutral loss in a positive
value.}

\item{neuloss_composition}{The chemical composition of a neutral loss.}

\item{record_id}{The index of a modification.}
}
\description{
Adds a new node of modification.
}
