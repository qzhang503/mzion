% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{prep_pepfdr_td}
\alias{prep_pepfdr_td}
\title{Prepares target-decoy data.}
\usage{
prep_pepfdr_td(
  td = NULL,
  out_path,
  enzyme = "trypsin_p",
  nes_fdr_group = "base",
  fdr_group = "base",
  is_notched = FALSE,
  only_notch_zero = TRUE
)
}
\arguments{
\item{td}{A data frame of targets and decoys (for Percolator).}

\item{out_path}{A file path of outputs.}

\item{enzyme}{A character string; the proteolytic specificity of the assumed
enzyme will be used to generate peptide sequences from protein entries. The
default is \code{Trypsin_P}. See also parameter \code{custom_enzyme}.}

\item{nes_fdr_group}{A character string in one of \code{c("all",
"all_cterm_tryptic", "all_cterm_nontryptic", "base", "base_cterm_tryptic",
"base_cterm_nontryptic")}. All peptides will be used in the classifications
of targets and decoys at \code{"all"}. Peptides with the chemistry of
C-terminal K or R will be used at \code{"all_cterm_tryptic"} (peptides from
protein C-terminals being excluded). Peptides without C-terminal K or R
will be used at \code{"all_cterm_nontryptic"}. The same applied to
\code{"base_cterm_tryptic"} and \code{"base_cterm_nontryptic"} with the
difference of only peptides from the \code{base} group being used. See also
parameter \code{fdr_group}.}

\item{fdr_group}{A character string; the modification group(s) for uses in
peptide FDR controls. The value is in one of \code{c("all", "base")}. The
\code{base} corresponds to the modification group with the largest number
of matches.}

\item{is_notched}{Logical; is a search with MS1 notches or not.}

\item{only_notch_zero}{Logical; if TRUE, use only data at the zero notch.}
}
\description{
Prepares target-decoy data.
}
