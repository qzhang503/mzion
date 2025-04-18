% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{split_fastaseqs}
\alias{split_fastaseqs}
\title{Helper of \link{calc_pepmasses2}.}
\usage{
split_fastaseqs(
  fasta = NULL,
  enzyme = "trypsin_p",
  custom_enzyme = c(Cterm = NULL, Nterm = NULL),
  acc_type = "uniprot_acc",
  acc_pattern = NULL,
  maxn_fasta_seqs = 200000L,
  max_miss = 2L
)
}
\arguments{
\item{fasta}{Character string(s) to the name(s) of fasta file(s) with
prepended directory path. The experimenter needs to supply the files.}

\item{enzyme}{A character string; the proteolytic specificity of the assumed
enzyme will be used to generate peptide sequences from proteins. The enzyme
is currently \code{trypsin}.}

\item{custom_enzyme}{Regular expression(s) for custom enzyme specificity. The
  default is NULL. Uses of custom enzyme specificity is probably rather
  infrequent. Should there be such need, the argument \code{enzyme} will be
  ignored and the following may be applied:

  \cr ## Examples \cr \cr # Equivalent to Trypsin \cr # at the Cterm of K or
  R but not followed by P \cr # (the quantifiers "\{1\}" can be skipped at a
  small cost of speed) \cr custom_enzyme = c(Cterm =
  "([KR]\{1\})([^P]\{1\})")

  \cr # GluN again \cr custom_enzyme = c(Nterm = "([A-Z]{1})([E]\{1\})")

  \cr # Trypsin_P + GluN \cr custom_enzyme = c(Cterm = "([KR]\{1\})", Nterm =
  "([E]\{1\})")

  \cr # Faked: Trypsin, proline not allowed on neither Nterm or Cterm \cr
  custom_enzyme = c(Cterm = "([KR]\{1\})([^P]\{1\})", Nterm =
  "([^P]\{1\})([KR]\{1\})")
  
  \cr # Faked: N-term trypsin, serine not allowed on Nterm \cr
  custom_enzyme = c(Nterm = "([^S]\{1\})([KR]\{1\})")}

\item{acc_type}{Character string(s); the types of protein accessions in one
of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
the corresponding regular expressions need to be supplied via argument
\code{acc_pattern}.}

\item{acc_pattern}{Regular expression(s) describing the patterns in
separating the header lines of fasta entries. At the \code{NULL} default,
the pattern will be automated when \code{acc_type} are among
c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). See also
\link{load_fasta2} for custom examples.}

\item{maxn_fasta_seqs}{Integer; the maximum number of protein sequences in
fasta files.}

\item{max_miss}{The maximum number of mis-cleavages per peptide sequence.}
}
\value{
Two named list of "fwds" and "revs". List "fwds" contains peptide
  sequences split from forward fasta and "revs" from reversed fasta.
}
\description{
Prior to the calculation of peptide masses; for the base with fixed
modification only.
}
\details{
Each entry in the outputs contains two lists: 1. terminal sequences without
N-term methionine; 2. all sequences without methionine clipping.
}
\examples{
if (FALSE) {
  x <- split_fastaseqs("~/mzion/dbs/fasta/uniprot/uniprot_mm_2020_11.fasta")
  x <- split_fastaseqs("~/mzion/dbs/fasta/crap/crap.fasta")
}
}
