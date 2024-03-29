% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastas.R
\name{find_acc_pattern}
\alias{find_acc_pattern}
\title{Helper for \link{load_fasta2}.}
\usage{
find_acc_pattern(acc_type)
}
\arguments{
\item{acc_type}{Character string(s); the types of protein accessions in one
of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
the corresponding regular expressions need to be supplied via argument
\code{acc_pattern}.}
}
\description{
Not used for custom acc_type, i.e. acc_type = "my_acctype".
}
