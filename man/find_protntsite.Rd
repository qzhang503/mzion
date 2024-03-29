% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispatch.R
\name{find_protntsite}
\alias{find_protntsite}
\alias{find_anyntsite}
\alias{find_anysite}
\alias{find_protctsite}
\alias{find_anyctsite}
\alias{contain_protntsite}
\alias{contain_anyntsite}
\alias{contain_anysite}
\alias{contain_protctsite}
\alias{contain_anyctsite}
\title{Finds the sites of amino-acid residues at the position of \code{Protein
N-term}.}
\usage{
find_protntsite(vmods, posns)

find_anyntsite(vmods, posns)

find_anysite(vmods, posns)

find_protctsite(vmods, posns)

find_anyctsite(vmods, posns)

contain_protntsite(vmods, posns, len)

contain_anyntsite(vmods, posns, len)

contain_anysite(vmods, posns, len)

contain_protctsite(vmods, posns, len)

contain_anyctsite(vmods, posns, len)
}
\arguments{
\item{vmods}{A named list of variable modifications.}

\item{posns}{The position (e.g., \code{Protein N-term}, \code{Anywhere},
etc.) of \code{vmods}. The argument can be obtained from \code{vmods} but
passed as a parameter for vectorization.}

\item{len}{The count of \code{vmods} (passed as a parameter for
vectorization).}
}
\description{
Flowchart (1-nt): Dimethyl (Protein N-term = P)

Flowchart (3-nt): Gln->pyro Glu (N-term = Q)

Flowchart (5): Oxidation (M)

Flowchart f(1-ct): Dehydrated (Protein C-term = N)

Flowchart f(3-ct): Oxidation (C-term = G)
}
\details{
In the less common cases of multiple modifications to the same
\code{Anywhere} site, e.g., "Oxidation (M)" and "Carbamyl (M)", an amino-acid
sequence should contain at least two "M"s. This is due to the design that
\code{vmods} in an \code{aa_masses} are realized. The additional subsetting
by the counts of residues is applied at this stage of the peptide
distributions. A downstream step in the calculations of MS2 ions will also
guard against this.

In the case of, e.g., "Oxidation (M)" and "Protein (N-term)", and the peptide
sequence has only one "M" on the N-term, the entry is still kept with
the "additive" effect of modifications.
}
\examples{
\donttest{
library(mzion)

## `Protein N-term = P`
sites <- list(`Dimethyl (Protein N-term = P)` = "P",
              `Acetyl (Protein N-term)` = "N-term", 
              `Oxidation (M)` = "M",
              `Deamidated (N)` = "N")
positions <- c("Protein N-term", "Protein N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

# (pretend `vmods` in a real workflow)
vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
posns <- names(vmods)
# stopifnot(identical(positions, posns))

stopifnot(mzion:::contain_protntsite(vmods, names(vmods), length(vmods)))

ans <- mzion:::find_protntsite(vmods, posns)
stopifnot(identical(ans, vmods[1]))
}
\donttest{
library(mzion)

## Gln->pyro Glu (N-term = Q)
sites <- list(`Gln->pyro Glu (N-term = Q)` = "Q",
              `Acetyl (N-term)` = "N-term",
              `Oxidation (M)` = "M",
              `Deamidated (N)` = "N")
positions <- c("Any N-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
posns <- names(vmods)

stopifnot(mzion:::contain_anyntsite(vmods, posns, length(vmods)))

ans <- mzion:::find_anyntsite(vmods, posns)
stopifnot(identical(ans, vmods[1])) # M and N
}
\donttest{
library(mzion)

## `Oxidation (M)` and `Deamidated (N)`
sites <- list(`Acetyl (N-term)` = "N-term",
              `Oxidation (M)` = "M",
              `Deamidated (N)` = "N")
positions <- c("Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
posns <- names(vmods)

stopifnot(mzion:::contain_anysite(vmods, posns, length(vmods)))

ans <- mzion:::find_anysite(vmods, posns)
stopifnot(length(ans) == 2L)


## `Oxidation (M)` and `Carbamyl (M)`
sites <- list(`Dimethyl (Protein N-term = P)` = "P",
              `Oxidation (M)` = "M", `Carbamyl (M)` = "M",
              `Deamidated (N)` = "N")
positions <- c("Any N-term", "Anywhere", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
vmods <- vmods[!duplicated.default(vmods)]

posns <- names(vmods)
stopifnot(length(posns) < length(positions))

ans <- mzion:::find_anysite(vmods, posns)
stopifnot(length(ans) == 2L)
}
\donttest{
library(mzion)

## `Dehydrated (Protein C-term = N)`
sites <- list(`Dehydrated (Protein C-term = N)` = "N", 
              `Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Protein C-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
posns <- names(vmods)

stopifnot(mzion:::contain_protctsite(vmods, posns, length(vmods)))

ans <- mzion:::find_protctsite(vmods, posns)
stopifnot(identical(ans, vmods[1]))
}
\donttest{
library(mzion)

## `Oxidation (C-term = G)`
sites <- list(`Oxidation (C-term = G)` = "G", 
              `Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Any C-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))

vmods <- unname(vmods)
vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
posns <- names(vmods)

stopifnot(mzion:::contain_anyctsite(vmods, posns, length(vmods)))

ans <- mzion:::find_anyctsite(vmods, posns)
stopifnot(identical(ans, vmods[1]))
}
}
\seealso{
contain_protntany
}
