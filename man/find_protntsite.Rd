% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ions.R
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
\title{Finds the sites of amino-acid residues at corresponding variable
modifications.}
\usage{
find_protntsite(vmods)

find_anyntsite(vmods)

find_anysite(vmods)

find_protctsite(vmods)

find_anyctsite(vmods)

contain_protntsite(vmods)

contain_anyntsite(vmods)

contain_anysite(vmods)

contain_protctsite(vmods)

contain_anyctsite(vmods)
}
\arguments{
\item{vmods}{A named list of variable modifications.}
}
\description{
Flowchart (1-nt): `Dimethyl (Protein N-term = P)`

Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`

Flowchart (5): `Oxidation (M)`

Flowchart f(1-ct): `Dehydrated (Protein C-term = N)`

Flowchart f(3-ct): `Oxidation (C-term = G)`
}
\examples{
# `Protein N-term = P`
sites <- list(`Dimethyl (Protein N-term = P)` = "P", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Protein N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
contain_protntsite(vmods)
find_protntsite(vmods)
# Gln->pyro Glu (N-term = Q)
sites <- list(`Gln->pyro Glu (N-term = Q)` = "Q", 
              `Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Any N-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
contain_anyntsite(vmods)
find_anyntsite(vmods)
# `Oxidation (M)` and `Deamidated (N)`
sites <- list(`Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
contain_anysite(vmods)
find_anysite(vmods)
# `Dehydrated (Protein C-term = N)`
sites <- list(`Dehydrated (Protein C-term = N)` = "N", 
              `Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Protein C-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
contain_protctsite(vmods)
find_protctsite(vmods)
# `Oxidation (C-term = G)`
sites <- list(`Oxidation (C-term = G)` = "G", 
              `Acetyl (N-term)` = "N-term", 
              `Oxidation (M)` = "M", 
              `Deamidated (N)` = "N")
positions <- c("Any C-term", "Any N-term", "Anywhere", "Anywhere")
vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
contain_anyctsite(vmods)
find_anyctsite(vmods)
}
\seealso{
contain_protntany
}