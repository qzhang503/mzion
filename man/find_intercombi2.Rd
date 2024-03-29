% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vmod_ms1_labels.R
\name{find_intercombi2}
\alias{find_intercombi2}
\title{Finds the combinations across residues (with attributes).}
\usage{
find_intercombi2(vmodsets = NULL, maxn_vmods_per_pep = 5L)
}
\arguments{
\item{vmodsets}{The results from \link{make_ms1_vmodsets}.}

\item{maxn_vmods_per_pep}{A non-negative integer; the maximum number of
\code{Anywhere} (non-terminal) variable modifications per peptide. The
default is 5.}
}
\description{
Used after \link{make_ms1_vmodsets} across residues.
}
\examples{
\donttest{
library(mzion)

## Simple
fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")

varmods <- c("Acetyl (Protein N-term)",
             "Oxidation (M)", "Deamidated (N)", 
             "Gln->pyro-Glu (N-term = Q)")

aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                               varmods = varmods,
                               maxn_vmods_setscombi = 64,
                               out_path = NULL)

maxn_vmods_per_pep <- 5L
maxn_sites_per_vmod <- 3L

# M and N
i <- 10L
aa_masses <- aa_masses_all[[i]]

vmodsets <- mzion:::make_ms1_vmodsets(aa_masses_all = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod)

ms1vmods <- mzion:::find_intercombi2(vmodsets)

## More complex
fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")

varmods <- c("Acetyl (Protein N-term)",
             "Oxidation (M)", "Carbamidomethyl (M)",
             "Deamidated (N)", "Carbamyl (M)",
             "Gln->pyro-Glu (N-term = Q)")

aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                               varmods = varmods,
                               maxn_vmods_setscombi = 64,
                               out_path = NULL)

maxn_vmods_per_pep <- 5L
maxn_sites_per_vmod <- 3L

# Only `Oxidation (M)` 
i <- 8L
aa_masses <- aa_masses_all[[i]]

vmodsets <- mzion:::make_ms1_vmodsets(aa_masses_all = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod)

ms1vmods <- mzion:::find_intercombi2(vmodsets)

## `Oxidation (M)`, "Carbamidomethyl (M)"
fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")

varmods <- c("Acetyl (Protein N-term)",
             "Oxidation (M)", "Carbamidomethyl (M)",
             "Deamidated (N)")

aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                               varmods = varmods,
                               maxn_vmods_setscombi = 64,
                               out_path = NULL)

maxn_vmods_per_pep <- 5L
maxn_sites_per_vmod <- 3L

i <- 16L
aa_masses <- aa_masses_all[[i]]

vmodsets <- mzion:::make_ms1_vmodsets(aa_masses_all = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod)

ms1vmods <- mzion:::find_intercombi2(vmodsets)
}
}
