% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms1_precursors.R
\name{calc_aamasses}
\alias{calc_aamasses}
\title{Calculates molecular weight of a polypeptide ([MH]+).}
\usage{
calc_aamasses(
  fixedmods = c("TMT6plex (K)", "Carbamidomethyl (. = C)"),
  varmods = c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
    "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)"),
  rm_dup_term_anywhere = TRUE,
  aa_masses = NULL,
  varlabs = NULL,
  mod_motifs = NULL,
  maxn_vmods_setscombi = 64L,
  out_path = NULL
)
}
\arguments{
\item{fixedmods}{A character vector of fixed modifications. See also
\link{parse_unimod} for grammars.}

\item{varmods}{A character vector of variable modifications.}

\item{rm_dup_term_anywhere}{Logical; if TRUE, removes combinations in
variable modifications with site(s) in positions of both terminal and
anywhere, e.g., "Gln->pyro-Glu (N-term = Q)" and "Deamidated (Q).}

\item{aa_masses}{An amino-acid mass lookup.}

\item{varlabs}{Character string(s) of variable isotopic labels. See examples
of SILAC for details. Can be but not typically used in standard alone
searches of labeled residues.}

\item{mod_motifs}{The motifs to restrict \code{Anywhere} variable
modification. For example, provided the \code{Anywhere} variable
modifications containing \code{c("Oxidation (M)", "Deamidated (N)")} and
\code{mod_motifs = list(`Deamidated (N)` = c("NG", "NM"), `Oxidation (M)` =
c("NM", "MP"))},
variable modifications will only be considered at sites that satisfy the
motifs.}

\item{maxn_vmods_setscombi}{Integer; the maximum number of combinatorial variable
modifications and neutral losses.}

\item{out_path}{An output path.}
}
\description{
Calculates molecular weight of a polypeptide ([MH]+).
}
\examples{
\donttest{
library(mzion)

x <- calc_aamasses()
x_att <- lapply(x, attributes)
names(x_att[[1]])
x_vmods <- lapply(x_att, `[`, c("vmods"))
x_fmods <- lapply(x_att, `[`, c("fmods"))

x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"), 
                   varmods   = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))

stopifnot(length(x) == 6L)

# Fixed N-term mod (no coercion to variable mod)
x <- calc_aamasses(fixedmods = "TMT6plex (N-term)", varmods = NULL)
x[[1]][["N-term"]]

# Fixed N-term mod (coerced to variable mod)
x <- calc_aamasses(fixedmods = "TMT6plex (N-term)", varmods = "Acetyl (Protein N-term)")
lapply(x, `[[`, "N-term")
x[[1]]["N-term"]
x[[2]][["Acetyl (Protein N-term)"]]

# No fixed mod
x <- calc_aamasses(fixedmods = NULL)
stopifnot(length(x) == 16L)

x <- calc_aamasses(fixedmods = NULL, varmods = NULL)
stopifnot(length(x) == 1L)

# Fixed mod, no NL
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
                                 "Carbamidomethyl (. = C)"), varmods = NULL)

stopifnot(length(x) == 1L)

# Fixed mod + NL
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
                                 "Carbamidomethyl (. = M)"), varmods = NULL)

stopifnot(length(x) == 1L)

# Fixed mod, no NL; var mod, no NL
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
                                 "Carbamidomethyl (. = C)"),
                   varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)"))

stopifnot(length(x) == 3L)

# Fixed mod + NL; var mod + NL
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
                     "Carbamidomethyl (. = M)",
                     "Deamidated (. = R)"),
                   varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)",
                     "Hex(5)HexNAc(2) (N)"))

stopifnot(length(x) == 6L)

x <- calc_aamasses(c(fixedmods = "TMT6plex (N-term)", "TMT6plex (K)",
                     "Carbamidomethyl (. = M)", "Deamidated (. = R)"),
                   varmods = c("Acetyl (N-term)", "Carbamyl (. = M)",
                     "Gln->pyro-Glu (N-term = Q)", "Hex(5)HexNAc(2) (N)"))

stopifnot(length(x) == 18L)

## Coercion       
# No fixed terminal or fixed anywhere coercion
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                     "Carbamidomethyl (C)"),
                   varmods = c("Carbamidomethyl (M)"))

stopifnot(length(x) == 2L)

x <- calc_aamasses(fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"), 
                   varmods = c("Acetyl (Protein N-term)", "TMT6plex (N-term)", 
                     "Oxidation (M)", "Carbamidomethyl (M)"))

stopifnot(length(x) == 12L)

# Fixed terminal coercion
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                     "Carbamidomethyl (C)"),
                   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)"))
                   
stopifnot(length(x) == 4L)

# Fixed anywhere coercion
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                     "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
                   varmods = c("Oxidation (M)"))
                   
stopifnot(length(x) == 3L)
                   
# Both fixed terminal and fixed anywhere coercion
x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                     "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
                   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)"))

stopifnot(length(x) == 6L)

}
\dontrun{
# conflicts
x <- calc_aamasses(fixedmods = c("Carbamidomethyl (N-term)", "TMT2plex (N-term)"), 
                   varmods = NULL)

# need separate S and T
x <- calc_aamasses(fixedmods = NULL, varmods = "Phospho (ST)")
}
}
