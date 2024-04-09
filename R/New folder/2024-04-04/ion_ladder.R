#' Helper: switches among ion types for calculating MS2 masses.
#'
#' @param aam A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values (note that argument \code{aas}
#'   corresponds to residues without masses).
#' @param ntmass The mass of a fixed or variable N-term modification.
#' @param ctmass The mass of a fixed or variable C-term modification.
#' @inheritParams matchMS
ms2ions_by_type <- function (aam, ntmass, ctmass, type_ms2ions = "by") 
{
  switch(type_ms2ions, 
         by = byions(ntmass = ntmass, ctmass = ctmass, aam = aam), 
         cz = czions(ntmass = ntmass, ctmass = ctmass, aam = aam), 
         ax = axions(ntmass = ntmass, ctmass = ctmass, aam = aam), 
         stop("Unknown type.", call. = FALSE))
}


#' Masses of singly-charged b- and y-ions.
#' 
#' b-ions first, then y-ions
#' 
#' @rdname bions_base
byions <- function (ntmass, ctmass, aam) 
  c(cumsum(c(ntmass, aam))[-1], cumsum(c(ctmass, aam[length(aam):1L]))[-1])


#' Masses of singly-charged c- and z-ions.
#'
#' @rdname bions_base
czions <- function (ntmass, ctmass, aam)
  c(cumsum(c(ntmass + 17.026549, aam))[-1], 
    cumsum(c(ctmass - 17.026549, aam[length(aam):1L]))[-1])


#' Masses of singly-charged a- and x-ions.
#'
#' @rdname bions_base
axions <- function (ntmass, ctmass, aam) 
  c(cumsum(c(ntmass - 27.9949146, aam))[-1], 
    cumsum(c(ctmass + 25.9792646, aam[length(aam):1L]))[-1])


###
# No direct uses of the followings.
###

#' B-ions.
#'
#' For (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-".
#'
#' @param aam A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values.
#'
#'   The masses reflects fixed/variable modifications, and/or fixed/variable
#'   neutral losses.
#'   
#' @param ntmass The mass of a fixed or variable N-term modification.
#'
#' @importFrom stringr str_split
#' @examples
#' \donttest{
#' library(mzion)
#' library(stringr)
#' 
#' ## (1) "amods- tmod- vnl- fnl-"
#' # (Fixed N-term mods; also for no N-term mod)
#'
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Oxidation (M)", "Deamidated (N)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses <- aa_masses_all[[1]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (!length(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (!length(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MAKEMASSPECFUN"
#' aas <- stringr::str_split(aa_seq, "", simplify = TRUE)
#' aam <- aa_masses[aas]
#'
#' b <- mzion:::bions_base(aam, ntmass)
#' y <- mzion:::yions_base(aam, ctmass)
#'
#'
#' ## (2) "amods- tmod+ vnl- fnl-"
#' # (2a, N-term)
#' fixedmods <- c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
#'              "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses <- aa_masses_all[[3]]
#'
#' # (Fixed or variable C-term mods +/- makes no difference on b-ions;
#' # and vice versa for y-ions)
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (!length(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (!length(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MAKEMASSPECFUN"
#' aas <- stringr::str_split(aa_seq, "", simplify = TRUE)
#' aam <- aa_masses[aas]
#'
#' b <- mzion:::bions_base(aam, ntmass)
#' y <- mzion:::yions_base(aam, ctmass)
#'
#'
#' # (2b, C-term)
#' fixedmods <- c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("TMT6plex (N-term)", "Amidated (Protein C-term)", "Oxidation (M)",
#'              "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' # `TMT6plex (N-term)`; `Amidated (Protein C-term)`
#' aa_masses <- aa_masses_all[[7]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (!length(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (!length(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' b <- mzion:::bions_base(aam, ntmass)
#' y <- mzion:::yions_base(aam, ctmass)
#' }
bions_base <- function (aam, ntmass) cumsum(c(ntmass, aam))[-1]


#' Y-ions.
#' 
#' # (1) OH (C-term), + H (neutralizes the N-term on a fragment) + H+ 
#' # (2) Other C-term (other than OH) + H + H+: X + 1.007825 + 1.00727647
#' 
#' @param ctmass The mass of a fixed or variable C-term modification.
#' @rdname bions_base
yions_base <- function (aam, ctmass) cumsum(c(ctmass, aam[length(aam):1L]))[-1]


#' C-ions.
#' 
#' \code{NH3 = 17.026549}
#' 
#' @rdname bions_base
cions_base <- function (aam, ntmass) cumsum(c(ntmass + 17.026549, aam))[-1]


#' Z-ions.
#' 
#' @rdname bions_base
zions_base <- function (aam, ctmass) 
  cumsum(c(ctmass - 17.026549, aam[length(aam):1L]))[-1]


#' C2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
c2ions <- function (aam, ntmass, n = 2L) (cions_base(aam, ntmass) + 1.00727647)/n


#' Z2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
z2ions <- function (aam, ctmass, n = 2L) (zions_base(aam, ctmass) + 1.00727647)/n


#' A-ions.
#' 
#' \code{CO = 27.9949146}

#' @rdname bions_base
aions_base <- function (aam, ntmass) cumsum(c(ntmass - 27.9949146, aam))[-1]


#' X-ions.
#' 
#' \code{+CO -H2 = 27.9949146 - 2 * 1.007825}
#' 
#' @rdname bions_base
xions_base <- function (aam, ctmass) 
  cumsum(c(ctmass + 25.9792646, aam[length(aam):1L]))[-1]


#' A2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
a2ions <- function (aam, ntmass, n = 2L) (aions_base(aam, ntmass) + 1.00727647)/n


#' A*-ions.
#' 
#' \code{-CO -NH3 = -(27.9949146 + 17.026549)}
#' 
#' @rdname bions_base
astarions <- function (aam, ntmass) cumsum(c(ntmass - 45.0214636, aam))[-1]


#' A*2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
astar2ions <- function (aam, ntmass, n = 2L) (astarions(aam, ntmass) + 1.00727647)/n


#' A0-ions.
#' 
#' \code{-CO -H2O = -(27.9949146 + 18.010565)}
#' 
#' @rdname bions_base
a0ions <- function (aam, ntmass) cumsum(c(ntmass - 46.0054796, aam))[-1]


#' A02-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
a02ions <- function (aam, ntmass, n = 2L) (a0ions(aam, ntmass) + 1.00727647)/n


#' X2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
x2ions <- function (aam, ctmass, n = 2L) (xions_base(aam, ctmass) + 1.00727647)/n


