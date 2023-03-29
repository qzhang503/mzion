#' Helper: switches among ion types for calculating MS2 masses.
#'
#' @param aam A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values (note that argument \code{aas}
#'   corresponds to residues without masses).
#' @param ntmass The mass of a fixed or variable N-term modification.
#' @param ctmass The mass of a fixed or variable C-term modification.
#' @inheritParams matchMS
ms2ions_by_type <- function (aam, ntmass, ctmass, type_ms2ions = "by", 
                             digits = 4L) 
{
  switch(type_ms2ions, 
         by = byions(ntmass, ctmass, aam, digits), 
         cz = czions(ntmod, ctmod, aam, digits), 
         ax = axions(ntmod, ctmod, aam, digits), 
         stop("Unknown type.", call. = FALSE))
}


#' Masses of singly-charged b- and y-ions.
#' 
#' @inheritParams ms2ions_by_type
#' @rdname bions_base
byions <- function (ntmass, ctmass, aam, digits = 4L) 
  c(bions_base(aam, ntmass, digits), yions_base(aam, ctmass, digits))


#' Masses of singly-charged c- and z-ions.
#'
#' @rdname bions_base
czions <- function (ntmass, ctmass, aam, digits = 4L) 
  c(cions_base(aam, ntmass, digits), zions_base(aam, ctmass, digits))

#' Masses of singly-charged a- and x-ions.
#'
#' @rdname bions_base
axions <- function (ntmass, ctmass, aam, digits = 4L) 
  c(aions_base(aam, ntmass, digits), xions_base(aam, ctmass, digits))


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
#' @param digits Integer; the number of decimal places to be used.
#' @param tmass The mass of a fixed or variable N-term or C-term modification.
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
#'
#' }
bions_base <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' Y-ions.
#' 
#' @rdname bions_base
yions_base <- function (aam, tmass, digits = 4L) 
{
  # (1) OH (C-term), + H (neutralizes the N-term on a fragment) + H+ 
  # (2) Other C-term (other than OH) + H + H+: X + 1.007825 + 1.00727647
  ions <- c(tmass, aam[length(aam):1L])
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' B2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
b2ions_base <- function (aam, tmass, digits = 4L, n = 2L) 
  (bions_base(aam, tmass, digits) + 1.00727647)/n


#' B*-ions.
#' 
#' @rdname bions_base
bstarions <- function (aam, tmass, digits = 4L) 
{
  # -NH3:17.026549
  ions <- c(tmass - 17.026549, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' B*2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
bstar2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (bstarions(aam, tmass, digits) + 1.00727647)/n


#' B0-ions.
#' 
#' \code{H2O = 18.010565}.
#' 
#' @rdname bions_base
b0ions <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass - 18.010565, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' B02-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
b02ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (b0ions(aam, tmass, digits) + 1.00727647)/n


#' Y2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
y2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (yions_base(aam, tmass, digits) + 1.00727647)/n


#' Y*-ions.
#' 
#' @rdname bions_base
ystarions <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass - 17.026549, aam[length(aam):1L])
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' Y*2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
ystar2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (ystarions(aam, tmass, digits) + 1.00727647)/n


#' Y0-ions.
#' 
#' @rdname bions_base
y0ions <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass - 18.010565, aam[length(aam):1L])
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' Y02-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
y02ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (y0ions(aam, tmass, digits) + 1.00727647)/n


#' C-ions.
#' 
#' @rdname bions_base
cions_base <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass + 17.026549, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' C2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
c2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (cions_base(aam, tmass, digits) + 1.00727647)/n


#' Z-ions.
#' 
#' @rdname bions_base
zions_base <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass - 17.026549, aam[length(aam):1L])
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' Z2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
z2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (zions_base(aam, tmass, digits) + 1.00727647)/n


#' A-ions.
#' 
#' @rdname bions_base
aions_base <- function (aam, tmass, digits = 4L) 
{
  ions <- c(tmass - 27.9949146, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' A2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
a2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (aions_base(aam, tmass, digits) + 1.00727647)/n


#' A*-ions.
#' 
#' @rdname bions_base
astarions <- function (aam, tmass, digits = 4L) 
{
  # -CO -NH3 = -(27.9949146 + 17.026549)
  ions <- c(tmass - 45.0214636, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' A*2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
astar2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (astarions(aam, tmass, digits) + 1.00727647)/n


#' A0-ions.
#' 
#' @rdname bions_base
a0ions <- function (aam, tmass, digits = 4L) 
{
  # -CO -H2O = -(27.9949146 + 18.010565)
  ions <- c(tmass - 46.0054796, aam)
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' A02-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
a02ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (a0ions(aam, tmass, digits) + 1.00727647)/n


#' X-ions.
#' 
#' @rdname bions_base
xions_base <- function (aam, tmass, digits = 4L) 
{
  # +CO -H2 = 27.9949146 - 2*1.007825
  ions <- c(tmass + 25.9792646, aam[length(aam):1L])
  ions <- cumsum(ions)
  ions <- ions[-1]
  
  round(ions, digits = digits)
}


#' X2-ions.
#' 
#' @param n The charge state.
#' @rdname bions_base
x2ions <- function (aam, tmass, digits = 4L, n = 2L) 
  (xions_base(aam, tmass, digits) + 1.00727647)/n

