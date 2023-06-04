#' Calculates the masses of MS2 ion series.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#' 
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param ms1vmods The set of all possible MS1 vmod labels at a given aa_masses.
#' @param ms2vmods Matrices of labels of variable modifications. Each
#'   permutation in a row for each matrix.
#' @param ntmass The mass of N-terminal.
#' @param ctmass The mass of C-terminal.
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @param vmods_nl The attribute of \code{vmods_nl} from an \code{aa_masses}.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' 
#' @seealso \link{bions_base}, \link{yions_base}.
#' 
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' # (2) "amods- tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
#'              "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses = aa_masses_all[[2]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (!length(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647 # + proton
#' }
#'
#' if (!length(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- mzion:::gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                         aa_masses = aa_masses, ntmod = NULL, ctmod = NULL, 
#'                         ntmass = ntmass, ctmass = ctmass, 
#'                         amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
#'                         mod_indexes = mod_indexes)
#'                         
#' # (1) "amods- tmod- vnl- fnl-"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Oxidation (M)", "Deamidated (N)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses <- aa_masses_all[[1]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (!length(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647
#' }
#'
#' if (!length(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- mzion:::gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                         aa_masses = aa_masses, ntmod = NULL, ctmod = NULL, 
#'                         ntmass = ntmass, ctmass = ctmass, 
#'                         amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
#'                         mod_indexes = mod_indexes)
#' }
gen_ms2ions_base <- function (aa_seq = NULL, ms1_mass = NULL, 
                              aa_masses = NULL, ms1vmods = NULL, ms2vmods = NULL, 
                              ntmod = NULL, ctmod = NULL, 
                              ntmass = NULL, ctmass = NULL, 
                              amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                              mod_indexes = NULL, 
                              type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                              maxn_sites_per_vmod = 3L, 
                              
                              # dummy
                              maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                              maxn_vmods_sitescombi_per_pep = 64L, 
                              digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  naa <- length(aas)
  
  nm <- .Internal(paste0(list(rep_len("0", naa)), collapse = "", recycle0 = FALSE))
  af <- ms2ions_by_type(aam, ntmass, ctmass, type_ms2ions, digits)
  
  av <- list(calc_rev_ms2(af, aas))
  names(av) <- NA_character_
  af <- list(af)
  names(af) <- nm
  c(af, av)
}


#' Calculates the masses of MS2 ion series.
#'
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @rdname gen_ms2ions_base
#' 
#' @examples 
#' \donttest{
#' library(mzion)
#' 
#' # (5) "amods- tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (N-term)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Acetyl (Protein N-term)")
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'   
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' aa_masses <- aa_masses_all[[2]]
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
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#' 
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' 
#' # variable `TMT6plex (N-term)` + `fixed Oxidation (M)`
#' # (additive varmod on top of fixedmod allowed)
#' 
#' out <- mzion:::gen_ms2ions_a0_vnl0_fnl1(
#'    aa_seq = aa_seq, ms1_mass = NULL, 
#'    aa_masses = aa_masses, ntmod = NULL, ctmod = NULL, 
#'    ntmass = ntmass, ctmass = ctmass, 
#'    amods = NULL, vmods_nl = NULL, fmods_nl = fmods_nl, 
#'    mod_indexes = mod_indexes)
#' }
gen_ms2ions_a0_vnl0_fnl1 <- function (aa_seq, ms1_mass = NULL, 
                                      aa_masses = NULL, ms1vmods = NULL, ms2vmods = NULL, 
                                      ntmod = NULL, ctmod = NULL, # not used
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, vmods_nl = NULL, # not used
                                      fmods_nl = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_fnl_per_seq = 3L, 
                                      
                                      # dummy
                                      maxn_vnl_per_seq = 3L, 
                                      
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      digits = 4L) 
{
  if (maxn_fnl_per_seq < 2L)
    return(
      gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                       aa_masses = aa_masses, ms1vmods = NULL, ms2vmods = NULL, 
                       ntmod = NULL, ctmod = NULL, 
                       ntmass = ntmass, ctmass = ctmass, 
                       amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                       mod_indexes = mod_indexes, 
                       type_ms2ions = type_ms2ions, 
                       maxn_vmods_per_pep = maxn_vmods_per_pep, 
                       maxn_sites_per_vmod = maxn_sites_per_vmod, 
                       maxn_vmods_sitescombi_per_pep = 
                         maxn_vmods_sitescombi_per_pep, 
                       digits = digits))

  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  # (no pep_seq dispatching by Anywhere fmod residues -> possible no matched sites)
  
  sites <- names(fmods_nl)
  
  pattern <- if (length(sites) > 1L)
    .Internal(paste0(list(sites), collapse = "|", recycle0 = FALSE))
  else
    sites

  if (!grepl(pattern, aa_seq)) 
    return(
      gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                       aa_masses = aa_masses, ms1vmods = NULL, ms2vmods = NULL, 
                       ntmod = NULL, ctmod = NULL, 
                       ntmass = ntmass, ctmass = ctmass, 
                       amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                       mod_indexes = mod_indexes, 
                       type_ms2ions = type_ms2ions, 
                       maxn_vmods_per_pep = maxn_vmods_per_pep, 
                       maxn_sites_per_vmod = maxn_sites_per_vmod, 
                       maxn_vmods_sitescombi_per_pep = 
                         maxn_vmods_sitescombi_per_pep, 
                       digits = digits))
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  naa <- length(aas)
  
  # At varmods "Oxidation (M)", pep_seq(s) must contain "M" 
  #   (with an additional entry of "Oxidation (M)" in aa_masses)
  # 
  # At fixedmods "Oxidation (M)", pep_seq(s) may not contain "M"; 
  #   (as `distri_peps` does not filter pep_seq by fixedmods)
  
  idxes <- .Internal(which(aas %fin% names(fmods_nl)))
  
  if (length(idxes) > maxn_vmods_per_pep)
    idxes <- idxes[1:maxn_vmods_per_pep]
  
  # ---
  fmods_combi <- aas[idxes]

  if (length(fmods_combi) == 1L) {
    fnls <- fmods_nl[[fmods_combi]]
    len <- length(fnls)
    ans <- vector("list", len)
    
    for (i in 1:len) {
      ans[[i]] <- fnls[[i]]
      names(ans[[i]]) <- fmods_combi
    }
  }
  else {
    ans <- expand_grid_rows(fmods_nl[fmods_combi], nmax = maxn_fnl_per_seq, 
                            use.names = FALSE)
    len <- length(ans)
  }

  av <- af <- vector("list", len)
  aam <- aa_masses[aas]
  af[[1]] <- af1 <- ms2ions_by_type(aam, ntmass, ctmass, type_ms2ions, digits)
  av[[1]] <- av1 <- calc_rev_ms2(af1, aas)
  
  if (len > 1L) {
    aami <- aam
    aamii <- aami[idxes]
    
    for (i in 2:len) {
      fnl_combi_i <- ans[[i]]
      aami[idxes] <- aamii - fnl_combi_i
      af[[i]] <- afi <- ms2ions_by_type(aami, ntmass, ctmass, type_ms2ions, digits)
      av[[i]] <- calc_rev_ms2(afi, aas)
    }
  }
  
  nm <- .Internal(paste0(list(rep_len("0", naa)), collapse = "", recycle0 = FALSE))
  nm <- .Internal(paste0(list(nm, " [", as.character(seq_len(len)), "]"), 
                         collapse = NULL, recycle0 = FALSE))
  names(af) <- nm
  names(av) <- NA_character_
  c(af, av)
}


#' Calculates the masses of MS2 ion series.
#'
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#' 
#' @rdname gen_ms2ions_base
#' 
#' @param amods The attribute \code{amods} from a \code{aa_masses}.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#'  @rdname gen_ms2ions_base
#'
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' # (8a) "amods+ tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("Deamidated (N)", "Carbamidomethyl (S)",
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'                        
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 8L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide("MHQGVMNVGMGQKMNS",
#'                                fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[8]][2] # 2077.9256
#' 
#' # 144 us
#' out <- mzion:::gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod,
#'                                 ntmass = ntmass, ctmass = ctmass, amods = amods, 
#'                                 vmods_nl = NULL, fmods_nl = NULL, 
#'                                 mod_indexes = mod_indexes)
#' 
#' # (8b)
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'             "Gln->pyro-Glu (N-term = Q)")
#' 
#' fixedmods <- sort(fixedmods)
#' varmods <- sort(varmods)
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'                        
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 4L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' 
#' aa_seq <- "EKNALVNEADSADVLQVANTDDEGGPENHRENFNNNNNNSVAVSSLNNGR"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[3]][1] # 5824.7551
#' 
#' out <- mzion:::gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod,
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = amods, 
#'                                 vmods_nl = NULL, fmods_nl = NULL, 
#'                                 mod_indexes = mod_indexes)
#' 
#' ## same site at multiple variable modifications
#' fixedmods <- c("TMT6plex (N-term)", "Carbamidomethyl (C)")
#' varmods   <- c("Acetyl (K)", "TMT6plex (K)", "Gln->pyro-Glu (N-term = Q)")
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep  <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 5L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#' ntmod <- attr(aa_masses, "ntmod")
#' ctmod <- attr(aa_masses, "ctmod")
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
#' 
#' aa_seq <- "HQGVMKVGMGQKMNK"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[4]][3]
#' 
#' out <- mzion:::gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                         aa_masses = aa_masses, 
#'                                         ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                         ntmod = ntmod, ctmod = ctmod, 
#'                                         ntmass = ntmass, ctmass = ctmass, 
#'                                         amods = amods, 
#'                                         vmods_nl = vmods_nl, fmods_nl = NULL, 
#'                                         mod_indexes = mod_indexes)
#' }
gen_ms2ions_a1_vnl0_fnl0 <- function (aa_seq, ms1_mass = NULL, aa_masses = NULL, 
                                      ms1vmods = NULL, ms2vmods = NULL, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, 
                                      vmods_nl = NULL, fmods_nl = NULL, # not used
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      
                                      # dummy
                                      maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                                      
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods[["inds"]]
  ms2vmods <- ms2vmods[oks]
  
  # filtered by `ms1_mass`; may be no match as
  # `idxes` may be beyond `maxn_vmods_sitescombi_per_pep`
  idxes <- check_ms1_mass_vmods(ms2vmods = ms2vmods, aam = aam, 
                                aa_masses = aa_masses, 
                                ntmod = ntmod, ctmod = ctmod, 
                                ms1_mass = ms1_mass)
  ms2vmods <- ms2vmods[idxes]
  
  if (!length(ms2vmods)) 
    return(NULL)
  
  # most likely a list-one
  # `[1]` in the rare case of multiple combinations have the same `ms1_mass`
  ms2vmods <- ms2vmods[[1]]

  if (attr(ms2vmods, "single")) {
    P <- find_vmodposU(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M <- attr(P, "mods", exact = TRUE)
    
    af <- calc_ms2ions_a1_vnl0_fnl0(
      M = M, P = P, aam = aam, aa_masses = aa_masses, 
      ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
      mod_indexes = mod_indexes, digits = digits)
  }
  else {
    P <- find_vmodposM(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M <- attr(P, "mods", exact = TRUE)

    af <- lapply(split_matrix(M, by = "row"), calc_ms2ions_a1_vnl0_fnl0, 
                 P = P, aam = aam, aa_masses = aa_masses, 
                 ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
                 mod_indexes = mod_indexes, digits = digits)
    af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
  }

  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(af, av)
}


#' Helper for the calculation of MS2 ion series.
#' 
#' @param M A modification matrix or vector.
#' @param P A matrix of position permutations.
#' @param mod_indexes Modification indexes.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_var_masses
calc_ms2ions_a1_vnl0_fnl0 <- function (M, P, aam, aa_masses, ntmass, ctmass, 
                                       type_ms2ions = "by", mod_indexes, 
                                       digits = 4L) 
{
  ds  <- aa_masses[M]
  nvm <- nrow(P)
  out <- vector("list", nvm)
  
  naa <- length(aam)
  hx0 <- rep_len("0", naa)

  for (i in 1:nvm) {
    vi <- P[i, ]
    aam_i <- aam
    aam_i[vi] <- aam_i[vi] + ds
    out[[i]] <- ms2ions_by_type(aam_i, ntmass, ctmass, type_ms2ions, digits)
    
    h <- hx0
    h[vi] <- mod_indexes[M]
    names(out)[i] <- .Internal(paste0(list(h), collapse = "", recycle0 = FALSE))
  }

  out
}


#' Checks the MS1 mass for proceeding with MS2 matches.
#'
#' Maybe missed if a "match" is beyond `maxn_vmods_sitescombi_per_pep`.
#'
#' 'bare + fixed terminals (aa_masses["N-term"] + aa_masses["C-term"]) +
#' variable terminals + anywhere'.
#' 
#' @param ms2vmods Lists of variable modifications. 
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param tol The tolerance in mass.
#' @inheritParams calc_ms2ions_a1_vnl0_fnl0
#' @return A logical vector.
check_ms1_mass_vmods <- function (ms2vmods, aam, aa_masses, ntmod, ctmod, 
                                  ms1_mass = NULL, tol = 1e-3) 
{
  if (is.null(ms1_mass))
    return(FALSE)

  bare <- sum(aam) + aa_masses["N-term"] + aa_masses["C-term"]
  ok_n <- length(ntmod)
  ok_c <- length(ctmod)

  # No need of is_empty(ntmod) && is_empty(ctmod)
  delta <- if (!(ok_n || ok_c))
    0
  else if (ok_n && ok_c)
    aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  else if (ok_n)
    aa_masses[names(ntmod)]
  else if (ok_c)
    aa_masses[names(ctmod)]
  
  bd <- bare + delta
  
  len <- length(ms2vmods)
  ans <- vector("logical", len)
  
  for (i in 1:len) {
    vi <- ms2vmods[[i]]
    
    ans[i] <- if (length(vi))
      (abs(bd + sum(aa_masses[vi]) - ms1_mass) <= tol)
    else
      FALSE
  }

  ans
}


#' Calculates the masses of MS2 ion series.
#'
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#' 
#' @rdname gen_ms2ions_base
#' 
#' @param aa_masses An amino-acid mass lookup.
#' @param amods The attribute \code{amods} from a \code{aa_masses}.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#'  @rdname gen_ms2ions_base
#'
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' # (12) "amods+ tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (K)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' 
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 2L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' ms1_mass <- ms1_masses$mass[[2]][2] # 2041.8958
#'
#' out <- mzion:::gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod,
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = amods, 
#'                                 vmods_nl = NULL, fmods_nl = fmods_nl, 
#'                                 mod_indexes = mod_indexes)
#' 
#' # No "M", no "S"
#' aa_seq <- "HQGVNVGGQKN"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[2]][2] # 1367.6996
#' 
#' ## same site at multiple variable modifications
#' fixedmods <- c("Oxidation (M)", "dHex (S)")
#' varmods <- c("Acetyl (K)", "TMT6plex (K)", "Deamidated (N)")
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' 
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 5L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#' 
#' aa_seq <- "HQGVMKVGMGQKMNSK"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' ms1_mass <- ms1_masses$mass[[5]][2] # 2041.8958
#' 
#' out <- mzion:::gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                         aa_masses = aa_masses, 
#'                                         ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                         ntmod = ntmod, ctmod = ctmod,
#'                                         ntmass = ntmass, ctmass = ctmass, 
#'                                         amods = amods, 
#'                                         vmods_nl = NULL, fmods_nl = fmods_nl, 
#'                                         mod_indexes = mod_indexes)
#' }
gen_ms2ions_a1_vnl0_fnl1 <- function (aa_seq = NULL, ms1_mass = NULL, 
                                      aa_masses = NULL, 
                                      ms1vmods = NULL, ms2vmods = NULL, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, 
                                      vmods_nl = NULL, # not used
                                      fmods_nl = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      maxn_fnl_per_seq = 3L, 
                                      
                                      # dummy
                                      maxn_vnl_per_seq = 3L, 
                                      digits = 4L) 
{
  if (maxn_fnl_per_seq < 2L)
    return(
      gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                               aa_masses = aa_masses, 
                               ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
                               ntmod = ntmod, ctmod = ctmod, 
                               ntmass = ntmass, ctmass = ctmass, 
                               amods = amods, mod_indexes = mod_indexes, 
                               type_ms2ions = type_ms2ions, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               maxn_vmods_sitescombi_per_pep = 
                                 maxn_vmods_sitescombi_per_pep, 
                               digits = digits))

  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  # (no pep_seq dispatching by fmod residues -> possible no matched sites)
  sites <- names(fmods_nl)
  # pattern <- .Internal(paste(list(sites), sep = " ", collapse = "|", recycle0 = FALSE))
  pattern <- if (length(sites) > 1L)
    .Internal(paste0(list(sites), collapse = "|", recycle0 = FALSE))
  else
    sites
  
  if (!grepl(pattern, aa_seq)) 
    return(
      gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                               aa_masses = aa_masses, 
                               ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
                               ntmod = ntmod, ctmod = ctmod, 
                               ntmass = ntmass, ctmass = ctmass, 
                               amods = amods, mod_indexes = mod_indexes, 
                               type_ms2ions = type_ms2ions, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               maxn_vmods_sitescombi_per_pep = 
                                 maxn_vmods_sitescombi_per_pep, 
                               digits = digits))
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods[["inds"]]
  ms2vmods <- ms2vmods[oks]
  
  idxes <- check_ms1_mass_vmods(ms2vmods = ms2vmods, aam = aam, 
                                aa_masses = aa_masses, 
                                ntmod = ntmod, ctmod = ctmod, 
                                ms1_mass = ms1_mass)
  ms2vmods <- ms2vmods[idxes]
  
  if (!length(ms2vmods)) 
    return(NULL)
  
  fnl_idxes <- .Internal(which(aas %fin% names(fmods_nl)))
  fmods_combi <- aas[fnl_idxes]
  
  if (length(fmods_combi) == 1L) {
    fnls <- fmods_nl[[fmods_combi]]
    len <- length(fnls)
    fnl_combi <- vector("list", len)
    
    for (i in 1:len) {
      fnl_combi[[i]] <- fnls[[i]]
      names(fnl_combi[[i]]) <- fmods_combi
    }
  }
  else {
    fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi], nmax = maxn_fnl_per_seq, 
                                  use.names = FALSE)
  }
  
  # most likely a list-one
  # `[1]` in case of multiple combinations have the same `ms1_mass`
  ms2vmods <- ms2vmods[[1]]

  if (attr(ms2vmods, "single")) {
    P <- find_vmodposU(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M <- attr(P, "mods", exact = TRUE)
    
    af <- calc_ms2ions_a1_vnl0_fnl1(
      M = M, P = P, fnl_combi = fnl_combi, 
      fnl_idxes = fnl_idxes, aam = aam, aa_masses = aa_masses, ntmass = ntmass, 
      ctmass = ctmass, type_ms2ions = type_ms2ions, mod_indexes = mod_indexes, 
      digits = digits)
  }
  else {
    P <- find_vmodposM(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M <- attr(P, "mods", exact = TRUE)

    af <- lapply(split_matrix(M, by = "row"), calc_ms2ions_a1_vnl0_fnl1, 
                 P = P, fnl_combi = fnl_combi, 
                 fnl_idxes = fnl_idxes, aam = aam, aa_masses = aa_masses, 
                 ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
                 mod_indexes = mod_indexes, digits = digits)
    
    af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
  }

  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(af, av)
}


#' Calculates
#'
#' @param M A modification matrix or vector.
#' @param fnl_idxes The position indexes of amino acids containing fixed neutral
#'   losses.
#' @param mod_indexes Modification indexes.
#' @inheritParams calc_ms2ions_a1_vnl0_fnl0
#' @inheritParams hms1_a0_vnl0_fnl1
calc_ms2ions_a1_vnl0_fnl1 <- function (M, P, fnl_combi, fnl_idxes, 
                                       aam, aa_masses, ntmass, ctmass, 
                                       type_ms2ions = "by", mod_indexes, 
                                       digits = 4L) 
{
  ds  <- aa_masses[M]
  nvm <- nrow(P)
  nnl <- length(fnl_combi)
  len <- nvm * nnl
  out <- vector("list", len)
  r <- 1L
  
  naa <- length(aam)
  hx0 <- rep_len("0", naa)
  
  for (i in 1:nvm) {
    vi <- P[i, ]
    aam_i <- aam
    aam_i[vi] <- aam_i[vi] + ds
    
    # the first fnl are all 0's
    out[[r]] <- ms2ions_by_type(aam_i, ntmass, ctmass, type_ms2ions, digits)
    r <- r + 1L
    
    if (nnl > 1L) {
      for (j in 2:nnl) {
        aam_j <- aam_i
        delta_nl <- .Internal(unlist(fnl_combi[[j]], recursive = FALSE, use.names = FALSE))
        aam_j[fnl_idxes] <- aam_j[fnl_idxes] - delta_nl
        out[[r]] <- ms2ions_by_type(aam_j, ntmass, ctmass, type_ms2ions, digits)
        r <- r + 1L
      }
    }

    h <- hx0
    h[vi] <- mod_indexes[M]
    h <- .Internal(paste0(list(h), collapse = "", recycle0 = FALSE))
    
    # Syntax: `(` for `vnl` and `[` for fnl
    names(out)[((i-1)*j+1L):(i*j)] <- 
      .Internal(paste0(list(h, 
                            " [", 
                            as.character(1:j), 
                            "]"), 
                       collapse = NULL, recycle0 = FALSE))
  }

  out
}


#' Calculates the masses of MS2 ion series.
#'
#' (9) "amods+ tmod- vnl+ fnl-", (10) "amods+ tmod+ vnl+ fnl-"
#' 
#' @rdname gen_ms2ions_base
#' 
#' 
#' @param amods The attribute \code{amods} from a \code{aa_masses}.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#'  @rdname gen_ms2ions_base
#'
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' # (10) "amods+ tmod+ vnl+ fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("dHex (S)", "Oxidation (M)", "Deamidated (N)", 
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#'
#' i <- 16L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide("HQGVMNVGMGQKMNS",
#'                                fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[16]][2] # 2197.9679
#'
#' out <- mzion:::gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod, 
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = amods, 
#'                                 vmods_nl = vmods_nl, fmods_nl = NULL, 
#'                                 mod_indexes = mod_indexes)
#' 
#' # Not in the category; 
#' # should be at least one `amod` with vnl+ 
#' # (aa_seq <- "HQGVVGGQK")
#' # (aa_seq <- "HQNGVVGGQK")
#' 
#' # Mismatches between `vmods_nl` and `aa_seq`
#' #  All of M, N, S should be present after pep_seq dispatching 
#' #    -> "correct" vmods_nl with all and only M, N, S (+/- tmods)
#' # (aa_seq <- "HQNGVVGGQKM") # no "S"
#' 
#' 
#' # (10) "amods+ tmod+ vnl+ fnl-"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                "Carbamidomethyl (C)")
#' varmods <- c("Acetyl (Protein N-term)", "Oxidation (M)", 
#'              "Carbamidomethyl (M)")
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#'
#' i <- 8L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
#' 
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' 
#' # No-matches at `ms1_mass = 2123.9424` as the matching index 
#' # exceeds "maxn_vmods_sitescombi_per_pep"
#' out <- mzion:::gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = 2123.9424, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod, 
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = amods, 
#'                                 vmods_nl = vmods_nl, fmods_nl = NULL,
#'                                 mod_indexes = mod_indexes)
#' 
#' ## same site at multiple variable modifications
#' fixedmods <- c("TMT6plex (K)", "dHex (S)")
#' varmods   <- c("Carbamidomethyl (M)", "Carbamyl (M)",
#'                "Deamidated (N)", "Acetyl (Protein N-term)")
#' mod_indexes <- seq_along(c(fixedmods, varmods)) |>
#'   as.hexmode() |>
#'   `names<-`(c(fixedmods, varmods))
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep  <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' 
#' i <- 12L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#' ntmod <- attr(aa_masses, "ntmod")
#' ctmod <- attr(aa_masses, "ctmod")
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
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
#' vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
#' 
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[i]][3] # 2135.9601
#' 
#' out <- mzion:::gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                         aa_masses = aa_masses, 
#'                                         ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                         ntmod = ntmod, ctmod = ctmod, 
#'                                         ntmass = ntmass, ctmass = ctmass, 
#'                                         amods = amods, 
#'                                         vmods_nl = vmods_nl, fmods_nl = NULL, 
#'                                         mod_indexes = mod_indexes)
#' }
gen_ms2ions_a1_vnl1_fnl0 <- function (aa_seq = NULL, ms1_mass = NULL, 
                                      aa_masses = NULL, 
                                      ms1vmods = NULL, ms2vmods = NULL, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, vmods_nl = NULL, 
                                      fmods_nl = NULL, # not used
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      
                                      # dummy
                                      maxn_fnl_per_seq = 3L, 
                                      maxn_vnl_per_seq = 3L, 
                                      digits = 4L) 
{
  if (maxn_vnl_per_seq < 2L)
    return(gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits))

  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods[["inds"]]
  ms2vmods <- ms2vmods[oks]
  
  idxes <- check_ms1_mass_vmods(ms2vmods = ms2vmods, aam = aam, 
                                aa_masses = aa_masses, 
                                ntmod = ntmod, ctmod = ctmod, 
                                ms1_mass = ms1_mass)
  ms2vmods <- ms2vmods[idxes]
  
  if (!length(ms2vmods)) 
    return(NULL)
  
  # most likely a list-one
  # `[1]` in case of multiple combinations have the same `ms1_mass`
  ms2vmods <- ms2vmods[[1]]

  if (attr(ms2vmods, "single")) {
    P  <- find_vmodposU(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M  <- attr(P, "mods", exact = TRUE)
    nP <- nrow(P)
    
    if (nP >= maxn_vmods_sitescombi_per_pep)
      af <- calc_ms2ions_a1_vnl0_fnl0(
        M = M, P = P, aam = aam, aa_masses = aa_masses, 
        ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
        mod_indexes = mod_indexes, digits = digits)
    else {
      nnl <- min(maxn_vmods_sitescombi_per_pep %/% nP, maxn_vnl_per_seq)
      
      if (nnl <= 1L)
        af <- calc_ms2ions_a1_vnl0_fnl0(
          M = M, P = P, aam = aam, aa_masses = aa_masses, 
          ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
          mod_indexes = mod_indexes, digits = digits)
      else
        af <- calc_ms2ions_a1_vnl1_fnl0(
          N = expand_grid_rows(vmods_nl[ms2vmods], nmax = nnl, use.names = FALSE), 
          M = M, 
          P = P, 
          aam = aam, 
          aa_masses = aa_masses, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          type_ms2ions = type_ms2ions, 
          mod_indexes = mod_indexes, 
          digits = digits)
    }
  }
  else {
    P  <- find_vmodposM(M = ms2vmods, aas = aas, nmax = maxn_vmods_sitescombi_per_pep)
    M  <- attr(P, "mods", exact = TRUE)
    nP <- nrow(P)
    nM <- nrow(M)
    n1 <- nP * nM

    if (n1 > maxn_vmods_sitescombi_per_pep) {
      l <- maxn_vmods_sitescombi_per_pep %/% nP
      M <- M[1:l, ] # l >= 1L

      if (l == 1L) {
        af <- calc_ms2ions_a1_vnl0_fnl0(
          M = M, P = P, aam = aam, aa_masses = aa_masses, 
          ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
          mod_indexes = mod_indexes, digits = digits)
      }
      else {
        af <- lapply(split_matrix(M, by = "row"), calc_ms2ions_a1_vnl0_fnl0, 
                     P = P, aam = aam, aa_masses = aa_masses, 
                     ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
                     mod_indexes = mod_indexes, digits = digits)
        af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
      }
    }
    else {
      M <- split_matrix(M, by = "row")
      l <- maxn_vmods_sitescombi_per_pep  %/% n1
      
      if (l <= 1L) {
        af <- lapply(M, calc_ms2ions_a1_vnl0_fnl0, 
                     P = P, aam = aam, aa_masses = aa_masses, 
                     ntmass = ntmass, ctmass = ctmass, type_ms2ions = type_ms2ions, 
                     mod_indexes = mod_indexes, digits = digits)
        af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
      }
      else {
        n2 <- n1 * prod(lengths(vmods_nl))
        
        N <- if (n2 > maxn_vmods_sitescombi_per_pep)
          lapply(M, function (x) 
            expand_grid_rows(vmods_nl[x], nmax = 2L, use.names = FALSE))
        else
          lapply(M, function (x) 
            expand_grid_rows(vmods_nl[x], nmax = maxn_vnl_per_seq, use.names = FALSE))
        
        af <- mapply(
          calc_ms2ions_a1_vnl1_fnl0, 
          N, M, 
          MoreArgs = list(
            P = P, 
            aam = aam, 
            aa_masses = aa_masses, 
            ntmass = ntmass, 
            ctmass = ctmass, 
            type_ms2ions = type_ms2ions, 
            mod_indexes = mod_indexes, 
            digits = digits
          ), 
          SIMPLIFY = FALSE, 
          USE.NAMES = FALSE)
        
        af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
      }
    }
  }

  # hexcodes of the reversed entries are not yet reversed; 
  #   they are not used and for time efficiency just leave them as they are
  # names are `pep_ivmod`; NA is the indicator for reversed entries
  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(af, av)
}


#' Calculates MS2 ions.
#'
#' @param M A vector of modifications.
#' @param P A matrix of positions.
#' @param N Lists of combinations of neutral losses for corresponding \code{P}.
#'   Each list contains a table where each column corresponds to a set of
#'   neutral loss. The first column corresponds to the combination without NLs.
#' @param mod_indexes Modification indexes.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_var_masses
#' @inheritParams matchMS
calc_ms2ions_a1_vnl1_fnl0 <- function (N, M, P, aam, aa_masses, 
                                       ntmass, ctmass, 
                                       type_ms2ions = "by", mod_indexes, 
                                       digits = 4L) 
{
  ds  <- aa_masses[M]
  nnl <- length(N)
  nvm <- nrow(P)
  len <- nvm * nnl
  out <- vector("list", len)
  naa <- length(aam)
  hx0 <- rep_len("0", naa)
  r <- 1L

  for (i in 1:nvm) {
    vi <- P[i, ]
    aam_i <- aam
    aam_i[vi] <- aam_i[vi] + ds
    
    for (j in 1:nnl) {
      aam_j <- aam_i
      delta_nl <- .Internal(unlist(N[[j]], recursive = FALSE, use.names = FALSE))
      aam_j[vi] <- aam_j[vi] - delta_nl
      out[[r]] <- ms2ions_by_type(aam_j, ntmass, ctmass, type_ms2ions, digits)
      r <- r + 1L
    }

    # both i and j must exist
    h <- hx0
    h[vi] <- mod_indexes[M]
    h <- .Internal(paste0(list(h), collapse = "", recycle0 = FALSE))
    
    # Syntax: `(` for `vnl` and `[` for fnl
    names(out)[((i-1)*j+1L):(i*j)] <- 
      .Internal(paste0(list(h, 
                            " (", 
                            as.character(1:j), 
                            ")"), 
                       collapse = NULL, recycle0 = FALSE))
  }
  
  out
}


