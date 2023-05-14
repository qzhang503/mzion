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
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
                              maxn_fnl_per_seq = 64L, maxn_vnl_per_seq = 64L, 
                              maxn_vmods_sitescombi_per_pep = 64L, 
                              
                              digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  l <- length(aas)
  nm <- .Internal(paste0(list(rep("0", l)), collapse = "", recycle0 = FALSE))
  # currently no subsetting by ms1_mass
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
#' library(magrittr)
#' 
#' # (5) "amods- tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (N-term)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Acetyl (Protein N-term)")
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
#' out <- mzion:::gen_ms2ions_a0_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = NULL, 
#'                                 aa_masses = aa_masses, ntmod = NULL, ctmod = NULL, 
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = NULL, vmods_nl = NULL, fmods_nl = fmods_nl, 
#'                                 mod_indexes = mod_indexes)
#' 
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
                                      maxn_fnl_per_seq = 8L, 
                                      
                                      # dummy
                                      maxn_vnl_per_seq = 8L, 
                                      
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      digits = 4L) 
{
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  # (no pep_seq dispatching by Anywhere fmod residues -> possible no matched sites)
  
  sites <- names(fmods_nl)
  pattern <- .Internal(paste0(list(sites), collapse = "|", recycle0 = FALSE))
  
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
  len_a <- length(aas)
  
  # At varmods "Oxidation (M)", pep_seq(s) must contain "M" 
  #   (with an additional entry of "Oxidation (M)" in aa_masses)
  # 
  # At fixedmods "Oxidation (M)", pep_seq(s) may not contain "M"; 
  #   (as `distri_peps` does not filter pep_seq by fixedmods)
  
  idxes <- .Internal(which(aas %in% names(fmods_nl)))
  
  if (length(idxes) > maxn_vmods_per_pep)
    idxes <- idxes[1:maxn_vmods_per_pep]
  
  # ---
  fmods_combi <- aas[idxes]
  names(fmods_combi) <- idxes
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi], nmax = maxn_fnl_per_seq, 
                                use.names = FALSE)
  len <- length(fnl_combi)
  av <- af <- vector("list", len)
  aam <- aa_masses[aas]
  af[[1]] <- af1 <- ms2ions_by_type(aam, ntmass, ctmass, type_ms2ions, digits)
  av[[1]] <- av1 <- calc_rev_ms2(af1, aas)
  
  if (len > 1L) {
    aami <- aam
    aamii <- aami[idxes]
    
    for (i in 2:len) {
      fnl_combi_i <- fnl_combi[[i]]
      
      aami[idxes] <- aamii - fnl_combi_i
      af[[i]] <- afi <- ms2ions_by_type(aami, ntmass, ctmass, type_ms2ions, digits)
      av[[i]] <- calc_rev_ms2(afi, aas)
    }
  }
  
  nm <- .Internal(paste0(list(rep("0", len_a)), collapse = "", recycle0 = FALSE))
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
#' library(magrittr)
#' 
#' # (8a) "amods+ tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("Deamidated (N)", "Carbamidomethyl (S)",
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
#' # 379 us
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
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
                                      maxn_fnl_per_seq = 64L, maxn_vnl_per_seq = 64L, 
                                      
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods$inds
  ms2vmods <- ms2vmods[oks]
  
  vmods_combi <- find_vmodscombi(aas = aas, ms2vmods = ms2vmods, 
                                 maxn_vmods_sitescombi_per_pep = 
                                   maxn_vmods_sitescombi_per_pep)
  
  # filtered by `ms1_mass`; may be no match as
  # `idxes` may be beyond `maxn_vmods_sitescombi_per_pep`
  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aam = aam, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  vmods_combi <- vmods_combi[idxes]
  
  if (!length(vmods_combi))
    return(NULL)
  
  # outputs
  af <- lapply(vmods_combi, calc_ms2ions_a1_vnl0_fnl0, 
               aam = aam, 
               aa_masses = aa_masses, 
               ntmass = ntmass, 
               ctmass = ctmass, 
               type_ms2ions = type_ms2ions, 
               digits = digits)
  
  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(add_hexcodes(af, vmods_combi, length(aas), mod_indexes), av)
}


#' Helper for the calculation of MS2 ion series.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_var_masses
calc_ms2ions_a1_vnl0_fnl0 <- function (vmods_combi, aam, aa_masses, 
                                       ntmass, ctmass, type_ms2ions = "by", 
                                       digits = 4L) 
{
  delta <- aa_masses[vmods_combi]
  idxes <- as.integer(names(vmods_combi))
  aam[idxes] <- aam[idxes] + delta
  
  ms2ions_by_type(aam, ntmass, ctmass, type_ms2ions, digits)
}


#' Checks the MS1 mass for proceeding with MS2 matches.
#'
#' Maybe missed if a "match" is beyond `maxn_vmods_sitescombi_per_pep`.
#'
#' 'bare + fixed terminals (aa_masses["N-term"] + aa_masses["C-term"]) +
#' variable terminals + anywhere'.
#'
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param tol The tolerance in mass.
#' @inheritParams calc_ms2ions_a1_vnl0_fnl0
#' @importFrom purrr is_empty
#' @return A logical vector.
check_ms1_mass_vmods2 <- function (vmods_combi, aam, aa_masses, ntmod, ctmod, 
                                   ms1_mass, tol = 1e-3) 
{
  len <- length(vmods_combi)
  
  if (len && !is.null(ms1_mass)) {
    # (direct addition may be faster than introducing a new argument:
    # ftmass <- unname(aa_masses["N-term"] + aa_masses["C-term"]))
    
    # bare <- sum(aam) + 18.010565
    bare <- sum(aam) + aa_masses["N-term"] + aa_masses["C-term"]
    
    len_n <- length(ntmod)
    len_c <- length(ctmod)
    
    # No need of is_empty(ntmod) && is_empty(ctmod)
    delta <- if (!(len_n || len_c))
      0
    else if (len_n && len_c)
      aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
    else if (len_n)
      aa_masses[names(ntmod)]
    else if (len_c)
      aa_masses[names(ctmod)]
    
    bd <- bare + delta
    
    idxes <- vector("logical", len)
    
    for (i in 1:len) {
      vmods_combi_i <- vmods_combi[[i]]
      vmass <- sum(aa_masses[vmods_combi_i])
      ok_mass <- bd + vmass
      
      idxes[i] <- (abs(ms1_mass - ok_mass) <= tol)
    }
  }
  else {
    idxes <- NULL
  }
  
  idxes
}


#' Adds hex codes (without NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @param ms2ions A series of MS2 ions with masses.
#' @param len The number of amino acid residues for the sequence indicated in
#'   \code{ms2ions}.
#' @inheritParams calc_aamasses
#' @inheritParams ms2match
add_hexcodes <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) 
{
  hexs <- rep("0", len)
  rows <- lapply(ms2ions, function (x) !is.null(x))
  rows <- .Internal(unlist(rows, recursive = FALSE, use.names = FALSE))
  
  vmods_combi <- vmods_combi[rows]
  ms2ions <- ms2ions[rows]
  
  hexs2 <- lapply(vmods_combi, function (x) {
    idxes <- .Internal(unlist(x, recursive = FALSE, use.names = FALSE))
    nms <- names(x)
    hexs[as.integer(nms)] <- mod_indexes[idxes]
    .Internal(paste0(list(hexs), collapse = "", recycle0 = FALSE))
  })
  
  names(ms2ions) <- hexs2
  
  invisible(ms2ions)
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
#' library(magrittr)
#' 
#' # (12) "amods+ tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (K)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
                                      maxn_fnl_per_seq = 64L, 
                                      
                                      # dummy
                                      maxn_vnl_per_seq = 64L, 
                                      
                                      digits = 4L) 
{
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  # (no pep_seq dispatching by fmod residues -> possible no matched sites)
  sites <- names(fmods_nl)
  pattern <- .Internal(paste(list(sites), sep = " ", collapse = "|", recycle0 = FALSE))
  
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
  oks <- ms1vmods$inds
  ms2vmods <- ms2vmods[oks]
  
  vmods_combi <- find_vmodscombi(aas = aas, ms2vmods = ms2vmods, 
                                 maxn_vmods_sitescombi_per_pep = 
                                   maxn_vmods_sitescombi_per_pep)
  
  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aam = aam, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  vmods_combi <- vmods_combi[idxes]
  
  if (!length(vmods_combi)) 
    return(NULL)
  
  # NLs of fixedmods
  fnl_idxes <- .Internal(which(aas %in% names(fmods_nl)))
  fmods_combi <- aas[fnl_idxes]
  names(fmods_combi) <- fnl_idxes
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi], nmax = maxn_fnl_per_seq)

  # go through each vmods_combi
  af <- lapply(vmods_combi, 
               calc_ms2ions_a1_vnl0_fnl1, 
               fnl_combi, fnl_idxes, aam, aa_masses, ntmass, ctmass, 
               type_ms2ions, digits = digits)
  
  af <- mapply(
    add_hexcodes_fnl2, 
    ms2ions = af, vmods_combi = vmods_combi, 
    MoreArgs = list(len = length(aas), mod_indexes  = mod_indexes), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE)
  
  af <- flatten_list(af)
  
  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(af, av)
}


#' Calculates
#'
#' @param fnl_idxes The position indexes of amino acids containing fixed neutral
#'   losses.
#' @inheritParams calc_ms2ions_a1_vnl0_fnl0
#' @inheritParams hms1_a0_vnl0_fnl1
calc_ms2ions_a1_vnl0_fnl1 <- function (vmods_combi, fnl_combi, fnl_idxes, 
                                       aam, aa_masses, 
                                       ntmass, ctmass, type_ms2ions = "by", 
                                       digits = 4L) 
{
  # updates amod masses
  delta_amod <- aa_masses[vmods_combi]
  amod_idxes <- as.numeric(names(vmods_combi))
  aam[amod_idxes] <- aam[amod_idxes] + delta_amod
  
  # updates fnl masses
  len <- length(fnl_combi)
  out <- vector("list", len)
  
  for (i in 1:len) {
    aam_i <- aam
    delta_nl <- .Internal(unlist(fnl_combi[[i]], recursive = FALSE, use.names = FALSE))
    aam_i[fnl_idxes] <- aam_i[fnl_idxes] - delta_nl
    out[[i]] <- ms2ions_by_type(aam_i, ntmass, ctmass, type_ms2ions, digits)
  }
  
  invisible(out)
}


#' Adds hex codes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @inheritParams add_hexcodes
#' @inheritParams calc_ms2ions_a1_vnl0_fnl1
#' @inheritParams ms2match
add_hexcodes_fnl2 <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) 
{
  idxes <- .Internal(unlist(vmods_combi, recursive = FALSE, use.names = FALSE))
  nms <- names(vmods_combi)
  
  hex_mods = rep("0", len)
  hex_mods[as.numeric(nms)] <- mod_indexes[idxes]
  hex_mods <- .Internal(paste0(list(hex_mods), collapse = "", recycle0 = FALSE))
  
  # Syntax: `(` for `vnl` and `[` for fnl
  names(ms2ions) <- .Internal(paste0(list(hex_mods, 
                                          " [", 
                                          as.character(seq_along(ms2ions)), 
                                          "]"), 
                                     collapse = NULL, recycle0 = FALSE))
  
  invisible(ms2ions)
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
#' library(magrittr)
#' 
#' # (10) "amods+ tmod+ vnl+ fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("dHex (S)", "Oxidation (M)", "Deamidated (N)", 
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
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
#' # exceeds `maxn_vmods_sitescombi_per_pep`
#' out <- mzion:::gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = 2123.9424, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod, 
#'                                 ntmass = ntmass, ctmass = ctmass, 
#'                                 amods = amods, 
#'                                 vmods_nl = vmods_nl, fmods_nl = NULL,
#'                                 mod_indexes = mod_indexes)
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
                                      maxn_fnl_per_seq = 64L, 
                                      
                                      maxn_vnl_per_seq = 64L, 
                                      digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aam <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods[["inds"]]
  ms2vmods <- ms2vmods[oks]
  
  vmods_combi <- find_vmodscombi(aas = aas, ms2vmods = ms2vmods, 
                                 maxn_vmods_sitescombi_per_pep = 
                                   maxn_vmods_sitescombi_per_pep)
  
  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aam = aam, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  vmods_combi <- vmods_combi[idxes]
  
  if (!length(vmods_combi)) 
    return(NULL)
  
  vnl_combi <- 
    lapply(vmods_combi, 
           function (x) expand_grid_rows(vmods_nl[x], nmax = maxn_vnl_per_seq))

  # theoretical MS2 of forward sequences
  af <- mapply(
    calc_ms2ions_a1_vnl1_fnl0, 
    vmods_combi = vmods_combi, 
    vnl_combi = vnl_combi, 
    MoreArgs = list(
      aam = aam, 
      aa_masses = aa_masses, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      type_ms2ions = type_ms2ions, 
      digits = digits
    ), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE)
  
  af <- mapply(
    add_hexcodes_vnl2, 
    ms2ions = af, 
    vmods_combi = vmods_combi, 
    MoreArgs = list(
      len = length(aas), 
      mod_indexes  = mod_indexes), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE)
  
  af <- .Internal(unlist(af, recursive = FALSE, use.names = TRUE))
  
  if (length(af) > maxn_vmods_sitescombi_per_pep) 
    af <- af[1:maxn_vmods_sitescombi_per_pep]
  
  # hexcodes of the reversed entries are not yet reversed; 
  #  they are not used and for time efficiency just leave them as are
  # names are `pep_ivmod`; NA is the indicator for reversed entries
  av <- lapply(af, calc_rev_ms2, aas)
  names(av) <- NA_character_
  c(af, av)
}


#' Calculates MS2 ions.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @param vnl_combi Lists of combinations of neutral losses for corresponding
#'   \code{vmods_combi}. Each list contains a table where each column
#'   corresponds to a set of neutral loss. The first column corresponds to the
#'   combination without NLs.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_var_masses
calc_ms2ions_a1_vnl1_fnl0 <- function (vmods_combi, vnl_combi, aam, aa_masses, 
                                       ntmass, ctmass, type_ms2ions = "by", 
                                       digits = 4L) 
{
  # updates vmod masses
  delta_amod <- aa_masses[vmods_combi]
  idxes <- as.integer(names(vmods_combi))
  aam[idxes] <- aam[idxes] + delta_amod
  
  # updates vnl masses
  len <- length(vnl_combi)
  out <- vector("list", len)
  
  # the first vnl masses are always all zeros
  out[[1]] <- ms2ions_by_type(aam, ntmass, ctmass, type_ms2ions, digits)
  
  if (len > 1L) {
    for (i in 2:len) {
      aam_i <- aam
      delta_nl <- .Internal(unlist(vnl_combi[[i]], recursive = FALSE, use.names = FALSE))
      aam_i[idxes] <- aam_i[idxes] - delta_nl
      out[[i]] <- ms2ions_by_type(aam_i, ntmass, ctmass, type_ms2ions, digits)
    }
  }
  
  invisible(out)
}


#' Adds hexcodes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @inheritParams add_hexcodes
#' @inheritParams ms2match
add_hexcodes_vnl2 <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) 
{
  nms <- names(vmods_combi)
  
  hexs <- rep("0", len)
  hexs[as.integer(nms)] <- mod_indexes[vmods_combi]
  hexs <- .Internal(paste0(list(hexs), collapse = "", recycle0 = FALSE))
  
  # Syntax: `(` for `vnl` and `[` for fnl
  names(ms2ions) <- .Internal(paste0(list(hexs, 
                                          " (", 
                                          as.character(seq_along(ms2ions)), 
                                          ")"), 
                                     collapse = NULL, recycle0 = FALSE))
  invisible(ms2ions)
}


