#' Matching MS2 ions.
#' 
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#' 
#' @param amods The attribute of \code{amods} from an \code{aa_masses}.
#' @rdname ms2match_base
#' @import purrr
#' @import parallel
#' @import dplyr
ms2match_a1_vnl0_fnl0 <- function (i, aa_masses, ms1vmods, ms2vmods, 
                                   ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, 
                                   mod_indexes, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 64L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 115L, df0 = NULL, digits = 4L) 
{
  tempdata <- purge_search_space(i, aa_masses, mgf_path, detect_cores(16L), ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  theopeps2 <- tempdata$theopeps2
  rm(list = c("tempdata"))
  
  if (!length(mgf_frames) || !length(theopeps)) {
    saveRDS(df0, file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")))
    return(df0)
  }
  
  n_cores <- detect_cores(32L)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), envir = environment(fastmatch::fmatch))
                          
  parallel::clusterExport(
    cl,
    c("frames_adv", 
      "frames_adv2", 
      "gen_ms2ions_a1_vnl0_fnl0", 
      "match_mvmods", 
      "expand_grid_rows", 
      "find_vmodscombi", 
      "combi_namesiteU", 
      "find_vmodposU", 
      "vec_to_list", 
      "sim_combn", 
      "combi_namesiteM", 
      "find_vmodposM", 
      "match_aas_indexes", 
      "check_ms1_mass_vmods2", 
      "calc_ms2ions_a1_vnl0_fnl0", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "bions_base", "yions_base",
      "cions_base", "zions_base", 
      "aions_base", "xions_base", 
      "add_hexcodes", 
      "search_mgf", 
      "find_ms2_bypep", 
      "fuzzy_match_one", 
      "fuzzy_match_one2", 
      "post_frame_adv"), 
    envir = environment(proteoM:::frames_adv)
  )

  if (is.null(theopeps2)) {
    out <- parallel::clusterMap(
      cl, frames_adv, 
      mgf_frames, theopeps, 
      MoreArgs = list(aa_masses = aa_masses, 
                      ms1vmods = ms1vmods, 
                      ms2vmods = ms2vmods, 
                      ntmod = ntmod, 
                      ctmod = ctmod, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = amods, 
                      vmods_nl = NULL, 
                      fmods_nl = NULL, 
                      mod_indexes = mod_indexes, 
                      type_ms2ions = type_ms2ions, 
                      maxn_vmods_per_pep = 
                        maxn_vmods_per_pep, 
                      maxn_sites_per_vmod = 
                        maxn_sites_per_vmod, 
                      maxn_vmods_sitescombi_per_pep = 
                        maxn_vmods_sitescombi_per_pep, 
                      minn_ms2 = minn_ms2, 
                      ppm_ms1 = ppm_ms1, 
                      ppm_ms2 = ppm_ms2, 
                      min_ms2mass = min_ms2mass, 
                      digits = digits, 
                      FUN = gen_ms2ions_a1_vnl0_fnl0), 
      .scheduling = "dynamic")
  }
  else {
    out <- parallel::clusterMap(
      cl, frames_adv2, 
      mgf_frames, theopeps, theopeps2, 
      MoreArgs = list(aa_masses = aa_masses, 
                      ms1vmods = ms1vmods, 
                      ms2vmods = ms2vmods, 
                      ntmod = ntmod, 
                      ctmod = ctmod, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = amods, 
                      vmods_nl = NULL, 
                      fmods_nl = NULL, 
                      mod_indexes = mod_indexes, 
                      type_ms2ions = type_ms2ions, 
                      maxn_vmods_per_pep = 
                        maxn_vmods_per_pep, 
                      maxn_sites_per_vmod = 
                        maxn_sites_per_vmod, 
                      maxn_vmods_sitescombi_per_pep = 
                        maxn_vmods_sitescombi_per_pep, 
                      minn_ms2 = minn_ms2, 
                      ppm_ms1 = ppm_ms1, 
                      ppm_ms2 = ppm_ms2, 
                      min_ms2mass = min_ms2mass, 
                      digits = digits), 
      .scheduling = "dynamic")
  }
  
  
  parallel::stopCluster(cl)
  
  out <- dplyr::bind_rows(out)
  
  out <- post_ms2match(out, i, aa_masses, out_path)
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
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'                        
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#' 
#' i <- 8L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (is_empty(ctmod)) {
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
#' out <- gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, 
#'                                 ms1vmods = ms1vmods, ms2vmods = ms2vmods, 
#'                                 ntmod = ntmod, ctmod = ctmod,
#'                                 ntmass = ntmass, ctmass = ctmass, amods = amods, 
#'                                 vmods_nl = NULL, fmods_nl = NULL, 
#'                                 mod_indexes = mod_indexes)
#' 
#' # (8b)
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
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
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'                        
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#' 
#' i <- 4L
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
#' 
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#' 
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#' 
#' if (is_empty(ctmod)) {
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
#' out <- gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
                                      maxn_vmods_sitescombi_per_pep = 64L, 
                                      digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aas2 <- aa_masses[aas]
  
  ms1vmods <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
  oks <- ms1vmods$inds
  ms2vmods <- ms2vmods[oks]
  
  vmods_combi <- find_vmodscombi(aas = aas, ms2vmods = ms2vmods, 
                                 maxn_vmods_sitescombi_per_pep = 
                                   maxn_vmods_sitescombi_per_pep)

  # filtered by `ms1_mass`; may be no match as
  # `idxes` may be beyond `maxn_vmods_sitescombi_per_pep`
  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aas2 = aas2, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  if (!any(idxes)) 
    return(NULL)
  
  vmods_combi <- vmods_combi[idxes]
  
  # outputs
  len <- length(vmods_combi)
  out <- vector("list", len)
  
  for (i in 1:len) 
    out[[i]] <- calc_ms2ions_a1_vnl0_fnl0(vmods_combi = vmods_combi[[i]], 
                                          aas2 = aas2, 
                                          aa_masses = aa_masses, 
                                          ntmass = ntmass, 
                                          ctmass = ctmass, 
                                          type_ms2ions = type_ms2ions, 
                                          digits = digits)
  
  out <- add_hexcodes(out, vmods_combi, length(aas), mod_indexes)

  invisible(out)
}


#' Helper for the calculation of MS2 ion series.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_fixvar_masses
#' @inheritParams ms2match_base
calc_ms2ions_a1_vnl0_fnl0 <- function (vmods_combi, aas2, aa_masses, 
                                       ntmass, ctmass, type_ms2ions = "by", 
                                       digits = 4L) 
{
  # mass delta
  delta <- aa_masses[vmods_combi]
  
  idxes <- as.numeric(names(vmods_combi))
  aas2[idxes] <- aas2[idxes] + delta
  
  ms2ions_by_type(aas2, ntmass, ctmass, type_ms2ions, digits)
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
check_ms1_mass_vmods2 <- function (vmods_combi, aas2, aa_masses, ntmod, ctmod, 
                                   ms1_mass, tol = 1e-3) 
{
  len <- length(vmods_combi)
  
  if (len && !is.null(ms1_mass)) {
    # (direct addition may be faster than introducing a new argument:
    # ftmass <- unname(aa_masses["N-term"] + aa_masses["C-term"]))
    
    # bare <- sum(aas2) + 18.010565
    bare <- sum(aas2) + aa_masses["N-term"] + aa_masses["C-term"]
    
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
  hex_mods = rep("0", len)
  rows <- lapply(ms2ions, function (x) !is.null(x))
  rows <- .Internal(unlist(rows, recursive = FALSE, use.names = FALSE))
  
  vmods_combi <- vmods_combi[rows]
  ms2ions <- ms2ions[rows]
  
  hex_mods2 <- lapply(vmods_combi, function (x) {
    idxes <- .Internal(unlist(x, recursive = FALSE, use.names = FALSE))
    nms <- names(x)
    
    hex_mods[as.numeric(nms)] <- mod_indexes[idxes]
    hex_mods <- .Internal(paste0(list(hex_mods), collapse = "", recycle0 = FALSE))
  })
  
  names(ms2ions) <- hex_mods2
  
  invisible(ms2ions)
}

