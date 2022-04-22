#' Matching MS2 ions.
#' 
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @inheritParams gen_ms2ions_a1_vnl0_fnl1
#' @rdname ms2match_base
ms2match_a0_vnl0_fnl1 <- function (i, aa_masses, ms1vmods, ms2vmods, 
                                   ntmass, ctmass, 
                                   fmods_nl, mod_indexes, mgf_path, out_path, 
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
  gc()
  
  if (!length(mgf_frames) || !length(theopeps)) {
    qs::qsave(df0, file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")))
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
      "gen_ms2ions_a0_vnl0_fnl1", 
      "expand_grid_rows", 
      "gen_ms2ions_base", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "bions_base", "yions_base",
      "cions_base", "zions_base", 
      "aions_base", "xions_base", 
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
                      ntmod = NULL, 
                      ctmod = NULL, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = NULL, vmods_nl = NULL, fmods_nl = fmods_nl, 
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
                      FUN = gen_ms2ions_a0_vnl0_fnl1), 
      .scheduling = "dynamic")
  }
  else {
    out <- parallel::clusterMap(
      cl, frames_adv2, 
      mgf_frames, theopeps, theopeps2, 
      MoreArgs = list(aa_masses = aa_masses, 
                      ms1vmods = ms1vmods, 
                      ms2vmods = ms2vmods, 
                      ntmod = NULL, 
                      ctmod = NULL, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = NULL, vmods_nl = NULL, fmods_nl = fmods_nl, 
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
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @rdname gen_ms2ions_base
#' 
#' @examples 
#' \donttest{
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
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#' 
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' 
#' # variable `TMT6plex (N-term)` + `fixed Oxidation (M)`
#' # (additive varmod on top of fixedmod allowed)
#' 
#' out <- gen_ms2ions_a0_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = NULL, 
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
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  
  idxes <- which(aas %in% names(fmods_nl))
  
  # At varmods "Oxidation (M)", pep_seq(s) must contain "M" 
  #   (with an additional entry of "Oxidation (M)" in aa_masses)
  # 
  # At fixedmods "Oxidation (M)", pep_seq(s) may not contain "M"; 
  #   (as `distri_peps` does not filter pep_seq by fixedmods)
  
  len_i <- length(idxes)
  if (len_i > maxn_vmods_per_pep) 
    idxes <- idxes[1:maxn_vmods_per_pep]

  # ---
  fmods_combi <- aas[idxes]
  names(fmods_combi) <- idxes
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi])
  
  len <- length(fnl_combi)
  
  if (len > maxn_vmods_sitescombi_per_pep) 
    fnl_combi <- fnl_combi[1:maxn_vmods_sitescombi_per_pep, ]

  out <- vector("list", len)
  
  aas2 <- aa_masses[aas]

  for (i in 1:len) {
    aas2_i <- aas2
    delta <- .Internal(unlist(fnl_combi[[i]], recursive = FALSE, use.names = FALSE))
    aas2_i[idxes] <- aas2_i[idxes] - delta
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }
  
  len_a <- length(aas)
  nm <- .Internal(paste0(list(rep("0", len_a)), collapse = "", recycle0 = FALSE))
  
  names(out) <- .Internal(paste0(list(nm, " [", as.character(seq_len(len)), "]"), 
                                 collapse = NULL, recycle0 = FALSE))

  invisible(out)
}

