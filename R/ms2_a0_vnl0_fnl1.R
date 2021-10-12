#' Matching MS2 ions.
#' 
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @rdname ms2match_base
ms2match_a0_vnl0_fnl1 <- function (i, aa_masses, ntmass, ctmass, fmods_nl, 
                                   mod_indexes, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 32L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 110L, digits = 4L) {
  
  n_cores <- detect_cores()
  
  tempdata <- purge_search_space(i, aa_masses, mgf_path, n_cores, ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  rm(list = c("tempdata"))
  gc()
  
  if (!length(mgf_frames) || !length(theopeps)) return(NULL)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("%>%"), 
                          envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), 
                          envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), 
                          envir = environment(fastmatch::fmatch))
  
  # ms2_a0_vnl0_fnl1.R: (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
  #   ms2match_a0_vnl0_fnl1 
  #     purge_search_space
  #     hms2_a0_vnl0_fnl1
  #       frames_adv
  #         gen_ms2ions_a0_vnl0_fnl1
  #           // early return
  #           gen_ms2ions_base (ms2base.R)
  #             ms2ions_by_type (ion_ladder.R)
  #               byions, czions, axions (ion_ladder.R)
  #           // regular return
  #           ms2ions_by_type (ion_ladder.R)
  #             byions, czions, axions (ion_ladder.R)
  #         search_mgf2 (ms2base.R)
  #           find_ms2_bypep (ms2base.R)
  #             fuzzy_match_one (ms2base.R)
  #       post_frame_adv (ms2base.R)
  #     post_ms2match (utils_engine.R)
  
  parallel::clusterExport(
    cl,
    c("frames_adv", 
      "gen_ms2ions_a0_vnl0_fnl1", 
      "gen_ms2ions_base", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "search_mgf2", 
      "find_ms2_bypep", 
      "fuzzy_match_one", 
      "fuzzy_match_one2", 
      "post_frame_adv"), 
    envir = environment(proteoM:::frames_adv)
  )

  out <- parallel::clusterMap(
    cl, hms2_a0_vnl0_fnl1, 
    mgf_frames, theopeps, 
    MoreArgs = list(aa_masses = aa_masses, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    fmods_nl = fmods_nl, 
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
    .scheduling = "dynamic") %>% 
    dplyr::bind_rows() %>% 
    post_ms2match(i, aa_masses, out_path)
  
  parallel::stopCluster(cl)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#'
#' `res[[i]]` contains results for multiple mgfs within a frame (the number of
#' entries equals to the number of mgf frames).
#'
#' @inheritParams ms2match_base
#' @rdname hms2_base
hms2_a0_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, ntmass, ctmass, 
                               fmods_nl, mod_indexes, type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                               min_ms2mass = 110L, digits = 4L) {
  
  res <- frames_adv(mgf_frames = mgf_frames, 
                    theopeps = theopeps, 
                    aa_masses = aa_masses, 
                    ntmod = NULL, ctmod = NULL, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    amods = NULL, vmods_nl = NULL, fmods_nl = fmods_nl, 
                    mod_indexes = mod_indexes, 
                    type_ms2ions = type_ms2ions, 
                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                    maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
                    minn_ms2 = minn_ms2, 
                    ppm_ms1 = ppm_ms1, 
                    ppm_ms2 = ppm_ms2, 
                    min_ms2mass = min_ms2mass, 
                    digits = digits, 
                    FUN = gen_ms2ions_a0_vnl0_fnl1)
  
  res <- post_frame_adv(res, mgf_frames)

  rm(list = "mgf_frames", "theopeps")
  
  invisible(res)
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, add_varmasses = FALSE, 
#'                                add_nlmasses = FALSE)
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
#' 
gen_ms2ions_a0_vnl0_fnl1 <- function (aa_seq, ms1_mass = NULL, aa_masses, 
                                      ntmod = NULL, ctmod = NULL, # not used
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, vmods_nl = NULL, # not used
                                      fmods_nl = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 32L, 
                                      digits = 4L) {
  
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  # (no pep_seq dispatching by Anywhere fmod residues -> possible no matched sites)
  
  sites <- names(fmods_nl)
  pattern <- paste(sites, collapse = "|")
  
  if (!grepl(pattern, aa_seq)) {
    out <- gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                            aa_masses = aa_masses, 
                            ntmass = ntmass, ctmass = ctmass, 
                            mod_indexes = mod_indexes, 
                            type_ms2ions = type_ms2ions, 
                            maxn_vmods_per_pep = maxn_vmods_per_pep, 
                            maxn_sites_per_vmod = maxn_sites_per_vmod, 
                            maxn_vmods_sitescombi_per_pep = 
                              maxn_vmods_sitescombi_per_pep, 
                            digits = digits)
    return(out)
  }

  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  aas <- .Internal(strsplit(aa_seq, "", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  
  idxes <- which(aas %in% names(fmods_nl))
  
  # At varmods "Oxidation (M)", pep_seq(s) must contain "M" 
  #   (with an additional entry of "Oxidation (M)" in aa_masses)
  # 
  # At fixedmods "Oxidation (M)", pep_seq(s) may not contain "M"; 
  #   (as `distri_peps` does not filter pep_seq by fixedmods)
  
  len_i <- length(idxes)

  if (len_i > maxn_vmods_per_pep) {
    len_i <- maxn_vmods_per_pep
    idxes <- idxes[1:len_i]
  }
  
  # a faster way to check maxn_sites_per_vmod (but avoid slow R table)?
  # ...

  # ---
  fmods_combi <- aas[idxes]
  names(fmods_combi) <- idxes
  fnl_combi <- expand.grid(fmods_nl[fmods_combi], KEEP.OUT.ATTRS = FALSE, 
                           stringsAsFactors = FALSE)
  
  len <- nrow(fnl_combi)
  
  if (len > maxn_vmods_sitescombi_per_pep) {
    len <- maxn_vmods_sitescombi_per_pep
    fnl_combi <- fnl_combi[1:len, ]
  }
  
  out <- vector("list", len)
  
  aas2 <- aa_masses[aas]

  for (i in 1:len) {
    aas2_i <- aas2
    delta <- .Internal(unlist(fnl_combi[i, ], recursive = FALSE, use.names = FALSE))
    aas2_i[idxes] <- aas2_i[idxes] - delta
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }
  
  len_a <- length(aas)
  nm <- .Internal(paste0(list(rep(0, len_a)), collapse = "", recycle0 = FALSE))
  
  names(out) <- .Internal(paste0(list(nm, " [", seq_len(nrow(fnl_combi)), "]"), 
                                 collapse = NULL, recycle0 = FALSE))

  invisible(out)
}

