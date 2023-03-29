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
                                   maxn_fnl_per_seq = 64L, 
                                   maxn_vmods_sitescombi_per_pep = 64L, 
                                   minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                                   min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                                   df0 = NULL, digits = 4L) 
{
  tempd <- purge_search_space(i, aa_masses = aa_masses, mgf_path = mgf_path, 
                                 n_cores = detect_cores(16L), ppm_ms1 = ppm_ms1)

  mgf_frames <- tempd$mgf_frames
  theopeps <- tempd$theopeps
  rm(list = c("tempd"))
  gc()
  
  if (!length(mgf_frames) || !length(theopeps)) {
    qs::qsave(df0, file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")))
    return(df0)
  }
  
  n_cores <- detect_cores(96L)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), envir = environment(fastmatch::fmatch))

  parallel::clusterExport(
    cl,
    c("frames_adv", 
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
    envir = environment(mzion:::frames_adv)
  )

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
                    maxn_fnl_per_seq = 
                      maxn_fnl_per_seq, 
                    maxn_vmods_sitescombi_per_pep = 
                      maxn_vmods_sitescombi_per_pep, 
                    minn_ms2 = minn_ms2, 
                    ppm_ms1 = ppm_ms1, 
                    ppm_ms2 = ppm_ms2, 
                    min_ms2mass = min_ms2mass, 
                    index_mgf_ms2 = index_mgf_ms2, 
                    digits = digits, 
                    FUN = gen_ms2ions_a0_vnl0_fnl1), 
    .scheduling = "dynamic")

  parallel::stopCluster(cl)
  
  out <- dplyr::bind_rows(out)
  post_ms2match(out, i, aa_masses, out_path)
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
                                      maxn_fnl_per_seq = 64L, 
                                      
                                      # dummy
                                      maxn_vnl_per_seq = 64L, 
                                      
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
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi], use.names = FALSE)
  
  if ((len <- length(fnl_combi)) > maxn_fnl_per_seq) {
    fnl_combi <- fnl_combi[1:maxn_fnl_per_seq]
    len <- length(fnl_combi)
  }
  
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


