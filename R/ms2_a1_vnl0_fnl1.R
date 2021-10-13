#' Matching MS2 ions.
#' 
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#' 
#' @param amods The attribute of \code{amods} from an \code{aa_masses}.
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' 
#' @rdname ms2match_base
#' @import purrr
#' @import parallel
#' @import dplyr
ms2match_a1_vnl0_fnl1 <- function (i, aa_masses, ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, fmods_nl, 
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
  
  if (!length(mgf_frames) || !length(theopeps)) return(NULL)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), envir = environment(fastmatch::fmatch))
  
  # ms2_a1_vnl0_fnl1.R: (11, 12) "amods+ tmod+ vnl- fnl+", "amods+ tmod- vnl- fnl+"
  #   ms2match_a1_vnl0_fnl1 
  #     purge_search_space
  #     hms2_a1_vnl0_fnl1
  #       frames_adv_a1_vnl0_fnl1
  #         gen_ms2ions_a1_vnl0_fnl1
  #           combi_mvmods2 (ms2_a1_vnl0_fnl0.R)
  #             combi_vmods2 (ms2_a1_vnl0_fnl0.R)
  #           find_intercombi_p2 (ms2_a1_vnl0_fnl0.R)
  #           check_ms1_mass_vmods2 (ms2_a1_vnl0_fnl0.R)
  #           calc_ms2ions_a1_vnl0_fnl1
  #             ms2ions_by_type
  #               byions, czions, axions
  #           add_hexcodes_fnl2
  #         search_mgf2 (ms2base.R)
  #           find_ms2_bypep (ms2base.R)
  #             fuzzy_match_one (ms2base.R)
  #       post_frame_adv (ms2base.R)
  #     post_ms2match (utils_engine.R)
  
  parallel::clusterExport(
    cl,
    c("frames_adv", 
      "gen_ms2ions_a1_vnl0_fnl1", 
      "combi_mvmods2", 
      "combi_vmods2", 
      "find_intercombi_p2", 
      "check_ms1_mass_vmods2", 
      "calc_ms2ions_a1_vnl0_fnl1", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "add_hexcodes_fnl2", 
      "search_mgf2", 
      "find_mass_error_range", 
      "find_ms2_bypep", 
      "find_ms1_interval", 
      "fuzzy_match_one", 
      "fuzzy_match_one2", 
      "post_frame_adv"), 
    envir = environment(proteoM:::frames_adv)
  )

  out <- parallel::clusterMap(
    cl, hms2_a1_vnl0_fnl1, 
    mgf_frames, theopeps, 
    MoreArgs = list(aa_masses = aa_masses, 
                    ntmod = ntmod, 
                    ctmod = ctmod, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    amods = amods, 
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
    dplyr::bind_rows() %>% # across nodes
    post_ms2match(i, aa_masses, out_path)
  
  parallel::stopCluster(cl)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#' 
#' @inheritParams ms2match_a1_vnl0_fnl0
#' @rdname hms2_base
hms2_a1_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, 
                               ntmod = NULL, ctmod = NULL, 
                               ntmass, ctmass, 
                               amods, fmods_nl, mod_indexes, type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                               min_ms2mass = 110L, digits = 4L) {
  
  res <- frames_adv(mgf_frames = mgf_frames, 
                    theopeps = theopeps, 
                    aa_masses = aa_masses, 
                    ntmod = ntmod, 
                    ctmod = ctmod, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    amods = amods, 
                    vmods_nl = NULL, # not used
                    fmods_nl = fmods_nl, 
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
                    FUN = gen_ms2ions_a1_vnl0_fnl1)
  
  res <- post_frame_adv(res, mgf_frames)

  rm(list = "mgf_frames", "theopeps")
  
  invisible(res)
}


#' Calculates the masses of MS2 ion series.
#'
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
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
#' # (12) "amods+ tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (K)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' ms1_mass <- ms1_masses$mass[[2]][2] # 2041.8958
#'
#' out <- gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
#'                                 aa_masses = aa_masses, ntmod = ntmod, ctmod = ctmod,
#'                                 ntmass = ntmass, ctmass = ctmass, amods = amods, 
#'                                 vmods_nl = NULL, fmods_nl = fmods_nl, 
#'                                 mod_indexes = mod_indexes)
#' 
#' # No "M", no "S"
#' aa_seq <- "HQGVNVGGQKN"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[2]][2] # 1367.6996
#' 
#' }
gen_ms2ions_a1_vnl0_fnl1 <- function (aa_seq = NULL, ms1_mass = NULL, aa_masses = NULL, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, 
                                      vmods_nl = NULL, # not used
                                      fmods_nl = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 32L, 
                                      digits = 4L) {
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  # (no pep_seq dispatching by fmod residues -> possible no matched sites)
  sites <- names(fmods_nl)
  pattern <- paste(sites, collapse = "|")

  if (!grepl(pattern, aa_seq)) {
    out <- gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    return(out)
  }
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  aas <- .Internal(strsplit(aa_seq, "", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aas2 <- aa_masses[aas]

  vmods_combi <- combi_mvmods2(amods = amods, 
                               aas = aas, 
                               aa_masses = aa_masses, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               maxn_vmods_sitescombi_per_pep = 
                                 maxn_vmods_sitescombi_per_pep, 
                               digits = digits)
  
  vmods_combi <- find_intercombi_p2(vmods_combi, maxn_vmods_sitescombi_per_pep)

  # filtered by MS1 masses
  if (length(vmods_combi) && !is.null(ms1_mass)) {
    idxes <- lapply(vmods_combi, check_ms1_mass_vmods2, aas2, aa_masses, 
                    ntmod, ctmod, ms1_mass)
    idxes <- .Internal(unlist(idxes, recursive = FALSE, use.names = FALSE))
    
    vmods_combi <- vmods_combi[idxes]
    rm(list = c("idxes"))
  }

  # NLs of fixedmods
  fnl_idxes <- which(aas %in% names(fmods_nl))

  fmods_combi <- aas[fnl_idxes]
  names(fmods_combi) <- fnl_idxes
  fnl_combi <- expand.grid(fmods_nl[fmods_combi], KEEP.OUT.ATTRS = FALSE, 
                           stringsAsFactors = FALSE)
  
  # go through each vmods_combi
  out <- lapply(vmods_combi, 
                calc_ms2ions_a1_vnl0_fnl1, 
                fnl_combi, fnl_idxes, aas2, aa_masses, ntmass, ctmass, 
                type_ms2ions, digits = digits)
  
  out <- mapply(
    add_hexcodes_fnl2, 
    ms2ions = out, vmods_combi = vmods_combi, 
    MoreArgs = list(len = length(aas), mod_indexes  = mod_indexes), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE
  )
  
  purrr::flatten(out)
}


#' Calculates
#'
#' @param fnl_idxes The position indexes of amino acids containing fixed neutral
#'   losses.
#' @inheritParams calc_ms2ions_a1_vnl0_fnl0
#' @inheritParams hms1_a0_vnl0_fnl1
calc_ms2ions_a1_vnl0_fnl1 <- function (vmods_combi, fnl_combi, fnl_idxes, 
                                       aas2, aa_masses, 
                                       ntmass, ctmass, type_ms2ions, digits) {

  # updates amod masses
  delta_amod <- aa_masses[vmods_combi]
  amod_idxes <- as.numeric(names(vmods_combi))
  aas2[amod_idxes] <- aas2[amod_idxes] + delta_amod
  
  # updates fnl masses
  len <- nrow(fnl_combi)
  
  out <- vector("list", len)

  for (i in 1:len) {
    aas2_i <- aas2
    delta_nl <- .Internal(unlist(fnl_combi[i, ], recursive = FALSE, use.names = FALSE))
    aas2_i[fnl_idxes] <- aas2_i[fnl_idxes] - delta_nl
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }

  invisible(out)
}


#' Adds hex codes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @inheritParams add_hexcodes
#' @inheritParams calc_ms2ions_a1_vnl0_fnl1
add_hexcodes_fnl2 <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) {

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


