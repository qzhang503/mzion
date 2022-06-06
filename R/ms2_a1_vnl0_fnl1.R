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
ms2match_a1_vnl0_fnl1 <- function (i, aa_masses, ms1vmods, ms2vmods, 
                                   ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, fmods_nl, 
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
  rm(list = c("tempdata"))
  
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
      "gen_ms2ions_a1_vnl0_fnl1", 
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
      "calc_ms2ions_a1_vnl0_fnl1", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "bions_base", "yions_base",
      "cions_base", "zions_base", 
      "aions_base", "xions_base", 
      "add_hexcodes_fnl2", 
      "search_mgf", 
      "find_ms2_bypep", 
      "fuzzy_match_one", 
      "fuzzy_match_one2", 
      "post_frame_adv"), 
    envir = environment(proteoM:::frames_adv)
  )
  
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
                    digits = digits, 
                    FUN = gen_ms2ions_a1_vnl0_fnl1), 
    .scheduling = "dynamic")

  parallel::stopCluster(cl)
  
  out <- dplyr::bind_rows(out)
  
  out <- post_ms2match(out, i, aa_masses, out_path)
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' 
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#' 
#' i <- 2L
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
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' ms1_mass <- ms1_masses$mass[[2]][2] # 2041.8958
#'
#' out <- gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
                                      digits = 4L) 
{
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  # (no pep_seq dispatching by fmod residues -> possible no matched sites)
  sites <- names(fmods_nl)
  
  pattern <- paste(sites, collapse = "|")
  
  if (!grepl(pattern, aa_seq)) 
    return(
      gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                               aa_masses = aa_masses, 
                               ms1vmods = NULL, ms2vmods = NULL, 
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

  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aas2 = aas2, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  if (!any(idxes)) return(NULL)
  
  vmods_combi <- vmods_combi[idxes]
  
  # NLs of fixedmods
  fnl_idxes <- which(aas %in% names(fmods_nl))

  fmods_combi <- aas[fnl_idxes]
  names(fmods_combi) <- fnl_idxes
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi])
  
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
                                       ntmass, ctmass, type_ms2ions = "by", 
                                       digits = 4L) 
{
  # updates amod masses
  delta_amod <- aa_masses[vmods_combi]
  amod_idxes <- as.numeric(names(vmods_combi))
  aas2[amod_idxes] <- aas2[amod_idxes] + delta_amod
  
  # updates fnl masses
  len <- length(fnl_combi)
  
  out <- vector("list", len)

  for (i in 1:len) {
    aas2_i <- aas2
    delta_nl <- .Internal(unlist(fnl_combi[[i]], recursive = FALSE, use.names = FALSE))
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


