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
                                   maxn_fnl_per_seq = 64L, 
                                   maxn_vmods_sitescombi_per_pep = 64L, 
                                   minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                                   min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                                   df0 = NULL, digits = 4L) 
{
  n_splits <- detect_cores(ceiling(16L*6L/minn_ms2))
  
  tempd <- purge_search_space(i, aa_masses = aa_masses, mgf_path = mgf_path, 
                              n_cores = n_splits, ppm_ms1 = ppm_ms1)

  mgf_frames <- tempd$mgf_frames
  theopeps <- tempd$theopeps
  rm(list = c("tempd", "n_splits"))
  
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
                    FUN = gen_ms2ions_a1_vnl0_fnl1), 
    .scheduling = "dynamic")

  parallel::stopCluster(cl)
  
  out <- dplyr::bind_rows(out)
  post_ms2match(out, i, aa_masses, out_path)
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
#' library(proteoM)
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
#' ms1vmods_all <- lapply(aa_masses_all, proteoM:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' 
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, proteoM:::make_ms2vmods))
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
#' out <- proteoM:::gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
  fnl_combi <- expand_grid_rows(fmods_nl[fmods_combi])
  
  if (length(fnl_combi) > maxn_fnl_per_seq)
    fnl_combi <- fnl_combi[1:maxn_fnl_per_seq]
  
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


