#' Matching MS2 ions.
#' 
#' (9) "amods+ tmod- vnl+ fnl-", (10) "amods+ tmod+ vnl+ fnl-"
#' 
#' @param amods The attribute of \code{amods} from an \code{aa_masses}.
#' @param vmods_nl The attribute of \code{vmods_nl} from an \code{aa_masses}.
#' 
#' @rdname ms2match_base
#' @import purrr
#' @import parallel
#' @import dplyr
ms2match_a1_vnl1_fnl0 <- function (i, aa_masses, ms1vmods, ms2vmods, 
                                   ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, vmods_nl, 
                                   mod_indexes, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 32L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 110L, digits = 4L) 
{
  tempdata <- purge_search_space(i, aa_masses, mgf_path, detect_cores(16L), ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  theopeps2 <- tempdata$theopeps2
  rm(list = c("tempdata"))
  
  if (!length(mgf_frames) || !length(theopeps)) 
    return(NULL)
  
  n_cores <- detect_cores(32L)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), envir = environment(fastmatch::fmatch))

  parallel::clusterExport(
    cl,
    c("frames_adv", 
      "frames_adv2", 
      "gen_ms2ions_a1_vnl1_fnl0", 
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
      "calc_ms2ions_a1_vnl1_fnl0", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "bions_base", "yions_base",
      "cions_base", "zions_base", 
      "aions_base", "xions_base", 
      "add_hexcodes_vnl2", 
      "search_mgf2", 
      "find_ms2_bypep", 
      # "fuzzy_match_one", 
      # "fuzzy_match_one2", 
      "post_frame_adv"), 
    envir = environment(proteoM:::frames_adv)
  )

  if (is.null(theopeps2)) {
    out <- parallel::clusterMap(
      cl, frames_adv, 
      mgf_frames = mgf_frames, theopeps = theopeps, 
      MoreArgs = list(aa_masses = aa_masses, 
                      ms1vmods = ms1vmods, 
                      ms2vmods = ms2vmods, 
                      ntmod = ntmod, 
                      ctmod = ctmod, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = amods, 
                      vmods_nl = vmods_nl, 
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
                      FUN = gen_ms2ions_a1_vnl1_fnl0), 
      .scheduling = "dynamic")
  }
  else {
    out <- parallel::clusterMap(
      cl, frames_adv2, 
      mgf_frames = mgf_frames, theopeps = theopeps, theopeps2 = theopeps2, 
      MoreArgs = list(aa_masses = aa_masses, 
                      ms1vmods = ms1vmods, 
                      ms2vmods = ms2vmods, 
                      ntmod = ntmod, 
                      ctmod = ctmod, 
                      ntmass = ntmass, 
                      ctmass = ctmass, 
                      amods = amods, 
                      vmods_nl = vmods_nl, 
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
#' # (9) "amods+ tmod+ vnl+ fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("dHex (S)", "Oxidation (M)", "Deamidated (N)", 
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
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
#' vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide("HQGVMNVGMGQKMNS",
#'                                fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[16]][2] # 2197.9679
#'
#' out <- gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, 
#'                                add_varmasses = FALSE, 
#'                                add_nlmasses = FALSE)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
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
#' out <- gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = 2123.9424, 
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
                                      maxn_vmods_sitescombi_per_pep = 32L, 
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

  idxes <- check_ms1_mass_vmods2(vmods_combi = vmods_combi, 
                                 aas2 = aas2, aa_masses = aa_masses, 
                                 ntmod = ntmod, ctmod = ctmod, 
                                 ms1_mass = ms1_mass)
  
  if (!any(idxes)) return(NULL)
  
  vmods_combi <- vmods_combi[idxes]
  vnl_combi <- lapply(vmods_combi, function (x) expand_grid_rows(vmods_nl[x]))

  ## --- (tentative) to restricts the total number of vnl_combi's
  # nrows <- lapply(vnl_combi, function (x) length(attributes(x)$row.names)) # faster than nrow
  # nrows <- .Internal(unlist(nrows, recursive = FALSE, use.names = FALSE))
  # counts <- cumsum(nrows)
  
  # oks <- which(counts <= maxn_vmods_sitescombi_per_pep)
  # vnl_combi <- vnl_combi[oks]
  # vmods_combi <- vmods_combi[oks]
  ## ---
  
  # 725 us
  out <- mapply(
    calc_ms2ions_a1_vnl1_fnl0, 
    vmods_combi = vmods_combi, 
    vnl_combi = vnl_combi, 
    MoreArgs = list(
      aas2 = aas2, 
      aa_masses = aa_masses, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      type_ms2ions = type_ms2ions, 
      digits = digits
    ), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE
  )
  
  # 360 us
  out <- mapply(
    add_hexcodes_vnl2, 
    ms2ions = out, 
    vmods_combi = vmods_combi, 
    MoreArgs = list(
      len = length(aas), 
      mod_indexes  = mod_indexes), 
    SIMPLIFY = FALSE, 
    USE.NAMES = FALSE
  )
  
  out <- .Internal(unlist(out, recursive = FALSE, use.names = TRUE))
  
  if (length(out) > maxn_vmods_sitescombi_per_pep) 
    out <- out[1:maxn_vmods_sitescombi_per_pep]
  
  invisible(out)
}


#' Calculates MS2 ions.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @param vnl_combi Lists of combinations of neutral losses for corresponding
#'   \code{vmods_combi}. Each list contains a table where each column
#'   corresponds to a set of neutral loss. The first column corresponds to the
#'   combination without NLs.
#' @inheritParams ms2ions_by_type
#' @inheritParams add_fixvar_masses
#' @inheritParams ms2match_base
calc_ms2ions_a1_vnl1_fnl0 <- function (vmods_combi, vnl_combi, aas2, aa_masses, 
                                       ntmass, ctmass, type_ms2ions, digits) 
{
  # updates vmod masses
  delta_amod <- aa_masses[vmods_combi]
  idxes <- as.numeric(names(vmods_combi))
  aas2[idxes] <- aas2[idxes] + delta_amod
  
  # updates vnl masses
  len <- length(vnl_combi)
  
  out <- vector("list", len)
  
  for (i in 1:len) {
    aas2_i <- aas2
    delta_nl <- .Internal(unlist(vnl_combi[[i]], recursive = FALSE, 
                                 use.names = FALSE))

    aas2_i[idxes] <- aas2_i[idxes] - delta_nl
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }

  invisible(out)
}


#' Adds hex codes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @param vmods_combi Lists of variable modifications.
#' @inheritParams add_hexcodes
add_hexcodes_vnl2 <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) 
{
  idxes <- .Internal(unlist(vmods_combi, recursive = FALSE, use.names = FALSE))
  nms <- names(vmods_combi)
  
  hex_mods = rep("0", len)
  hex_mods[as.numeric(nms)] <- mod_indexes[idxes]
  hex_mods <- .Internal(paste0(list(hex_mods), collapse = "", recycle0 = FALSE))

  # Syntax: `(` for `vnl` and `[` for fnl
  names(ms2ions) <- .Internal(paste0(list(hex_mods, 
                                          " (", 
                                          as.character(seq_along(ms2ions)), 
                                          ")"), 
                                     collapse = NULL, recycle0 = FALSE))

  invisible(ms2ions)
}


