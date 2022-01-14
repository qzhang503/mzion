#' Matching MS2 ions.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#'
#' @param i Integer; the index for a set of corresponding aa_masses and
#'   theoretical peptides.
#' @param ms1vmods The set of all possible MS1 vmod labels at a given aa_masses.
#' @param ms2vmods Matrices of labels of variable modifications. Each
#'   permutation in a row for each matrix.
#' @param ntmass The mass of N-terminal.
#' @param ctmass The mass of C-terminal.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams add_fixvar_masses
#' @import purrr
#' @import parallel
ms2match_base <- function (i, aa_masses, ms1vmods, ms2vmods, ntmass, ctmass, 
                           mod_indexes, mgf_path, out_path, type_ms2ions = "by", 
                           maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L, 
                           maxn_vmods_sitescombi_per_pep = 32L, 
                           minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                           min_ms2mass = 110L, digits = 4L) 
{
  # note: split into 16^2 lists
  tempdata <- purge_search_space(i, aa_masses, mgf_path, detect_cores(16L), ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  theopeps2 <- tempdata$theopeps2
  rm(list = c("tempdata"))
  gc()
  
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
      "gen_ms2ions_base", 
      "ms2ions_by_type", 
      "byions", "czions", "axions", 
      "bions_base", "yions_base",
      "cions_base", "zions_base", 
      "aions_base", "xions_base", 
      "search_mgf2", 
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
                      amods = NULL, 
                      vmods_nl = NULL, 
                      fmods_nl = NULL, 
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
                      FUN = gen_ms2ions_base), 
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
                      amods = NULL, 
                      vmods_nl = NULL, 
                      fmods_nl = NULL, 
                      mod_indexes = mod_indexes, 
                      type_ms2ions = type_ms2ions, 
                      maxn_vmods_per_pep = maxn_vmods_per_pep, 
                      maxn_sites_per_vmod = maxn_sites_per_vmod, 
                      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
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


#' Frames advancement.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#'
#' @param mgf_frames A group of mgf frames. Each frame contains one to multiple
#'   MGFs whose MS1 masses are in the same interval.
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @param FUN A function pointer to, e.g., \link{gen_ms2ions_base}.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
frames_adv <- function (mgf_frames = NULL, theopeps = NULL, 
                        aa_masses = NULL, ms1vmods = NULL, ms2vmods = NULL, 
                        ntmod = NULL, ctmod = NULL, 
                        ntmass = NULL, ctmass = NULL, 
                        amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                        mod_indexes = NULL, 
                        type_ms2ions = "by", 
                        maxn_vmods_per_pep = 5L, 
                        maxn_sites_per_vmod = 3L, 
                        maxn_vmods_sitescombi_per_pep = 32L, 
                        minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                        min_ms2mass = 110L, digits = 4L, FUN) 
{
  len <- length(mgf_frames)
  out <- vector("list", len) 
  
  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr$frame[1]
  
  bf_idx <- 1L
  theos_bf_ms1 <- theopeps[[bf_idx]] 
  theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
  theomasses_bf_ms1 <- theos_bf_ms1$mass
  
  cr_idx <- bf_idx + 1L
  theos_cr_ms1 <- theopeps[[cr_idx]]
  theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
  theomasses_cr_ms1 <- theos_cr_ms1$mass
  
  theos_bf_ms2 <- mapply(
    FUN, 
    aa_seq = theopeps_bf_ms1, 
    ms1_mass = theomasses_bf_ms1, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmod = ntmod, ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  names(theos_bf_ms2) <- theopeps_bf_ms1
  
  theos_cr_ms2 <- mapply(
    FUN, 
    aa_seq = theopeps_cr_ms1, 
    ms1_mass = theomasses_cr_ms1, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmod = ntmod, ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  names(theos_cr_ms2) <- theopeps_cr_ms1
  
  ## --- iteration ---
  for (i in seq_len(len)) {
    exptmasses_ms1 <- mgfs_cr$ms1_mass
    exptmoverzs_ms2 <- mgfs_cr$ms2_moverz
    
    af_idx <- cr_idx + 1L
    
    theos_af_ms1 <- theopeps[[af_idx]]
    theopeps_af_ms1 <- theos_af_ms1$pep_seq
    theomasses_af_ms1 <- theos_af_ms1$mass
    
    theos_af_ms2 <- mapply(
      FUN, 
      aa_seq = theopeps_af_ms1, 
      ms1_mass = theomasses_af_ms1, 
      MoreArgs = list(
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmod = ntmod, ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    names(theos_af_ms2) <- theopeps_af_ms1
    
    # each `out` for the results of multiple mgfs in one frame
    
    out[[i]] <- mapply(
      search_mgf2, 
      expt_mass_ms1 = exptmasses_ms1, 
      expt_moverz_ms2 = exptmoverzs_ms2, 
      MoreArgs = list(
        theomasses_bf_ms1 = theomasses_bf_ms1, 
        theomasses_cr_ms1 = theomasses_cr_ms1, 
        theomasses_af_ms1 = theomasses_af_ms1, 
        theos_bf_ms2 = theos_bf_ms2, 
        theos_cr_ms2 = theos_cr_ms2, 
        theos_af_ms2 = theos_af_ms2, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    
    if (i == len) break
    
    # advance to the next frame
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr$frame[1]
    
    if (isTRUE(new_frame == (frame + 1L))) {
      cr_idx <- cr_idx + 1L
      
      theos_bf_ms1 <- theos_cr_ms1
      theomasses_bf_ms1 <- theomasses_cr_ms1
      theos_bf_ms2 <- theos_cr_ms2
      
      theos_cr_ms1 <- theos_af_ms1
      theomasses_cr_ms1 <- theomasses_af_ms1
      theos_cr_ms2 <- theos_af_ms2
    } 
    else if (isTRUE(new_frame == (frame + 2L))) {
      cr_idx <- cr_idx + 2L
      
      theos_bf_ms1 <- theos_af_ms1
      theomasses_bf_ms1 <- theomasses_af_ms1
      theos_bf_ms2 <- theos_af_ms2
      
      theos_cr_ms1 <- theopeps[[cr_idx]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_cr_ms2 <- mapply(
        FUN, 
        aa_seq = theopeps_cr_ms1, 
        ms1_mass = theomasses_cr_ms1, 
        MoreArgs = list(
          aa_masses = aa_masses, 
          ms1vmods = ms1vmods, 
          ms2vmods = ms2vmods, 
          ntmod = ntmod, ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(theos_cr_ms2) <- theopeps_cr_ms1
    } 
    else {
      cr_idx <- cr_idx + 3L
      bf_idx <- cr_idx - 1L
      
      theos_bf_ms1 <- theopeps[[bf_idx]]
      theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
      theomasses_bf_ms1 <- theos_bf_ms1$mass
      
      theos_cr_ms1 <- theopeps[[cr_idx]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- mapply(
        FUN, 
        aa_seq = theopeps_bf_ms1, 
        ms1_mass = theomasses_bf_ms1, 
        MoreArgs = list(
          aa_masses = aa_masses, 
          ms1vmods = ms1vmods, 
          ms2vmods = ms2vmods, 
          ntmod = ntmod, ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(theos_bf_ms2) <- theopeps_bf_ms1
      
      theos_cr_ms2 <- mapply(
        FUN, 
        aa_seq = theopeps_cr_ms1, 
        ms1_mass = theomasses_cr_ms1, 
        MoreArgs = list(
          aa_masses = aa_masses, 
          ms1vmods = ms1vmods, 
          ms2vmods = ms2vmods, 
          ntmod = ntmod, ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(theos_cr_ms2) <- theopeps_cr_ms1
    }
    
    frame <- new_frame
  }
  
  out <- post_frame_adv(out, mgf_frames)
}


#' Frames advancement (with cached \code{ms2masses}.
#'
#' Not yet used. readRDS takes about 30 secs and clusterApply(cl, x, readRDS)
#' takes about 16 seconds for one module and couple minutes for all twelve
#' modules. For small data sets, the gain in using pre-calcucated ms2masses did
#' not outweigh the loss in the loading of the large files. Needs a faster
#' parallel readRDS.
#'
#' @param mgf_frames A group of mgf frames. Each frame contains one to multiple
#'   MGFs whose MS1 masses are in the same interval.
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param theopeps2 Binned theoretical MS2 ions corresponding to an i-th
#'   \code{aa_masses}.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
frames_adv2 <- function (mgf_frames = NULL, theopeps = NULL, theopeps2 = NULL,
                         aa_masses = NULL, ms1vmods = NULL, ms2vmods = NULL, 
                         ntmod = NULL, ctmod = NULL, 
                         ntmass = NULL, ctmass = NULL, 
                         amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                         mod_indexes = NULL, 
                         type_ms2ions = "by", 
                         maxn_vmods_per_pep = 5L, 
                         maxn_sites_per_vmod = 3L, 
                         maxn_vmods_sitescombi_per_pep = 32L, 
                         minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                         min_ms2mass = 110L, digits = 4L) 
{
  len <- length(mgf_frames)
  out <- vector("list", len) 
  
  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr$frame[1]

  bf_idx <- 1L
  theos_bf_ms1 <- theopeps[[bf_idx]] 
  theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
  theomasses_bf_ms1 <- theos_bf_ms1$mass
  
  cr_idx <- bf_idx + 1L
  theos_cr_ms1 <- theopeps[[cr_idx]]
  theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
  theomasses_cr_ms1 <- theos_cr_ms1$mass
  
  # name order are different
  theos_bf_ms2 <- theopeps2[[as.character(frame - 1L)]]
  theos_cr_ms2 <- theopeps2[[as.character(frame)]]

  ## --- iteration ---
  for (i in seq_len(len)) {
    exptmasses_ms1 <- mgfs_cr$ms1_mass
    exptmoverzs_ms2 <- mgfs_cr$ms2_moverz
    
    af_idx <- cr_idx + 1L
    frame_af <- frame + 1L
    
    theos_af_ms1 <- theopeps[[af_idx]]
    theopeps_af_ms1 <- theos_af_ms1$pep_seq
    theomasses_af_ms1 <- theos_af_ms1$mass
    
    theos_af_ms2 <- theopeps2[[as.character(frame_af)]]

    out[[i]] <- mapply(
      search_mgf2, 
      expt_mass_ms1 = exptmasses_ms1, 
      expt_moverz_ms2 = exptmoverzs_ms2, 
      MoreArgs = list(
        theomasses_bf_ms1 = theomasses_bf_ms1, 
        theomasses_cr_ms1 = theomasses_cr_ms1, 
        theomasses_af_ms1 = theomasses_af_ms1, 
        theos_bf_ms2 = theos_bf_ms2, 
        theos_cr_ms2 = theos_cr_ms2, 
        theos_af_ms2 = theos_af_ms2, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    
    if (i == len) break
    
    # advance to the next frame
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr$frame[1]

    if (isTRUE(new_frame == (frame + 1L))) {
      cr_idx <- cr_idx + 1L
      
      theos_bf_ms1 <- theos_cr_ms1
      theomasses_bf_ms1 <- theomasses_cr_ms1
      theos_bf_ms2 <- theos_cr_ms2
      
      theos_cr_ms1 <- theos_af_ms1
      theomasses_cr_ms1 <- theomasses_af_ms1
      theos_cr_ms2 <- theos_af_ms2
    } 
    else if (isTRUE(new_frame == (frame + 2L))) {
      cr_idx <- cr_idx + 2L
      
      theos_bf_ms1 <- theos_af_ms1
      theomasses_bf_ms1 <- theomasses_af_ms1
      theos_bf_ms2 <- theos_af_ms2
      
      theos_cr_ms1 <- theopeps[[cr_idx]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_cr_ms2 <- theopeps2[[as.character(new_frame)]]
    } 
    else {
      cr_idx <- cr_idx + 3L
      bf_idx <- cr_idx - 1L
      
      theos_bf_ms1 <- theopeps[[bf_idx]]
      theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
      theomasses_bf_ms1 <- theos_bf_ms1$mass
      
      theos_cr_ms1 <- theopeps[[cr_idx]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- theopeps2[[as.character(new_frame - 1L)]]
      theos_cr_ms2 <- theopeps2[[as.character(new_frame)]]
    }
    
    frame <- new_frame
  }
  
  out <- post_frame_adv(out, mgf_frames)
}


#' Calculates the masses of MS2 ion series.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#' 
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' 
#' @seealso \link{bions_base}, \link{yions_base}.
#' 
#' @examples
#' \donttest{
#' # (2) "amods- tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
#'              "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses = aa_masses_all[[2]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647 # + proton
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses <- aa_masses_all[[1]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
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
                              maxn_vmods_sitescombi_per_pep = 32L, 
                              digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aas2 <- aa_masses[aas]
  
  out <- ms2ions_by_type(aas2, ntmass, ctmass, type_ms2ions, digits)
  
  len_a <- length(aas)
  nm <- .Internal(paste0(list(rep("0", len_a)), collapse = "", recycle0 = FALSE))

  out <- list(out)
  names(out) <- nm
  
  invisible(out)
}


#' Matches by indexes
#'
#' @param expt_mass_ms1 Numeric; the experimental MS1 mass.
#' @param expt_moverz_ms2 A numeric list; the experimental MS2 m/z's.
#' @param theomasses_bf_ms1 Numeric vector; the theoretical MS1 masses at the
#'   preceding \code{-1} frame.
#' @param theomasses_cr_ms1 Numeric vector; the theoretical MS1 masses at the
#'   current frame.
#' @param theomasses_af_ms1 Numeric vector; the theoretical MS1 masses at the
#'   following \code{+1} frame.
#' @param theos_bf_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   preceding \code{-1} frame.
#' @param theos_cr_ms2 Numeric vector; the theoretical MS2 m/z's at the current
#'   frame.
#' @param theos_af_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   following \code{+1} frame.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @examples
#' \donttest{
#' expt_ms2 <-
#'   c(1628,3179,7677,9129,13950,14640,18571,19201,19205,19830,19833,20454,
#'     20457,21030,21073,21077,21687,24644,25232,37146,42042,43910,43920,44811,
#'     44824,45298,47494,55080,55901,56677,59014,66693,72396,72402,72720,73043,
#'     82411,91067,91838,93101,95572,98301,98665,100270,102081,102305,102744,106013,
#'     107998,108102,113713,113898,115045,115140,117669,119131,120730,123859,124029,124200,
#'     126199,126208,126610,126693,126775,128157,129447,129603,132396,135402,135475,138158,
#'     140397,141566,141634,141702,142183,142580,144189,147799,147926,148678,148860,149911,
#'     149973,153047,155607,158520,158631,162612,162717,163346,169537,170401,171249,171344,
#'     178012,178620,181980,188455)
#'
#' theo_ms2 <-
#'   c(-26231,62754,105787,129278,151731,161552,174924,184489,196534,204867,212917,219771,
#'     236270,106013,129447,148679,163242,178619,187776,197630,203310,212976,219825,227451,
#'     234026,237018)
#' cr <- which(expt_ms2 %fin% theo_ms2)
#' pr <- which((expt_ms2-1) %fin% theo_ms2)
#' af <- which((expt_ms2+1) %fin% theo_ms2)
#' c(cr, pr, af)
#'
#' }
search_mgf2 <- function (expt_mass_ms1, expt_moverz_ms2, 
                         theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                         theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                         minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                         min_ms2mass = 110L) 
{
  # --- subsets from the `before` and the `after` by MS1 mass tolerance 
  d <- expt_mass_ms1 * ppm_ms1/1E6
  bf_allowed <- which(theomasses_bf_ms1 >= (expt_mass_ms1 - d)) # 2 us
  af_allowed <- which(theomasses_af_ms1 <= (expt_mass_ms1 + d)) # 1.9 us

  theomasses_bf_ms1 <- theomasses_bf_ms1[bf_allowed]
  theomasses_af_ms1 <- theomasses_af_ms1[af_allowed]
  
  theos_bf_ms2 <- theos_bf_ms2[bf_allowed]
  theos_af_ms2 <- theos_af_ms2[af_allowed]
  
  # --- find MS2 matches ---
  ans_bf <- if (length(theos_bf_ms2)) 
    lapply(theos_bf_ms2, find_ms2_bypep, 
           expts = expt_moverz_ms2, 
           ppm_ms2 = ppm_ms2, 
           min_ms2mass = min_ms2mass, 
           minn_ms2 = minn_ms2)
  else 
    theos_bf_ms2
  
  ans_cr <- if (length(theos_cr_ms2)) 
    lapply(theos_cr_ms2, find_ms2_bypep, 
           expts = expt_moverz_ms2, 
           ppm_ms2 = ppm_ms2, 
           min_ms2mass = min_ms2mass, 
           minn_ms2 = minn_ms2)
  else 
    theos_cr_ms2
  
  ans_af <- if (length(theos_af_ms2)) 
    lapply(theos_af_ms2, find_ms2_bypep, 
           expts = expt_moverz_ms2, 
           ppm_ms2 = ppm_ms2, 
           min_ms2mass = min_ms2mass, 
           minn_ms2 = minn_ms2)
  else 
    theos_af_ms2
  
  ans <- c(ans_bf, ans_cr, ans_af)
  
  ## cleans up
  # (1) within a list: removes vmods+ positions that are NULL (< minn_ms2)
  # (no effects on vmods-; need `type` info if to limit to vmods+)
  oks <- lapply(ans, function (this) {
    oks <- lapply(this, function (x) !is.null(x$theo))
    .Internal(unlist(oks, recursive = FALSE, use.names = FALSE))
  })
  
  # USE.NAMES = TRUE 
  # (lapply loses names by `[[` whereas map2 reserves names when available)
  ans <- mapply(function (x, y) x[y], ans, oks, 
                SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  # (2)  removes empty lists
  oks2 <- lapply(ans, function(x) length(x) > 0L)
  oks2 <- .Internal(unlist(oks2, recursive = FALSE, use.names = FALSE))
  ans <- ans[oks2]
  
  theomasses_ms1 <- c(theomasses_bf_ms1, 
                      theomasses_cr_ms1, 
                      theomasses_af_ms1)
  
  theomasses_ms1 <- theomasses_ms1[oks2]

  # ---
  ans <- mapply(
    function (x, y) {
      attr(x, "theo_ms1") <- y
      x
    }, 
    ans, theomasses_ms1, 
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE)

  # ---
  # `length(ans) == N(theos_peps)` within the ppm window
  # 
  # ATIPIFFDMMLCEYQR
  # (1) ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # (2) ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # $KADEQMESMTYSTER
  # ...
  
  ## No evidence of M
  # 
  # $ATIPIFFDMMLCEYQR
  # $ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  # 
  # $ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  
  invisible(ans)
}


#' Helper: matches between theoretical and experimental MS2 ions.
#'
#' Run time about 45 us.
#'
#' @param expts Numeric vector; one series of experimental MS2s.
#' @param theos Numeric vector; one to multiple series of theoretical MS2s.
#' @inheritParams frames_adv
#' @inheritParams load_mgfs
#' @importFrom purrr map
#' @importFrom fastmatch fmatch %fin% 
#' @examples 
#' \donttest{
#' # experimental `322.18704` fit to both b- and y- ions
#' expts <- c(101.07140,102.05540,107.04956,110.07165,111.07500,
#'            112.08729,113.07130,115.08693,116.07092,120.08105,
#'            121.08428,126.12794,127.12494,127.13121,128.12825,
#'            128.13452,129.13164,129.13789,130.06493,130.09772,
#'            130.13495,130.14124,131.13831,132.14160,134.07635,
#'            136.07573,157.10826,158.09238,159.09117,170.12067,
#'            173.12846,173.14980,175.11896,176.12245,176.15956,
#'            184.11806,186.15303,188.15977,190.09743,193.63914,
#'            207.11292,210.12289,229.16670,230.17030,231.17410,
#'            235.10779,240.14383,248.18073,262.15305,265.67874,
#'            273.21210,285.13416,301.20700,305.16055,312.17740,
#'            314.69138,322.18704,365.23856,369.24496,371.70316,
#'            374.18283,376.27573,392.19308,393.23337,394.23743,
#'            399.20920,400.27567,400.72501,401.22650,401.27139,
#'            401.58405,401.72778,409.21918,410.22241,423.24564,
#'            433.29709,452.27066,462.25473,480.26578,481.26828,
#'            498.27670,530.35126,572.28278,573.28625,599.33899,
#'            600.34039,609.32202,626.29218,627.29303,627.33948,
#'            628.33417,629.30463,630.30219,643.31946,644.31763,
#'            644.35913,645.32520,646.32880,647.32825,648.32892)
#' 
#' theos <- c(114.0913,251.1503,322.1874,423.2350,480.2565,
#'            627.3249,783.4260,175.1190,322.1874,379.2088,
#'            480.2565,551.2936,688.3525,801.4366)
#' 
#' pep <- "ABCDEFG"
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' names(theos) <- c(nms, rev(nms))
#' theos <- list(`0000000` = theos)
#' 
#' x <- find_ms2_bypep(theos, expts)
#' 
#' ## Both expts 74953 and 74955 fit to theos 74954
#' pep <- "DIAVEEDLSSTPLFKDLLALMR"
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' theos <- c(158.0448,271.1288,342.1660,441.2344,570.2770,699.3196,
#'            814.3465,927.4306,1094.4289,1181.4610,1282.5086,1379.5614,
#'            1492.6455,1639.7139,1996.9718,2111.9987,2225.0828,2338.1668,
#'            2409.2040,2522.2880,2653.3285,2809.4296,175.1190,306.1594,
#'            419.2435,490.2806,603.3647,716.4487,831.4757,1188.7336,
#'            1335.8020,1448.8861,1545.9388,1646.9865,1734.0185,1901.0169,
#'            2014.1010,2129.1279,2258.1705,2387.2131,2486.2815,2557.3186,
#'            2670.4027,2785.4296)
#' names(theos) <- c(nms, rev(nms))
#' theos <- list(`0000000070000000000000 (1)` = theos)
#' 
#' expts <- c(126.12768, 127.13107, 128.12868, 128.13484, 130.13541, 
#'            130.14117, 136.07610, 167.08173, 178.27960, 228.13425, 
#'            230.17053, 238.11922, 248.18092, 256.12946, 257.12473, 
#'            276.10172, 278.15823, 283.14044, 321.21228, 327.13010, 
#'            358.71371, 368.22955, 376.27615, 394.73138, 396.15076, 
#'            396.22488, 400.27673, 407.24142, 414.16251, 426.76105, 
#'            445.29724, 447.31360, 458.76022, 486.28998, 490.79019, 
#'            491.29218, 500.26715, 514.27716, 516.33514, 522.78906, 
#'            523.28961, 528.33496, 535.27863, 543.27417, 549.31561, 
#'            572.32416, 572.39752, 572.82666, 576.35638, 603.31360, 
#'            603.36700, 604.36945, 619.36670, 631.30743, 632.30969, 
#'            641.41919, 642.42230, 646.36780, 647.36938, 658.36029, 
#'            675.03418, 675.36835, 681.03772, 681.37146, 690.27856, 
#'            701.34332, 707.45813, 707.69342, 708.02728, 716.41888, 
#'            716.45227, 726.35089, 726.85413, 760.40997, 787.48865, 
#'            788.45386, 789.45624, 796.44757, 813.47363, 814.47614, 
#'            815.48041, 874.51715, 906.91064, 907.41144, 916.51337, 
#'            918.38928, 919.39429, 926.55829, 953.56543, 971.57239, 
#'            972.57562, 1031.47534, 1044.57398, 1069.55029, 1085.45569, 
#'            1143.64014, 1144.64404, 1214.67590, 1321.77869, 1322.78186)
#' 
#' x <- find_ms2_bypep(theos, expts)
#' }
find_ms2_bypep <- function (theos = NULL, expts = NULL, ppm_ms2 = 25L, 
                            min_ms2mass = 110L, minn_ms2 = 6L) 
{
  ##############################################################################
  # `theos` may be empty: 
  #   e.g. the matched one is after `maxn_vmods_sitescombi_per_pep` 
  #   and never get matched.
  # 
  # length(theos) may be greater than one with site permutation and/or NLs.
  #  
  # ex: `expts` in integers
  # th_i: the i-th `theos` in integers
  # 
  # ex has no duplicated entries; 
  # th_i might.
  # 
  # Forward matching: match(theos, expts)
  # (i) allowed, e.g., b4- and y5 theo ions matched to the same ex value: 
  #   match(c(2,2,3,4), c(1:2, 5:10))
  # (ii) multiple ex' value's to the same th_i value not allowed; 
  #   otherwise longer length `c(expts[bps], expts[yps])` than lhs.
  #   e.g. ex's 74953, 74955 both fit to th_i 74954 and the best one is applied.
  #   (after a th_i is matched, it will be removed from further matching)
  # 
  # Backward matching: match(expts, theos)
  # (i) %in% and %fin% only shows the first match for duplicated entries th_i:
  #   match(1:4, c(1, 2, 2, 5))
  #   (so no worry about th_i duplication)
  ##############################################################################
  
  len <- length(theos)

  if (!len) 
    return(list(theo = NULL, expt = NULL))
  
  d <- ppm_ms2/1e6
  ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
  
  # ---
  out <- vector("list", len)
  
  for (i in 1:len) {
    theos_i <- theos[[i]]
    th_i <- ceiling(log(theos_i/min_ms2mass)/log(1+d))
    
    # forward matches
    mi <- th_i %fin% ex
    bf <- (th_i - 1L) %fin% ex
    af <- (th_i + 1L) %fin% ex
    ps <- mi | bf | af

    # backward matches
    n_ps <- sum(ps)
    
    if(n_ps >= minn_ms2) {
      
      # (i) place held by theoretical values to indicate matches
      # (ii) replaced the matches by experimental values
      
      es <- theos_i
      es[!ps] <- NA
      
      # separated b and y matches (to handled double-dipping between b and y)
      
      ## bps <- fuzzy_match_one2(ex, th_i[1:mid])
      lth <- length(ps)
      mid <- lth/2L
      
      ex_bf <- ex - 1L
      ex_af <- ex + 1L
      
      y_1 <- th_i[1:mid]
      mi_1 <- ex %fin% y_1
      bf_1 <- ex_bf %fin% y_1
      af_1 <- ex_af %fin% y_1
      ps_1 <- mi_1 | bf_1 | af_1
      
      ## yps <- fuzzy_match_one2(ex, th_i[(mid+1L):lth])
      y_2 <- th_i[(mid+1L):lth]
      mi_2 <- ex %fin% y_2
      bf_2 <- ex_bf %fin% y_2
      af_2 <- ex_af %fin% y_2
      ps_2 <- mi_2 | bf_2 | af_2
      
      expt_12 <- c(expts[ps_1], expts[ps_2])
      
      # (occur rarely; OK to recalculate freshly `expt_12`)
      if (n_ps != length(expt_12)) {
        ps_12 <- fuzzy_match_one2(ex, th_i)
        expt_12 <- expts[ps_12]
      }
      
      es[ps] <- expt_12
      
      out[[i]] <- list(theo = theos_i, expt = es)
    } 
    else
      out[[i]] <- list(theo = NULL, expt = NULL)
  }

  names(out) <- names(theos)
  
  out
}


#' Fuzzy matches with a +/-1 window.
#' 
#' Not used but called the codes inside directly.
#' 
#' @param x A vector to be matched.
#' @param y A vector to be matched against.
#' @importFrom fastmatch fmatch %fin% 
#' @examples 
#' ans1 <- fuzzy_match_one(c(74953, 74955), rep(74954, 2))
#' ans2 <- fuzzy_match_one(c(74953, 74955), 74954)
#' 
#' stopifnot(identical(ans1, ans2))
#' stopifnot(ans1 == c(TRUE, TRUE))
fuzzy_match_one <- function (x, y) 
{
  mi <- x %fin% y
  bf <- (x - 1L) %fin% y
  af <- (x + 1L) %fin% y
  
  ps <- mi | bf | af
}


#' Fuzzy matches with a +/-1 window.
#'
#' No multiple dipping of \code{y} matches. A \code{y} value will be removed (or
#' became 0) if matched,
#'
#' @param x A vector to be matched.
#' @param y A vector to be matched against.
#' @importFrom fastmatch fmatch %fin%
#' @examples
#' ans1 <- fuzzy_match_one2(c(74953, 74955), rep(74954, 2))
#' ans2 <- fuzzy_match_one2(c(74953, 74955), 74954)
#'
#' stopifnot(identical(ans1, ans2))
#' stopifnot(ans1 == c(FALSE, TRUE))
#'
#' ans3 <- fuzzy_match_one2(c(74953, 74955, 80000), c(74955, 80000))
fuzzy_match_one2 <- function (x, y) 
{
  mi <- x %fin% y
  if (any(mi)) y[y %fin% x[mi]] <- 0L
  
  x2 <- x - 1L
  bf <- x2 %fin% y
  if (any(bf)) y[y %fin% x2[bf]] <- 0L

  af <- (x + 1L) %fin% y
  
  ps <- mi | bf | af
}


