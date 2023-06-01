#' Pairs MGF queries to theoretical MS1 masses and peptide sequences.
#'
#' @param mgf_path The path to MGF files
#' @param n_modules The number of modules (\code{length(aa_masses_all)}) or one
#' @param .path_bin The path to binned theoretical masses
#' @param ppm_ms1_bin The tolerance in precursor mass error after mass binning.
#' @param by_modules Logical; if TRUE, results are saved with one mgf to one
#'   theo module. At FALSE, results are saved with one mgf paired to all theo
#'   modules
#' @inheritParams ms2match
#' @inheritParams matchMS
pair_mgftheo <- function (mgf_path, n_modules, .path_bin, by_modules = TRUE, 
                          reframe_mgfs = FALSE, min_mass = 200L, 
                          ppm_ms1_bin = 10L, first_search = FALSE)
{
  message("Pairing experimental and theoretical data.")
  
  tempfiles <- if (by_modules)
    list.files(mgf_path, pattern = "^expttheo_", full.names = TRUE)
  else
    list.files(mgf_path, pattern = "^mgftheo_",  full.names = TRUE)
  
  if (length(tempfiles))
    unlink(tempfiles)

  # MGFs (in data frame) split by frame indexes
  mgf_files <- list.files(mgf_path, pattern = "^mgf_queries_\\d+\\.rds$", 
                          full.names = TRUE)
  mgf_frames <- lapply(mgf_files, qs::qread)

  # for MGF calibrations
  if (first_search) {
    mgf_frames <- lapply(mgf_frames, function (x) {
      min_mgfmass <- min(x$ms1_mass, na.rm = TRUE)
      max_mgfmass <- max(x$ms1_mass, na.rm = TRUE)
      oks_min <- with(x, ms1_mass <= min_mgfmass + 10L)
      oks_max <- with(x, ms1_mass >= max_mgfmass - 10L)
      
      mgfa <- x[oks_max, ]
      mgfb <- x[oks_min, ]
      mgfc <- x[!(oks_max | oks_min), ]
      rows <- (1:nrow(mgfc)) %% 10L == 1L
      dplyr::bind_rows(mgfa, mgfc[rows, ], mgfb)
    })
  }
  
  mgf_frames <- dplyr::bind_rows(mgf_frames)
  
  if (reframe_mgfs) {
    mgf_frames[["frame"]] <- 
      find_ms1_interval(mgf_frames[["ms1_mass"]], from = min_mass, 
                        ppm = ppm_ms1_bin)
  }

  mgf_frames <- dplyr::group_by(mgf_frames, frame)
  mgf_frames <- dplyr::group_split(mgf_frames)
  fr_names   <- lapply(mgf_frames, function (x) x[["frame"]][1])
  names(mgf_frames) <- unlist(fr_names, recursive = FALSE, use.names = FALSE)
  
  # -> chunks: each chunk has multiple frames: each frame multiple precursors
  ranges <- seq_along(mgf_frames)
  
  n_chunks <- if (n_modules == 1L || by_modules)
    min(detect_cores(96L)^2, 1024L)
  else if (n_modules >= 96L)
    min(length(ranges), length(mgf_files) * n_modules * 2L)
  else
    min(length(ranges), length(mgf_files) * n_modules)
  
  labs   <- levels(cut(ranges, n_chunks))
  lower  <- floor(as.numeric( sub("\\((.+),.*", "\\1", labs)))
  grps   <- findInterval(ranges, lower)
  mgf_frames <- split(mgf_frames, grps)
  rm(list = c("fr_names", "ranges", "labs", "lower", "grps"))
  
  # (1) splits `theos` in accordance to `mgf_frames` with
  #     preceding and following frames: (o)|range of mgf_frames[[1]]|(o)
  mfrs <- lapply(mgf_frames, function (x) as.integer(names(x)))
  mins <- lapply(mfrs, function (x) if (length(x)) min(x, na.rm = TRUE) else 0L)
  mins <- .Internal(unlist(mins, recursive = FALSE, use.names = FALSE))
  maxs <- lapply(mfrs, function (x) if (length(x)) max(x, na.rm = TRUE) else 0L)
  maxs <- .Internal(unlist(maxs, recursive = FALSE, use.names = FALSE))
  
  anstheo <- vector("list", n_modules)
  
  for (i in seq_len(n_modules)) {
    theos <- 
      qs::qread(file.path(.path_bin, paste0("binned_theopeps_", i, ".rds")))
    
    if (is.null(theos))
      next
    
    theos <- lapply(theos, function (x) x[, c("pep_seq", "mass")])
    thfrs <- as.integer(names(theos))
    
    # separates into intervals (intersecting mgf_frames)
    anstheo[[i]] <- mapply(function (x, y) {
      theos[which(thfrs >= (x - 1L) & thfrs <= (y + 1L))]
    }, mins, maxs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  rm(list = c("theos", "mfrs", "thfrs", "mins", "maxs"))
  
  # (2) removes mgf_frames[[i]] not be found in anstheo[[i]]
  for (i in seq_along(mgf_frames)) {
    mi  <- mgf_frames[[i]]
    fmi <- names(mi)
    ti  <- lapply(anstheo, `[[`, i) # theos at chunk[[i]] (for all modules)
    
    # expt frames found in (any) theo module
    oks <- lapply(ti, function (x) fmi %in% names(x))
    oks <- Reduce(`|`, oks)
    mgf_frames[[i]] <- mi[oks]
  }
  rm(list = c("mi", "ti", "fmi", "oks"))
  
  # (3) removes unused frames of `anstheo`
  # (mgf_frames determines the length of each anstheo[[i]]; 
  #  more effective to first generate all bracketed mgf_frames indexes and apply
  #  the same set of indexes to each anstheo[[i]])
  for (i in seq_along(anstheo))
    anstheo[[i]] <- mapply(subset_theoframes, mgf_frames, anstheo[[i]], 
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  # (4) reverses the order (longer/heavier peptides towards the beginning)
  #     do the difficult ones first when paralleling with LB
  ord <- rev(seq_along(anstheo[[1]]))
  mgf_frames <- mgf_frames[ord]
  anstheo <- lapply(anstheo, function (x) x[ord])
  
  # (5) outputs (in chunks)
  if (by_modules) {
    for (i in seq_len(n_modules))
      qs::qsave(list(mgf_frames = mgf_frames, theopeps = anstheo[[i]]), 
                file.path(mgf_path, paste0("expttheo_", i, ".rds")), 
                preset = "fast")
  }
  else {
    for (i in seq_along(mgf_frames))
      qs::qsave(list(mgf_frames = mgf_frames[[i]], 
                     theopeps = lapply(anstheo, `[[`, i)), 
                file.path(mgf_path, paste0("mgftheo_", i, ".rds")), 
                preset = "fast")
  }
  
  invisible(NULL)
}


#' Help of \link{ms2match_all}
#'
#' By MGF chunks
#'
#' @param aa_masses_all All look-ups of \code{aa_masses}
#' @param funs_ms2 The functions in correspondence to \code{aa_masses_all} for
#'   calculating MS2 ion series.
#' @param ms1vmods_all All possible labels of MS1 variable modifications in
#'   correspondence to \code{aa_masses_all}
#' @param ms2vmods_all All possible labels of MS2 variable modifications in
#'   correspondence to \code{aa_masses_all}
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications
#' @param df0 An output template with zero rows
#' @inheritParams matchMS
#' @inheritParams pair_mgftheo
hms2match <- function (aa_masses_all, funs_ms2, ms1vmods_all, ms2vmods_all, 
                       mod_indexes, mgf_path, out_path, 
                       type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                       maxn_sites_per_vmod = 3L, maxn_fnl_per_seq = 3L, 
                       maxn_vnl_per_seq = 3L, 
                       maxn_vmods_sitescombi_per_pep = 64L, 
                       minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                       by_modules = FALSE, df0 = NULL, digits = 4L)
{
  pat   <- if (by_modules) "^expttheo_" else "^mgftheo_"
  mgths <- list.files(mgf_path, pattern = paste0(pat, "[0-9]+"))
  ord   <- order(as.integer(gsub(paste0(pat, "([0-9]+)\\.rds$"), "\\1", mgths)))
  mgths <- mgths[ord]
  
  message("\n===  MS2 ion searches started at ", Sys.time(), ". ===\n")
  
  n_cores <- detect_cores(96L) - 1L
  logs <- file.path(out_path, "temp", "log.txt")
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores), outfile = logs)
  parallel::clusterExport(cl, list("fmatch", "%fin%"), 
                          envir = environment(fastmatch::fmatch))
  parallel::clusterExport(
    cl,
    c("ms2match_all", "gen_ms2ions_base", 
      "gen_ms2ions_a0_vnl0_fnl1", "gen_ms2ions_a1_vnl0_fnl1", 
      "gen_ms2ions_a1_vnl0_fnl0", "gen_ms2ions_a1_vnl1_fnl0", 
      "calc_rev_ms2", "match_mvmods", "expand_grid_rows", "expand_gr", 
      "find_vmodposU", "find_vmodposM", 
      "vec_to_list", "split_matrix", "check_ms1_mass_vmods", 
      "calc_ms2ions_a1_vnl0_fnl0", "calc_ms2ions_a1_vnl0_fnl1", 
      "calc_ms2ions_a1_vnl1_fnl0", "ms2ions_by_type", 
      "byions", "czions", "axions", "bions_base", "yions_base", 
      "cions_base", "zions_base", "aions_base", "xions_base", 
      "find_ms2_bypep", "fuzzy_match_one", 
      "fuzzy_match_one2", "post_frame_adv"), 
    envir = environment(mzion::matchMS))
  
  if (by_modules) {
    parallel::clusterExport(cl, c("frames_adv"), envir = environment(mzion::matchMS))
    
    for (i in seq_along(aa_masses_all))
      ms2match_one(
        pep_mod_group = i, 
        aa_masses = aa_masses_all[[i]], 
        FUN = funs_ms2[[i]],  
        ms1vmods = ms1vmods_all[[i]], 
        ms2vmods = ms2vmods_all[[i]], 
        cl = cl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_fnl_per_seq = maxn_fnl_per_seq, 
        maxn_vnl_per_seq = maxn_vnl_per_seq, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, index_mgf_ms2 = index_mgf_ms2, 
        df0 = df0, digits = digits)
  }
  else {
    message("Check search progress at: ", logs)
    parallel::clusterExport(cl, c("mframes_adv"), envir = environment(mzion::matchMS))
    
    df <- parallel::clusterApplyLB(
      cl, mgths, ms2match_all, 
      aa_masses_all = aa_masses_all, 
      funs_ms2 = funs_ms2, 
      ms1vmods_all = ms1vmods_all, 
      ms2vmods_all = ms2vmods_all, 
      mod_indexes = mod_indexes, 
      mgf_path = mgf_path, 
      out_path = out_path, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_fnl_per_seq = maxn_fnl_per_seq, 
      maxn_vnl_per_seq = maxn_vnl_per_seq, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      minn_ms2 = minn_ms2, 
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2, 
      min_ms2mass = min_ms2mass, 
      index_mgf_ms2 = index_mgf_ms2, 
      df0 = df0, 
      digits = digits)
  }
  
  parallel::stopCluster(cl)
  
  message("\n===  MS2 ion searches completed at ", Sys.time(), ". ===\n")
  
  invisible(df)
}


#' Matches experimentals and theoreticals
#' 
#' One MGF chunk and all modules
#' 
#' @param mgth Experimental MGFs and paired theoreticals (from all modules)
#' @param df0 An output template
#' @inheritParams hms2match
ms2match_all <- function (mgth, aa_masses_all, funs_ms2, ms1vmods_all, 
                          ms2vmods_all, mod_indexes, mgf_path, out_path, 
                          type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                          maxn_sites_per_vmod = 3L, maxn_fnl_per_seq = 3L, 
                          maxn_vnl_per_seq = 3L, 
                          maxn_vmods_sitescombi_per_pep = 64L, 
                          minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                          min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                          df0 = NULL, digits = 4L)
{
  msg <- paste0("Matching expt-theo pair: ", mgth)
  write(msg, stdout())
  
  out_name   <- gsub("^mgftheo", "ion_matches", mgth)
  mgftheo    <- qs::qread(file.path(mgf_path, mgth))
  mgf_frames <- mgftheo[["mgf_frames"]]
  theopeps   <- mgftheo[["theopeps"]]
  rm("mgftheo")
  
  if (!length(mgf_frames)) {
    qs::qsave(df0, file.path(out_path, "temp", out_name))
    return(df0)
  }
  
  df <- mframes_adv(
    mgf_frames = mgf_frames, 
    theopeps = theopeps, 
    aa_masses_all = aa_masses_all, 
    funs_ms2 = funs_ms2, 
    ms1vmods_all = ms1vmods_all, 
    ms2vmods_all = ms2vmods_all, 
    mod_indexes = mod_indexes, 
    type_ms2ions = type_ms2ions, 
    maxn_vmods_per_pep = maxn_vmods_per_pep, 
    maxn_sites_per_vmod = maxn_sites_per_vmod, 
    maxn_fnl_per_seq = maxn_fnl_per_seq, 
    maxn_vnl_per_seq = maxn_vnl_per_seq, 
    maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
    minn_ms2 = minn_ms2, 
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
    min_ms2mass = min_ms2mass, 
    index_mgf_ms2 = index_mgf_ms2, 
    digits = digits)

  if (!dir.exists(tempdir <- file.path(out_path, "temp")))
    create_dir(tempdir)

  if (is.null(df)) {
    qs::qsave(df, file.path(tempdir, out_name), preset = "fast")
    return(NULL)
  }
  
  df <- dplyr::rename(df, 
                      pep_ret_range = ret_time, 
                      pep_scan_title = scan_title,
                      pep_exp_mz = ms1_moverz, 
                      pep_n_ms2 = ms2_n, 
                      pep_exp_mr = ms1_mass, 
                      pep_tot_int = ms1_int, 
                      pep_scan_num = scan_num, 
                      pep_exp_z = ms1_charge, 
                      pep_ms2_moverzs = ms2_moverz, 
                      pep_ms2_ints = ms2_int, 
                      pep_frame = frame)
  df[["pep_scan_num"]] <- as.character(df[["pep_scan_num"]])
  df <- reloc_col_after(df, "raw_file", "pep_scan_num")
  qs::qsave(df, file.path(tempdir, out_name), preset = "fast")
  
  message("Completed: ", mgth)
  
  invisible(df)
}


#' Frame advancing
#' 
#' For all modules of fixed and variable modifications
#' 
#' @param mgf_frames A group of mgf frames (from chunk splitting). A
#'   \code{mgf_frames[[i]]} contains one to multiple MGFs whose MS1 masses are
#'   in the same interval. The \code{mgf_frames} are ordered by increasing
#'   values in \code{frame} for progressive searches.
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams hms2match
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
mframes_adv <- function (mgf_frames = NULL, theopeps = NULL, 
                         aa_masses_all, funs_ms2, ms1vmods_all, ms2vmods_all, 
                         mod_indexes = NULL, 
                         type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                         maxn_sites_per_vmod = 3L, maxn_fnl_per_seq = 3L, 
                         maxn_vnl_per_seq = 3L, 
                         maxn_vmods_sitescombi_per_pep = 64L, 
                         minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                         min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                         digits = 4L) 
{
  lenm <- length(mgf_frames)
  
  if (!lenm)
    return(NULL)
  
  FUNs <- lapply(funs_ms2, as.symbol)
  lenf <- length(funs_ms2)
  out  <- vector("list", lenm) 
  list_null <- list(NULL)
  
  ntmod_all    <- lapply(aa_masses_all, attr, "ntmod", exact = TRUE)
  ctmod_all    <- lapply(aa_masses_all, attr, "ctmod", exact = TRUE)
  ntmass_all   <- lapply(aa_masses_all, find_nterm_mass)
  ctmass_all   <- lapply(aa_masses_all, find_cterm_mass)
  fmods_nl_all <- lapply(aa_masses_all, attr, "fmods_nl", exact = TRUE)
  vmods_nl_all <- lapply(aa_masses_all, attr, "vmods_nl", exact = TRUE)
  amods_all    <- lapply(aa_masses_all, attr, "amods", exact = TRUE)
  
  fmods_nl_all <- lapply(fmods_nl_all, function (x) if (length(x)) x else NULL)
  vmods_nl_all <- lapply(vmods_nl_all, function (x) if (length(x)) x else NULL)
  amods_all    <- lapply(amods_all,    function (x) if (length(x)) x else NULL)
  ms1vmods_all <- lapply(ms1vmods_all, function (x) if (length(x)) x else NULL)
  ms2vmods_all <- lapply(ms2vmods_all, function (x) if (length(x)) x else NULL)

  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame   <- mgfs_cr[["frame"]][1]
  
  thaf_ms2s <- thcr_ms2s <- thbf_ms2s <- 
    thaf_peps <- thcr_peps <- thbf_peps <- 
    thaf_masses <- thcr_masses <- thbf_masses <- 
    thaf <- thcr <- thbf <- vector("list", lenf)
  
  bfi <- 1L
  cri <- bfi + 1L
  
  for (j in 1:lenf) {
    thbf_ms1s_j <- theopeps[[j]][[bfi]]
    
    if (is.null(thbf_ms1s_j))
      thbf_ms2s[j] <- thbf_masses[j] <- thbf_peps[j] <- thbf[j] <- list_null
    else {
      thbf[[j]] <- thbf_ms1s_j
      thbf_peps_j <- thbf_peps[[j]] <- thbf_ms1s_j[["pep_seq"]]
      thbf_masses[[j]] <- thbf_ms1s_j[["mass"]]
      
      thbf_ms2s[[j]] <- mapply(
        FUNs[[j]], 
        aa_seq = thbf_peps_j, 
        ms1_mass = thbf_masses[[j]], 
        MoreArgs = list(
          aa_masses = aa_masses_all[[j]], 
          ms1vmods = ms1vmods_all[[j]], 
          ms2vmods = ms2vmods_all[[j]], 
          ntmod = ntmod_all[[j]], 
          ctmod = ctmod_all[[j]], 
          ntmass = ntmass_all[[j]], 
          ctmass = ctmass_all[[j]], 
          amods = amods_all[[j]], 
          vmods_nl = vmods_nl_all[[j]], 
          fmods_nl = fmods_nl_all[[j]], 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_fnl_per_seq = maxn_fnl_per_seq, 
          maxn_vnl_per_seq = maxn_vnl_per_seq, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      # temporarily share peptide names between targets and decoys; 
      # later is.na(pep_ivmod) -> decoys -> add "-" to prot_acc -> reverse sequence
      names(thbf_ms2s[[j]]) <- thbf_peps_j
    }
    
    thcr_ms1s_j <- theopeps[[j]][[cri]]
    
    if (is.null(thcr_ms1s_j))
      thcr_ms2s[j] <- thcr_masses[j] <- thcr_peps[j] <- thcr[j] <- list_null
    else {
      thcr[[j]] <- thcr_ms1s_j
      thcr_peps_j <- thcr_peps[[j]] <- thcr_ms1s_j[["pep_seq"]]
      thcr_masses[[j]] <- thcr_ms1s_j[["mass"]]
      
      thcr_ms2s[[j]] <- mapply(
        FUNs[[j]], 
        aa_seq = thcr_peps_j, 
        ms1_mass = thcr_masses[[j]], 
        MoreArgs = list(
          aa_masses = aa_masses_all[[j]], 
          ms1vmods = ms1vmods_all[[j]], 
          ms2vmods = ms2vmods_all[[j]], 
          ntmod = ntmod_all[[j]], 
          ctmod = ctmod_all[[j]], 
          ntmass = ntmass_all[[j]], 
          ctmass = ctmass_all[[j]], 
          amods = amods_all[[j]], 
          vmods_nl = vmods_nl_all[[j]], 
          fmods_nl = fmods_nl_all[[j]], 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_fnl_per_seq = maxn_fnl_per_seq, 
          maxn_vnl_per_seq = maxn_vnl_per_seq, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(thcr_ms2s[[j]]) <- thcr_peps_j
    }
  }
  
  ## --- iteration ---
  for (i in seq_len(lenm)) {
    ms1_exptmasses  <- mgfs_cr[["ms1_mass"]]
    ms2_exptmoverzs <- mgfs_cr[["ms2_moverz"]]
    
    afi <- cri + 1L
    
    for (j in 1:lenf) {
      thaf_ms1s_j <- theopeps[[j]][[afi]]
      
      if (is.null(thaf_ms1s_j))
        thaf_ms2s[j] <- thaf_masses[j] <- thaf_peps[j] <- thaf[j] <- list_null
      else {
        thaf[[j]] <- thaf_ms1s_j
        thaf_peps_j <- thaf_peps[[j]] <- thaf_ms1s_j[["pep_seq"]]
        thaf_masses[[j]] <- thaf_ms1s_j[["mass"]]
        
        thaf_ms2s[[j]] <- mapply(
          FUNs[[j]], 
          aa_seq = thaf_peps_j, 
          ms1_mass = thaf_masses[[j]], 
          MoreArgs = list(
            aa_masses = aa_masses_all[[j]], 
            ms1vmods = ms1vmods_all[[j]], 
            ms2vmods = ms2vmods_all[[j]], 
            ntmod = ntmod_all[[j]], 
            ctmod = ctmod_all[[j]], 
            ntmass = ntmass_all[[j]], 
            ctmass = ctmass_all[[j]], 
            amods = amods_all[[j]], 
            vmods_nl = vmods_nl_all[[j]], 
            fmods_nl = fmods_nl_all[[j]], 
            mod_indexes = mod_indexes, 
            type_ms2ions = type_ms2ions, 
            maxn_vmods_per_pep = maxn_vmods_per_pep, 
            maxn_sites_per_vmod = maxn_sites_per_vmod, 
            maxn_fnl_per_seq = maxn_fnl_per_seq, 
            maxn_vnl_per_seq = maxn_vnl_per_seq, 
            maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
            digits = digits
          ), 
          SIMPLIFY = FALSE,
          USE.NAMES = FALSE
        )
        names(thaf_ms2s[[j]]) <- thaf_peps_j
      }
    }
    
    thbca_masses <- mapply(`c`, thbf_masses, thcr_masses, thaf_masses)
    thbca_ms2s   <- mapply(`c`, thbf_ms2s, thcr_ms2s, thaf_ms2s)
    lens_bca <- .Internal(unlist(lapply(thbca_masses, length), 
                                 recursive = FALSE, use.names = FALSE))
    mod_grps <- rep.int(1:lenf, times = lens_bca)
    
    # each `out` for the results of multiple mgfs in one frame
    out[[i]] <- mapply(
      search_mgf, 
      expt_mass_ms1 = ms1_exptmasses, 
      expt_moverz_ms2 = ms2_exptmoverzs, 
      MoreArgs = list(
        theomasses_ms1 = flatten_list(thbca_masses), 
        theomasses_ms2 = flatten_list(thbca_ms2s), 
        pep_mod_groups = mod_grps, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        index_mgf_ms2 = index_mgf_ms2, 
        by_modules = FALSE
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    
    # out[[i]]
    #   MGF_1
    #     theo: "ADCLVPSEIRKLK" (grp_1)
    #     theo: "ASRDLLKEFPQPK" (grp_1)
    #     ...
    #     theo: "TYLNLMGKSKK"   (grp_10)
    #     ...
    #   MGF_2
    #     theo: 
    
    if (i == lenm) 
      break
    
    # advance to the next frame
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr[["frame"]][1]
    
    if (isTRUE(new_frame == (frame + 1L))) {
      cri <- cri + 1L
      
      thbf <- thcr
      thbf_masses <- thcr_masses
      thbf_peps <- thcr_peps
      thbf_ms2s <- thcr_ms2s
      
      thcr <- thaf
      thcr_masses <- thaf_masses
      thcr_peps <- thaf_peps
      thcr_ms2s <- thaf_ms2s
    } 
    else if (isTRUE(new_frame == (frame + 2L))) {
      cri <- cri + 2L
      
      thbf <- thaf
      thbf_masses <- thaf_masses
      thbf_peps <- thaf_peps
      thbf_ms2s <- thaf_ms2s
      
      for (j in 1:lenf) {
        thcr_ms1s_j <- theopeps[[j]][[cri]]
        
        if (is.null(thcr_ms1s_j))
          thcr_ms2s[j] <- thcr_masses[j] <- thcr_peps[j] <- thcr[j] <- list_null
        else {
          thcr[[j]] <- thcr_ms1s_j
          thcr_peps_j <- thcr_peps[[j]] <- thcr_ms1s_j[["pep_seq"]]
          thcr_masses[[j]] <- thcr_ms1s_j[["mass"]]
          
          thcr_ms2s[[j]] <- mapply(
            FUNs[[j]], 
            aa_seq = thcr_peps_j, 
            ms1_mass = thcr_masses[[j]], 
            MoreArgs = list(
              aa_masses = aa_masses_all[[j]], 
              ms1vmods = ms1vmods_all[[j]], 
              ms2vmods = ms2vmods_all[[j]], 
              ntmod = ntmod_all[[j]], 
              ctmod = ctmod_all[[j]], 
              ntmass = ntmass_all[[j]], 
              ctmass = ctmass_all[[j]], 
              amods = amods_all[[j]], 
              vmods_nl = vmods_nl_all[[j]], 
              fmods_nl = fmods_nl_all[[j]], 
              mod_indexes = mod_indexes, 
              type_ms2ions = type_ms2ions, 
              maxn_vmods_per_pep = maxn_vmods_per_pep, 
              maxn_sites_per_vmod = maxn_sites_per_vmod, 
              maxn_fnl_per_seq = maxn_fnl_per_seq, 
              maxn_vnl_per_seq = maxn_vnl_per_seq, 
              maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
              digits = digits
            ), 
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
          )
          names(thcr_ms2s[[j]]) <- thcr_peps_j
        }
      }
    } 
    else {
      cri <- cri + 3L
      bfi <- cri - 1L
      
      for (j in 1:lenf) {
        thbf_ms1s_j <- theopeps[[j]][[bfi]]
        
        if (is.null(thbf_ms1s_j))
          thbf_ms2s[j] <- thbf_masses[j] <- thbf_peps[j] <- thbf[j] <- list_null
        else {
          thbf[[j]] <- thbf_ms1s_j
          thbf_peps_j <- thbf_peps[[j]] <- thbf_ms1s_j[["pep_seq"]]
          thbf_masses[[j]] <- thbf_ms1s_j[["mass"]]
          
          thbf_ms2s[[j]] <- mapply(
            FUNs[[j]], 
            aa_seq = thbf_peps_j, 
            ms1_mass = thbf_masses[[j]], 
            MoreArgs = list(
              aa_masses = aa_masses_all[[j]], 
              ms1vmods = ms1vmods_all[[j]], 
              ms2vmods = ms2vmods_all[[j]], 
              ntmod = ntmod_all[[j]], 
              ctmod = ctmod_all[[j]], 
              ntmass = ntmass_all[[j]], 
              ctmass = ctmass_all[[j]], 
              amods = amods_all[[j]], 
              vmods_nl = vmods_nl_all[[j]], 
              fmods_nl = fmods_nl_all[[j]], 
              mod_indexes = mod_indexes, 
              type_ms2ions = type_ms2ions, 
              maxn_vmods_per_pep = maxn_vmods_per_pep, 
              maxn_sites_per_vmod = maxn_sites_per_vmod, 
              maxn_fnl_per_seq = maxn_fnl_per_seq, 
              maxn_vnl_per_seq = maxn_vnl_per_seq, 
              maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
              digits = digits
            ), 
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
          )
          names(thbf_ms2s[[j]]) <- thbf_peps_j
        }
        
        thcr_ms1s_j <- theopeps[[j]][[cri]]
        
        if (is.null(thcr_ms1s_j))
          thcr_ms2s[j] <- thcr_masses[j] <- thcr_peps[j] <- thcr[j] <- list_null
        else {
          thcr[[j]] <- thcr_ms1s_j
          thcr_peps_j <- thcr_peps[[j]] <- thcr_ms1s_j[["pep_seq"]]
          thcr_masses[[j]] <- thcr_ms1s_j[["mass"]]
          
          thcr_ms2s[[j]] <- mapply(
            FUNs[[j]], 
            aa_seq = thcr_peps_j, 
            ms1_mass = thcr_masses[[j]], 
            MoreArgs = list(
              aa_masses = aa_masses_all[[j]], 
              ms1vmods = ms1vmods_all[[j]], 
              ms2vmods = ms2vmods_all[[j]], 
              ntmod = ntmod_all[[j]], 
              ctmod = ctmod_all[[j]], 
              ntmass = ntmass_all[[j]], 
              ctmass = ctmass_all[[j]], 
              amods = amods_all[[j]], 
              vmods_nl = vmods_nl_all[[j]], 
              fmods_nl = fmods_nl_all[[j]], 
              mod_indexes = mod_indexes, 
              type_ms2ions = type_ms2ions, 
              maxn_vmods_per_pep = maxn_vmods_per_pep, 
              maxn_sites_per_vmod = maxn_sites_per_vmod, 
              maxn_fnl_per_seq = maxn_fnl_per_seq, 
              maxn_vnl_per_seq = maxn_vnl_per_seq, 
              maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
              digits = digits
            ), 
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
          )
          names(thcr_ms2s[[j]]) <- thcr_peps_j
        }
      }
    }
    
    frame <- new_frame
  }
  
  out <- post_frame_adv(out, mgf_frames)
}


#' Fuzzy matches with a +/-1 window.
#' 
#' Not used but called the codes inside directly.
#' 
#' @param x A vector to be matched.
#' @param y A vector to be matched against.
#' @importFrom fastmatch fmatch %fin% 
#' @examples 
#' library(mzion)
#' 
#' ans1 <- mzion:::fuzzy_match_one(c(74953, 74955), rep(74954, 2))
#' ans2 <- mzion:::fuzzy_match_one(c(74953, 74955), 74954)
#' 
#' stopifnot(identical(ans1, ans2))
#' stopifnot(ans1 == c(TRUE, TRUE))
fuzzy_match_one <- function (x, y) 
{
  mi <- x %fin% y
  bf <- (x - 1L) %fin% y
  af <- (x + 1L) %fin% y
  
  mi | bf | af
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
#' library(mzion)
#' 
#' ans1 <- mzion:::fuzzy_match_one2(c(74953, 74955), rep(74954, 2))
#' ans2 <- mzion:::fuzzy_match_one2(c(74953, 74955), 74954)
#'
#' stopifnot(identical(ans1, ans2))
#' stopifnot(ans1 == c(FALSE, TRUE))
#'
#' ans3 <- mzion:::fuzzy_match_one2(c(74953, 74955, 80000), c(74955, 80000))
#' 
#' ## The x3 example from "find_ms2_bypep"
#' x <- c(-9185, -3369, -1973, -626, 59, 714, 3326, 7106, 7711, 7715, 8316, 8320, 
#'        8916, 8920, 9511, 9515, 10102, 10688, 11211, 12945, 16807, 24001, 24481, 
#'        31480, 32350, 32805, 37050, 37875, 42986, 53028, 53377, 53711, 56940, 58542, 
#'        59172, 61310, 62482, 70941, 73801, 77575, 78046, 78047, 84120, 85881, 89313, 
#'        91185, 96328, 101503, 102916, 104302, 113257, 113411, 116563, 118593, 
#'        121336, 121405, 121474, 123450, 123841, 125826, 127823, 130750, 131786, 
#'        131842, 131903, 134568, 135267, 135956, 139090, 139200, 146310, 146801, 
#'        146902, 149442, 152081, 152174, 153544, 153635, 160913, 160995, 161078, 
#'        162794, 162875, 163036, 163117, 163191, 163271, 168686, 169869, 169943, 
#'        173741, 173812, 173951, 174856, 174922, 174990, 175059, 175128, 175197, 
#'        175266)
#' 
#' aas <- unlist(strsplit("SLAAEEEAAR", ""))
#' 
#' y <- c(317.2022, 430.2863, 501.3234, 572.3605, 701.4031, 
#'        830.4457, 959.4883, 1030.5254, 1101.5625, 1257.6636, 
#'        175.1190, 246.1561, 317.1932, 446.2358, 575.2784, 
#'        704.3210, 775.3581, 846.3952, 959.4793, 1046.5113)
#' 
#' names(y) <- c(aas, rev(aas))
#' 
#' ppm_ms2 <- 13L
#' min_ms2mass <- 115L
#' d <- ppm_ms2/1E6
#' y <- ceiling(log(y/min_ms2mass)/log(1+d))
#' 
#' ans <- mzion:::fuzzy_match_one2(x, y)
fuzzy_match_one2 <- function (x, y) 
{
  mi <- x %fin% y
  if (any(mi)) y[y %fin% x[mi]] <- 0L
  
  x2 <- x - 1L
  bf <- x2 %fin% y
  if (any(bf)) y[y %fin% x2[bf]] <- 0L
  
  af <- (x + 1L) %fin% y
  
  mi | bf | af
}


#' Helper: matches between theoretical and experimental MS2 ions.
#'
#' @param expts Numeric vector; one series of experimental MS2s.
#' @param theos Numeric vector; one to multiple series of theoretical MS2s.
#' @param ex Converted expts as integers.
#' @param d ppm_ms2 divided by 1E6.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @importFrom purrr map
#' @importFrom fastmatch fmatch %fin%
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' ## Experimental 322.18704 fit to both b- and y- ions
#' #  (one expt to multiple theos)
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
#' ppm_ms2 <- 25L
#' min_ms2mass <- 115L
#' 
#' d <- ppm_ms2/1E6
#' ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
#' 
#' pep <- "PEPTIDE"
#' nms <- unlist(stringr::str_split(pep, ""))
#'
#' names(theos) <- c(nms, rev(nms))
#' theos <- list(`0000000` = theos)
#'
#' x1 <- mzion:::find_ms2_bypep(theos, expts, ex, d)
#'
#' ## Both expts 74953 and 74955 fit to theos 74954
#' #  (multiple expts to one theo)
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
#' ppm_ms2 <- 25L
#' min_ms2mass <- 115L
#' 
#' d <- ppm_ms2/1E6
#' ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
#' 
#' x2 <- mzion:::find_ms2_bypep(theos, expts, ex, d)
#' 
#' ## Experimental 317.20001 & 317.19315 match to theoreitcal 317.2022 & 317.1932;
#' #  experimental 959.48468 is also multiple dipping
#' #  (multiple to multiple)
#' aas <- unlist(strsplit("SLAAEEEAAR", ""))
#' theos <- c(317.2022, 430.2863, 501.3234, 572.3605, 701.4031, 
#'            830.4457, 959.4883, 1030.5254, 1101.5625, 1257.6636, 
#'            175.1190, 246.1561, 317.1932, 446.2358, 575.2784, 
#'            704.3210, 775.3581, 846.3952, 959.4793, 1046.5113)
#' names(theos) <- c(aas, rev(aas))
#' theos <- list(`0000000000` = theos)
#' 
#' expts <- c(102.05550, 110.07173, 112.08743, 114.06662, 115.08708, 116.07106, 
#'            120.08112, 126.12806, 127.12510, 127.13136, 128.12843, 128.13467, 
#'            129.13181, 129.13805, 130.13509, 130.14134, 131.13843, 132.14159, 
#'            133.04318, 136.07600, 143.08151, 157.10843, 158.09259, 173.14983, 
#'            175.11917, 176.15997, 186.15306, 188.15982, 201.08725, 229.12970, 
#'            230.17041, 231.17366, 241.08228, 246.15643, 248.18086, 255.17380, 
#'            259.09235, 289.20801, 300.16696, 315.25952, 317.19315, 317.20001, 
#'            343.25476, 351.20572, 367.22961, 376.27621, 402.29141, 430.28644,
#'            438.26654, 446.23587, 501.32376, 502.32788, 523.34491, 537.33649, 
#'            556.83990, 557.34125, 557.84283, 572.36127, 575.27863, 590.31671, 
#'            605.84045, 629.33453, 637.86853, 638.33740, 638.83984, 661.35687, 
#'            667.39923, 673.40417, 701.40350, 702.40936, 770.42517, 775.35840, 
#'            776.37720, 802.44659, 830.44556, 831.44958, 846.39349, 847.39465, 
#'            931.49158, 932.48187, 933.48218, 954.54474, 955.54791, 957.55493, 
#'            958.55646, 959.48468, 960.48773, 1030.52563, 1046.50830, 1047.50684, 
#'            1100.53101, 1101.54858, 1103.53040, 1116.59692, 1117.54749, 
#'            1118.54334, 1119.54565, 1120.55347, 1121.55408, 1122.55737)
#' 
#' ppm_ms2 <- 13L
#' min_ms2mass <- 115L
#' d <- ppm_ms2/1E6
#' ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
#' 
#' x3 <- mzion:::find_ms2_bypep(theos, expts, ex, d)
#' 
#' ## 7 b-ions and 9 y-bios 
#' #  (an even total, n_ps matched but identities off)
#' theos <- c(344.2131,504.2438,617.3278,732.3548,845.4389,946.4865,
#'            1003.5080,1102.5764,1258.6775,175.1190,274.1874,331.2088,
#'            432.2565,545.3406,660.3675,773.4516,933.4822,1047.5252)
#' names(theos) <- c("N","C","I","D","I","T","G","V","R",
#'                   "R","V","G","T","I","D","I","C","N")
#' theos <- list(`000000000` = theos)
#' 
#' expts <- c(115.08701,116.07098,120.08125,126.12807,127.12508,127.13135,
#'            128.12840,128.13467,129.10252,129.13177,129.13803,130.09766,
#'            130.13507,130.14136,131.13843,136.07590,142.09769,157.09749,
#'            157.10825,158.09258,173.14984,175.11916,175.15634,176.15977,
#'            186.15321,188.15988,215.13939,230.17043,231.17412,247.19698,
#'            248.18083,254.16139,255.14516,257.16092,272.17194,273.21265,
#'            295.17020,314.18243,316.21832,331.21051,331.21661,344.21338,
#'            345.19727,376.27597,377.27975,389.14813,400.27567,415.22980,
#'            432.25732,445.25986,458.28168,462.25909,473.74744,476.24893,
#'            491.30301,497.73514,501.31610,504.24396,506.24860,520.24512,
#'            545.34082,546.34387,584.76740,588.35773,589.32990,606.37183,
#'            616.83417,617.32935,629.88574,630.33142,638.32935,638.38837,
#'            638.84747,638.89111,639.34772,639.39099,660.36957,672.38916,
#'            715.32959,717.40179,732.35510,733.35760,773.45404,845.43958,
#'            900.49872,901.49738,902.42206,933.48560,946.48560,1003.50885,
#'            1047.52295,1102.53992,1104.54358,1118.53979,1118.65076,
#'            1119.55566,1120.54919,1121.56397,1122.56702,1123.57056)
#' 
#' ppm_ms2 <- 13L
#' min_ms2mass <- 115L
#' d <- ppm_ms2/1E6
#' ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
#' 
#' x4 <- mzion:::find_ms2_bypep(theos, expts, ex, d)
#' 
#' 
#' ## fewer matches with "find_ppm_outer_bycombi" and check minn_ms2 again
#' #  (doesn't really return NULL since now with "ppm_ms2 * 2")
#' pep <- "EFINSLRLYR"
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' theos <- c(359.2128,506.2812,619.3653,734.3922,821.4243,934.5083,1090.6094,
#'            1203.6935,1366.7568,1522.8579,175.1190,338.1823,451.2663,607.3675,
#'            720.4515,807.4835,922.5105,1035.5946,1182.6630,1311.7056)
#' 
#' names(theos) <- c(nms, rev(nms))
#' theos <- list(`0005000000` = theos)
#' 
#' expts <- c(126.12811,127.12556,127.13139,128.12862,128.13455,129.13194,
#'            129.13786,130.13542,130.14139,131.13852,136.07597,173.14980,
#'            175.11916,175.15663,176.15979,186.15309,188.15987,215.13940,
#'            227.10303,230.17044,231.17381,247.19708,248.18105,249.18443,
#'            273.21262,316.21869,344.21350,345.19760,345.21689,353.19775,
#'            358.22925,364.14948,376.27615,377.27942,397.20880,479.28204,
#'            480.28650,507.27704,507.31665,508.28040,508.31964,550.31934,
#'            578.31403,579.31818,601.33063,602.33466,620.40070,621.40363,
#'            679.36194,680.40875,680.90851,690.32965,707.35663,708.35950,
#'            721.44830,722.45142,736.94965,761.93817,762.43274,762.93304,
#'            770.79028,770.94440,771.44598,786.48285,786.98523,792.44543,
#'            820.44061,821.44318,834.53247,835.53638,877.45514,893.49561,
#'            903.47803,904.48071,921.48846,922.49133,945.56488,963.57489,
#'            964.57990,1006.57715,1016.55975,1017.56122,1034.60925,1035.57776,
#'            1035.61169,1165.61340,1166.61731,1197.67346,1198.67700,1311.71948,
#'            1312.72559,1366.73694,1368.74927,1369.74463,1382.75659,1383.75610,
#'            1384.75915,1385.76355,1386.76611,1387.76599)
#' 
#' ppm_ms2 <- 13L
#' min_ms2mass <- 115L
#' d <- ppm_ms2/1E6
#' ex <- ceiling(log(expts/min_ms2mass)/log(1+d))
#' 
#' x5 <- mzion:::find_ms2_bypep(theos, expts, ex, d, ppm_ms2)
#' 
#' }
#' 
#' @return Lists of (1) theo, (2) expt, (3) ith, (4) iex and (5) m.
find_ms2_bypep <- function (theos = NULL, expts = NULL, ex = NULL, d = NULL, 
                            ppm_ms2 = 10L, min_ms2mass = 115L, minn_ms2 = 6L, 
                            index_mgf_ms2 = FALSE) 
{
  ##############################################################################
  # `theos`
  #   the same pep_seq at different applicable ivmods and NLs
  # 
  # length(theos) may be greater than one with site permutation and/or NLs.
  # 
  # `theos` may be empty: 
  #   e.g. the matched one is after `maxn_vmods_sitescombi_per_pep` 
  #   and never get matched.
  # 
  # ex: `expts` in integers
  # th_i: the i-th `theos` in integers
  # 
  # ex has no duplicated entries; 
  # th_i can.
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
    return(list(theo = NULL, expt = NULL, ith = NULL, iex = NULL, m = NULL))
  
  # ---
  out <- vector("list", len)
  
  ## forward matches
  if (len > 3L) {
    mths  <- index_mz(.Internal(unlist(theos, recursive = FALSE, use.names = FALSE)), 
                      min_ms2mass, d)
    pss  <- mths %fin% ex | (mths - 1L) %fin% ex | (mths + 1L) %fin% ex
    pss  <- fold_vec(pss, len)
    mths <- fold_vec(mths, len)
    ipss <- lapply(pss, function (x) .Internal(which(x)))
    
    for (i in 1:len) {
      theos_i <- theos[[i]]
      th_i <- mths[[i]]
      ps <- pss[[i]]
      ips <- ipss[[i]]
      
      ### the remaining are the same as those under "else" ###
      
      ## backward matches
      #  expts are in ascending orders, but theos in b1, b2, ... , y1, y2, ...
      #  separate matches to theos_b and theos_y, each are in ascending order
      n_ps <- length(ips)
      
      if(n_ps >= minn_ms2) {
        # separated b and y matches (to handled double-dipping between b and y)
        # (adj: bps <- fuzzy_match_one2(ex, th_i[1:mid]))
        lth <- length(ps)
        mid <- lth/2L
        
        # experimental es initially filled by NA and matched theoretical values, 
        # and at the end replaced with matched experimental values
        es <- theos_i
        es[!ps] <- NA_real_
        
        ex_bf <- ex - 1L
        ex_af <- ex + 1L
        
        # b-ions
        y_1 <- th_i[1:mid]
        ps_1 <- ex %fin% y_1 | ex_bf %fin% y_1 | ex_af %fin% y_1
        ips_1 <- .Internal(which(ps_1))
        
        # y-ions
        y_2 <- th_i[(mid+1L):lth]
        ps_2 <- ex %fin% y_2 | ex_bf %fin% y_2 | ex_af %fin% y_2
        ips_2 <- .Internal(which(ps_2))
        
        # b- and y-ions
        expt_1 <- expts[ips_1]
        expt_2 <- expts[ips_2]
        expt_12 <- c(expt_1, expt_2)
        ips_12 <- c(ips_1, ips_2)
        len_12 <- length(expt_12)
        
        # (occur rarely; OK to recalculate freshly `expt_12`)
        if (n_ps != len_12) {
          # "* 2" for three-frame searches
          # also ensure that "ith = ips" in ascending order, not "iex = ips_12"
          out_i <- find_ppm_outer_bycombi(expts, theos_i, ppm_ms2 * 2L)
          
          if (sum(!is.na(out_i[["expt"]])) < minn_ms2) {
            out[[i]] <- list(theo = NULL, expt = NULL, ith = NULL, iex = NULL, m = NULL)
            next
          }
          
          out[[i]] <- out_i
          next
        }
        
        es[ps] <- expt_12
        
        out[[i]] <- list(theo = theos_i, expt = es, ith = ips, iex = ips_12, m = len_12)
      } 
      else
        out[[i]] <- list(theo = NULL, expt = NULL, ith = NULL, iex = NULL, m = NULL)
    }
  }
  else {
    # mths <- lapply(theos, index_mz, min_ms2mass, d)
    # pss <- lapply(mths, function (x) x %fin% ex | (x - 1L) %fin% ex | (x + 1L) %fin% ex)
    # ipss <- lapply(pss, function (x) .Internal(which(x)))
    
    for (i in 1:len) {
      ## forward matches
      theos_i <- theos[[i]]
      th_i <- index_mz(theos_i, min_ms2mass, d)
      ps <- th_i %fin% ex | (th_i - 1L) %fin% ex | (th_i + 1L) %fin% ex
      ips <- .Internal(which(ps))
      
      ## "ith = ips" in ascending order, not "iex = ips_12"
      
      ### the remaining are the same under "if" ###
      
      ## backward matches
      #  expts are in ascending orders, but theos in b1, b2, ... , y1, y2, ...
      #  separate matches to theos_b and theos_y, each are in ascending order
      n_ps <- length(ips)
      
      if(n_ps >= minn_ms2) {
        # separated b and y matches (to handled double-dipping between b and y)
        # (adj: bps <- fuzzy_match_one2(ex, th_i[1:mid]))
        lth <- length(ps)
        mid <- lth/2L
        
        # experimental es initially filled by NA and matched theoretical values, 
        # and at the end replaced with matched experimental values
        es <- theos_i
        es[!ps] <- NA_real_
        
        ex_bf <- ex - 1L
        ex_af <- ex + 1L
        
        # b-ions
        y_1 <- th_i[1:mid]
        ps_1 <- ex %fin% y_1 | ex_bf %fin% y_1 | ex_af %fin% y_1
        ips_1 <- .Internal(which(ps_1))
        
        # y-ions
        y_2 <- th_i[(mid+1L):lth]
        ps_2 <- ex %fin% y_2 | ex_bf %fin% y_2 | ex_af %fin% y_2
        ips_2 <- .Internal(which(ps_2))
        
        # b- and y-ions
        expt_1 <- expts[ips_1]
        expt_2 <- expts[ips_2]
        expt_12 <- c(expt_1, expt_2)
        ips_12 <- c(ips_1, ips_2)
        len_12 <- length(expt_12)
        
        # (occur rarely; OK to recalculate freshly `expt_12`)
        if (n_ps != len_12) {
          # "* 2" for three-frame searches
          # also ensure that "ith = ips" in ascending order, not "iex = ips_12"
          out_i <- find_ppm_outer_bycombi(expts, theos_i, ppm_ms2 * 2L)
          
          if (sum(!is.na(out_i[["expt"]])) < minn_ms2) {
            out[[i]] <- list(theo = NULL, expt = NULL, ith = NULL, iex = NULL, m = NULL)
            next
          }
          
          out[[i]] <- out_i
          next
        }
        
        es[ps] <- expt_12
        
        out[[i]] <- list(theo = theos_i, expt = es, ith = ips, iex = ips_12, m = len_12)
      } 
      else
        out[[i]] <- list(theo = NULL, expt = NULL, ith = NULL, iex = NULL, m = NULL)
    }
  }
  
  names(out) <- names(theos)
  
  out
}


#' Matches an MGF query
#'
#' @param expt_mass_ms1 Numeric; the experimental MS1 mass
#' @param expt_moverz_ms2 A numeric list; the experimental MS2 m/z's
#' @param theomasses_ms1 Numeric vector; the theoretical MS1 masses at the
#'   preceding \code{-1}, the current and the following \code{+1} frames
#' @param theomasses_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   preceding \code{-1}, the current and the following \code{+1} frames
#' @param pep_mod_groups The index(es) of peptide modification groups; single
#'   value at \code{by_modules = TRUE}
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @inheritParams pair_mgftheo
#' @examples
#' \donttest{
#' library(mzion)
#' library(fastmatch)
#'
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
search_mgf <- function (expt_mass_ms1 = NULL, expt_moverz_ms2 = NULL, 
                        theomasses_ms1 = NULL, theomasses_ms2 = NULL, 
                        pep_mod_groups = NULL, 
                        minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                        min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                        by_modules = FALSE) 
{
  # --- find MS2 matches ---
  d2 <- ppm_ms2/1E6
  
  ex <- if (index_mgf_ms2) # already indexed
    expt_moverz_ms2
  else
    index_mz(expt_moverz_ms2, min_ms2mass, d2)
  
  # lapply by the same pep_seq at different ivmods and/or NLs
  ans <- if (length(theomasses_ms2)) 
    lapply(theomasses_ms2, find_ms2_bypep, 
           expts = expt_moverz_ms2, 
           ex = ex, 
           d = d2, 
           ppm_ms2 = ppm_ms2, 
           min_ms2mass = min_ms2mass, 
           minn_ms2 = minn_ms2, 
           index_mgf_ms2 = index_mgf_ms2)
  else 
    theomasses_ms2
  
  ## Not faster
  # if (is.null(.Internal(unlist(ans, recursive = TRUE, use.names = FALSE)))) return(list())
  
  ## cleans up
  # (1) within a list: removes vmods+ positions that are NULL (< minn_ms2)
  # (no effects on vmods-; need `type` info if to limit to vmods+)
  oks <- lapply(ans, function (this) {
    oks <- lapply(this, function (x) !is.null(x$theo))
    .Internal(unlist(oks, recursive = FALSE, use.names = FALSE))
  })
  
  # USE.NAMES = TRUE 
  # (lapply loses names by `[[` whereas map2 reserves names when available)
  ans <- mapply(function (x, y) x[y], ans, oks, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  # (2)  removes empty lists
  oks2 <- lapply(ans, function(x) length(x) > 0L)
  oks2 <- .Internal(unlist(oks2, recursive = FALSE, use.names = FALSE))
  ans <- ans[oks2]
  theomasses_ms1 <- theomasses_ms1[oks2]
  
  if (!by_modules)
    pep_mod_groups <- pep_mod_groups[oks2]

  ans <- mapply(
    function (x, y, g) {
      attr(x, "theo_ms1") <- y
      attr(x, "pep_mod_group") <- g
      x
    }, 
    ans, theomasses_ms1, pep_mod_groups, 
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


#' Matches experimentals and theoreticals
#'
#' For a single module
#'
#' @param pep_mod_group The index of peptide modification groups
#' @param aa_masses An amino-acid look-up
#' @param FUN A function, e.g., \link{gen_ms2ions_base}, with an i-th module of
#'   \code{aa_masses}
#' @param ms1vmods All possible labels of MS1 variable modifications with
#'   an i-th \code{aa_masses}
#' @param ms2vmods All possible labels of MS2 variable modifications with
#'   an i-t \code{aa_masses}
#' @param cl The value of clusters for parallel processes
#' @param df0 An output template
#' @inheritParams hms2match
ms2match_one <- function (pep_mod_group, aa_masses, FUN, 
                          ms1vmods, ms2vmods, cl, 
                          mod_indexes, mgf_path, out_path, type_ms2ions = "by", 
                          maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L, 
                          maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                          maxn_vmods_sitescombi_per_pep = 64L, 
                          minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                          min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                          df0 = NULL, digits = 4L) 
{
  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
  
  message("Matching against: ", 
          if (nchar(nm_vmods) == 0L) nm_fmods else paste0(nm_fmods, " | ", nm_vmods))
  
  mgth       <- paste0("expttheo_", pep_mod_group, ".rds")
  out_name   <- gsub("^expttheo", "ion_matches", mgth)
  mgftheo    <- qs::qread(file.path(mgf_path, mgth))
  mgf_frames <- mgftheo[["mgf_frames"]]
  theopeps   <- mgftheo[["theopeps"]]
  rm("mgftheo")
  
  if (!length(mgf_frames)) {
    qs::qsave(df0, file.path(out_path, "temp", out_name))
    return(df0)
  }
  
  ntmod    <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod    <- attr(aa_masses, "ctmod", exact = TRUE)
  ntmass   <- find_nterm_mass(aa_masses)
  ctmass   <- find_cterm_mass(aa_masses)
  fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
  vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
  amods    <- attr(aa_masses, "amods", exact = TRUE)
  
  fmods_nl <- if (length(fmods_nl)) fmods_nl else NULL
  vmods_nl <- if (length(vmods_nl)) vmods_nl else NULL
  amods    <- if (length(amods)) amods else NULL
  ms1vmods <- if (length(ms1vmods)) ms1vmods else NULL
  ms2vmods <- if (length(ms2vmods)) ms2vmods else NULL
  
  FUN  <- as.symbol(FUN)
  
  df <- parallel::clusterMap(
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
                    vmods_nl = vmods_nl, 
                    fmods_nl = fmods_nl, 
                    pep_mod_group = pep_mod_group, 
                    mod_indexes = mod_indexes, 
                    type_ms2ions = type_ms2ions, 
                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                    maxn_fnl_per_seq = maxn_fnl_per_seq, 
                    maxn_vnl_per_seq = maxn_vnl_per_seq, 
                    maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
                    minn_ms2 = minn_ms2, 
                    ppm_ms1 = ppm_ms1, 
                    ppm_ms2 = ppm_ms2, 
                    min_ms2mass = min_ms2mass, 
                    index_mgf_ms2 = index_mgf_ms2, 
                    digits = digits, 
                    FUN = FUN), 
    .scheduling = "dynamic")

  df <- dplyr::bind_rows(df)
  
  out_nm <- file.path(out_path, "temp", paste0("ion_matches_", pep_mod_group, ".rds"))
  
  if (is.null(df)) {
    qs::qsave(df, out_nm, preset = "fast")
    return(NULL)
  }
  
  # fields not yet available with `ms2match_all`
  df[["pep_fmod"]] <- nm_fmods
  df[["pep_vmod"]] <- nm_vmods
  df[["pep_mod_group"]] <- pep_mod_group
  
  df <- dplyr::rename(df, 
                      pep_ret_range = ret_time, 
                      pep_scan_title = scan_title,
                      pep_exp_mz = ms1_moverz, 
                      pep_n_ms2 = ms2_n, 
                      pep_exp_mr = ms1_mass, 
                      pep_tot_int = ms1_int, 
                      pep_scan_num = scan_num, 
                      pep_exp_z = ms1_charge, 
                      pep_ms2_moverzs = ms2_moverz, 
                      pep_ms2_ints = ms2_int, 
                      pep_frame = frame)
  df[["pep_scan_num"]] <- as.character(df[["pep_scan_num"]])
  
  df <- reloc_col_after(df, "raw_file", "scan_num")
  df <- reloc_col_after(df, "pep_mod_group", "raw_file")
  
  qs::qsave(df, out_nm, preset = "fast")
}


#' Frames advancement.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#'
#' @param mgf_frames A group of mgf frames (from chunk splitting). A
#'   \code{mgf_frames[[i]]} contains one to multiple MGFs whose MS1 masses are
#'   in the same interval. The \code{mgf_frames} are ordered by increasing
#'   values in \code{frame} for progressive searches.
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses}.
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses}.
#' @param ntmass The mass of a fixed or variable N-term modification.
#' @param ctmass The mass of a fixed or variable C-term modification.
#' @param amods \code{Anywhere} variable modifications.
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @param vmods_nl The attribute of \code{vmods_nl} from an \code{aa_masses}.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @param FUN A function pointer to, e.g., \link{gen_ms2ions_base}.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_one
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
frames_adv <- function (mgf_frames = NULL, theopeps = NULL, 
                        aa_masses = NULL, ms1vmods = NULL, ms2vmods = NULL, 
                        ntmod = NULL, ctmod = NULL, 
                        ntmass = NULL, ctmass = NULL, 
                        amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                        pep_mod_group = NULL, mod_indexes = NULL, 
                        type_ms2ions = "by", 
                        maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L, 
                        maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                        maxn_vmods_sitescombi_per_pep = 64L, 
                        minn_ms2 = 6L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                        min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                        digits = 4L, FUN) 
{
  len <- length(mgf_frames)
  
  if (!len)
    return(NULL)
  
  out <- vector("list", len) 
  
  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr[["frame"]][1]
  
  bfi <- 1L
  thbf <- theopeps[[bfi]] 
  thbf_peps <- thbf[["pep_seq"]]
  thbf_masses <- thbf[["mass"]]
  
  cri <- bfi + 1L
  thcr <- theopeps[[cri]]
  thcr_peps <- thcr[["pep_seq"]]
  thcr_masses <- thcr[["mass"]]
  
  # generate both target and decoy MS2
  thbf_ms2s <- mapply(
    FUN, 
    aa_seq = thbf_peps, 
    ms1_mass = thbf_masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_fnl_per_seq = maxn_fnl_per_seq, 
      maxn_vnl_per_seq = maxn_vnl_per_seq, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  # temporarily share peptide names between targets and decoys; 
  # later is.na(pep_ivmod) -> decoys -> add "-" to prot_acc -> reverse sequence
  names(thbf_ms2s) <- thbf_peps
  
  thcr_ms2s <- mapply(
    FUN, 
    aa_seq = thcr_peps, 
    ms1_mass = thcr_masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_fnl_per_seq = maxn_fnl_per_seq, 
      maxn_vnl_per_seq = maxn_vnl_per_seq, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  names(thcr_ms2s) <- thcr_peps
  
  ## --- iteration ---
  for (i in seq_len(len)) {
    exptmasses_ms1  <- mgfs_cr$ms1_mass
    exptmoverzs_ms2 <- mgfs_cr$ms2_moverz
    
    ### Slower to subset + passed as argument 
    #   compared to direct calculation at ~ 4us
    # 
    # exptimoverzs_ms2 <- mgfs_cr$ms2_imoverzs
    ###
    
    afi <- cri + 1L
    
    thaf <- theopeps[[afi]]
    thaf_peps <- thaf[["pep_seq"]]
    thaf_masses <- thaf[["mass"]]
    
    thaf_ms2s <- mapply(
      FUN, 
      aa_seq = thaf_peps, 
      ms1_mass = thaf_masses, 
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
        maxn_fnl_per_seq = maxn_fnl_per_seq, 
        maxn_vnl_per_seq = maxn_vnl_per_seq, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    names(thaf_ms2s) <- thaf_peps
    
    # each `out` for the results of multiple mgfs in one frame
    out[[i]] <- mapply(
      search_mgf, 
      expt_mass_ms1 = exptmasses_ms1, 
      expt_moverz_ms2 = exptmoverzs_ms2, 
      MoreArgs = list(
        pep_mod_groups = pep_mod_group, 
        theomasses_ms1 = c(thbf_masses, thcr_masses, thaf_masses), 
        theomasses_ms2 = c(thbf_ms2s, thcr_ms2s, thaf_ms2s), 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        index_mgf_ms2 = index_mgf_ms2, 
        by_modules = TRUE
      ), 
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    
    if (i == len) break
    
    # advance to the next frame
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr[["frame"]][1]
    
    if (isTRUE(new_frame == (frame + 1L))) {
      cri <- cri + 1L
      
      thbf <- thcr
      thbf_masses <- thcr_masses
      thbf_ms2s <- thcr_ms2s
      
      thcr <- thaf
      thcr_masses <- thaf_masses
      thcr_ms2s <- thaf_ms2s
    } 
    else if (isTRUE(new_frame == (frame + 2L))) {
      cri <- cri + 2L
      
      thbf <- thaf
      thbf_masses <- thaf_masses
      thbf_ms2s <- thaf_ms2s
      
      thcr <- theopeps[[cri]]
      thcr_peps <- thcr[["pep_seq"]]
      thcr_masses <- thcr[["mass"]]
      
      thcr_ms2s <- mapply(
        FUN, 
        aa_seq = thcr_peps, 
        ms1_mass = thcr_masses, 
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
          maxn_fnl_per_seq = maxn_fnl_per_seq, 
          maxn_vnl_per_seq = maxn_vnl_per_seq, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(thcr_ms2s) <- thcr_peps
    } 
    else {
      cri <- cri + 3L
      bfi <- cri - 1L
      
      thbf <- theopeps[[bfi]]
      thbf_peps <- thbf[["pep_seq"]]
      thbf_masses <- thbf[["mass"]]
      
      thcr <- theopeps[[cri]]
      thcr_peps <- thcr[["pep_seq"]]
      thcr_masses <- thcr[["mass"]]
      
      thbf_ms2s <- mapply(
        FUN, 
        aa_seq = thbf_peps, 
        ms1_mass = thbf_masses, 
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
          maxn_fnl_per_seq = maxn_fnl_per_seq, 
          maxn_vnl_per_seq = maxn_vnl_per_seq, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(thbf_ms2s) <- thbf_peps
      
      thcr_ms2s <- mapply(
        FUN, 
        aa_seq = thcr_peps, 
        ms1_mass = thcr_masses, 
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
          maxn_fnl_per_seq = maxn_fnl_per_seq, 
          maxn_vnl_per_seq = maxn_vnl_per_seq, 
          maxn_vmods_sitescombi_per_pep = 
            maxn_vmods_sitescombi_per_pep, 
          digits = digits
        ), 
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      names(thcr_ms2s) <- thcr_peps
    }
    
    frame <- new_frame
  }
  
  out <- post_frame_adv(out, mgf_frames)
}


