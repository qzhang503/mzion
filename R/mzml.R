#' Helper in preparing peaklists
#' 
#' @param filelist A list of peaklist files.
#' @param data_type A data type of either RAW or mzML.
#' @param ppm_ms1_bin The three-bin MS1 ppm.
#' @param ppm_ms2_bin The three-bin MS2 ppm.
#' @param tol_ms1 The (original) MS1 tolerance (\code{20 * 1E-6}).
#' @inheritParams load_mgfs
readmzML <- function (filelist = NULL, out_path = NULL, mgf_path = NULL, 
                      data_type = "mzml", topn_ms2ions = 150L, 
                      topn_dia_ms2ions = 2400L, delayed_diams2_tracing = FALSE, 
                      maxn_dia_precurs = 1000L, 
                      n_dia_ms2bins = 1L, n_dia_scans = 6L, 
                      min_mass = 200L, max_mass = 4500L, 
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      min_ms1_charge = 2L, max_ms1_charge = 4L, 
                      min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                      min_ret_time = 0, max_ret_time = Inf, 
                      tol_ms1 = 2E-5, ppm_ms1_bin = 10L, ppm_ms2_bin = 10L, 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, 
                      is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE, 
                      mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                      enzyme = "trypsin_p", 
                      deisotope_ms2 = TRUE, grad_isotope = 1.6, fct_iso2 = 3.0,
                      max_ms2_charge = 3L, use_defpeaks = FALSE, 
                      maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                      ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                      quant = "none", 
                      use_lfq_intensity = TRUE, ppm_ms1trace = 5L, 
                      bypass_rawexe = FALSE, digits = 4L)
{
  # - hloadMZML (helper)
  #   - loadMZML
  #     - extrDIA
  #     - extrDDA
  # 
  # [Y] is_dia?
  # - hdeisoDIA (helper)
  #   - deisoDIA
  # - hsubDIAMS1 (helper)
  #   - subDIAMS1
  # - htraceDIA (helper)
  #   - traceDIA
  #     -traceLCMS
  # 
  # [Y] is_dda?
  # - hdeisoDDA (helper)
  #   - deisoDDA
  #     - getMS1xyz
  #     - getMS2xyz
  
  
  # traceLCMS
  #   - collapse_xyz
  #     - mapcoll_xyz
  #     - find_gates
  #       - find_gate_edges
  #   - find_lc_gates
  #     - fill_lc_gaps
  #     - find_gate_edges
  # 
  # - find_mdda_mms1s
  #   - collapse_mms1ints
  #     - find_gates
  #       - find_gate_edges
  #   - find_ms1stat
  # 
  
  options(warn = 1L)
  
  # mzML may contains NULL entries and need additional handling
  if (!dir.exists(path_ms1 <- file.path(out_path, "ms1data"))) {
    path_ms1 <- create_dir(path_ms1) # for MBR in proteoQ
  }
  qs::qsave(data_type, file.path(mgf_path, "data_type.rds"), preset = "fast")
  qs::qsave(data_type, file.path(path_ms1, "data_type.rds"), preset = "fast") # for proteoQ
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))

  if (data_type == "raw") {
    message("Processing RAW files.")
    peakfiles <- readRAW(mgf_path = mgf_path, filelist = filelist)
    gc()
  }
  else if (data_type == "pasef") {
    message("Processing RAW PASEF files.")
    peakfiles <- readPASEF(mgf_path = mgf_path, filelist = filelist, 
                           topn_ms2ions = topn_ms2ions, 
                           bypass_rawexe = bypass_rawexe)
    gc()
  }
  else if (data_type == "mzml") {
    message("Processing mzML files.")
    peakfiles <- hloadMZML(filelist, mgf_path, temp_dir)
    gc() # free up xml pointers
  }
  
  if (TRUE) {
    peakfile1 <- peakfiles[[1]]
    is_dia    <- attr(peakfile1, "is_dia", exact = TRUE)
    mzml_type <- attr(peakfile1, "mzml_type", exact = TRUE)
    peakfiles <- unlist(peakfiles)
  }
  else {
    if (data_type == "pasef") {
      peakfiles <- list.files(file.path(mgf_path, "temp_dir"), pattern = "\\.d\\.rds$")
    }
    else if (data_type == "raw") {
      peakfiles <- list.files(file.path(mgf_path, "temp_dir"), pattern = "\\.raw\\.rds$")
    }

    peakfiles <- peakfiles[!grepl("^predeisoDDA_", peakfiles)]
    peakfiles <- peakfiles[!grepl("^pasefms1_", peakfiles)]
    
    is_dia    <- FALSE
    mzml_type <- "raw"
  }
  
  # Coercion and handling after knowing `mzml_type`
  # (mzml_type == "raw" for both Thermo and timsTOF)
  if (mzml_type == "raw" && !maxn_mdda_precurs) {
    warning("Require a minimum `maxn_mdda_precurs = 1`.")
    maxn_mdda_precurs <- 1L
  }
  
  # maxn_mdda_precurs == 0 only to test mzML default deisotoping
  if (grepl("^pwiz", mzml_type)) {
    if (!(maxn_mdda_precurs || use_defpeaks)) {
      warning("Coerce to `use_defpeaks = TRUE`.")
      use_defpeaks <- TRUE
    }
  }
  
  # temporary for bypassing Mzion deisotoping against PASEF???
  data_format <- local({
    file <- file.path(mgf_path, "info_format.rds")
    if (file.exists(file)) qs::qread(file)$data_format else NULL
  })
  is_pasef <- isTRUE(data_format == "Bruker-RAW")
  
  lenf <- length(peakfiles)
  rams <- find_free_mem() / 1024
  n_pcs <- detect_cores(64L) - 1L
  n_cores <- max(min(n_pcs, ceiling(rams / if (is_pasef) 30 else 5), lenf), 1L)
  n_para  <- max(min(n_pcs, round(n_pcs / n_cores)), 1L)
  
  if (is_pasef) {
    if (n_mdda_flanks) {
      n_mdda_flanks <- 0L
      warning("Coerce to `n_mdda_flanks = 0` for PASEF data.")
    }
    
    file_indexes <- seq_along(peakfiles)
    multi_pasef  <- lenf > n_cores
    # n_para <- 1L
    
    if (multi_pasef) {
      pasef_n_secs <- ceiling(lenf / n_cores) # lenf2 <- ceiling(lenf / n_cores) * n_cores
      pasef_cuts   <- ceiling(lenf / pasef_n_secs) * seq_len(pasef_n_secs - 1L)
      pasef_brs    <- findInterval(seq_len(lenf), pasef_cuts, left.open = TRUE)
      
      pasef_fis <- split(peakfiles, pasef_brs)
      pasef_ids <- split(seq_len(lenf), pasef_brs)
      pasef_n_cores <- lengths(pasef_fis, use.names = FALSE)
      pasef_n_paras <-  floor(min(n_pcs, 25L) / pasef_n_cores)
      raws <- vector("list", lenf)
    }
    # else {
    #   pasef_fis  <- peakfiles
    #   pasef_ids <- seq_len(lenf)
    #   pasef_n_cores <- n_cores
    #   pasef_n_paras <-  floor(min(n_pcs, 25L) / pasef_n_cores)
    # }
  }
  else { # Thermo's
    file_indexes <- seq_along(peakfiles)
    multi_pasef  <- FALSE
  }
  
  if (isTRUE(is_dia)) {
    message("Deisotoping DIA-MS.")
    
    # `else` is only slightly faster than `if` 
    if (n_cores <= 1L) {
      
      # need to average flanking spectra prior to deisotoping...
      # need to include `ppm_ms1trace`
      
      if (TRUE) {
        dia_files <- lapply(
          peakfiles, hdeisoDIA,
          temp_dir = temp_dir, 
          min_mass = min_mass, max_mass = max_mass,
          min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
          maxn_dia_precurs = maxn_dia_precurs, 
          topn_dia_ms2ions = topn_dia_ms2ions, 
          min_ms1_charge = min_ms1_charge, 
          max_ms1_charge = max_ms1_charge, 
          min_ret_time = min_ret_time, 
          max_ret_time = max_ret_time, 
          ppm_ms1_deisotope = ppm_ms1_deisotope, 
          ppm_ms2_deisotope = ppm_ms2_deisotope, 
          deisotope_ms2 = deisotope_ms2, 
          max_ms2_charge = max_ms2_charge, 
          grad_isotope = grad_isotope, 
          fct_iso2 = fct_iso2, 
          quant = quant, 
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          n_para = n_para)
      }
      else {
        dia_files <- "dia_23aug2017_hela_serum_timecourse_wide_1a.raw.rds"
      }
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      dia_files <- parallel::clusterApply(
        cl, peakfiles, hdeisoDIA,
        temp_dir = temp_dir, 
        min_mass = min_mass, max_mass = max_mass,
        min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
        maxn_dia_precurs = maxn_dia_precurs, 
        topn_dia_ms2ions = topn_dia_ms2ions, 
        min_ms1_charge = min_ms1_charge, 
        max_ms1_charge = max_ms1_charge, 
        min_ret_time = min_ret_time, 
        max_ret_time = max_ret_time, 
        ppm_ms1_deisotope = ppm_ms1_deisotope, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        deisotope_ms2 = deisotope_ms2, 
        max_ms2_charge = max_ms2_charge, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, 
        quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region, 
        n_para = n_para)
      parallel::stopCluster(cl)
    }
    
    message("Subset DIA-MS1 by isolation windows.")
    n_para2 <- min(n_para, 8L) # 8L: not rate-limiting
    
    if (n_cores <= 1L) {
      if (TRUE) {
        subfiles <- lapply(dia_files, hsubDIAMS1, temp_dir = temp_dir, 
                           n_para = n_para2)
      }
      else {
        subfiles <- qs::qread(file.path(temp_dir, "subfiles.rds"))
      }
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      subfiles <- parallel::clusterApply(
        cl, dia_files, hsubDIAMS1, temp_dir = temp_dir, n_para = n_para2)
      parallel::stopCluster(cl)
    }
    
    message("Tracing DIA-MS.")
    n_para3 <- max(floor(min(n_para, 26L/lenf)), 1L)
    if (lenf > n_para3) n_para3 <- find_min_ncores(lenf, n_para3)
    
    if (n_cores <= 1L) {
      raws <- mapply(
        htraceDIA, 
        subfiles, seq_along(subfiles), 
        MoreArgs = list(
          temp_dir = temp_dir, 
          mgf_path = mgf_path, 
          min_ret_time = min_ret_time, 
          max_ret_time = max_ret_time, 
          min_mass = min_mass, 
          max_mass = max_mass, 
          min_ms2mass = min_ms2mass, 
          max_ms2mass = max_ms2mass, 
          ppm_ms1 = ppm_ms1_deisotope, 
          ppm_ms2 = ppm_ms2_deisotope, 
          n_dia_ms2bins = n_dia_ms2bins, 
          n_dia_scans = n_dia_scans, 
          topn_dia_ms2ions = topn_dia_ms2ions, 
          delayed_diams2_tracing = delayed_diams2_tracing, 
          mgf_cutmzs = mgf_cutmzs, 
          mgf_cutpercs = mgf_cutpercs, 
          quant = quant, 
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          n_para = n_para3), 
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      raws <- parallel::clusterMap(
        cl, htraceDIA, 
        subfiles, seq_along(subfiles), 
        MoreArgs = list(
          temp_dir = temp_dir, 
          mgf_path = mgf_path, 
          min_ret_time = min_ret_time, 
          max_ret_time = max_ret_time, 
          min_mass = min_mass, 
          max_mass = max_mass, 
          min_ms2mass = min_ms2mass, 
          max_ms2mass = max_ms2mass, 
          ppm_ms1 = ppm_ms1_deisotope, 
          ppm_ms2 = ppm_ms2_deisotope, 
          n_dia_ms2bins = n_dia_ms2bins, 
          n_dia_scans = n_dia_scans,
          topn_dia_ms2ions = topn_dia_ms2ions, 
          delayed_diams2_tracing = delayed_diams2_tracing, 
          mgf_cutmzs = mgf_cutmzs, 
          mgf_cutpercs = mgf_cutpercs, 
          quant = quant, 
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          n_para = n_para3), 
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
      parallel::stopCluster(cl)
    }
  }
  else {
    message("Deisotoping DDA-MS at: ", Sys.time())
    ## revisit this, maybe no longer an issue...
    # yco cannot be greater than the co in PASEF extraction
    # too large of yco may cause disparity btw. X and Y values, esp. w/ PASEF
    # since values < yco replaced with NA and can become all NA
    # the current yco = 10 for PASEF actually has no effect... 
    yco <- if (is_pasef) { 10 } else { 100 }
    
    if (multi_pasef) {
      for (i in seq_along(pasef_fis)) {
        pasef_subids <- pasef_ids[[i]]
        
        cl <- parallel::makeCluster(
          getOption("cl.cores", pasef_n_cores[[i]]), 
          outfile = file.path(mgf_path, "temp_dir", "log.txt"))
        raws[pasef_subids] <- parallel::clusterMap(
          cl, hdeisoDDA, 
          pasef_fis[[i]], pasef_subids, 
          MoreArgs = list(
            out_path = out_path, 
            mgf_path = mgf_path, 
            temp_dir = temp_dir, 
            mzml_type = mzml_type, 
            tol_ms1 = tol_ms1, 
            ppm_ms1_bin = ppm_ms1_bin, ppm_ms2_bin = ppm_ms2_bin, 
            maxn_mdda_precurs = maxn_mdda_precurs, 
            topn_ms2ions = topn_ms2ions, 
            n_mdda_flanks = n_mdda_flanks, 
            n_dia_scans = n_dia_scans, 
            min_mass = min_mass, max_mass = max_mass, 
            min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
            min_ms1_charge = min_ms1_charge, max_ms1_charge = max_ms1_charge, 
            min_ret_time = min_ret_time, max_ret_time = max_ret_time, 
            min_scan_num = min_scan_num, max_scan_num = max_scan_num, 
            deisotope_ms2 = deisotope_ms2, max_ms2_charge = max_ms2_charge, 
            ppm_ms1_deisotope = ppm_ms1_deisotope, 
            ppm_ms2_deisotope = ppm_ms2_deisotope, 
            grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
            mgf_cutmzs = mgf_cutmzs, 
            mgf_cutpercs = mgf_cutpercs, 
            quant = quant, 
            use_lfq_intensity = use_lfq_intensity, 
            ppm_ms1trace = ppm_ms1trace, 
            tmt_reporter_lower = tmt_reporter_lower, 
            tmt_reporter_upper = tmt_reporter_upper, 
            exclude_reporter_region = exclude_reporter_region, 
            use_defpeaks = use_defpeaks, 
            is_pasef = is_pasef,
            n_para = pasef_n_paras[[i]], 
            y_perc = .01,
            yco = yco
          ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
        parallel::stopCluster(cl)
      }
    }
    else {
      if (n_cores <= 1L) {
        raws <- mapply(
          hdeisoDDA, 
          peakfiles, file_indexes, 
          MoreArgs = list(
            out_path = out_path, 
            mgf_path = mgf_path, 
            temp_dir = temp_dir, 
            mzml_type = mzml_type, 
            tol_ms1 = tol_ms1, 
            ppm_ms1_bin = ppm_ms1_bin, ppm_ms2_bin = ppm_ms2_bin, 
            maxn_mdda_precurs = maxn_mdda_precurs, 
            topn_ms2ions = topn_ms2ions, 
            n_mdda_flanks = n_mdda_flanks, 
            n_dia_scans = n_dia_scans, 
            min_mass = min_mass, max_mass = max_mass, 
            min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
            min_ms1_charge = min_ms1_charge, max_ms1_charge = max_ms1_charge, 
            min_ret_time = min_ret_time, max_ret_time = max_ret_time, 
            min_scan_num = min_scan_num, max_scan_num = max_scan_num, 
            deisotope_ms2 = deisotope_ms2, max_ms2_charge = max_ms2_charge, 
            ppm_ms1_deisotope = ppm_ms1_deisotope, 
            ppm_ms2_deisotope = ppm_ms2_deisotope, 
            grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
            mgf_cutmzs = mgf_cutmzs, 
            mgf_cutpercs = mgf_cutpercs, 
            quant = quant, 
            use_lfq_intensity = use_lfq_intensity, 
            ppm_ms1trace = ppm_ms1trace, 
            tmt_reporter_lower = tmt_reporter_lower, 
            tmt_reporter_upper = tmt_reporter_upper, 
            exclude_reporter_region = exclude_reporter_region, 
            use_defpeaks = use_defpeaks, 
            is_pasef = is_pasef,
            n_para = n_para, 
            y_perc = .01, 
            yco = yco
          ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      }
      else {
        cl <- parallel::makeCluster(
          getOption("cl.cores", n_cores), 
          outfile = file.path(mgf_path, "temp_dir", "log.txt"))
        raws <- parallel::clusterMap(
          cl, hdeisoDDA, 
          peakfiles, file_indexes, 
          MoreArgs = list(
            out_path = out_path, 
            mgf_path = mgf_path, 
            temp_dir = temp_dir, 
            mzml_type = mzml_type, 
            tol_ms1 = tol_ms1, 
            ppm_ms1_bin = ppm_ms1_bin, ppm_ms2_bin = ppm_ms2_bin, 
            maxn_mdda_precurs = maxn_mdda_precurs, 
            topn_ms2ions = topn_ms2ions, 
            n_mdda_flanks = n_mdda_flanks, 
            n_dia_scans = n_dia_scans, 
            min_mass = min_mass, max_mass = max_mass, 
            min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
            min_ms1_charge = min_ms1_charge, max_ms1_charge = max_ms1_charge, 
            min_ret_time = min_ret_time, max_ret_time = max_ret_time, 
            min_scan_num = min_scan_num, max_scan_num = max_scan_num, 
            deisotope_ms2 = deisotope_ms2, max_ms2_charge = max_ms2_charge, 
            ppm_ms1_deisotope = ppm_ms1_deisotope, 
            ppm_ms2_deisotope = ppm_ms2_deisotope, 
            grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
            mgf_cutmzs = mgf_cutmzs, 
            mgf_cutpercs = mgf_cutpercs, 
            quant = quant, 
            use_lfq_intensity = use_lfq_intensity, 
            ppm_ms1trace = ppm_ms1trace, 
            tmt_reporter_lower = tmt_reporter_lower, 
            tmt_reporter_upper = tmt_reporter_upper, 
            exclude_reporter_region = exclude_reporter_region, 
            use_defpeaks = use_defpeaks, 
            is_pasef = is_pasef,
            n_para = n_para, 
            y_perc = .01,
            yco = yco
          ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
        parallel::stopCluster(cl)
      }
    }
  }
  
  message("Completed deisotoping at: ", Sys.time())
  # NULL raws (e.g., bad data witout MS2) removed
  raws <- unlist(raws, recursive = FALSE, use.names = TRUE)
  if (!length(raws)) {
    stop("All RAW files are without MS2.")
  }
  qs::qsave(raws, file.path(mgf_path, "raw_indexes.rds"), preset = "fast")
  
  type_acqu <- if (isTRUE(is_dia)) "dia" else "dda"
}


#' Helper of \link{loadMZML}.
#'
#' @param temp_dir A temporary file folder.
#' @inheritParams readMGF
#' @inheritParams matchMS
hloadMZML <- function (filelist = NULL, mgf_path = NULL, temp_dir = NULL)
{
  if (FALSE) {
    if (maxn_mdda_precurs > 1L && use_defpeaks) {
      warning("Default peaks not used at maxn_mdda_precurs > 1;", 
              "\nCoerce to use_defpeaks = FALSE.")
      use_defpeaks <- FALSE
    }
  }
  
  len <- length(filelist)
  files <- file.path(mgf_path, filelist)
  sizes <- max(unlist(lapply(files, file.size)))/1024^3
  n_cores <- min(detect_cores(32L), 
                 ceiling(find_free_mem()/1024/sizes/7.5 + 1), 
                 len)
  n_cores <- find_min_ncores(len, n_cores) 
  n_cores <- max(1L, n_cores)
  
  if (n_cores <= 1L) {
    out_names <- vector("list", len)
    
    for (i in 1:len) {
      out_names[[i]] <- loadMZML(files[[i]], temp_dir = temp_dir)
    }
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out_names <- parallel::clusterApply(cl, files, loadMZML, temp_dir = temp_dir)
    parallel::stopCluster(cl)
  }
  
  invisible(out_names)
}


#' Loads mzML from MSConvert.
#' 
#' For single mzML file.
#'
#' @param xml_file A file name of mzML.
#' @param temp_dir A temporary file folder.
loadMZML <- function (xml_file = NULL, temp_dir = NULL)
{
  ## spectrum
  xml_root <- xml2::read_xml(xml_file)
  mzML <- xml2::xml_child(xml_root)
  mzC <- xml2::xml_children(mzML)
  
  idx_sw <- which(xml2::xml_name(mzC) == "softwareList")
  softwares <- mzC[[idx_sw]]
  softwaresC <- xml2::xml_children(softwares)
  
  ###
  software_ids <- unlist(lapply(softwaresC, xml2::xml_attr, "id"))
  idx_pwiz <- which(software_ids == "pwiz")
  
  if (length(idx_pwiz)) {
    is_pwiz <- TRUE
    is_3d <- TRUE
  } else {
    idx_pwiz <- which(software_ids == "pwiz_Reader_Bruker") # 3L
    
    if (length(idx_pwiz)) {
      is_pwiz <- TRUE
      is_3d <- FALSE # 4-D ion mobility
    } else {
      is_pwiz <- TRUE
      is_3d <- FALSE
      idx_pwiz <- 3L
    }
  }
  
  if (is_pwiz) {
    if (is_3d) {
      mzml_type <- "pwiz_3d"
    } else {
      stop("MSConvert timsTOF supports by Mzion not yet available.")
    }
  }

  pwiz <- softwaresC[[idx_pwiz]]
  pwiz_ver <- xml2::xml_attr(pwiz, "version")
  pwiz_ver <- strsplit(pwiz_ver, ".", fixed = TRUE)[[1]]
  
  if ((len_pwiz <- length(pwiz_ver)) >= 2L) {
    pwiz_ver_major <- as.numeric(paste0(pwiz_ver[[1]], ".", pwiz_ver[[2]]))
  } else if (len_pwiz == 1L) {
    pwiz_ver_major <- as.numeric(pwiz_ver)
  } else {
    warning("Unknown MSConvert version.")
    pwiz_ver_major <- 3.0
  }
  
  if (pwiz_ver_major < 3) {
    stop("Use MSConvert version >= 3.0.")
  }
  
  raw_file <- local({
    idx_file <- which(xml2::xml_name(mzC) == "fileDescription")
    file_des <- mzC[[idx_file]]
    idx_srcl <- which(xml2::xml_name(xml2::xml_children(file_des)) == "sourceFileList")
    info_raw <- xml2::xml_children(file_des)[[idx_srcl]]
    idx_srcf <- which(xml2::xml_name(xml2::xml_children(info_raw)) == "sourceFile")
    
    if (length(idx_srcf) == 1L) {
      info_fi  <- xml2::xml_children(info_raw)[[idx_srcf]]
      raw_file <- xml2::xml_attr(info_fi, "name")
    } else { # timsTOF
      idx_srcf <- idx_srcf[[1]]
      info_fi  <- xml2::xml_children(info_raw)[[idx_srcf]]
      raw_file <- xml2::xml_attr(info_fi, "id")
      raw_file <- gsub("(^.*\\.d)_.*\\.tdf$", "\\1", raw_file)
    }
  })
  
  idx_run <- which(xml2::xml_name(mzC) == "run")
  run <- mzC[[idx_run]]
  idx_specs <- which(xml2::xml_name(xml2::xml_children(run)) == "spectrumList")
  spec <- xml2::xml_children(xml2::xml_children(run)[[idx_specs]])
  rm(list = c("mzML", "mzC", "idx_run", "run", "idx_specs", "idx_sw", 
              "softwares", "softwaresC", "idx_pwiz", "pwiz", "pwiz_ver", 
              "pwiz_ver_major", "len_pwiz"))
  
  # - mzML
  #   - ...
  #   - fileDescription
  #     - ...
  #     - sourceFileList
  #       - sourceFile
  
  # - run
  #  - spectrumList
  #   - spectrum # x
  #     - ...
  #     - <cvParam ... name="ms level" value="2"/>
  #     - <cvParam ... name="total ion current" value="2578.8672"/>
  #     - <cvParam ... name="spectrum title" value="23aug2017_hela_serum_timecourse_4mz_narrow_1.1.1. 
  #                    File:"FILENAME.raw", ...originalScan=1 demux=0 scan=1/>
  #     - scanList # xc
  #       - ...
  #       - <scan ... originalScan=3 demux=1 scan=4> # originalScan=3 demux=1 scan=4 only only for MS2
  #         - <cvParam ... name="scan start time" value="0.012728725" .../>
  #         - <cvParam ... name="ion injection time" value="54.999999701977" .../>
  #     - precursorList # xc
  #       - <precursor ... originalScan=3 demux=1 scan=4"> # redundant?
  #         - <cvParam ... name="isolation window target m/z" value="401.432342529297" .../>
  #         - selectedIonList
  #           - selectedIon
  #             -  <cvParam ... name="selected ion m/z" value="401.432342529297" .../>
  #             -  <cvParam ... name="charge state" value="3" .../> # No charge state if DIA
  #             -  <cvParam ... name="peak intensity" value="0" .../> # 0 if DIA
  #     - binaryDataArrayList # xc
  #       - binaryDataArray
  #         - binary
  #       - binaryDataArray
  #         - binary
  #  - chromatogramList
  #   - chromatogram
  #    - binaryDataArrayList
  #     - binaryDataArray # name="time array"
  #     - binaryDataArray # name="intensity array"
  #     - binaryDataArray # name="non-standard data array"
  

  ## the first scan
  x <- spec[[1]]
  ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))[[1]]
  ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))
  id_nms <- lapply(ids, `[[`, 1)
  idx_demux <- which(id_nms == "demux")
  is_demux <- if (length(idx_demux)) TRUE else FALSE
  
  if (is_3d) {
    idx_sc <- which(id_nms == "scan")
    if (!length(idx_osc <- which(id_nms == "originalScan"))) idx_osc <- idx_sc
    idx_sc_start <- idx_sc_end <- idx_osc
  } else {
    idx_sc_start <- which(id_nms == "scanStart") # 3
    idx_sc_end <- which(id_nms == "scanEnd") # 4
    idx_osc <- idx_sc <- idx_sc_start
    # c(idx_sc_start, idx_sc_end) # length(2)
  }
  
  xc <- xml2::xml_children(x)
  xcp_attrs <- xml2::xml_attrs(xc[which(xml2::xml_name(xc) == "cvParam")])
  xcp_names <- lapply(xcp_attrs, `[[`, "name")
  xcp_vals <- lapply(xcp_attrs, `[[`, "value")
  idx_mslev <- which(xcp_names == "ms level")
  idx_title <- which(xcp_names == "spectrum title")
  rm(list = c("x", "ids", "xc", "xcp_attrs", "xcp_names", "xcp_vals", "id_nms"))
  
  # MS2 indexes (no MS1 with DDA -> MS2 goes first)
  len <- length(spec)
  rng <- max(floor(len/2L), 1L):len
  allowance <- min(100L, len)
  count <- 0L
  
  if (!len) {
    stop("No spectrum data found.")
  }
  
  for (i in rng) {
    # i=46
    # i = 6000
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (xml2::xml_attr(xc[[idx_mslev]], "value") == "2") {
      if (is_3d) {
        idx_scan_lwr_2 <- grep("lowest observed m/z", xc)
        if (!length(idx_scan_lwr_2)) {
          warning("Fields of `lowest observed m/z` not found.")
          idx_scan_lwr_2 <- 8L
        }
      } else {
        idx_scan_lwr_2 <- grep("ion mobility lower limit", xc)
        if (!length(idx_scan_lwr_2)) {
          warning("Fields of `ion mobility lower limit` not found.")
          idx_scan_lwr_2 <- 8L
        }
      }
      
      if (is_3d) {
        idx_scan_upr_2 <- grep("highest observed m/z", xc)
        if (!length(idx_scan_upr_2)) {
          warning("Fields of `highest observed m/z` not found.")
          idx_scan_upr_2 <- 9L
        }
      } else {
        idx_scan_upr_2 <- grep("ion mobility upper limit", xc)
        if (!length(idx_scan_upr_2)) {
          warning("Fields of `ion mobility upper limit` not found.")
          idx_scan_upr_2 <- 9L
        }
      }
      
      idx_precursor_2 <- grep("precursorList", xc)
      if (!length(idx_precursor_2)) {
        warning("Fields of `precursorList` not found.")
        idx_precursor_2 <- if (is_3d) 12L else 11L
      }
      
      idx_scanList_2 <- grep("scanList", xc)
      if (!length(idx_scanList_2)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_2 <- if (is_3d) 11L else 10L
      }
      
      idx_bin_2 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_2)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_2 <- if (is_3d) 13L else 12L
      }
      
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      
      # !!! a range if timsTOF
      idx_rt_2 <- which(xml2::xml_name(scanList) == "scan")
      if (length(idx_rt_2)) {
        idx_rt_2 <- idx_rt_2[[1]]
      } else {
        warning("Fields of `scan` not found.")
        idx_rt_2 <- 2L
      }
      
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      scanList_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      
      idx_scan_start_2 <- which(scanList_ret_attrs == "scan start time") # 1
      if (!length(idx_scan_start_2)) {
        warning("Fields of `scan start time` not found.")
        idx_scan_start_2 <- 1L
      }
      
      if (is_3d) {
        if (is_pwiz) {
          idx_ms2_reso <- which(scanList_ret_attrs == "mass resolving power")
          if (length(idx_ms2_reso)) {
            ms2_reso <- scanList_ret[[idx_ms2_reso]]
            ms2_reso <- as.integer(xml2::xml_attr(ms2_reso, "value"))
          } else {
            warning("Fields of `mass resolving power` not found.")
            ms2_reso <- 120000L
          }
        } else { # MSFragger timsTOF
          idx_ms2_reso <- which(scanList_ret_attrs == "inverse reduced ion mobility")
          if (length(idx_ms2_reso)) {
            ms2_reso <- scanList_ret[[idx_ms2_reso]]
            ms2_reso <- as.numeric(xml2::xml_attr(ms2_reso, "value"))
          } else {
            warning("Fields of `inverse reduced ion mobility` not found.")
            ms2_reso <- 1L
          }
        }
      } else {
        idx_ms2_reso <- which(scanList_ret_attrs == "inverse reduced ion mobility")
        if (length(idx_ms2_reso)) {
          ms2_reso <- scanList_ret[[idx_ms2_reso]]
          ms2_reso <- as.numeric(xml2::xml_attr(ms2_reso, "value"))
        } else {
          warning("Fields of `inverse reduced ion mobility` not found.")
          ms2_reso <- 1
        }
      }
      
      # entire MS2 is empty
      if (length(xc) < idx_precursor_2) {
        next
      }
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # assume one precursor
      precursorc <- xml2::xml_children(precursor)
      
      idx_selectedIonList <- grep("selectedIonList", precursorc)
      if (!length(idx_selectedIonList)) {
        warning("Fields of `selectedIonList` not found.")
        idx_selectedIonList <- 2L
      }
      
      idx_isolationWindow <- grep("isolationWindow", precursorc)
      if (!length(idx_isolationWindow)) {
        warning("Fields of `isolationWindow` not found.")
        idx_isolationWindow <- 1L
      }
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_nms <- lapply(isolationWindowc, function (x) xml2::xml_attr(x, "name"))
      
      idx_ctrmz <- which(iso_nms == "isolation window target m/z")
      if (!length(idx_ctrmz)) {
        warning("Fields of `isolation window target m/z` not found.")
        idx_ctrmz <- 1L
      }
      
      idx_lwrmz <- which(iso_nms == "isolation window lower offset")
      if (!length(idx_lwrmz)) {
        warning("Fields of `isolation window lower offset` not found.")
        idx_lwrmz <- 2L
      }
      
      idx_uprmz <- which(iso_nms == "isolation window upper offset")
      if (!length(idx_uprmz)) {
        warning("Fields of `isolation window upper offset` not found.")
        idx_uprmz <- 3L
      }
      
      lwr_offset <- 
        as.numeric(xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value"))
      upr_off_set <- 
        as.numeric(xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value"))
      iso_width <- upr_off_set + lwr_offset
      
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      selion_nms <- lapply(selectedIonc, function (x) xml2::xml_attr(x, "name"))
      
      idx_moverz <- which(selion_nms == "selected ion m/z")
      if (!length(idx_moverz)) {
        warning("Fields of `selected ion m/z` not found.")
        idx_moverz <- 1L
      }
      
      # not with timsTOF?
      idx_ms1int <- which(selion_nms == "peak intensity") # DIA: zero intensity; timsTOF: length(0)
      if (!length(idx_ms1int)) {
        if (count <= allowance) {
          count <- count + 1L
          next
        }
        else {
          warning("The field of `peak intensity` not found.")
          idx_ms1int <- 3L
        }
      }
      
      idx_charge <- which(selion_nms == "charge state") # DIA: no "charge state"
      is_dia <- if (length(idx_charge)) FALSE else TRUE
      
      # timsTOF
      if (is_dia && !is_3d) is_dia <- FALSE
      
      rm(list = c("scanList", "scanList_ret", "scanList_ret_attrs", 
                  "precursorList", "precursor", "precursorc", 
                  "selectedIon", "selectedIonc", "selion_nms"))
      break
    }
  }
  
  # MS1 indexes
  for (i in rng) {
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (as.integer(xml2::xml_attr(xc[[idx_mslev]], "value")) == 1) {
      if (is_3d) {
        idx_scan_lwr_1 <- grep("lowest observed m/z", xc) # "ion mobility lower limit"
        if (!length(idx_scan_lwr_1)) {
          warning("Fields of `lowest observed m/z` not found.")
          idx_scan_lwr_1 <- 8L
        }
        
        idx_scan_upr_1 <- grep("highest observed m/z", xc) # "ion mobility upper limit"
        if (!length(idx_scan_upr_1)) {
          warning("Fields of `highest observed m/z` not found.")
          idx_scan_upr_1 <- 9L
        }
      } else {
        idx_scan_lwr_1 <- grep("ion mobility lower limit", xc) # "ion mobility lower limit"
        if (!length(idx_scan_lwr_1)) {
          warning("Fields of `ion mobility lower limit` not found.")
          idx_scan_lwr_1 <- 8L
        }
        
        idx_scan_upr_1 <- grep("ion mobility upper limit", xc) # "ion mobility upper limit"
        if (!length(idx_scan_upr_1)) {
          warning("Fields of `ion mobility upper limit` not found.")
          idx_scan_upr_1 <- 9L
        }
      }
      
      idx_scanList_1 <- grep("scanList", xc)
      if (!length(idx_scanList_1)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_1 <- if (is_3d) 11L else 10L
      }
      
      idx_bin_1 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_1)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_1 <- if (is_3d) 12L else 11L
      }
      
      local({
        # [[1]] name="m/z array"; [[2]] name="intensity array"; [[3]] name="mean inverse reduced ion mobility array"
        binData <- xml2::xml_children(xc[[idx_bin_1]])[[1]]
        binDataC <- xml2::xml_children(binData)
        binDataC <- binDataC[which(xml2::xml_name(binDataC) == "cvParam")]
        
        lapply(binDataC, function (v) {
          attrs <- xml2::xml_attrs(v)
          oks <- grepl("compression", attrs)
          
          if (any(oks))
            if (attrs[oks][[1]] == "zlib compression")
              stop("Please uncheck zlib compression in  mzML generations.")
        })
      })
      
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      idx_rt_1 <- which(xml2::xml_name(scanList) == "scan")
      if (!length(idx_rt_1)) {
        warning("Fields of `scan` not found.")
        idx_rt_1 <- 2L
      }
      
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      
      scan_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      ## for later removal of the assumption of a fixed resolution with Astral
      # ms1_reso <- scanList_ret[[which(scan_ret_attrs == "mass resolving power")]]
      # ms1_reso <- as.integer(xml2::xml_attr(ms1_reso, "value"))
      idx_scan_start_1 <- which(scan_ret_attrs == "scan start time")
      if (!length(idx_scan_start_1)) {
        warning("Fields of `scan start time` not found.")
        idx_scan_start_1 <- 1L
      }
      
      # idx_scan_start_im <- which(scan_ret_attrs == "inverse reduced ion mobility")
      
      rm(list = c("scanList", "scanList_ret", "scan_ret_attrs"))
      break
    }
  }
  
  rm(list = c("x", "xc", "i"))
  
  if (is_dia) {
    filename <- extrDIA(
      spec = spec, raw_file = raw_file, 
      temp_dir = temp_dir, is_demux = is_demux, 
      idx_sc = idx_sc, idx_osc = idx_osc, idx_mslev = idx_mslev, 
      idx_title = idx_title, idx_scanList_1 = idx_scanList_1, 
      idx_scanList_2 = idx_scanList_2, idx_rt_1 = idx_rt_1, 
      idx_rt_2 = idx_rt_2, idx_scan_start_1 = idx_scan_start_1, 
      idx_scan_start_2 = idx_scan_start_2, 
      idx_precursor_2 = idx_precursor_2, 
      idx_isolationWindow = idx_isolationWindow, 
      idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz, 
      idx_uprmz = idx_uprmz, 
      idx_selectedIonList = idx_selectedIonList, 
      idx_demux = idx_demux, idx_bin_1 = idx_bin_1, 
      idx_bin_2 = idx_bin_2, 
      idx_scan_lwr_1 = idx_scan_lwr_1, 
      idx_scan_upr_1 = idx_scan_upr_1, 
      idx_scan_lwr_2 = idx_scan_lwr_2, 
      idx_scan_upr_2 = idx_scan_upr_2)
  }
  else {
    filename <- extrDDA(
      spec = spec, raw_file = raw_file, temp_dir = temp_dir, 
      idx_sc = idx_sc, idx_osc = idx_osc, 
      idx_mslev = idx_mslev, idx_title = idx_title, 
      idx_scanList_2 = idx_scanList_2, idx_rt_2 = idx_rt_2, 
      idx_scan_start_2 = idx_scan_start_2, 
      idx_precursor_2 = idx_precursor_2, 
      idx_isolationWindow = idx_isolationWindow, 
      idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz,
      idx_uprmz = idx_uprmz, 
      idx_selectedIonList = idx_selectedIonList, 
      idx_moverz = idx_moverz, idx_charge = idx_charge, 
      idx_ms1int = idx_ms1int, idx_bin_2 = idx_bin_2, 
      idx_scanList_1 = idx_scanList_1, 
      idx_rt_1 = idx_rt_1, idx_scan_start_1 = idx_scan_start_1, 
      idx_bin_1 = idx_bin_1, 
      idx_scan_lwr_1 = idx_scan_lwr_1, idx_scan_upr_1 = idx_scan_upr_1, 
      idx_scan_lwr_2 = idx_scan_lwr_2, idx_scan_upr_2 = idx_scan_upr_2)
  }
  
  attr(filename, "is_dia") <- is_dia
  attr(filename, "iso_width") <- iso_width
  attr(filename, "mzml_type") <- mzml_type

  invisible(filename)
}


#' Proc mzML data for mDDA workflows.
#'
#' @param spec Spectrum list.
#' @param raw_file A raw file name.
#' @param temp_dir A temporary directory.
#' @param idx_sc Index of scan numbers.
#' @param idx_osc Index of original scan numbers.
#' @param idx_mslev Index of MS levels.
#' @param idx_title Index of scan titles.
#' @param idx_scanList_2 Index of MS2 scanList.
#' @param idx_rt_2 Index of MS2 retention times.
#' @param idx_scan_start_2 Index of MS2 scan start.
#' @param idx_precursor_2 Index of MS2 precursor.
#' @param idx_isolationWindow Index of an isolationWindow.
#' @param idx_ctrmz Index of center m-over-z of an isolationWindow.
#' @param idx_lwrmz Index of lower m-over-z an isolationWindow.
#' @param idx_uprmz Index of upper m-over-z an isolationWindow.
#' @param idx_selectedIonList Index of selectedIonList.
#' @param idx_moverz Index of m-over-z sequences.
#' @param idx_charge Index of charge states.
#' @param idx_ms1int Index of MS1 intensities.
#' @param idx_bin_2 Index of MS2 binary data.
#' @param idx_scan_lwr_1 Index of lowest MS1 mass.
#' @param idx_scan_upr_1 Index of highest MS1 mass.
#' @param idx_scan_lwr_2 Index of lowest MS2 mass.
#' @param idx_scan_upr_2 Index of highest MS2 mass.
#' @param idx_scanList_1 Index of MS1 scanList.
#' @param idx_rt_1 Index of precursor retention times.
#' @param idx_scan_start_1 Index of MS1 scan starts.
#' @param idx_bin_1 Index of MS1 binary data.
#' @param raw_file The RAW file name of an mzML.
#' @param is_regular Logical; TRUE at Thermo's and FALSE at timsTOF's .d.
extrDDA <- function (spec = NULL, raw_file = NULL, temp_dir = NULL, 
                     idx_sc = 3L, idx_osc = 3L, idx_mslev = 2L, 
                     idx_title = 10L, idx_scanList_2 = 11L, idx_rt_2 = 2L, 
                     idx_scan_start_2 = 1L, idx_precursor_2 = 12L, 
                     idx_isolationWindow = 1L, idx_ctrmz = 1L, idx_lwrmz = 2L,
                     idx_uprmz = 3L, idx_selectedIonList = 2L, idx_moverz = 1L, 
                     idx_charge = 2L, idx_ms1int = 3L, idx_bin_2 = 13L, 
                     idx_scanList_1 = 11L, idx_rt_1 = 2L, 
                     idx_scan_start_1 = 1L, idx_bin_1 = 12L, 
                     idx_scan_lwr_1 = 8L, idx_scan_upr_1 = 9L, 
                     idx_scan_lwr_2 = 8L, idx_scan_upr_2 = 9L, is_regular = TRUE)
{
  len <- length(spec)
  
  ret_times <- orig_scans <- scan_nums <- scan_titles <- 
    iso_ctr <- iso_lwr <- iso_upr <- character(len)
  
  # msx_: both MS1 and MS2, differentiated by ms_lev
  # ms0_: precursor info by other peak-pickings, e.g., MSConvert
  
  ms0_moverzs <- ms0_ints <- ms0_charges <- 
    msx_moverzs <- msx_ints <- msx_charges <- vector("list", len)
  ms_levs <- msx_ns <- integer(len)
  
  if (!is_regular) {
    msx_mobils <- vector("list", len)
  }

  for (i in 1:len) {
    x <- spec[[i]]
    ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))[[1]]
    ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))
    scan_nums[[i]]  <- ids[[idx_sc]][[2]]
    orig_scans[[i]] <- ids[[idx_osc]][[2]]
    
    xc <- xml2::xml_children(x)
    ms_levs[[i]] <- ms_lev <- as.integer(xml2::xml_attr(xc[[idx_mslev]], "value"))
    scan_titles[[i]] <- xml2::xml_attr(xc[[idx_title]], "value")
    
    if (ms_lev == 2L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_2]], "value")
      
      # entire MS2 is empty
      if (length(xc) < idx_precursor_2) {
        next
      }

      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # (assume one precursor by MSConvert)
      precursorc <- xml2::xml_children(precursor)
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_ctr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_ctrmz]], "value")
      iso_lwr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value")
      iso_upr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value")
      
      selectedIon  <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      
      # precursor(s) of MS2
      ms0_moverzs[[i]] <- as.numeric(xml2::xml_attr(selectedIonc[[idx_moverz]], "value"))
      ms0_charges[[i]] <- as.numeric(xml2::xml_attr(selectedIonc[[idx_charge]], "value"))
      # no precursor intensity: 01CPTAC3_Benchmarking_P_BI_20170523_BL_f12.mzML; scan=54837
      ms0_ints[[i]] <- if (length(selectedIonc) >= idx_ms1int)
        as.numeric(xml2::xml_attr(selectedIonc[[idx_ms1int]], "value"))
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_2]]))
    }
    else if (ms_lev == 1L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_1]], "value")
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_1]]))
    }
    
    # full MS1 or MS2 spectrum
    msData <- xml2::xml_contents(binData)
    len_ms <- length(msData)
    
    if (len_ms == 2L) {
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_moverzs[[i]] <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ints[[i]] <- readBin(r2, "double", n = msx_n, size = 8L)
    } else if (len_ms == 3L) { # timsTOF
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      r3 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[3]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_moverzs[[i]] <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ints[[i]] <- readBin(r2, "double", n = msx_n, size = 8L)
      msx_mobils[[i]] <- readBin(r3, "double", n = msx_n, size = 8L)
    }
  }
  
  out <- list(
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    ms1_moverzs = ms0_moverzs, 
    ms1_charges = ms0_charges, 
    ms1_ints = ms0_ints, 
    
    scan_title = scan_titles,
    raw_file = raw_file, # single
    ms_level = ms_levs, 
    # mzML: ret_times in minutes; MGF: in seconds
    ret_time = as.numeric(ret_times) * 60, 
    scan_num = as.integer(scan_nums), 
    orig_scan = orig_scans,
    iso_ctr = as.numeric(iso_ctr), 
    iso_lwr = as.numeric(iso_lwr), 
    iso_upr = as.numeric(iso_upr)
  )
  
  ###
  out$iso_lwr <- out$iso_ctr - out$iso_lwr
  out$iso_upr <- out$iso_ctr + out$iso_upr
  ###
  
  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Helper of \link{deisoDDA}.
#'
#' @param raw_id A raw file id.
#' @param n_para The allowance of parallel processing.
#' @param mgf_cutmzs The cut points in peak lists.
#' @param mgf_cutpercs The percentage of cuts.
#' @param tol_ms1 The (original) MS1 tolerance (\code{20 * 1E-6}).
#' @param ppm_ms1_bin The three-bin MS1 ppm.
#' @param ppm_ms2_bin The three-bin MS2 ppm.
#' @inheritParams deisoDDA
#' @inheritParams matchMS
hdeisoDDA <- function (filename, raw_id = 1L, 
                       out_path = NULL, mgf_path = NULL, temp_dir = NULL, 
                       mzml_type = "raw", ppm_ms1_bin = 10L, ppm_ms2_bin = 10L, 
                       tol_ms1 = 2E-5, maxn_mdda_precurs = 5L, 
                       topn_ms2ions = 150L, n_mdda_flanks = 6L, n_dia_scans = 6L, 
                       min_mass = 200L, max_mass = 4500L,
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       quant = "none", 
                       use_lfq_intensity = TRUE, ppm_ms1trace = 6L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, use_defpeaks = FALSE, 
                       is_pasef = FALSE, n_para = 1L, y_perc = .01, yco = 100)
{
  df <- deisoDDA(
    filename, 
    out_path = out_path, 
    mgf_path = mgf_path, 
    temp_dir = temp_dir, 
    mzml_type = mzml_type, 
    tol_ms1 = tol_ms1, 
    ppm_ms1_bin = ppm_ms1_bin, 
    ppm_ms2_bin = ppm_ms2_bin, 
    maxn_mdda_precurs = maxn_mdda_precurs, 
    topn_ms2ions = topn_ms2ions, 
    n_mdda_flanks = n_mdda_flanks, 
    n_dia_scans = n_dia_scans, 
    min_mass = min_mass, 
    max_mass = max_mass, 
    min_ms2mass = min_ms2mass, 
    max_ms2mass = max_ms2mass, 
    min_ms1_charge = min_ms1_charge, 
    max_ms1_charge = max_ms1_charge, 
    min_ret_time = min_ret_time, 
    max_ret_time = max_ret_time, 
    min_scan_num = min_scan_num, 
    max_scan_num = max_scan_num, 
    deisotope_ms2 = deisotope_ms2, 
    max_ms2_charge = max_ms2_charge, 
    ppm_ms1_deisotope = ppm_ms1_deisotope, 
    ppm_ms2_deisotope = ppm_ms2_deisotope, 
    grad_isotope = grad_isotope, 
    fct_iso2 = fct_iso2, 
    quant = quant, 
    use_lfq_intensity = use_lfq_intensity, 
    ppm_ms1trace = ppm_ms1trace, 
    tmt_reporter_lower = tmt_reporter_lower, 
    tmt_reporter_upper = tmt_reporter_upper, 
    exclude_reporter_region = exclude_reporter_region, 
    use_defpeaks = use_defpeaks, 
    is_pasef = is_pasef, 
    n_para = n_para, 
    y_perc = y_perc,
    yco = yco)
  
  # e.g., bad data without MS2
  if (is.null(df)) {
    return(NULL)
  }
  
  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = df[["ms2_moverzs"]], 
                          ms2_ints = df[["ms2_ints"]], 
                          ms2_charges = df[["ms2_charges"]], 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  df[["ms2_moverzs"]] <- mz_n_int[["ms2_moverzs"]]
  df[["ms2_ints"]] <- mz_n_int[["ms2_ints"]]
  df[["ms2_charges"]] <- mz_n_int[["ms2_charges"]]
  df[["ms2_n"]] <- mz_n_int[["lens"]]
  
  # fuzzy precursor masses and the corresponding scan_nums are duplicated 
  #  In `calcpepsc`, uniq_id: scan_num + raw_file + ms1_offset
  #  In `calc_pepfdr`, uniq_id: scan_num + raw_file
  #  In `post_pepfdr`, `calc_peploc`, `hcalc_tmtint`, `add_rptrs`
  
  post_readmgf(df, raw_id = raw_id, mgf_path = mgf_path)
}



#' Deisotopes DDA.
#' 
#' @param filename A peaklist filename.
#' @param temp_dir A temp_dir to the filename.
#' @param mzml_type The type of peak lists. To trigger MS1 deisotoping at
#'   \code{mzml_type == "raw"} (MS1 X, y and Z are NA with Mzion).
#' @param is_pasef Logical; is TIMS TOF data or not.
#' @param n_para The allowance of parallel processing.
#' @param y_perc The cut-off in intensity values in relative to the base peak.
#' @param yco The absolute cut-off in intensity values.
#' @param tol_ms1 The (original) MS1 tolerance (\code{20 * 1E-6}).
#' @param ppm_ms1_bin The three-bin MS1 ppm.
#' @param ppm_ms2_bin The three-bin MS2 ppm.
#' @inheritParams matchMS
deisoDDA <- function (filename = NULL, 
                      out_path = NULL, mgf_path = NULL, temp_dir = NULL, 
                      mzml_type = "raw", tol_ms1 = 2E-5, 
                      ppm_ms1_bin = 10L, ppm_ms2_bin = 10L, 
                      maxn_mdda_precurs = 5L, topn_ms2ions = 150L, 
                      n_mdda_flanks = 6L, n_dia_scans = 6L, 
                      min_mass = 200L, max_mass = 4500L,
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      min_ms1_charge = 2L, max_ms1_charge = 4L, 
                      min_ret_time = 0, max_ret_time = Inf, 
                      min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                      deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                      ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                      grad_isotope = 1.6, fct_iso2 = 3.0, 
                      quant = "none", 
                      use_lfq_intensity = TRUE, ppm_ms1trace = 6L, 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, use_defpeaks = FALSE, 
                      is_pasef = FALSE, n_para = 1L, y_perc = .01, yco = 100)
{
  ###
  # msx_: full spectra of ms1 and ms2, differentiated by ms_lev
  ###
  
  if (TRUE) {
    df <- predeisoDDA(
      filename, 
      temp_dir = temp_dir, 
      mzml_type = mzml_type, 
      ppm_ms1_bin = ppm_ms1_bin, 
      ppm_ms2_bin = ppm_ms2_bin, 
      maxn_mdda_precurs = maxn_mdda_precurs, 
      topn_ms2ions = topn_ms2ions, 
      n_mdda_flanks = n_mdda_flanks, 
      n_dia_scans = n_dia_scans, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass, 
      max_ms2mass = max_ms2mass, 
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      deisotope_ms2 = deisotope_ms2, 
      max_ms2_charge = max_ms2_charge, 
      ppm_ms1_deisotope = ppm_ms1_deisotope, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2, 
      quant = quant, 
      use_lfq_intensity = use_lfq_intensity, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region, 
      use_defpeaks = use_defpeaks, 
      is_pasef = is_pasef, 
      n_para = n_para)
    
    qs::qsave(df, file.path(temp_dir, paste0("predeisoDDA_", filename)), 
              preset = "fast")
  }
  else {
    df <- qs::qread(file.path(temp_dir, paste0("predeisoDDA_", filename)))
  }
  
  # e.g., bad data without MS2
  if (is.null(df)) {
    return(df)
  }
  
  ## LFQ: replaces intensities with apex values
  # No MS1 info at maxn_mdda_precurs == 0L
  if (maxn_mdda_precurs) {
    if (!dir.exists(path_ms1 <- file.path(out_path, "ms1data"))) {
      path_ms1 <- create_dir(path_ms1) # for MBR in proteoQ
    }
    
    qs::qsave(
      df[with(df, ms_level == 1L), 
         c("ret_time", "scan_num", "orig_scan", "msx_moverzs", "msx_ints")], 
      file.path(path_ms1, paste0("ms1full_", filename)), preset = "fast")
  }
  
  if (use_lfq_intensity && maxn_mdda_precurs) {
    df     <- get_ms1xs_space(df)
    rows1  <- which(df$ms_level == 1L)
    rt_gap <- estimate_rtgap(df$ret_time[rows1], d = 180)
    rt_gap <- as.integer(rt_gap * 1.5) # for early RTs
    
    # exclude MS2 scans before the first MS1 (timsTOF)
    if (row1 <- rows1[[1]] - 1L) {
      df <- df[-seq_len(row1), ]
    }
    rm(list = c("rows1", "row1"))
    
    step_tr <- ppm_ms1trace * 1e-6
    
    if (TRUE) {
      ans_prep <- pretraceXY(
        df[, c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", "ms_level", 
               "msx_moverzs", "msx_ints", "msx_charges", "orig_scan", "ret_time")], 
        from = min_mass, step = step_tr, 
        # n_chunks = ceiling(sum(df$ms_level == 1L) / 512L), 
        rt_gap = rt_gap, rt_tol = 180)
      # qs::qsave(ans_prep, file.path(path_ms1, paste0("msxspace_", filename)), preset = "fast")
    }
    else {
      ans_prep <- qs::qread(file.path(path_ms1, paste0("msxspace_", filename)))
    }
    
    dfs  <- ans_prep$dfs
    df1s <- ans_prep$df1s
    lenv <- length(df1s)
    # gaps <- unlist(ans_prep$gaps, use.names = FALSE, recursive = FALSE)
    # gaps_bf <- ans_prep$gaps_bf <- c(0L, gaps[1:(lenv - 1L)])
    # gaps_af <- ans_prep$gaps_af <- c(gaps[2:lenv], 0L)
    gaps_bf <- ans_prep$gaps_bf
    gaps_af <- ans_prep$gaps_af
    rm(list = c("ans_prep", "rt_gap"))
    
    fn_traces <- gsub("\\.rds$", "", filename)
    fn_traces <- paste0("ms1apexes_", fn_traces, "_", seq_along(dfs), ".rds")
    
    cols  <- c("ms_level", "ms1_moverz", "ms1_int", "orig_scan")
    
    if (is_pasef) {
      min_y   <- 500
      min_n1  <- 10L
      min_n2  <- 20L
      min_n3  <- 15L
      ytot_co <- 5000
    }
    else {
      min_y   <- 2e6
      min_n1  <- 10L
      min_n2  <- 20L
      min_n3  <- 15L
      ytot_co <- 2e5
    }

    if (TRUE) {
      cl <- parallel::makeCluster(getOption("cl.cores", 2L))
      out <- parallel::clusterMap(
        cl, htraceXY, 
        xs = lapply(df1s, `[[`, "msx_moverzs"), 
        ys = lapply(df1s, `[[`, "msx_ints"), 
        ss = lapply(df1s, `[[`, "orig_scan"), 
        ts = lapply(df1s, `[[`, "ret_time"), 
        df = lapply(dfs, `[`, cols), 
        gap_bf = gaps_bf, 
        gap_af = gaps_af, 
        out_name = fn_traces, 
        MoreArgs = list(
          n_dia_scans = n_dia_scans, from = min_mass, step_tr = step_tr, 
          tol_ms1 = tol_ms1, y_perc = y_perc, yco = yco, ytot_co = ytot_co, 
          min_y = min_y, min_n1 = min_n1, min_n2 = min_n2, min_n3 = min_n3, 
          path_ms1 = path_ms1
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
      parallel::stopCluster(cl)
      out <- dplyr::bind_rows(out)
      # qs::qsave(out, file.path(path_ms1, paste0("htraceXY_", filename)), preset = "fast")
    }
    else if (FALSE) {
      vxs <- lapply(df1s, `[[`, "msx_moverzs")
      vys <- lapply(df1s, `[[`, "msx_ints")
      vss <- lapply(df1s, `[[`, "orig_scan") # for troubleshooting
      vts <- lapply(df1s, `[[`, "ret_time")
      vdf <- lapply(dfs,  `[`,   cols)
      out <- vector("list", lenv)
      # sapply(vss, `[[`, 1L)
      
      for (i in 1:lenv) {
        out[[i]] <- htraceXY(
          xs = vxs[[i]], ys = vys[[i]], ss = vss[[i]], ts = vts[[i]], 
          df = vdf[[i]], gap_bf = gaps_bf[[i]], gap_af = gaps_af[[i]], 
          out_name = fn_traces[[i]], n_dia_scans = n_dia_scans, 
          from = min_mass, step_tr = step_tr, tol_ms1 = tol_ms1, 
          y_perc = y_perc, yco = yco, min_y = min_y, ytot_co = ytot_co, 
          min_n1 = min_n1, min_n2 = min_n2, min_n3 = min_n3, 
          path_ms1 = path_ms1)
      }
      rm(list = c("vxs", "vys", "vdf", "lenv"))
      
      out <- dplyr::bind_rows(out)
    }
    else {
      out <- qs::qread(file.path(path_ms1, paste0("htraceXY_", filename)))
    }
    
    qs::qsave(
      out[, c("ms_level", "orig_scan", "apex_ps", "apex_xs", "apex_ys", 
              "apex_ts")] |> 
        dplyr::filter(ms_level > 1L), 
      file.path(path_ms1, paste0("ms1apexes_", filename)), preset = "fast")
    
    if (nrow(df) == nrow(out)) {
      df[, cols] <- out[, cols]
    }
    else {
      stop("Developer: checks for row dropping in tracing MS1.")
    }
    
    # is_list <- class(out$apex_scan_num) == "list" # must be
    df$apex_scan_num <- out$apex_scan_num
    rm(list = ("out"))
    
    df <- local({
      rows2 <- df$ms_level != 1L
      df2 <- df[rows2, ]
      empties <- lengths(df2$apex_scan_num) == 0L # list(character(0))
      df2$apex_scan_num[empties] <- as.list(NA_integer_) # or as.list(0L)
      
      # some apexes remain 0 or c(0, 0)...

      # if (use_ms2scan_for_ms1apex <- FALSE) {
      #   df2$apex_scan_num[empties] <- as.list(df2$scan_num[empties])
      # }
      # else {
      #   
      # }

      # chimeric MS1 entries
      lens2m <-lengths(df2$ms1_moverz)
      lens2a <- lengths(df2$apex_scan_num)
      bads   <- which(lens2m & (lens2m != lens2a))
      
      if (length(bads)) {
        dfx <- df2[bads, ]
        dfx$apex_scan_num <- 
          mapply(function (x, n) rep_len(x, n), 
                 dfx$apex_scan_num, lens2m[bads],
                 USE.NAMES = FALSE, SIMPLIFY = FALSE)
        df2[bads, ] <- dfx
      }

      # some undetermined is.null(df2$ms1_moverz) can have apex_scan_num
      df2$apex_scan_num[lengths(df2$ms1_moverz) == 0L] <- list(NULL)
      
      df[rows2, ] <- df2
      
      df
    })
  }
  else {
    df$apex_scan_num <- df$orig_scan
  }
  
  df <- reloc_col_after(df, "apex_scan_num", "orig_scan")
  
  ## cleans up
  rows <- df$ms_level == 1L
  df1 <- df[rows, ]
  df1$ms_level <- df1$msx_charges <- df1$apex_scan_num <- 
    df1$rptr_moverzs <- df1$rptr_ints <- NULL
  # for adding back apex_rts
  qs::qsave(df1, file.path(mgf_path, paste0("ms1_", filename)), preset = "fast")
  rm(list = "df1")

  df <- df[!rows, ]
  bads <- unlist(lapply(df$ms1_moverz, is.null))
  df <- df[!bads, ]
  
  ans_rts <- reset_rettimes(
    ret_times = df$ret_time, 
    min_ret_time = min_ret_time, 
    max_ret_time = max_ret_time)
  
  min_ret_time <- ans_rts$min_ret_time
  max_ret_time <- ans_rts$max_ret_time
  
  df <- dplyr::filter(df, ret_time >= min_ret_time, ret_time <= max_ret_time)
  
  if (maxn_mdda_precurs >= 0L) {
    bads <- lapply(df$ms1_mass, function (x) length(x) == 1L && is.na(x))
    bads <- unlist(bads)
    
    if (any(bads)) {
      df   <- df[!bads, ]
    }

    oks <- lapply(df$ms1_mass, function (x) 
      .Internal(which(x >= min_mass & x <= max_mass)))
    # oks <- lapply(df$ms1_mass, function (x) .Internal(which(x >= 230 & x <= max_mass)))
    
    # better check all MS1 columns with "list" properties...
    # also check all equal: lengths(df$ms1_mass), lengths(apex_scan_num)...
    df$ms1_mass      <- mapply(function (x, y) x[y], df$ms1_mass, oks, 
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_moverz    <- mapply(function (x, y) x[y], df$ms1_moverz, oks, 
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_charge    <- mapply(function (x, y) x[y], df$ms1_charge, oks, 
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_int       <- mapply(function (x, y) x[y], df$ms1_int, oks, 
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$apex_scan_num <- mapply(function (x, y) x[y], df$apex_scan_num, oks, 
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    # ms1_mass may again contain numeric(0) following the above filtration
    df <- df[lengths(df$ms1_mass) > 0L, ]
  }
  # may be deleted; ms1_charge, ms1_mass are list not scalar anymore
  else {
    df <- df[with(df, !is.na(ms1_mass)), ]
    df <- dplyr::filter(
      df, 
      ms1_charge >= min_ms1_charge, ms1_charge <= max_ms1_charge, 
      ms1_mass >= min_mass, ms1_mass <= max_mass)
    df <- dplyr::arrange(df, ms1_mass)
  }
  
  df <- dplyr::rename(
    df, 
    ms2_moverzs = msx_moverzs, 
    ms2_ints = msx_ints, 
    ms2_charges = msx_charges, 
    ms2_n = msx_n)
}




#' Obtain the space of monoisotopic MS1
#'
#' Followed by de-isotoping and temporarily puts MS1-XYZ of MS2 scans to the
#' preceding MS1 scans.
#'
#' Note that the Y-values were averaged from franking MS1 scans, not the real
#' precursor intensities for for ascribing X-values.
#'
#' @param df A data frame of MS1 and MS2.
get_ms1xs_space <- function (df)
{
  pos_levs <- getMSrowIndexes(df$ms_level)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  
  ms1_moverzs <- df[["ms1_moverz"]]
  ms1_ints    <- df[["ms1_int"]]
  ms1_charges <- df[["ms1_charge"]]
  ms1_masses  <- df[["ms1_mass"]]
  
  outx <- outy <- outz <- outm <- vector("list", length(ms1_stas))
  
  for (i in seq_along(ms1_stas)) {
    rng1 <- ms1_stas[[i]]
    rng2 <- ms2_stas[[i]]:ms2_ends[[i]]
    
    # isolation windows can have overlaps -> 
    #   the same precursor at multiple windows -> duplicated MS1 entries
    xs2 <- .Internal(unlist(ms1_moverzs[rng2], recursive = FALSE, use.names = FALSE))
    ys2 <- .Internal(unlist(ms1_ints[rng2],    recursive = FALSE, use.names = FALSE))
    zs2 <- .Internal(unlist(ms1_charges[rng2], recursive = FALSE, use.names = FALSE))
    ms2 <- .Internal(unlist(ms1_masses[rng2],  recursive = FALSE, use.names = FALSE))
    nx2 <- length(xs2)
    
    if (!nx2) { next }

    if (nx2 > 1L) {
      ord <- 
        .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, xs2))
      xs2 <- xs2[ord]
      ys2 <- ys2[ord]
      zs2 <- zs2[ord]
      ms2 <- ms2[ord]
      
      oks <- !duplicated(xs2)
      xs2 <- xs2[oks]
      ys2 <- ys2[oks]
      zs2 <- zs2[oks]
      ms2 <- ms2[oks]
    }
    
    outx[[i]] <- xs2
    outy[[i]] <- ys2
    outz[[i]] <- zs2
    outm[[i]] <- ms2
  }
  
  df[["ms1_moverz"]][ms1_stas] <- outx
  df[["ms1_int"]][ms1_stas]    <- outy
  df[["ms1_charge"]][ms1_stas] <- outz
  df[["ms1_mass"]][ms1_stas]   <- outm
  
  df
}


#' Estimates the number of MS1 scans at a 3-min range of retention times
#' 
#' @param rts A vector of retention times.
#' @param d The allowance of retention time distance in seconds.
#' @return The number of MS1 scans in a 3-min time frame.
estimate_rtgap <- function (rts, d = 180)
{
  rts <- sort(rts)
  len <- length(rts)
  mid <- max(as.integer(len / 2), 1L)
  rtm <- rts[[mid]]
  grs <- which(rts > rtm + d)
  upr <- if (length(grs)) grs[[1]] else len
  
  # rts[[upr]] - rtm
  
  max(upr - mid, 1L)
}


#' Helper of \link{deisoDDA}.
#' 
#' @param debug Logical; in a debug mode or not.
#' @param ppm_ms1_bin The three-bin MS1 ppm.
#' @param ppm_ms2_bin The three-bin MS2 ppm.
#' @inheritParams deisoDDA
predeisoDDA <- function (filename = NULL, temp_dir = NULL, mzml_type = "raw", 
                         ppm_ms1_bin = 10L, ppm_ms2_bin = 10L, 
                         maxn_mdda_precurs = 3L, topn_ms2ions = 150L, 
                         n_mdda_flanks = 6L, n_dia_scans = 6L, 
                         min_mass = 200L, max_mass = 4500L,
                         min_ms2mass = 115L, max_ms2mass = 4500L, 
                         min_ms1_charge = 2L, max_ms1_charge = 4L, 
                         min_ret_time = 0, max_ret_time = Inf, 
                         min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                         deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                         ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                         grad_isotope = 1.6, fct_iso2 = 3.0, 
                         quant = "none", use_lfq_intensity = TRUE, 
                         tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                         exclude_reporter_region = FALSE, use_defpeaks = FALSE, 
                         is_pasef = FALSE, n_para = 1L, debug = FALSE)
{
  ## reads parsed peak lists
  # ms1_moverzs, ms1_ints and ms1_charges: NA vectors
  # raw_file: scalar
  # mobility <- ans$mobility # NULL for Thermo's data
  
  ans <- qs::qread(file.path(temp_dir, filename))
  msx_moverzs <- ans$msx_moverzs
  msx_ints <- ans$msx_ints
  msx_ns <- ans$msx_ns
  ms1_fr <- ans$ms1_fr
  
  # by MSConvert or NA by Mzion
  # 0 if undetermined by Bruker metadata
  ms1_moverzs <- ans$ms1_moverzs
  ms1_ints    <- ans$ms1_ints
  ms1_charges <- ans$ms1_charges
  
  scan_title <- ans$scan_title
  ms_level   <- ans$ms_level
  ret_time   <- ans$ret_time
  scan_num   <- ans$scan_num
  orig_scan  <- ans$orig_scan
  iso_ctr    <- ans$iso_ctr
  iso_lwr    <- ans$iso_lwr
  iso_upr    <- ans$iso_upr
  raw_file   <- ans$raw_file # scalar
  
  ## NULL with Thermo's data
  mobility    <- ans$mobility
  slice_start <- ans$slice_start
  slice_end   <- ans$slice_end
  
  if (!is.list(ms1_moverzs)) {
    ms1_moverzs <- as.list(ms1_moverzs)
  }
  if (!is.list(ms1_ints)) {
    ms1_ints <- as.list(ms1_ints)
  }
  if (!is.list(ms1_charges)) {
    ms1_charges <- as.list(ms1_charges)
  }
  
  ###
  # May confirm that each vector of msx_moverzs are in ascending order
  ###
  
  len0 <- length(msx_moverzs) # safeguard against row drops
  msx_charges <- vector("list", len0)
  
  ###
  # try later to subset early by min_ and max_ masses...
  # even at unknown z, not very likely to have z >= 2 at small m/z...
  ### 
  
  # on the lookout for row drops that can cause mismatches between reporter data
  # and others
  restmt <- extract_mgf_rptrs(
    msx_moverzs, 
    msx_ints, 
    quant = quant, 
    tmt_reporter_lower = tmt_reporter_lower, 
    tmt_reporter_upper = tmt_reporter_upper, 
    exclude_reporter_region = exclude_reporter_region)
  
  msx_moverzs  <- restmt[["xvals"]]
  msx_ints     <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints    <- restmt[["rptr_ints"]]
  rm(list = "restmt")
  
  if (deisotope_ms2) {
    message("Deisotoping MS2: ", filename)
    oks2 <- ms_level == 2L
    
    # Bad samples without MS2
    if (!any(oks2)) {
      warning("No MS2 data: ", filename)
      return(NULL)
    }
    
    ans2 <- deisoDDAMS2(
      ms2_moverzs = msx_moverzs[oks2], ms2_ints = msx_ints[oks2], 
      topn_ms2ions = topn_ms2ions, ppm_ms2_deisotope = ppm_ms2_deisotope, 
      max_ms2_charge = max_ms2_charge, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, quant = quant, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region, 
      n_fwd = 10L, n_bwd = 10L, offset_upr = 30L, offset_lwr = 30L, 
      n_para = n_para)
    msx_moverzs[oks2] <- ans2[["msx_moverzs"]]
    msx_ints[oks2]    <- ans2[["msx_ints"]]
    msx_charges[oks2] <- ans2[["msx_charges"]]
    rm(list = c("ans2", "oks2"))
    
    message("Completed MS2 deisotoping: ", filename)
  }
  
  # MS1 X, Y and Z are NA from Mzion::readRAW and require deisotoping
  if (maxn_mdda_precurs) {
    message("Deisotoping MS1: ", filename)
    
    if (is_pasef) {
      if (n_mdda_flanks) {
        n_mdda_flanks <- 0L
        warning("Coerce to `n_mdda_flanks = 0`.")
      }
      
      if ((n_cores <- min(n_para, 32L)) > 1L) {
        grps <- find_ms2ends(ms_level, n_cores * 8L)
        
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        ans <- parallel::clusterMap(
          cl, pasefMS1xyz, 
          msx_moverzs = split(msx_moverzs, grps), msx_ints = split(msx_ints, grps), 
          ms1_fr = split(ms1_fr, grps), ms_level = split(ms_level, grps), 
          iso_ctr = split(iso_ctr, grps), iso_lwr = split(iso_lwr, grps), 
          iso_upr = split(iso_upr, grps), ms1_moverzs = split(ms1_moverzs, grps), 
          ms1_charges = split(ms1_charges, grps), ms1_ints = split(ms1_ints, grps), 
          MoreArgs = list(
            maxn_mdda_precurs = maxn_mdda_precurs, n_mdda_flanks = n_mdda_flanks,  
            topn_ms2ions = topn_ms2ions, max_ms1_charge = max_ms1_charge, 
            ppm_ms1_deisotope = ppm_ms1_deisotope, grad_isotope = grad_isotope, 
            fct_iso2 = fct_iso2, use_defpeaks = use_defpeaks, 
            min_mass = min_mass, filename = filename
          ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
        parallel::stopCluster(cl)
        
        for (nm in names(ans[[1]])) {
          assign(nm, unlist(lapply(ans, function (x) x[[nm]]), 
                            recursive = FALSE, use.names = FALSE))
        }
        
        rm(list = c("ans", "grps", "cl", "nm"))
      }
      else {
        ans <- pasefMS1xyz(
          msx_moverzs = msx_moverzs, msx_ints = msx_ints, ms1_fr = ms1_fr, 
          ms_level = ms_level, iso_ctr = iso_ctr, iso_lwr = iso_lwr, 
          iso_upr = iso_upr, ms1_moverzs = ms1_moverzs, 
          ms1_charges = ms1_charges, ms1_ints = ms1_ints, 
          maxn_mdda_precurs = maxn_mdda_precurs, n_mdda_flanks = n_mdda_flanks, 
          topn_ms2ions = topn_ms2ions, 
          max_ms1_charge = max_ms1_charge, ppm_ms1_deisotope = ppm_ms1_deisotope, 
          grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
          use_defpeaks = use_defpeaks, min_mass = min_mass, filename = filename)
        
        ms1_moverzs <- ans$ms1_moverzs
        ms1_masses <- ans$ms1_masses
        ms1_ints <- ans$ms1_ints
        ms1_charges <- ans$ms1_charges
        rm(list = "ans")
      }
      
      df <- tibble::tibble(
        msx_moverzs = msx_moverzs, 
        msx_ints = msx_ints, 
        msx_charges = msx_charges, # not in `df1`
        msx_ns = msx_ns, 
        
        ms1_fr = ms1_fr, 
        ms1_moverzs = ms1_moverzs, 
        ms1_ints = ms1_ints, 
        ms1_charges = ms1_charges, 
        ms1_masses = ms1_masses, # not in `df1`
        scan_title = scan_title, 
        ms_level = ms_level, 
        ret_time = ret_time, 
        scan_num = scan_num, 
        orig_scan = orig_scan, 
        iso_ctr = iso_ctr, 
        iso_lwr = iso_lwr, 
        iso_upr = iso_upr, 
        
        slice_start = slice_start,
        slice_end = slice_end, 
        mobility = mobility, )
      
      ## bring back PASEF full_ms1
      df  <- df[with(df, ms_level > 1L), ]
      
      fi_pasefms1 <- file.path(temp_dir, paste0("pasefms1_", filename))
      if (!file.exists(fi_pasefms1)) {
        stop("Full MS1 not found ", fi_pasefms1)
      }
      df1 <- qs::qread(fi_pasefms1)
      df1$ms1_moverzs <- df1$ms1_ints <- df1$ms1_charges <- df1$ms1_masses <- 
        df1$msx_charges <- vector("list", nrow(df1))
      
      nms <- names(df)
      if (length(setdiff(names(df1), nms))) {
        stop("Developer: check for inconsistency in column names")
      }
      
      df1 <- df1[, nms]
      
      if (!identical(lapply(df1, class), lapply(df, class))) {
        stop("Developer: inconsistent classes between MS1 and MS2 data.")
      }
      
      df <- dplyr::bind_rows(df1, df)
      df <- df |> dplyr::arrange(ms1_fr, ms_level)
      
      for (nm in nms) {
        assign(nm, df[[nm]])
      }
      
      rm(list = c("df", "df1", "nm", "nms", "fi_pasefms1"))
      # gc()
      message("Completed PASEF MS1 deisotoping at: ", Sys.time())
    }
    else {
      ans <- getMS1xyz(
        msx_moverzs = msx_moverzs, msx_ints = msx_ints, ms_level = ms_level, 
        iso_ctr = iso_ctr, iso_lwr = iso_lwr, iso_upr = iso_upr, 
        ms1_moverzs = ms1_moverzs, ms1_charges = ms1_charges, ms1_ints = ms1_ints, 
        maxn_mdda_precurs = maxn_mdda_precurs, n_mdda_flanks = n_mdda_flanks,  
        topn_ms2ions = topn_ms2ions, quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region, 
        max_ms1_charge = max_ms1_charge, ppm_ms1_deisotope = ppm_ms1_deisotope, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
        use_defpeaks = use_defpeaks, min_mass = min_mass, filename = filename)
      
      # at ms_level == 2L, `ms1_ints` are the sum over flanking scans,
      #   not the peak areas
      ms1_moverzs <- ans$ms1_moverzs
      ms1_masses  <- ans$ms1_masses
      ms1_ints    <- ans$ms1_ints
      ms1_charges <- ans$ms1_charges
      rm(list = "ans")
      message("Completed MS1 deisotoping: ", filename)
    }
    
    # look up MS2 for undetermined precursor charge states
    if (deisotope_ms2) {
      rows1 <- lapply(ms1_moverzs, is.null) # all MS1 and some MS2
      rows1 <- .Internal(unlist(rows1, recursive = FALSE, use.names = FALSE))
      rows2 <- ms_level == 2L
      rows  <- .Internal(which(rows2 & rows1)) # MS2 without precursor info
      nrows <- length(rows)
      
      if (nrows) {
        # msx_moverzs[rows] and msx_ints[rows] not NULL; msx_charges[rows]: NULL
        # msx_charges[rows] <- NA_real_
        ans2 <- mapply(
          find_ms1byms2, 
          moverzs = msx_moverzs[rows], msxints = msx_ints[rows], 
          charges = msx_charges[rows], center = iso_ctr[rows], 
          iso_lwr = iso_lwr[rows], iso_upr = iso_upr[rows], 
          SIMPLIFY = FALSE, USE.NAMES = FALSE)
        
        if (length(ans2) != length(rows)) {
          stop("Develeoper: checks for row drops.")
        }
        
        ans2 <- dplyr::bind_rows(ans2)
        ans2$ms1_moverzs <- as.list(ans2$ms1_moverzs)
        ans2$ms1_masses  <- as.list(ans2$ms1_masses)
        ans2$ms1_ints    <- as.list(ans2$ms1_ints)
        ans2$ms1_charges <- as.list(ans2$ms1_charges)
        
        # outputs contain NA and revert them back to list(NULL)
        nas <- which(is.na(ans2$ms1_moverzs))
        
        if (length(nas)) {
          ans2$ms1_ints[nas] <- ans2$ms1_charges[nas] <- 
            ans2$ms1_masses[nas] <- ans2$ms1_moverzs[nas] <- list(NULL)
        }
        
        if (length(ms1_moverzs[rows]) != length(ans2$ms1_moverzs)) {
          stop("Developer: check for entries dropping.")
        }
        
        ms1_moverzs[rows] <- ans2$ms1_moverzs
        ms1_masses[rows]  <- ans2$ms1_masses
        ms1_ints[rows]    <- ans2$ms1_ints
        ms1_charges[rows] <- ans2$ms1_charges
        rm(list = c("ans2", "nas"))
      }
      
      message("Completed auxillary MS1 deisotoping: ", filename)
      rm(list = c("rows1", "rows2", "rows"))
    }
    
    ###
    # Up to this point, ms1_moverzs are list(NULL) at ms_level == 1L;
    # Not to trace all MS1 features but only the monoisotopic 
    #  that have been assigned to MS2 scans.
    ###
  }
  else {
    # only for benchmarking default pwiz-mzML deisotoping; no LFQ-MS1 tracing
    ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                         ms1_moverzs, ms1_charges, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  # final check
  if (!is_pasef) {
    if (length(msx_moverzs) != len0) {
      stop("Developer: check for row dropping.")
    }
  }
  
  # msx_moverzs at ms_level == 1L correspond to full-spectrum ms1_moverzs
  # msx_charges at ms_level == 1L are list(NULL)
  df <- tibble::tibble(
    scan_title = scan_title,
    raw_file = raw_file,
    ms_level = ms_level, 
    ms1_moverz = ms1_moverzs, 
    ms1_mass = ms1_masses,
    ms1_int = ms1_ints, 
    ms1_charge = ms1_charges, 
    ret_time = ret_time, 
    scan_num = scan_num, 
    orig_scan = orig_scan,
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_charges = msx_charges, 
    msx_n = msx_ns, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints)
}


#' Deisotopes DDA-MS2
#' 
#' @param ms2_moverzs MS2 m-over-z values.
#' @param ms2_ints MS2 intensities.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries. The default is
#'   20 for MS1 and 10 for MS2.
#' @param n_bwd Backward looking up to \code{n_bwd} mass entries.
#' @param offset_upr A cardinal number of upper mass off-sets.
#' @param offset_lwr A cardinal number of lower mass off-sets.
#' @inheritParams deisoDDA
deisoDDAMS2 <- function (ms2_moverzs, ms2_ints, topn_ms2ions = 150L, 
                         ppm_ms2_deisotope = 8L, max_ms2_charge = 3L, 
                         n_fwd = 10L, n_bwd = 10L, offset_upr = 30L, 
                         offset_lwr = 30L, grad_isotope = 1.6, fct_iso2 = 3.0, 
                         quant = "none", 
                         tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                         exclude_reporter_region = FALSE, n_para = 1L)
  
{
  # for ensure the same number of entries before and after the process
  lenx <- length(ms2_moverzs)
  
  if (n_para > 1L) {
    n_cores <- n_para
    n_chunks <- n_cores * 4L
    grps <- sep_vec(ms2_moverzs, n_chunks)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out <- parallel::clusterMap(
      cl, getMS2xyz, 
      msx_moverzs = split(ms2_moverzs, grps), 
      msx_ints = split(ms2_ints, grps), 
      MoreArgs = list(
        topn_ms2ions = topn_ms2ions, 
        max_ms2_charge = max_ms2_charge, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, 
        quant = quant, 
        n_fwd = n_fwd, 
        n_bwd = n_bwd, 
        offset_upr = offset_upr, 
        offset_lwr = offset_lwr, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region
      ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
    parallel::stopCluster(cl)
    
    outx <- unlist(lapply(out, function (x) x[[1]]), recursive = FALSE, 
                   use.names = FALSE)
    outy <- unlist(lapply(out, function (x) x[[2]]), recursive = FALSE, 
                   use.names = FALSE)
    outz <- unlist(lapply(out, function (x) x[[3]]), recursive = FALSE, 
                   use.names = FALSE)
    
    if (length(outx) != lenx) {
      stop("Developer: check for entries dropping.")
    }
    
    out <- list(msx_moverzs = outx, msx_ints = outy, msx_charges = outz)
  }
  else {
    out <- getMS2xyz(
      msx_moverzs = ms2_moverzs, 
      msx_ints = ms2_ints, 
      topn_ms2ions = topn_ms2ions, 
      max_ms2_charge = max_ms2_charge, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2, 
      quant = quant, 
      n_fwd = n_fwd, 
      n_bwd = n_bwd, 
      offset_upr = offset_upr, 
      offset_lwr = offset_lwr, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region)
    
    if (length(out[[1]]) != lenx) {
      stop("Developer: check for entries dropping.")
    }
    
    names(out) <- c("msx_moverzs", "msx_ints", "msx_charges")
  }
  
  # no need to check rptr_moverzs, rptr_ints alignment since equal length
  
  out
}


#' De-isotoping DDA-MS1
#' 
#' @param msx_moverzs Full-spectrum MS1 or MS2 m-over-z values.
#' @param msx_ints Full-spectrum MS1 or MS2 intensities.
#' @param ms1_fr MS1 frame numbers.
#' @param ms_level Vectors of MS levels.
#' @param iso_ctr A vector of isolation centers.
#' @param iso_lwr A vector of isolation lowers.
#' @param iso_upr A vector of isolation uppers.
#' @param ms1_moverzs MS1 moverz values.
#' @param ms1_ints MS1 intensity values.
#' @param ms1_charges MS1 charge states.
#' @param filename A filename for logging.
#' @inheritParams matchMS
pasefMS1xyz <- function (msx_moverzs, msx_ints, ms1_fr, ms_level, 
                         iso_ctr = NULL, iso_lwr = NULL, iso_upr = NULL, 
                         # slice_start = NULL, slice_end = NULL, mobility = NULL, 
                         ms1_moverzs = NULL, ms1_charges = NULL, ms1_ints = NULL, 
                         maxn_mdda_precurs = 1L, n_mdda_flanks = 0L, 
                         topn_ms2ions = 150L, max_ms1_charge = 4L, 
                         ppm_ms1_deisotope = 8L, grad_isotope = 1.6, 
                         fct_iso2 = 3.0, use_defpeaks = FALSE, min_mass = 115L, 
                         filename = NULL)
{
  if (debug <- FALSE) {
    df0 <- tibble::tibble(
      msx_moverzs = msx_moverzs, msx_ints = msx_ints, 
      x0 = ms1_moverzs, y0 = ms1_ints, z0 = ms1_charges,
      ms1_fr = ms1_fr, ms_level = ms_level, 
      iso_ctr = iso_ctr, iso_lwr = iso_lwr, iso_upr = iso_upr, )
  }
  
  if (!use_defpeaks) {
    ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", length(msx_ints))
  }
  
  df <- tibble::tibble(
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    ms_level = ms_level, 
    iso_ctr = iso_ctr, 
    iso_lwr = iso_lwr, 
    iso_upr = iso_upr, 
    ms1_moverzs = ms1_moverzs,
    ms1_ints = ms1_ints,
    ms1_charges = ms1_charges, )
  
  dfs <- split(df, ms1_fr)
  frs <- unique(ms1_fr) # identical(as.integer(names(dfs)), frs)
  rm(list = "df")
  
  # for `n_mdda_flanks > 0L`, but mobilities seem to vary from frame to frame
  #   may still implement the integration by `mobility`, not by `slice_start`
  # 
  # i = 3000L
  # cols1 <- c("msx_moverzs", "msx_ints", "slice_start", "mobility")
  # rng1 <- max(1, i - n_mdda_flanks):min(frmax, i + n_mdda_flanks)
  # df1s <- dfs[rng1]
  # df1s <- lapply(df1s, function (x) x[x$ms_level == 1L, cols1])
  # 
  # compare mobility to flanking MS1 frames...
  
  cols <- c("ms1_moverzs", "ms1_ints", "ms1_charges")
  
  for (i in seq_along(frs)) {
    dfi  <- dfs[[i]]
    xsi  <- dfi$msx_moverzs
    ysi  <- dfi$msx_ints
    ctri <- dfi$iso_ctr
    lwri <- dfi$iso_lwr
    upri <- dfi$iso_upr
    
    pos_levs <- getMSrowIndexes(dfi$ms_level) # MS1 without following MS2s dropped
    ms1_stas <- pos_levs$ms1_stas
    ms2_stas <- pos_levs$ms2_stas
    ms2_ends <- pos_levs$ms2_ends
    
    for (j in seq_along(ms1_stas)) {
      sta1 <- ms1_stas[j]
      rng2 <- ms2_stas[[j]]:ms2_ends[[j]] # must have length
      
      ans <- find_mdda_mms1s(
        msx_moverzs = xsi[[sta1]], 
        msx_ints = ysi[[sta1]], 
        iso_ctr = ctri[rng2], iso_lwr = lwri[rng2], iso_upr = upri[rng2], 
        ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
        max_ms1_charge = max_ms1_charge, n_fwd = 15L, n_bwd = 20L, 
        offset_upr = 20L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
        use_defpeaks = use_defpeaks, min_mass = min_mass, margin = 0, 
        width_left = 0.26, width_right = 0.25, 
        is_pasef = TRUE, min_y = 10, ms1_min_ratio = 0.05)
      xs <- ans[["x"]]
      ys <- ans[["y"]] # as.integer; don't ys is list
      zs <- ans[["z"]]
      
      if (length(xs) != length(rng2)) { stop("Check for entry-dropping.") }
      
      if (debug) {
        dfx <- dplyr::bind_cols(tibble::tibble(x = xs, y = ys, z = zs), df0[[i]]) |>
          dplyr::select(-c("msx_moverzs", "msx_ints", "ms1_fr", "ms_level"))
      }
      
      # updates corresponding MS1 X, Y and Z for each MS2
      for (k in seq_along(rng2)) {
        xsk <- xs[[k]]
        if (length(xsk)) {
          rk <- rng2[[k]]
          dfi$ms1_moverzs[[rk]] <- xsk
          dfi$ms1_ints[[rk]] <- ys[[k]]
          dfi$ms1_charges[[rk]] <- zs[[k]]
        }
      }
    }
    
    dfs[[i]][cols] <- dfi[cols]
  }
  
  df <- dplyr::bind_rows(dfs)
  ms1_moverzs <- df$ms1_moverzs
  ms1_ints    <- df$ms1_ints
  ms1_charges <- df$ms1_charges
  ms1_masses  <- mapply(
    function (x, y) (x - 1.00727647) * y, 
    ms1_moverzs, ms1_charges, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  list(ms1_moverzs = ms1_moverzs, ms1_masses = ms1_masses, 
       ms1_ints = ms1_ints, ms1_charges = ms1_charges)
}


#' Obtains the indexes of MS1 and MS2 starts and ends.
#'
#' @param ms_level A vector of MS levels.
#' @param pad_nas Logical; if TRUE, pad NA values for \code{ms1_stas} without
#'   following \code{ms2_stas} and \code{ms2_ends}. Otherwise, the
#'   \code{ms1_stas} will drop.
getMSrowIndexes <- function (ms_level, pad_nas = FALSE)
{
  # add a trailing `1`
  len <- length(ms_level)
  
  if (ms_level[len] == 2L) {
    ms_level <- c(ms_level, 1L)
    len <- len + 1L
    add_trailing <- TRUE
  }
  else {
    add_trailing <- FALSE
  }
  
  idxes_ms1 <- which(ms_level == 1L)
  diff_ms1 <- c(0L, diff(idxes_ms1))
  oks <- which(diff_ms1 > 1L) # non-consecutive MS1s
  ms1_stas <- idxes_ms1[oks - 1L] # excludes the last non-consecutive MS1 index
  ms2_stas <- ms1_stas + 1L
  ms2_ends <- idxes_ms1[oks] - 1L
  
  if (pad_nas) {
    if (add_trailing) {
      idxes_ms1 <- idxes_ms1[-length(idxes_ms1)]
    }
    
    bads <- idxes_ms1[!idxes_ms1 %in% ms1_stas]
    nas <- rep_len(NA_integer_, length(bads))
    ms1_stas <- c(bads, ms1_stas)
    ms2_stas <- c(nas, ms2_stas)
    ms2_ends <- c(nas, ms2_ends)
    
    ord <- order(ms1_stas)
    ms1_stas <- ms1_stas[ord]
    ms2_stas <- ms2_stas[ord]
    ms2_ends <- ms2_ends[ord]
  }

  list(ms1_stas = ms1_stas, ms2_stas = ms2_stas, ms2_ends = ms2_ends)
}


#' Finds group indexes with end points at \code{vals = 2L}.
#' 
#' For chunk-splitting with each chuck ending at value 2.
#' 
#' @param vals An integer vector of (MS levels of) 1 or 2.
#' @param n_chunks The number of chunks.
#' 
#' @examples
#' vals1 <- c(rep(1, 2), rep(2, 5), 1, rep(2, 4), rep(1, 3), rep(2, 4), 1)
#' grps1 <- mzion:::find_ms2ends(vals1, 3)
#' 
#' vals2 <- c(rep(1, 2), rep(2, 5), 1, rep(2, 4), rep(1, 3), rep(2, 4))
#' grps2 <- mzion:::find_ms2ends(vals2, 3)
find_ms2ends <- function (vals, n_chunks = 3L)
{
  if (n_chunks <= 1L) {
    return(rep_len(1L, length(vals)))
  }

  ms2_ends <- getMSrowIndexes(vals)$ms2_ends
  brs <- floor(length(ms2_ends)/n_chunks) * 1:(n_chunks - 1L)
  findInterval(seq_along(vals), ms2_ends[brs] + 1L)
}


#' Deisotoping DDA-MS1.
#'
#' Inputs are full-spectrum X and Y values. MS levels are differentiated by
#' \code{ms_level}.
#' 
#' @param ms_level Vectors of MS levels
#' @param iso_ctr A vector of isolation centers.
#' @param iso_lwr A vector of isolation lowers.
#' @param iso_lwr A vector of isolation uppers.
#' @param ms1_moverzs MS1 m-over-z values.
#' @param ms1_ints MS1 intensity values.
#' @param ms1_charges MS1 charge states.
#' @param filename A file name (for troubleshooting).
#' @param debug Logical; debug mode or not
#' @inheritParams find_mdda_mms1s
#' @inheritParams matchMS
getMS1xyz <- function (msx_moverzs = NULL, msx_ints = NULL, ms_level = NULL, 
                       iso_ctr = NULL, iso_lwr = NULL, iso_upr = NULL, 
                       ms1_moverzs = NULL, ms1_charges = NULL, ms1_ints = NULL, 
                       maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                       topn_ms2ions = 150L, quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       max_ms1_charge = 4L, ppm_ms1_deisotope = 8L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, 
                       use_defpeaks = FALSE, min_mass = 200L, filename = NULL, 
                       debug = FALSE)
{
  if (debug) {
    fun <- "getMS1xyz"
  }
  else {
    fun <- as.character(match.call()[[1]]) # error if calling from do.call...
    fun <- fun[length(fun)]
  }

  ## Low priority: no data filtration by scan_nums; 
  #  plus, scan_nums with PASEF were altered in `proc_pasefs`
  #  better filter data by retention times
  
  len0 <- length(msx_moverzs) # to guard against row drops
  if (!len0) {
    return(NULL)
  }

  if (!use_defpeaks) {
    ms1_moverzs <- ms1_charges <- ms1_ints <- vector("list", len0)
  }
  
  pos_levs <- getMSrowIndexes(ms_level)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  len1     <- length(ms1_stas)

  # go from z = min_ms1_charge:max_ms1_charge first,  
  # then if (max_ms1_charge < 6) max_ms1_charge:6

  # obtain precursor XYZ values from multiple adjacent MS1 scans
  for (i in 1:len1) {
    rng1 <- ms1_stas[max(1L, i - n_mdda_flanks):min(len1, i + n_mdda_flanks)]
    rng2 <- ms2_stas[[i]]:ms2_ends[[i]]

    # ys2: summed over flanking MS1 scans, not the precursor peak area
    ans2 <- find_mdda_mms1s(
      msx_moverzs = msx_moverzs[rng1], msx_ints = msx_ints[rng1], 
      iso_ctr = iso_ctr[rng2], iso_lwr = iso_lwr[rng2], iso_upr = iso_upr[rng2],
      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
      max_ms1_charge = max_ms1_charge, n_fwd = 15L, n_bwd = 20L, 
      offset_upr = 20L, offset_lwr = 30L, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks, min_mass = min_mass)
    xs2 <- ans2[["x"]]
    ys2 <- ans2[["y"]]
    zs2 <- ans2[["z"]]

    if (length(xs2) != length(rng2)) {
      stop("Check for entries drop in ", fun, " for ", filename) 
    }
    
    # updates corresponding MS1 X, Y and Z for each MS2
    if (length(oks <- .Internal(which(lengths(xs2) > 0L)))) {
      ms1_moverzs[rng2][oks] <- xs2[oks]
      ms1_charges[rng2][oks] <- zs2[oks]
      ms1_ints[rng2][oks]    <- ys2[oks]
    }
  }

  ms1_masses <- mapply(
    function (x, y) (x - 1.00727647) * y, ms1_moverzs, ms1_charges, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  if (length(ms1_moverzs) != len0) {
    stop("Developer: check for entries dropping.")
  }

  list(ms1_moverzs = ms1_moverzs, ms1_masses = ms1_masses, ms1_ints = ms1_ints, 
       ms1_charges = ms1_charges)
}


#' De-isotoping DDA-MS2.
#' 
#' @param msx_moverzs Lists of MS2-X values.
#' @param msx_ints Lists of MS2-Y values.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries. The default is
#'   20 for MS1 and 10 for MS2.
#' @param n_bwd Backward looking up to \code{n_bwd} mass entries.
#' @param offset_upr A cardinal number of upper mass off-sets.
#' @param offset_lwr A cardinal number of lower mass off-sets.
#' @inheritParams matchMS
getMS2xyz <- function (msx_moverzs = NULL, msx_ints = NULL, topn_ms2ions = 150L, 
                       quant = "none", n_fwd = 10L, n_bwd = 10L, 
                       offset_upr = 20L, offset_lwr = 30L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       max_ms2_charge = 3L, ppm_ms2_deisotope = 10L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0)
{
  len <- length(msx_moverzs)
  msx_charges <- vector("list", len)
  is_tmt <- if (isTRUE(grepl("^tmt.*\\d+", quant))) TRUE else FALSE
  
  for (i in 1:len) {
    xi <- msx_moverzs[[i]]
    yi <- msx_ints[[i]]
    ni <- rep_len(1L, length(xi))
    
    mic <- find_ms1stat(
      moverzs = xi, msxints = yi, n_ms1s = ni, center = 0, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      # no deisotoping for TMT-MS2 reporter ions
      exclude_reporter_region = is_tmt, 
      ppm = ppm_ms2_deisotope, ms_lev = 2L, 
      maxn_feats = topn_ms2ions, 
      max_charge = max_ms2_charge,
      # smaller values for DDA-MS2
      n_fwd = n_fwd, n_bwd = n_bwd, offset_upr = offset_upr, 
      offset_lwr = offset_lwr, grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2)
    
    if (length(mic[["masses"]])) {
      msx_moverzs[[i]] <- mic[["masses"]]
      msx_ints[[i]] <- mic[["intensities"]]
      msx_charges[[i]] <- mic[["charges"]]
    }
  }
  
  list(msx_moverzs = msx_moverzs, msx_ints = msx_ints, msx_charges = msx_charges)
}


#' Extracts DIA mzML.
#' 
#' @param is_demux Is demux DIA or not.
#' @param idx_demux Index of demux.
#' @inheritParams hdeisoDIA
#' @inheritParams extrDDA
extrDIA <- function (spec = NULL, raw_file = NULL, temp_dir = NULL, 
                     is_demux = FALSE, 
                     idx_sc = 5L, idx_osc = 3L, 
                     idx_mslev = 2L, idx_title = 10L, idx_scanList_1 = 11L, 
                     idx_scanList_2 = 11L, idx_rt_1 = 2L, idx_rt_2 = 2L, 
                     idx_scan_start_1 = 1L, idx_scan_start_2 = 1L, 
                     idx_precursor_2 = 12L, idx_isolationWindow = 1L, 
                     idx_ctrmz = 1L, idx_lwrmz = 2L, idx_uprmz = 3L, 
                     idx_selectedIonList = 2L, idx_demux = 4L, 
                     idx_bin_1 = 12L, idx_bin_2 = 13L, 
                     idx_scan_lwr_1 = 8L, idx_scan_upr_1 = 9L, 
                     idx_scan_lwr_2 = 8L, idx_scan_upr_2 = 9L)
{
  if (!(len <- length(spec)))
    return(NULL)
  
  ret_times <- orig_scans <- scan_nums <- scan_titles <- 
    iso_ctr <- iso_lwr <- iso_upr <- character(len)
  ms_levs <- msx_ns <- integer(len)
  msx_moverzs <- msx_ints <- vector("list", len)
  demux <- if (is_demux) rep_len("0", len) else NULL
  
  for (i in 1:len) {
    x <- spec[[i]]
    ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))[[1]]
    ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))
    scan_nums[[i]] <- ids[[idx_sc]][[2]]
    orig_scans[[i]] <- ids[[idx_osc]][[2]]
    
    xc <- xml2::xml_children(x)
    ms_levs[[i]] <- ms_lev <- as.integer(xml2::xml_attr(xc[[idx_mslev]], "value"))
    scan_titles[[i]] <- xml2::xml_attr(xc[[idx_title]], "value")
    
    scanList <- xml2::xml_children(xc[[idx_scanList_1]])
    scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
    ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_1]], "value")
    
    if (ms_lev == 2L) {
      # entire MS2 is empty
      if (length(xc) < idx_precursor_2) 
        next
      
      if (is_demux) 
        demux[[i]] <- ids[[idx_demux]][[2]]
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # (assume always one precursor?)
      precursorc <- xml2::xml_children(precursor)
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_ctr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_ctrmz]], "value")
      iso_lwr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value")
      iso_upr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value")
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_2]]))
    }
    else if (ms_lev == 1L) {
      # iso_lwr[[i]] <- as.numeric(xml2::xml_attr(xc[[idx_scan_lwr_1]], "value"))
      # iso_upr[[i]] <- as.numeric(xml2::xml_attr(xc[[idx_scan_upr_1]], "value"))
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_1]]))
    }
    
    msData <- xml2::xml_contents(binData)
    
    if (length(msData) == 2L) {
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_moverzs[[i]] <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ints[[i]] <- readBin(r2, "double", n = msx_n, size = 8L)
    }
  }
  
  # mzML: ret_times in minutes; MGF: in seconds
  ret_times <- as.numeric(ret_times) * 60
  scan_nums <- as.integer(scan_nums)
  iso_ctr <- as.numeric(iso_ctr)
  iso_lwr <- as.numeric(iso_lwr) 
  iso_upr <- as.numeric(iso_upr) 
  demux <- as.integer(demux)
  
  out <- list(
    # either MS1 or MS2
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    scan_title = scan_titles,
    raw_file = raw_file, # single
    ms_level = ms_levs, # NA for MS1
    ret_time = ret_times, 
    scan_num = scan_nums, 
    orig_scan = orig_scans,
    iso_ctr = iso_ctr, 
    iso_lwr = iso_lwr, 
    iso_upr = iso_upr, 
    demux = demux
  )

  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Helper of \link{deisoDIA}.
#' 
#' @param filename A peaklist filename.
#' @param temp_dir A temp_dir to the filename.
#' @param n_para The allowance of parallel processing.
#' @inheritParams extrDDA
#' @inheritParams matchMS
hdeisoDIA <- function (filename = NULL, temp_dir = NULL, 
                       min_mass = 200L, max_mass = 4500L,
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       maxn_dia_precurs = 500L, topn_dia_ms2ions = 5000L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, 
                       quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, n_para = 1L)
{
  ans <- qs::qread(file.path(temp_dir, filename))
  msx_ns <- ans$msx_ns
  scan_title <- ans$scan_title
  raw_file <- ans$raw_file # scalar
  ms_level <- ans$ms_level
  ret_time <- ans$ret_time 
  scan_num <- ans$scan_num
  orig_scan <- ans$orig_scan
  iso_ctr <- ans$iso_ctr
  iso_lwr <- ans$iso_lwr
  iso_upr <- ans$iso_upr
  demux <- ans$demux
  msx_moverzs <- ans$msx_moverzs
  msx_ints <- ans$msx_ints
  rm(list = "ans")
  
  # need first disotope to learn charge states before data cleaning by 
  #  min_ms2mass and max_ms2mass
  if (n_para <= 1L) {
    out <- deisoDIA(
      msx_moverzs, msx_ints, ms_level, 
      min_mass = min_mass, max_mass = max_mass,
      min_ms2mass = min_ms2mass, max_ms2mass = max_ms2mass, 
      maxn_dia_precurs = maxn_dia_precurs, 
      topn_dia_ms2ions = topn_dia_ms2ions, 
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge, 
      ppm_ms1_deisotope = ppm_ms1_deisotope, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      deisotope_ms2 = deisotope_ms2, 
      max_ms2_charge = max_ms2_charge, 
      grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2, 
      quant = quant, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region)
    
    msx_moverzs <- out[[1]]
    msx_ints <- out[[2]]
    msx_charges <- out[[3]]
    msx_masses <- out[[4]] # new
    rm(list = "out")
  }
  else {
    grps <- sep_vec(msx_moverzs, n_para * 4L)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_para))
    out <- parallel::clusterMap(
      cl, deisoDIA, 
      split(msx_moverzs, grps), split(msx_ints, grps), split(ms_level, grps), 
      MoreArgs = list(
        min_mass = min_mass, 
        max_mass = max_mass,
        min_ms2mass = min_ms2mass, 
        max_ms2mass = max_ms2mass, 
        maxn_dia_precurs = maxn_dia_precurs, 
        topn_dia_ms2ions = topn_dia_ms2ions, 
        min_ms1_charge = min_ms1_charge, 
        max_ms1_charge = max_ms1_charge, 
        ppm_ms1_deisotope = ppm_ms1_deisotope, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        deisotope_ms2 = deisotope_ms2, 
        max_ms2_charge = max_ms2_charge, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, 
        quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region
      ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
    parallel::stopCluster(cl)
    
    msx_moverzs <- lapply(out, `[[`, 1L)
    msx_ints <- lapply(out, `[[`, 2L)
    msx_charges <- lapply(out, `[[`, 3L)
    msx_masses <<- lapply(out, `[[`, 4L) # new
    rm(list = "out")
    
    msx_moverzs <- unlist(msx_moverzs, recursive = FALSE, use.names = FALSE)
    msx_ints <- unlist(msx_ints, recursive = FALSE, use.names = FALSE)
    msx_charges <- unlist(msx_charges, recursive = FALSE, use.names = FALSE)
    msx_masses <- unlist(msx_masses, recursive = FALSE, use.names = FALSE) # new
  }
  
  ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", length(msx_moverzs))
  
  df <- tibble::tibble(
    scan_title = scan_title,
    raw_file = raw_file,
    ms_level = ms_level, 
    ret_time = ret_time, 
    scan_num = scan_num, 
    orig_scan = orig_scan,
    # full MS1 or MS2
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_charges = msx_charges, 
    msx_masses = msx_masses, # new
    msx_ns = msx_ns, 
    
    # MS1 placeholder
    ms1_moverzs = ms1_moverzs, 
    ms1_ints = ms1_ints, 
    # MS1: vector or empty vector; MS2: NULL
    ms1_charges = ms1_charges, 
    
    iso_ctr = iso_ctr, 
    iso_lwr = iso_lwr, 
    iso_upr = iso_upr, 
    demux = demux)
  
  out_name <- paste0("dia_", filename)
  qs::qsave(df, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Deisotopes DIA-MS from MSConvert peaklists.
#' 
#' @param msx_moverzs Vectors of moverzs.
#' @param msx_ints Vectors of intensities.
#' @param ms_level Vectors of MS levels.
#' @inheritParams hdeisoDIA
#' @inheritParams matchMS
deisoDIA <- function (msx_moverzs = NULL, msx_ints = NULL, ms_level = NULL, 
                      min_mass = 200L, max_mass = 4500L,
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      maxn_dia_precurs = 500L, topn_dia_ms2ions = 5000L, 
                      min_ms1_charge = 2L, max_ms1_charge = 4L, 
                      ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                      deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                      grad_isotope = 1.6, fct_iso2 = 3.0, quant = "none", 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE)
{
  if (!(len <- length(msx_moverzs)))
    return(NULL)
  
  msx_masses <- msx_charges <- vector("list", len)
  is_tmt <- if (isTRUE(grepl("^tmt.*\\d+", quant))) TRUE else FALSE
  # if (is_tmt) stop("TMT not yet supported with DIA workflows.")
  
  ### need to collapse MS1 or DIA-MS2
  # ans <- collMS(xs = msx_moverzs, ys = msx_ints, lwr = 115L, step = 1e-5, coll = TRUE)
  ###
  

  for (i in 1:len) {
    ms_lev <- ms_level[[i]]
    
    if (ms_lev == 2L && deisotope_ms2) {
      # to collapse MS2 if is DIA...
      
      xi <- msx_moverzs[[i]]
      yi <- msx_ints[[i]]
      ni <- rep_len(1L, length(xi))
      
      mic <- find_ms1stat(
        moverzs = xi, msxints = yi, n_ms1s = ni, 
        center = 0, exclude_reporter_region = is_tmt, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        is_dda = FALSE, ppm = ppm_ms2_deisotope, 
        ms_lev = 2L, maxn_feats = topn_dia_ms2ions, 
        max_charge = max_ms2_charge, n_fwd = 10L, n_bwd = 10L, 
        offset_upr = 30L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2)
      
      xs <- mic[["masses"]]
      ys <- mic[["intensities"]]
      zs <- mic[["charges"]]
      ms <- (xs - 1.00727647) * zs
      oks <- is.na(zs) | (ms >= min_ms2mass & ms <= max_ms2mass)

      msx_moverzs[[i]] <- xs[oks]
      msx_ints[[i]] <- ys[oks]
      msx_charges[[i]] <- zs[oks]
      # no need at ms_lev == 2L
      msx_masses[[i]] <- ms[oks]
    }
    else if (ms_lev == 1L) {
      xi <- msx_moverzs[[i]]
      yi <- msx_ints[[i]]
      ni <- rep_len(1L, length(xi))
      
      mic <- find_ms1stat(
        moverzs = xi, msxints = yi, n_ms1s = ni, 
        center = 0, exclude_reporter_region = FALSE, 
        is_dda = FALSE, ppm = ppm_ms1_deisotope, 
        ms_lev = 1L, 
        # by the width of MS1: 395 - 1005???
        maxn_feats = maxn_dia_precurs, 
        max_charge = max_ms1_charge, 
        n_fwd = 20L, n_bwd = 20L, offset_upr = 30L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2)
      
      xs <- mic[["masses"]]
      ys <- mic[["intensities"]]
      zs <- mic[["charges"]] # MS2: NULL; MS1: integer vectors
      ms <- (xs - 1.00727647) * zs
      oks <- is.na(zs) | (ms >= min_mass & ms <= max_mass)

      msx_moverzs[[i]] <- xs[oks]
      msx_ints[[i]] <- ys[oks]
      msx_charges[[i]] <- zs[oks]
      msx_masses[[i]] <- ms[oks]
    }
  }
  
  # msx_charges and msx_masses can be NA for MS2
  list(msx_moverzs = msx_moverzs, msx_ints = msx_ints, 
       msx_charges = msx_charges, msx_masses = msx_masses)
}


#' Helper of \link{subDIAMS1} for parallel processing.
#' 
#' @param filename A filename.
#' @param temp_dir A temporary file folder.
#' @param n_para The allowance of parallel processing.
hsubDIAMS1 <- function (filename = NULL, temp_dir = NULL, n_para = 1L)
{
  df <- qs::qread(file.path(temp_dir, filename))
  
  if (n_para <= 1L)
    return(subDIAMS1(df))
  
  idxes <- findInterval(1:nrow(df), which(df$ms_level == 1L))
  ends2 <- find_group_breaks(idxes, n_para, by_rngs = FALSE)
  stas1 <- c(1L, ends2[1:(n_para - 1L)] + 1L)
  
  dfs <- mapply(function (x, y) df[x:y, ], stas1, ends2, 
                SIMPLIFY = FALSE, USE.NAMES = FALSE)
  rm(list = c("df", "idxes", "ends2", "stas1"))
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_para))
  ans <- parallel::clusterApply(cl, dfs, subDIAMS1)
  parallel::stopCluster(cl)
  rm(list = "dfs")
  
  df <- dplyr::bind_rows(ans)
  out_name <- paste0("sub", filename)
  qs::qsave(df, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Gets precursor mass candidates for each DIA window.
#' 
#' Subsets by MS2 isolation windows.
#' 
#' @param df A data frame.
subDIAMS1 <- function (df)
{
  nr <- nrow(df)
  idxes_ms1 <- which(df$ms_level == 1L)
  len <- length(idxes_ms1)

  # adds a non-consecutive trailing index
  idxes_ms1x <- c(idxes_ms1, idxes_ms1[len] + 2L)
  
  # Keep the last MS1 for each segment of consecutive MS1 scans
  diff_ms1 <- c(0L, diff(idxes_ms1x))
  oks <- which(diff_ms1 > 1L)
  ms1_stas <- idxes_ms1x[oks - 1L]
  ms2_stas <- ms1_stas + 1L # the last value can be arbitrary if ends with MS1
  ms2_ends <- idxes_ms1x[oks] - 1L # the last is arbitrary
  
  if (nr > idxes_ms1[len]) { # with MS2s after the last MS1
    ms2_ends[len] <- nr
  }
  else { # no MS2s after the last MS1
    len1 <- length(ms1_stas)
    ms1_stas <- ms1_stas[-len1]
    ms2_stas <- ms2_stas[-len1]
    ms2_ends <- ms2_ends[-len1]
  }
  rm(list = c("oks", "diff_ms1", "idxes_ms1", "idxes_ms1x"))
  
  for (i in seq_along(ms1_stas)) {
    idx <- ms1_stas[i]
    df1 <- df[idx, ]
    xs <- df1$msx_moverzs[[1]]
    
    if ((nx <- length(xs)) == 0L || (nx == 1L && is.na(xs)))
      next
    
    ys  <- df1$msx_ints[[1]]
    css <- df1$msx_charges[[1]]
    
    sta <- ms2_stas[i]
    end <- ms2_ends[i]
    rows2 <- sta:end
    df2 <- df[rows2, ]
    
    half  <- df2$iso_upr[[1]]
    cents <- df2$iso_ctr
    
    # the first center not the smallest: up first then down
    if ((imin <- which.min(cents)) > 1L) {
      edg_b <- c(cents[[1]] - half, cents[1:(imin-1L)] + half)
      len_b <- length(edg_b)
      # edg_b[[1]] <- edg_b[[1]] - .5
      # edg_b[[len_b]] <- edg_b[[len_b]] + .5
      
      cuts_b <- findInterval(xs, edg_b)
      xs_b <- split(xs, cuts_b)
      ys_b <- split(ys, cuts_b)
      css_b <- split(css, cuts_b)
      
      bins_b <- as.integer(names(xs_b))
      oks_b <- bins_b > 0L & bins_b < len_b
      # bins_b <- bins_b[oks_b]
      xs_b <- xs_b[oks_b]
      ys_b <- ys_b[oks_b]
      css_b <- css_b[oks_b]
      
      edg_a <- c(cents[[imin]] - half, cents[imin:length(cents)] + half)
      len_a <- length(edg_a)
      # edg_a[[1]] <- edg_a[[1]] - .5
      # edg_a[[len_a]] <- edg_a[[len_a]] + .5
      
      cuts_a <- findInterval(xs, edg_a)
      xs_a <- split(xs, cuts_a)
      ys_a <- split(ys, cuts_a)
      css_a <- split(css, cuts_a)
      
      bins_a <- as.integer(names(xs_a))
      oks_a <- bins_a > 0L & bins_a < len_a
      # bins_a <- bins_a[oks_a]
      xs_a <- xs_a[oks_a]
      ys_a <- ys_a[oks_a]
      css_a <- css_a[oks_a]
      
      bins <- c(as.integer(names(xs_b)), as.integer(names(xs_a)) + imin - 1L)
      df$ms1_moverzs[rows2][bins] <- c(xs_b, xs_a)
      df$ms1_ints[rows2][bins] <- c(ys_b, ys_a)
      df$ms1_charges[rows2][bins] <- c(css_b, css_a)
    }
    else {
      edges <- c(cents[[1]] - half, cents + half)
      len_ab <- length(edges)
      # edges[[1]] <- edges[[1]] - .5
      # edges[[len_ab]] <- edges[[len_ab]] + .5
      
      cuts <- findInterval(xs, edges)
      xs_cuts <- split(xs, cuts)
      ys_cuts <- split(ys, cuts)
      css_cuts <- split(css, cuts)
      
      # excludes precursors outside the isolation window
      bins <- as.integer(names(xs_cuts))
      oks <- bins > 0L & bins < len_ab
      bins <- bins[oks]
      df$ms1_moverzs[rows2][bins] <- xs_cuts[oks]
      df$ms1_ints[rows2][bins] <- ys_cuts[oks]
      df$ms1_charges[rows2][bins] <- css_cuts[oks]
    }
  }
  
  df
}


#' Helper of tracing DIA.
#' 
#' @param filename A filename.
#' @param raw_id A raw file id.
#' @param temp_dir A temporary directory. 
#' @param n_para The allowance of parallel processing.
#' @inheritParams load_mgfs
htraceDIA <- function (filename, raw_id, temp_dir, mgf_path, 
                       min_ret_time = 0L, max_ret_time = Inf, 
                       min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       ppm_ms1 = 8L, ppm_ms2 = 8L, 
                       n_dia_ms2bins = 1L, n_dia_scans = 4L, 
                       delayed_diams2_tracing = FALSE, topn_dia_ms2ions = 2400L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, n_para = 1L)
{
  ### cleanDIAMS
  df <- qs::qread(file.path(temp_dir, filename))
  df$msx_masses <- NULL
  
  # NULL: no MS1 in the bins
  df <- df[with(df, ms_level != 1L), ]
  bads <- unlist(lapply(df$ms1_moverzs, is.null), recursive = FALSE, 
                 use.names = FALSE)
  df <- df[!bads, ]
  rm(list = "bads")
  
  df <- dplyr::rename(df, ms2_moverzs = msx_moverzs, ms2_ints = msx_ints, 
                      ms2_charges = msx_charges, ms2_n = msx_ns, 
                      ms1_moverz = ms1_moverzs, ms1_int = ms1_ints, 
                      ms1_charge = ms1_charges)
  rts <- reset_rettimes(ret_times = df$ret_time, min_ret_time = min_ret_time, 
                        max_ret_time = max_ret_time)
  min_ret_time <- rts$min_ret_time
  max_ret_time <- rts$max_ret_time
  df <- dplyr::filter(df, ret_time >= min_ret_time, ret_time <= max_ret_time)
  rm(list = "rts")
  
  ## traces MS2
  ws <- dnorm(-n_dia_scans:n_dia_scans, mean = 0, sd = 2)
  ws <- ws/ws[n_dia_scans+1L]
  
  cols <- c("ms1_moverz", "ms1_int", "ms1_charge", 
            "ms2_moverzs", "ms2_ints", "ms2_charges")
  
  if (!all(cols %in% names(df)))
    stop("Not all required columns found for tracing MS features.")
  
  dfs <- split(df, df$iso_ctr)
  rm(list = "df")
  
  subdir <- create_dir(file.path(temp_dir, paste0("sub_", raw_id)))
  
  if (n_para <= 1L) {
    dfs <- mapply(
      traceDIA, dfs, seq_along(dfs), 
      MoreArgs = list(
        ws = ws, n_dia_scans = n_dia_scans, n_dia_ms2bins = n_dia_ms2bins, 
        temp_dir = subdir, min_mass = min_mass, min_ms2mass = min_ms2mass, 
        step1 = ppm_ms1/1e6, step2 = ppm_ms2/1e6, 
        delayed_diams2_tracing = delayed_diams2_tracing
      ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_para))
    
    dfs <- parallel::clusterMap(
      cl, traceDIA, dfs, seq_along(dfs),
      MoreArgs = list(
        ws = ws, n_dia_scans = n_dia_scans, n_dia_ms2bins = n_dia_ms2bins, 
        temp_dir = subdir, 
        min_mass = min_mass, 
        min_ms2mass = min_ms2mass, 
        step1 = ppm_ms1/1e6, step2 = ppm_ms2/1e6, 
        delayed_diams2_tracing = delayed_diams2_tracing
      ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    parallel::stopCluster(cl)
  }
  
  df <- dplyr::bind_rows(dfs)
  rm(list = "dfs")
  unlink(subdir, recursive = TRUE)
  
  ## Clean up
  lens <- lengths(df$ms1_moverz)
  df <- df[lens > 0L, ]
  # use singular names for consistency with DDA
  df$ms1_mass <- mapply(function (x, y) (x - 1.00727647) * y, 
                        df$ms1_moverz, df$ms1_charge, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df <- reloc_col_before(df, "ms1_mass", "ms1_moverz")
  
  restmt <- extract_mgf_rptrs(df$ms2_moverzs, 
                              df$ms2_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  df$ms2_moverzs <- restmt[["xvals"]]
  df$ms2_ints <- restmt[["yvals"]]
  df$rptr_moverzs <- restmt[["rptr_moverzs"]]
  df$rptr_ints <- restmt[["rptr_ints"]]
  df$iso_lwr <- df$iso_upr <- df$iso_ctr <- NULL
  
  mz_n_int <- sub_mgftopn(ms2_moverzs = df[["ms2_moverzs"]], 
                          ms2_ints = df[["ms2_ints"]], 
                          ms2_charges = df[["ms2_charges"]], 
                          
                          topn_ms2ions = topn_dia_ms2ions, 
                          
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  df[["ms2_moverzs"]] <- mz_n_int[["ms2_moverzs"]]
  df[["ms2_ints"]] <- mz_n_int[["ms2_ints"]]
  df[["ms2_charges"]] <- mz_n_int[["ms2_charges"]]
  df[["ms2_n"]] <- mz_n_int[["lens"]]
  
  # fuzzy precursor masses and the corresponding scan_nums are duplicated 
  #  In `calcpepsc`, uniq_id: scan_num + raw_file + ms1_offset
  #  In `calc_pepfdr`, uniq_id: scan_num + raw_file
  #  In `post_pepfdr`, `calc_peploc`, `hcalc_tmtint`, `add_rptrs`
  
  post_readmgf(df, raw_id = raw_id, mgf_path = mgf_path)
}


#' Traces DIA-MS2 against MS1.
#'
#' For a single MS1 isolation range, e.g., 592.5 to 604.5.
#'
#' @param df A data frame.
#' @param icenter The index of an isolation center.
#' @param ws Weights.
#' @param temp_dir A temporary directory.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param step1 A step size for MS1.
#' @param step2 A step size for MS2.
#' @param join_ms Logical; if TRUE, combine adjacent entries.
#' @param spread_ohw logical; if TRUE, spread one-hit wonders to the nearest
#'   neighbors.
#' @inheritParams load_mgfs
traceDIA <- function (df = NULL, icenter = 1L, ws = NULL, n_dia_scans = 4L, 
                      n_dia_ms2bins = 1L, min_mass = 200L, min_ms2mass = 115L, 
                      step1 = 1E-5, step2 = 1E-5, temp_dir = NULL, 
                      join_ms = FALSE, spread_ohw = FALSE, 
                      delayed_diams2_tracing = FALSE)
{
  len <- nrow(df)

  mat1 <- traceLCMS(
    xs = df[["ms1_moverz"]], 
    ys = df[["ms1_int"]], 
    zs = df[["ms1_charge"]], 
    n_dia_scans = n_dia_scans, 
    from = min_mass, 
    step = step1, 
    icenter = icenter, 
    ms_lev = 1L, 
    temp_dir = temp_dir)
  
  mat1x <- mat1[["x"]] # columns: masses; rows: scans
  mat1y <- mat1[["y"]]
  mat1z <- mat1[["z"]]
  ns1 <- mat1[["n"]]
  ps1 <- mat1[["p"]] # e.g., matx[ps[[4]], 4] to bypass which(!is.na(matx[, 4]))
  rm(list = "mat1")
  gc()
  
  # may only borrow ps1 >= 2L
  
  ans1 <- flattenMSxyz(matx = mat1x, maty = mat1y, matz = mat1z, 
                       join_ms = join_ms)
  ansx1 <- ans1[["x"]]
  ansy1 <- ans1[["y"]]
  ansz1 <- ans1[["z"]]
  rm(list = "ans1")
  
  if (spread_ohw) {
    ans1_1 <- spreadMSohw(matx = mat1x, maty = mat1y, matz = mat1z, ns = ns1, 
                          ps = ps1, gap = 1L, join_ms = join_ms)
    ansx1_1 <- ans1_1[["x"]]
    ansy1_1 <- ans1_1[["y"]]
    ansz1_1 <- ans1_1[["z"]]
    rm(list = "ans1_1")
    
    for (i in seq_along(ansx1)) {
      ansx1[[i]] <- c(ansx1[[i]], ansx1_1[[i]])
      ansy1[[i]] <- c(ansy1[[i]], ansy1_1[[i]])
      ansz1[[i]] <- c(ansz1[[i]], ansz1_1[[i]])
    }
    rm(list = c("ansx1_1", "ansy1_1", "ansz1_1"))
  }

  if (delayed_diams2_tracing) {
    df$ms1_moverz <- ansx1
    df$ms1_int <- ansy1
    df$ms1_charge <- ansz1
    
    return(df)
  }

  mat2 <- traceLCMS(
    xs = df[["ms2_moverzs"]], 
    ys = df[["ms2_ints"]], 
    zs = df[["ms2_charges"]], 
    n_dia_scans = n_dia_scans, 
    from = min_ms2mass, 
    step = step2, 
    icenter = icenter, 
    ms_lev = 2L, 
    temp_dir = temp_dir)
  
  mat2x <- mat2[["x"]]
  mat2y <- mat2[["y"]]
  mat2z <- mat2[["z"]]
  ns2 <- mat2[["n"]]
  ps2 <- mat2[["p"]]
  rm(list = "mat2")
  gc()
  
  # need weights?
  # need to remove duplicated entries? should not occur
  
  # exclude high-confidence peak apexes, only join_ms from neighbors at ps2 <= 2L...
  
  # may only borrow ps1 >= 2L
  
  ans2 <- flattenMSxyz(matx = mat2x, maty = mat2y, matz = mat2z, 
                       join_ms = join_ms)
  ansx2 <- ans2[["x"]]
  ansy2 <- ans2[["y"]]
  ansz2 <- ans2[["z"]]
  rm(list = "ans2")
  gc()
  
  if (spread_ohw) {
    ans2_1 <- spreadMSohw(matx = mat2x, maty = mat2y, matz = mat2z, ns = ns2, 
                          ps = ps2, gap = 1L, join_ms = join_ms)
    ansx2_1 <- ans2_1[["x"]]
    ansy2_1 <- ans2_1[["y"]]
    ansz2_1 <- ans2_1[["z"]]
    rm(list = "ans2_1", "mat2x", "mat2y", "mat2z")
    gc()
    
    for (i in seq_along(ansx2)) {
      ansx2[[i]] <- c(ansx2[[i]], ansx2_1[[i]])
      ansy2[[i]] <- c(ansy2[[i]], ansy2_1[[i]])
      ansz2[[i]] <- c(ansz2[[i]], ansz2_1[[i]])
    }
    rm(list = c("ansx2_1", "ansy2_1", "ansz2_1"))
  }

  df$ms2_moverzs <- ansx2
  df$ms2_ints <- ansy2
  df$ms2_charges <- ansz2
  rm(list = c("ansx2", "ansy2", "ansz2"))

  # need to order or not???
  ords <- lapply(df$ms2_moverzs, order)
  df$ms2_moverzs <- mapply(function (x, y) x[y], df$ms2_moverzs, ords,
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df$ms2_ints <- mapply(function (x, y) x[y], df$ms2_ints, ords,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df$ms2_charges <- mapply(function (x, y) x[y], df$ms2_charges, ords,
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  df
}


#' Flattens MS data.
#' 
#' Joins data across adjacent scans af \code{combine = TRUE}.
#' 
#' @param matx The data matrix of moverzs.
#' @param maty The data matrix of intensities.
#' @param matz The data matrix of charge states.
#' @param join_ms Logical; if TRUE, combine adjacent entries.
flattenMSxyz <- function (matx, maty, matz, join_ms = FALSE)
{
  n_scans <- nrow(matx)
  n_masses <- ncol(matx)
  
  zs <- ys <- xs <- vector("list", n_scans)

  for (i in seq_len(n_scans)) {
    xi <- matx[i, ]
    yi <- maty[i, ]
    zi <- matz[i, ]
    oks <- .Internal(which(!is.na(xi)))
    xs[[i]] <- xi[oks]
    ys[[i]] <- yi[oks]
    zs[[i]] <- zi[oks]
  }
  rm(list = c("xi", "yi", "zi"))

  if (n_scans <= 3L || !join_ms)
    return(list(x = xs, y = ys, z = zs))

  ansz <- ansy <- ansx <- vector("list", n_scans)
  
  xbf <- xs[[1]]
  ybf <- ys[[1]]
  zbf <- zs[[1]]
  ansx[[1]] <- xbf
  ansy[[1]] <- ybf
  ansz[[1]] <- zbf
  
  for (i in 2:(n_scans-1)) {
    ix <- i+1
    xcr <- xs[[i]]
    ycr <- ys[[i]]
    zcr <- zs[[i]]
    xaf <- xs[[ix]]
    yaf <- ys[[ix]]
    zaf <- zs[[ix]]
    ansx[[i]] <- c(xbf, xcr, xaf)
    ansy[[i]] <- c(ybf, ycr, yaf)
    ansz[[i]] <- c(zbf, zcr, zaf)
    
    xbf <- xcr
    ybf <- ycr
    zbf <- zcr
  }
  
  list(x = ansx, y = ansy, z = ansz)
}


#' Spreads one-hit wonders to adjacent scans.
#'
#' @param matx The data matrix of moverzs.
#' @param maty The data matrix of intensities.
#' @param matz The data matrix of charge states.
#' @param ns The number of observing spectra that have contributed to an MS
#'   feature. One-hit wonders correspond to \code{ns == 1L}.
#' @param ps The row positions (along LC scans) of features.
#' @param gap The gap width for one-hit wonders.
#' @param join_ms Logical; if TRUE, combine adjacent entries.
spreadMSohw <- function (matx, maty, matz, ns, ps, gap = 2L, join_ms = FALSE)
{
  # stopifnot((gap >= 1L))
  
  n_scans <- nrow(matx)
  n_masses <- ncol(matx)
  # 2L: since the first neighbor already jointed in flattenMSxyz
  n_nearest <- if (join_ms) 2L else 1L
  
  zout <- yout <- xout <- vector("list", n_scans)

  # across columns of masses
  for (i in 1:n_masses) {
    ni <- ns[[i]]
    pi <- ps[[i]]
    rows1 <- pi[ni == 1L]
    
    # goes through one-hit wonders
    for (j in seq_along(rows1)) {
      rj1 <- rows1[[j]]
      x <- matx[rj1, i]
      y <- maty[rj1, i]
      z <- matz[rj1, i]

      rng1 <- c(max((rj1 - gap), 1):max((rj1 - n_nearest), 1), 
                min((rj1 + n_nearest), n_scans):min((rj1 + gap), n_scans))
      
      # go through the k neighbors of a one-hit wonder
      for (k in seq_along(rng1)) {
        rng1k <- rng1[[k]]
        xout[[rng1k]] <- c(xout[[rng1k]], x)
        yout[[rng1k]] <- c(yout[[rng1k]], y)
        zout[[rng1k]] <- c(zout[[rng1k]], z)
      }
    }
  }
  
  list(x = xout, y = yout, z = zout)
}


#' Finds the positions of logical gates.
#' 
#' The maximum width of a gate is 2L.
#' 
#' @param vals An integer vector.
#' 
#' @examples
#' library(mzion)
#' 
#' # starts at low and ends at low
#' # find_gates(c(100, 102, 104:108, 110, 112, 114:117, 119))
#' 
#' # low-high
#' # find_gates(c(100, 102, 104:108, 110, 112, 114:117))
#' 
#' # high-low
#' # find_gates(c(95:100, 102, 104:108, 110:111, 114:117, 119))
#' 
#' # high-high
#' # find_gates(c(95:100, 102, 104:108, 110, 112:113, 116:119))
#' 
#' # all-zeros (discrete)
#' # find_gates(c(6L, 12L, 18L, 164L))
#' 
#' # all-ones
#' # find_gates(1:5)
#' 
#' # up-width == 1L
#' # find_gates(c(143L, 159L, 310L, 311L, 316L))
#' 
#' # all discrete
#' # find_gates(c(143L, 159L, 310L, 316L))
#' 
#' # single TRUE
#' # find_gates(c(310L, 311L))
#' 
#' #' # single value
#' # find_gates(c(310L))
#' 
#' # find_gates(c(281, 326, 335, 336, 337, 437, 447, 448, 449, 450, 553, 554, 557))
find_gates <- function (vals)
{
  # all discrete
  if (is.null(edges <- find_gate_edges(vals))) {
    return(NULL)
  }

  ups  <- edges[["ups"]]
  dns  <- edges[["dns"]]
  ps   <- mapply(function (x, y) x:y, ups, dns, SIMPLIFY = FALSE, 
                 USE.NAMES = FALSE)
  lens  <- dns - ups + 1L # lengths(ps)
  longs <- lens > 2L
  
  if (any(longs)) {
    longs <- .Internal(which(longs))
    ps0   <- ps[-longs]
    ps1   <- ps[longs]
    lens1 <- lens[longs]
    
    for (i in seq_along(ps1)) {
      psi <- ps1[[i]]
      
      # adds a trailing 0L to maintain an even length
      if ((li <- lens1[[i]]) %% 2L) {
        psi <- c(ps1[[i]], 0L)
        li  <- li + 1L
      }
      
      ps1[[i]] <- split(psi, rep(1:(li / 2L), each = 2L))
    }
    
    ps1 <- .Internal(unlist(ps1, recursive = FALSE, use.names = FALSE))
    ps  <- c(ps0, ps1)
    
    ord <- lapply(ps, `[[`, 1)
    ord <- .Internal(unlist(ord, recursive = FALSE, use.names = FALSE))
    ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, 
                               TRUE, ord))
    ps  <- ps[ord]
  }
  
  ps
}


#' Helper to find the positions of signal edges
#'
#' @param vals A vector of integers or logical values.
#' @examples
#' vals <- c(1:3, 5:7, 10, 14:16)
#' find_gate_edges(vals)
#'
#' @return The indexes of ups and downs along \code{vals}. Note that discrete
#'   values without neighbors are excluded.
find_gate_edges <- function (vals)
{
  nvs <- length(vals)
  
  if (nvs <= 1L) {
    return(NULL)
  }

  dvs  <- diff(vals, 1L)
  vec  <- dvs == 1L
  lenv <- nvs - 1L
  
  if (vec[lenv]) {
    vec  <- c(vec, FALSE)
    lenv <- lenv + 1L
  }
  
  # all discrete
  if (!any(vec)) {
    return(NULL)
  }

  ds   <- diff(vec)
  ups  <- which(ds == 1L) + 1L
  dns  <- which(ds == -1L) + 1L
  lenu <- length(ups)
  lend <- length(dns)
  
  if (lenu == lend) {
    if (ups[[1]] > dns[[1]]) {
      ups <- c(1L, ups)
      dns <- c(dns, lenv)
    }
  }
  else if (lenu < lend) {
    ups <- c(1L, ups)
  }
  # should not occur with ensured trailing FALSE
  else if (lenu > lend) {
    dns <- c(dns, ups[lenu])
  }

  list(ups = ups, dns = dns)
}


#' Helper to find the positions of signal edges
#'
#' Results includes discrete values (one-hit-wonders).
#'
#' @param vals A vector of integers or logical values.
#' @examples
#' find_gate_edges2(c(1:3, 5:7, 10, 14:16))
#' find_gate_edges2(c(3, 5)); find_gate_edges2(c(3, 5:6, 8))
#' find_gate_edges2(2:5)
#' find_gate_edges2(c(2:3, 5:6, 8:9))
#' find_gate_edges2(c(3:5, 7, 12:14, 18))
#' vals <- c(3, 5)
#' vals <- c(3:5, 7:10, 12:14); vals <- c(vals, 18)
#' vals <- c(3:5, 7, 9)
#' vals <- c(3, 5:7, 10:11, 13)
#' vals <- c(3, 5, 7, 8, 10, 12:13)
#' 
#' @return The indexes of ups and downs along \code{vals}. Indexes of discrete
#'   values are also included.
find_gate_edges2 <- function (vals)
{
  nvs <- length(vals)
  if (!nvs) { return(list(ups = NULL, dns = NULL, ohw = NULL)) }
  if (nvs == 1L) { return(list(ups = 1L, dns = 1L, ohw = NULL)) }
  
  dvs  <- diff(vals, 1L)
  lenv <- nvs - 1L
  oks  <- dvs > 1L
  ioks <- .Internal(which(oks))
  noks <- length(ioks)
  
  # all consecutive; vals <- c(2:5)
  if (!noks) { return(list(ups = 1L, dns = nvs, ohw = NULL)) }
  
  ix  <- ioks + 1L
  ixn <- ix[[noks]]
  
  if (ixn == nvs) { # the last is discrete (ohw)
    # okx <- .Internal(which(dvs[ix[-noks]] > 1L)) # okx is one shorter than ix
    # iohs <- c(ix[okx], ixn)
    okx <- dvs[ix] > 1L # the last entry out of bound -> okx: c(TRUE, FALSE, NA)
    okx[[noks]] <- TRUE # NA -> TRUE
    iohs <- ix[okx]
  }
  else {
    okx  <- dvs[ix] > 1L # can be all FALSE: vals <- c(2:3, 5:6, 8:9)
    iohs <- ioks[okx] + 1L
    
    if (!length(iohs)) {
      iohs <- NULL
    }
  }
  
  if (ioks[[1]] == 1L) { # the first is ohw
    iohs <- c(1L, iohs)
  }
  
  ## this point on: essentially the same in `find_gate_edges`
  vec <- dvs == 1L
  
  if (vec[lenv]) {
    vec  <- c(vec, FALSE)
    lenv <- lenv + 1L
  }
  
  # all discrete
  if (!any(vec)) {
    # return(NULL)
    return(list(ups = NULL, dns = NULL, ohw = iohs))
  }
  
  ds   <- diff(vec)
  ups  <- which(ds == 1L) + 1L
  dns  <- which(ds == -1L) + 1L
  lenu <- length(ups)
  lend <- length(dns)
  
  if (lenu == lend) {
    if (ups[[1]] > dns[[1]]) {
      ups <- c(1L, ups)
      dns <- c(dns, lend)
    }
  }
  else if (lenu < lend) {
    ups <- c(1L, ups)
  }
  # should not occur with ensured trailing FALSE
  else if (lenu > lend) {
    dns <- c(dns, ups[lenu])
  }
  
  list(ups = ups, dns = dns, ohw = iohs)
}


#' Traces LCMS peaks.
#'
#' For both moverzs and intensities.
#'
#' @param xs Vectors of moverzs at a given isolation center (\code{icenter}).
#'   Each entry corresponds to a full vector of MS1 or MS2.
#' @param ys Vectors of intensities.
#' @param zs Vectors of charge states.
#' @param temp_dir A temporary directory.
#' @param icenter The index of an isolation center.
#' @param ms_lev The level of MS.
#' @param from The starting value for mass binning.
#' @param step A step size for mass binning.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param reord Logical; re-order data or not.
#' @param cleanup Logical; cleans up xs, ys and zs or not. Set the value to
#'   FALSE to maintain one-to-one correspondence between input (data frame) and
#'   the outputs. This will help, e.g., keep track of scan numbers in the input.
#' @param replace_ms1_by_apex Logical; if TRUE, fill all entries within a gate
#'   by its apex values.
#' @param direct_out Logical; if TRUE, return outputs directly.
#' @param filename A peaklist filename.
#' @param temp_dir A temp_dir to the filename.
traceLCMS <- function (xs, ys, zs, n_dia_scans = 4L, from = 115L, step = 8E-6, 
                       icenter = 1L, ms_lev = 2L, temp_dir = NULL, filename = NULL, 
                       reord = TRUE, cleanup = TRUE, replace_ms1_by_apex = FALSE, 
                       direct_out = FALSE)
{
  lens <- lengths(xs)
  
  if (reord) {
    ords <- lapply(xs, order)

    for (i in seq_along(xs)) {
      if (lens[[i]]) {
        ordi <- ords[[i]]
        xs[[i]] <- xs[[i]][ordi]
        ys[[i]] <- ys[[i]][ordi]
        zs[[i]] <- zs[[i]][ordi]
      }
    }
    rm(list = "ords", "ordi")
  }
  
  ## collapses MS data by the indexes of mass bins; 
  # matrix outputs; rows: scans; columns: masses at "x", intensities at "y" ...
  ans <- collapse_xyz(
    xs = xs, ys = ys, zs = zs, temp_dir = temp_dir, icenter = icenter, 
    ms_lev = ms_lev, lwr = from, step = step, cleanup = cleanup, 
    direct_out = direct_out)
  
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  ansz <- ans[["z"]]
  rm(list = c("ans"))
  gc()
  
  ## traces MS data matrices across LC scans; rows: scans; columns: masses
  nrc <- dim(ansy)
  nr <- nrc[[1]]
  nc <- nrc[[2]]

  xout <- yout <- matrix(rep_len(NA_real_, nc * nr), ncol = nc)
  zout <- matrix(rep_len(NA_integer_, nc * nr), ncol = nc)
  ranges <- apexes <- ns <- vector("list", nc)
  
  if (replace_ms1_by_apex) {
    for (i in 1:nc) {
      # temporary ts = NULL
      gates <- find_lc_gates(xs = NULL, ys = ansy[, i], ts = NULL, 
                             n_dia_scans = n_dia_scans)
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]
      
      for (j in seq_along(rows)) {
        rgj <- rngs[[j]]
        rwj <- rows[[j]]
        xout[rgj, i] <- ansx[rwj, i]
        yout[rgj, i] <- ansy[rwj, i]
        zout[rgj, i] <- ansz[rwj, i]
      }
    }
  }
  else {
    for (i in 1:nc) {
      # temporary ts = NULL
      gates <- find_lc_gates(xs = NULL, ys = ansy[, i], ts = NULL, 
                             n_dia_scans = n_dia_scans)
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]

      xout[rows, i] <- ansx[rows, i]
      yout[rows, i] <- ansy[rows, i]
      zout[rows, i] <- ansz[rows, i]
    }
  }

  rm(list = c("ansx", "ansy", "ansz"))
  gc()
  
  list(x = xout, y = yout, z = zout, n = ns, p = apexes, range = ranges)
}


#' Gather MS data
#'
#' Similar to \link{collapse_mms1ints} with the addition of zs.
#'
#' Values of \code{xs} can be smaller than \code{lwr} since no ms2_mass cut-offs
#' (some z undetermined).
#'
#' @param xs Vectors of m-over-z values within a given isolation center
#'   (\code{icenter}). Each entry corresponds to a full vector of MS1 or MS2.
#' @param ys Vectors of intensities.
#' @param zs Vectors of charge states.
#' @param temp_dir A temporary directory.
#' @param icenter The index of an isolation center.
#' @param ms_lev The level of MS.
#' @param lwr A lower bound as the starting point in mass binning.
#' @param step A step size for mass binning.
#' @param cleanup Logical; cleans up xs, ys and zs or not. Set the value to
#'   FALSE to maintain one-to-one correspondence between input (data frame) and
#'   the outputs. This will help, e.g., keep track of scan numbers in the input.
#' @param direct_out Logical; if TRUE, return outputs directly.
#' @importFrom fastmatch %fin%
#' @return Matrices of moverzs, intensities and charge states. Columns: masses;
#'   rows: scans.
collapse_xyz <- function (xs = NULL, ys = NULL, zs = NULL, temp_dir = NULL, 
                          icenter = 1L, ms_lev = 2L, lwr = 115L, step = 1e-5, 
                          cleanup = TRUE, direct_out = FALSE)
{
  if (cleanup) {
    null_out <- list(x = NULL, y = NULL, z = NULL)
    
    # all xs are NULL
    if (!any(oksx <- lengths(xs) > 0L)) {
      return(null_out)
    }

    oksx <- .Internal(which(oksx))
    xs <- xs[oksx]
    ys <- ys[oksx]
    zs <- zs[oksx]

    # remove zero intensities
    for (i in seq_along(oky <- lapply(ys, `>`, 0))) {
      oki <- .Internal(which(oky[[i]]))
      xs[[i]] <- xs[[i]][oki]
      ys[[i]] <- ys[[i]][oki]
      zs[[i]] <- zs[[i]][oki]
    }

    # does this again after ys removals
    if (!any(oks <- lengths(xs) > 0L)) {
      return(null_out)
    }

    oks <- .Internal(which(oks))
    xs <- xs[oks]
    ys <- ys[oks]
    zs <- zs[oks]
  }

  ixs <- lapply(xs, index_mz, lwr, step)
  
  # remove duplicated ixs
  for (i in seq_along(ixs)) {
    ix <- ixs[[i]]
    oks <- .Internal(which(!duplicated(ix)))
    ixs[[i]] <- ix[oks]
    xs[[i]]  <- xs[[i]][oks]
    ys[[i]]  <- ys[[i]][oks]
    zs[[i]]  <- zs[[i]][oks]
  }

  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups <- lapply(ixs, function (x) unv %in% x)
  
  xout <- mapcoll_xyz(vals = xs, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "xs", direct_out = direct_out)
  yout <- mapcoll_xyz(vals = ys, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "ys", direct_out = direct_out)
  zout <- mapcoll_xyz(vals = zs, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "zs", direct_out = direct_out)

  if (!direct_out) {
    x_nm <- paste0("xs", ms_lev, "_diauniv_", icenter, ".rds")
    y_nm <- paste0("ys", ms_lev, "_diauniv_", icenter, ".rds")
    z_nm <- paste0("zs", ms_lev, "_diauniv_", icenter, ".rds")
    xout <- qs::qread(file.path(temp_dir, x_nm))
    yout <- qs::qread(file.path(temp_dir, y_nm))
    zout <- qs::qread(file.path(temp_dir, z_nm))
  }

  ### collapses adjacent entries
  ps   <- find_gates(unv)
  lenp <- length(ps)
  
  # all discrete values
  if (is.null(ps)) {
    return(list(x = xout, y = yout, z = zout))
  }

  ps2 <- lapply(ps, `[[`, 2)
  ps2 <- .Internal(unlist(ps2, recursive = FALSE, use.names = FALSE))
  
  for (i in 1:lenp) {
    c12 <- ps[[i]]
    c2 <- c12[[2]]
    
    # with values in both columns, simply overwrite: 1 <- 2; 
    # c2 can be 0
    if (!c2) {
      next
    }

    c1 <- c12[[1]]
    oks <- .Internal(which(!is.na(xout[, c2])))
    xout[oks, c1] <- xout[oks, c2]
    yout[oks, c1] <- yout[oks, c2]
    zout[oks, c1] <- zout[oks, c2]
  }
  
  xout <- xout[, -ps2, drop = FALSE]
  yout <- yout[, -ps2, drop = FALSE]
  zout <- zout[, -ps2, drop = FALSE]
  
  list(x = xout, y = yout, z = zout)
}


#' Maps data onto the universe.
#'
#' For RAM efficiency.
#'
#' @param vals Data in one of moverzs, intensities or charge states.
#' @param ups Positions of \code{vals} in universe. Should be the same for among
#'   moverzs, intensities and charge states.
#' @param lenx The length of \code{vals} (number of rows in a matrix).
#' @param lenu The number of entries in the universe (number of columns in a
#'   matrix).
#' @param temp_dir A temporary directory.
#' @param icenter The index of an isolation center.
#' @param ms_lev The level of MS.
#' @param type The type of data in one of xs, ys and zs.
#' @param direct_out Logical; if TRUE, return outputs directly.
#' @return A matrix where rows correspond to \code{vals} and columns to
#'   \code{ups}.
mapcoll_xyz <- function (vals, ups, lenx, lenu, temp_dir, icenter = 1L, 
                         ms_lev = 2L, type = "xs", direct_out = FALSE)
{
  if (!length(vals)) {
    return(NULL)
  }

  out <- if (type == "zs") {
    rep_len(list(rep_len(NA_integer_, lenu)), lenx)
  }
  else if (type %in% c("xs", "ys")) {
    rep_len(list(rep_len(NA_real_, lenu)), lenx)
  }
  else {
    rep_len(list(rep_len(NA_real_, lenu)), lenx)
  }

  for (i in 1:lenx) {
    cols <- ups[[i]]
    out[[i]][cols] <- vals[[i]]
  }
  
  out <- do.call(rbind, out)
  
  if (direct_out) {
    return(out)
  }

  out_name <- paste0(type, ms_lev, "_diauniv_", icenter, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  
  invisible(NULL)
}


#' Finds the gates of retention times.
#'
#' @param xs A vector of m-over-z value for the same (approximate) mass along
#'   LC.
#' @param ys A vector of intensity value for the same (approximate) mass along
#'   LC.
#' @param ts The retention times corrresponding to \code{ys}.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param yco A cut-off in intensity (Y) values.
#' @param y_rpl A replacement value of intensities.
#' @param ytot_co A more permissive cut-off in peak area for finding all peaks.
#' @param max_perc The maximum percentage of baseline levels.
#' @param min_n The minimum number of points across a peak for consideration.
#' @param max_n The maximum number of points across a peak for consideration.
#' @examples
#' mzion:::find_lc_gates(ys = c(10,0,0,0,11,15,15,0,0,12,0,10,0,0,10), xs = rep_len(100, 15), ts = seq_len(15), n_dia_scans = 2)
#' mzion:::find_lc_gates(ys = c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 1), 20, 50), xs = rep_len(100, 18), ts = seq_len(18))
#' mzion:::find_lc_gates(ys = c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 4), 20, 50), xs = rep_len(100, 21), ts = seq_len(21))
#'
#' # all discrete
#' mzion:::find_lc_gates(ys = c(rep(0, 5), 100, rep(0, 6), 200, rep(0, 4), 50), xs = rep_len(100, 18), ts = seq_len(18))
#' 
#' # two passes
#' ys <- c(rep_len(0, 10), rep_len(100, 50), rep_len(0, 10), rep_len(100, 50), rep_len(3, 5), rep_len(50, 80), rep_len(0, 10))
#' xs <- rep(5, length(ys))
#' ts <- seq_along(ys)
#' mzion:::find_lc_gates(ys = ys, xs = xs, ts = ts)
#' 
#' # a tail gate on one end (right end) has Y values greater the valley between the first two "real" peaks
#' ys <- c(rep_len(NA_real_,15),196358.0,rep_len(NA_real_,4),72939.6,rep_len(NA_real_,2),
#'   200039.5,198492.1,251862.1,230444.3,87445.9,492812.3,331446.7,361864.8,475956.8,599853.3,
#'   964008.5,1023564.2,1388958.9,1629472.4,1475251.1,2075193.1,2352548.2,2846973.2,4058301.8,3947633.8,5008576.0,
#'   5398789.0,6200687.0,7283479.0,6824848.0,7587668.0,10106542.0,11173275.0,11035185.0,12388256.0,16368471.0,13374066.0,
#'   13641590.0,14868098.0,15204119.0,16079883.0,18423862.0,18902162.0,18423320.0,20154212.0,18065664.0,17772998.0,18746912.0,
#'   20148930.0,17781916.0,20552498.0,18365672.0,15672191.0,14854731.0,15693521.0,13779323.0,16107755.0,15195194.0,12831090.0,
#'   14584699.0,13121290.0,12413923.0,11027703.0,10901884.0,13319185.0,10706831.0,8747522.0,10857365.0,8368990.0,8653836.0,
#'   9097412.0,9122940.0,7533711.5,8070598.0,5867927.0,7320747.5,7169463.0,5733574.5,6551820.0,6305446.5,5464279.0,
#'   5223931.5,6137411.5,5216691.5,4426827.0,4302836.5,4035890.2,3729244.5,4129466.5,3168772.0,4331294.0,4424808.5,
#'   3833815.2,3709572.2,3032854.5,3455906.0,2951570.2,2819246.2,2591228.0,3529138.0,2646039.5,2864500.8,3268786.8,
#'   2451983.5,2314034.5,3149912.0,3287747.5,2938702.5,2763761.5,1903592.5,2861269.5,2467254.8,1946265.2,2613601.8,
#'   2208809.8,2101898.2,2171315.8,2446137.8,1964351.9,2023059.0,2979205.8,2976068.0,3081084.0,2944902.0,3687758.5,
#'   3594447.2,3483423.2,5176589.5,4857718.5,6138461.0,6738957.5,6161186.0,6903845.5,8240505.5,8328314.0,9211307.0,
#'   10451003.0,12589421.0,11609010.0,14125617.0,13773539.0,14815689.0,18450236.0,15245151.0,17122230.0,17680004.0,16380083.0,
#'   19416518.0,22584712.0,20886734.0,20594192.0,20921544.0,19620582.0,19714208.0,21572454.0,22365472.0,21324292.0,20963878.0,
#'   21934052.0,21021322.0,21450586.0,18756072.0,18683146.0,17521552.0,14446918.0,16657946.0,14589204.0,14723747.0,14555312.0,
#'   13498088.0,13474049.0,12276403.0,14057293.0,11991678.0,11141724.0,10432426.0,11457503.0,9032645.0,9328507.0,8703553.0,
#'   8696164.0,8235523.5,7634871.0,7186748.0,6313092.0,7423992.0,5722295.0,6303344.0,4935099.0,5156393.0,4614185.0,
#'   5217295.0,4403636.0,5018784.5,6287227.0,3746022.2,3344260.2,3794109.0,3256452.8,3016607.0,2951111.0,3258742.2,
#'   3226547.8,3049823.5,2950860.0,2663367.0,3134546.2,4605861.0,2545474.5,2621648.0,2636332.5,2311735.0,2633117.8,
#'   2127756.8,1885285.1,2323312.2,2767798.5,1994432.6,1745769.6,2981810.8,1478732.9,2110425.8,1489411.2,1631056.1,
#'   1818856.6,2130620.2,1759935.5,1224843.9,1711871.8,1441695.5,1556760.6,1915420.8,1760751.4,1566991.6,1577649.4,
#'   1843721.1,1457992.1,1948946.2,1435371.4,1499531.4,1370897.2,1157172.9,1998821.9,751608.8,1163820.2,1285921.9,
#'   1259784.9,984264.4,1167820.6,1032068.7,1146866.4,1104102.2,1045812.1,921146.6,1294800.5,1257817.6,1184590.6,
#'   1247893.4,1093933.8,573998.8,851412.1,1047905.9,1096414.8,1652794.8,1345635.9,806337.9,1378611.4,2103232.2,
#'   1747740.4,1875224.8,1178567.4,1043248.6,1302969.1,1364105.5,1417062.9,3013454.5,1512214.8,1474206.9,2211758.0,
#'   1434366.6,2015054.1,3293890.2,2786769.0,2825213.2,2017830.1,1835824.5,3046673.2,2046926.4,2337327.8,3682744.5)
#' 
#' ts <- c(2084.422,2084.803,2085.038,2085.634,2085.905,2086.525,2087.121,2087.500,2087.806,2088.077,2088.480,2088.858,2089.236,2089.650,
#'   2089.885,2090.156,2090.561,2091.013,2091.465,2091.699,2092.152,2092.457,2092.837,2093.180,2093.602,2094.092,2094.547,2094.783,
#'   2095.341,2095.941,2096.177,2096.484,2096.828,2097.428,2097.809,2098.045,2098.318,2098.698,2099.188,2099.641,2100.059,2100.442,
#'   2100.679,2101.207,2101.444,2101.717,2102.097,2102.551,2102.788,2103.314,2103.551,2103.932,2104.277,2104.513,2104.969,2105.348,
#'   2105.695,2106.111,2106.420,2106.763,2107.108,2107.622,2107.965,2108.200,2108.689,2109.034,2109.377,2109.902,2110.502,2110.701,
#'   2111.008,2111.497,2111.805,2112.221,2112.458,2112.838,2113.432,2113.669,2114.193,2114.538,2115.123,2115.616,2116.031,2116.376,
#'   2116.683,2117.208,2117.515,2117.966,2118.310,2118.607,2119.099,2119.443,2119.752,2120.131,2120.512,2121.001,2121.237,2121.618,
#'   2121.854,2122.126,2122.362,2122.815,2123.268,2123.651,2124.102,2124.590,2124.973,2125.352,2125.802,2126.037,2126.458,2127.021,
#'   2127.471,2127.960,2128.523,2128.903,2129.282,2129.518,2130.006,2130.352,2130.720,2131.100,2131.481,2131.826,2132.425,2133.020,
#'   2133.510,2134.110,2134.418,2134.799,2135.033,2135.378,2135.613,2136.028,2136.519,2136.898,2137.205,2137.513,2137.928,2138.238,
#'   2138.571,2138.805,2139.284,2139.592,2139.924,2140.400,2140.889,2141.199,2141.641,2141.912,2142.219,2142.660,2143.038,2143.273,
#'   2143.692,2144.071,2144.305,2144.900,2145.317,2145.661,2146.115,2146.640,2147.055,2147.361,2147.779,2148.122,2148.611,2148.845,
#'   2149.334,2149.898,2150.133,2150.478,2150.821,2151.201,2151.763,2152.148,2152.564,2152.909,2153.144,2153.743,2154.159,2154.652,
#'   2155.070,2155.380,2155.687,2156.068,2156.304,2156.721,2157.175,2157.413,2157.902,2158.138,2158.736,2159.044,2159.423,2159.985,
#'   2160.222,2160.711,2161.200,2161.436,2161.780,2162.305,2162.902,2163.137,2163.737,2164.262,2164.679,2165.279,2165.804,2166.330,
#'   2166.819,2167.055,2167.654,2168.253,2168.634,2169.015,2169.251,2169.777,2170.302,2170.537,2171.027,2171.626,2171.969,2172.492,
#'   2172.835,2173.071,2173.306,2173.614,2174.066,2174.554,2174.970,2175.565,2176.158,2176.645,2176.987,2177.222,2177.817,2178.050,
#'   2178.571,2178.877,2179.184,2179.550,2179.891,2180.327,2180.671,2181.050,2181.319,2181.607,2182.046,2182.280,2182.549,2182.961,
#'   2183.193,2183.426,2183.659,2184.072,2184.486,2184.862,2185.202,2185.545,2185.922,2186.227,2186.569,2186.874,2187.181,2187.486,
#'   2187.791,2188.023,2188.364,2188.742,2189.339,2189.536,2189.950,2190.183,2190.524,2190.867,2191.235,2191.611,2191.954,2192.259,
#'   2192.636,2193.014,2193.246,2193.586,2194.000,2194.233,2194.465,2194.878,2195.255,2195.741,2195.973,2196.460,2196.909,2197.467,
#'   2197.700,2198.225,2198.819,2199.268,2199.681,2200.023,2200.255,2200.815,2201.340,2201.935,2202.530,2203.125,2203.722,2203.955)
#' 
#' mzion:::find_lc_gates(ys = ys, ts = ts)
#' @importFrom fastmatch %fin%
#' @return Scan indexes of LC peaks.
find_lc_gates <- function (ys = NULL, xs = NULL, ts = NULL, n_dia_scans = 6L, 
                           yco = 100, y_rpl = 2.0, max_perc = .05, min_n = 10L, 
                           ytot_co = 2E5, max_n = 200L)
{
  # if (min_n <= 1L) { stop("Choose a value of `min_n` greater one.") }
  # if (n_dia_scans <= 0L) return(.Internal(which(ys > 0))) # should not occur
  
  ys[ys < yco] <- NA_real_
  
  if (all(is.na(ys))) {
    return(NULL) 
  }

  ys   <- fill_lc_gaps(ys = ys, n_dia_scans = n_dia_scans, y_rpl = y_rpl)
  ioks <- .Internal(which(ys >= y_rpl))
  nx   <- length(ioks)
  
  # one `one-hit wonder` across LC
  if (nx == 1L) {
    return(NULL)
  }
  
  edges <- find_gate_edges(ioks)
  
  # all discrete `one-hit wonders`
  if (is.null(edges)) {
    return(NULL)
  }
  
  ups  <- edges[["ups"]]
  dns  <- edges[["dns"]]
  
  rngs <- if (length(ups) == 1L) {
    list(ups:dns)
  }
  else {
    mapply(`:`, ups, dns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  # include or exclude counts of `y_rpl` replacement
  lens <- lengths(rngs)
  
  if (!(nps  <- length(ups))) {
    return(NULL)
  }
  
  stas <- ends <- vector("list", nps)
  
  for (i in seq_len(nps)) {
    rngi <- rngs[[i]]
    ri   <- ioks[rngi] # indexes in relative to ys
    lenr <- lens[[i]]
    ysi  <- ys[ri]
    ymax <- max(ysi, na.rm = TRUE)
    ny   <- sum(ysi > ymax * .05) # was .02
    ny2  <- sum(ysi > y_rpl)
    
    if (ny < min_n || ny2 < .7 * lenr || 
        is_partial_peak(ysi, width = 200L, min_p = 10L)
        ) {
      next
    }
    
    # an peak-specific baseline
    ymin <- find_baseline(
      vals = ysi, vmax = ymax, perc = .02, max_perc = max_perc, lwr = yco)
    # no "good" baseline
    if (is.null(ymin)) { next }

    if (FALSE) {
      data.frame(x = ts[ri]/60, y = ysi) |>
        ggplot2::ggplot() + 
        ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                              color = "gray", linewidth = .1)
    }
    
    if (lenr <= 60L) {
      stas[[i]] <- ri[[1]] # the starting index of gate-i in relative to ys
      ends[[i]] <- ri[lenr]
    }
    else {
      r_first <- ri[[1]]
      # z <- (seq_along(ysi) + r_first); plot(ysi ~ z)
      perc <- min(ymin / ymax, max_perc)
      ans  <- lapply(c(perc, .06 * 1:10 + perc), find_lc_gates2, 
                     ys = ysi, sta = r_first, n_dia_scans = n_dia_scans, 
                     y_rpl = y_rpl)
      ans_stas <- lapply(ans, `[[`, "stas")
      n_stas   <- lengths(ans_stas) # can be all zeros
      n_max    <- .Internal(which.max(n_stas))
      n_stas   <- n_stas[[n_max]] # can be zero
      
      if (n_stas <= 1L) {
        stas[[i]] <- r_first # also check left bending
        ends[[i]] <- ri[lenr] # also check right bending
        
        next
      }
      
      ans <- ans[[n_max]]
      ans_stas <- ans_stas[[n_max]]
      ans_ends <- ans[["ends"]]
      
      if (first_bad <- ans$first_bad) {
        first_dn <- ans$dn1
      }
      else {
        first_dn <- NULL
      }
      
      if (last_bad <- ans$last_bad) {
        last_up <- ans$upn
      }
      else {
        last_up <- NULL
      }
      
      stax <- endx <- vector("integer", n_stas)
      # stax[[1]] <-  ri[[1]]
      # endx[[n_stas]] <- ri[lenr]
      stax[[1]]      <- if (first_bad) first_dn else ri[[1]]
      endx[[n_stas]] <- if (last_bad)  last_up  else ri[lenr]
      
      
      for (j in 2:n_stas) {
        hend <- ans_ends[[j - 1]]
        hsta <- ans_stas[[j]]
        pmid <- as.integer((hend + hsta) / 2L)
        endx[[j-1]] <- max(pmid - 1L, hend)
        stax[[j]]   <- min(pmid + 1L, hsta)
      }
      
      stas[[i]] <- stax
      ends[[i]] <- endx
    }
  }
  
  # if (any(lengths(stas) > 1L)) {}
  # vector or list; NULL stas and ends also removed
  stas <- unlist(stas, recursive = FALSE, use.names = FALSE)
  ends <- unlist(ends, recursive = FALSE, use.names = FALSE)
  nps2 <- length(stas)
  
  if (!nps2) {
    return(NULL)
  }
  
  yints  <- fwhm <- vector("numeric", nps2)
  apexs  <- xstas <- vector("integer", nps2)
  ranges <- vector("list", nps2)
  
  for (i in 1:nps2) {
    ###
    # i <- 8L
    ###
    
    si  <- stas[[i]]
    ixi <- si:ends[[i]]
    tsi <- ts[ixi]
    ysi <- ys[ixi]
    yts <- calcAUC(ys = ysi, ts = tsi, rng = ixi, yco = yco, ytot_co = ytot_co, 
                   min_n = 15L, err = 2.0)
    
    if (is.null(yts)) { next }
    
    mi  <- yts[["idx"]]
    rgx <- yts[["rng"]]
    rg1 <- rgx[[1]]

    yints[[i]]  <- yts[["area"]]
    apexs[[i]]  <- rg1 + mi - 1L
    xstas[[i]]  <- rg1
    ranges[[i]] <- rgx
    fwhm[[i]]   <- yts[["fwhm"]]
  }
  
  # if (sum(yints) / sum(ys, na.rm = TRUE) < .01) { return(NULL) }
  
  ###
  # check the baseline against to guard against malformed peaks...
  # filter by 1% area
  ###
  
  ## outputs
  ioks <- 
    .Internal(which(apexs > 0L & yints >= max(yints, na.rm = TRUE) * .01))
  noks <- length(ioks)

  if (!noks) { # noks > 4L || 
    return(NULL)
  }
  
  ns <- lengths(ranges)
  if (noks < nps2) {
    apexs  <- apexs[ioks]
    yints  <- yints[ioks]
    ns     <- ns[ioks]
    ranges <- ranges[ioks]
    xstas  <- xstas[ioks]
    fwhm   <- fwhm[ioks]
  }
  
  list(apex = apexs, yints = yints, ns = ns, ranges = ranges, xstas = xstas, 
       fwhm = fwhm)
}


#' Check the convexness of a peak
#'
#' @param ys A vector of intensity values.
#' @param i_sta The left-start index of a peak along \code{ys}.
#' @param i_midsta The middle-start index of a peak along \code{ys}.
#' @param i_midend The middle-end index of a peak along \code{ys}.
#' @param i_end The right-end index of a peak along \code{ys}.
#' @param min_n The minimum points across the peak profile of \code{ys} for
#'   consideration.
check_peak_convex <- function (ys, i_sta, i_midsta, i_midend, i_end, 
                               min_n = 15L)
{
  len <- length(ys)
  
  if (len < min_n) {
    return(FALSE)
  }

  if (i_sta + 3L >= i_midsta) { # i_sta == i_midsta
    return(FALSE)
    # m1 <- 0.0
  }
  else {
    m1 <- median(ys[i_sta:(i_midsta - 1L)])
  }
  
  if (i_midend + 3L >= i_end) { # i_midend == i_end
    return(FALSE)
    # m3 <- 0.0
  }
  else {
    m3  <- median(ys[(i_midend + 1L):i_end])
  }
  
  m2 <- median(ys[i_midsta:i_midend])
  
  m2 > m1 * 1.25 && m2 > m3 * 1.25
}


#' Check the convexness of a peak
#'
#' Use a sub-sequence of intensity values with a pre-calculated left and right
#' margins (not much faster than \link{check_peak_convex}).
#'
#' @param ys A vector of intensity values.
#' @param i_midsta The middle-start index of a peak along \code{ys}.
#' @param i_midend The middle-end index of a peak along \code{ys}.
#' @param min_n The minimum points across the peak profile of \code{ys} for
#'   consideration.
check_peak_convex_sub <- function (ys, i_midsta, i_midend, min_n = 15L)
{
  # ys = ys0, i_midsta = ista - iok1 + 1L, i_midend = iend -iok1 + 1L
  
  len <- length(ys)
  
  if (len < min_n) {
    return(FALSE)
  }
  
  rg1 <- 1:(i_midsta - 1L)
  rg2 <- i_midsta:i_midend
  rg3 <- (i_midend + 1L):len

  m1 <- median(ys[rg1])
  m2 <- median(ys[rg2])
  m3 <- median(ys[rg3])
  
  m2 > m1 * 1.25 && m2 > m3 * 1.25
}



#' Calculate FWHM
#' 
#' Do not alter \code{length(ys)} since need \code{ista} and \code{iend}.
#' 
#' @param ts A vector of T values.
#' @param ys A vector of Y values.
#' @param yco The cut-off in Y values.
#' @param min_n The minimum number of points in \code{ys} for consideration.
calcFWHM <- function (ys, ts, yco = 100, min_n = 15L)
{
  # len0 <- length(ys)
  ymin <- min(ys, na.rm = TRUE) # baseline later...
  len  <- length(ys)
  if (len < min_n) { return(NULL) }
  imax <- .Internal(which.max(ys))

  if (imax == 1L || imax == len) { # <= 3L; imax + 3L >= len
    return(NULL)
  }

  ymax <- ys[[imax]]
  hmax <- (ymax + ymin) / 2
  ioks <- .Internal(which(ys >= hmax))
  noks <- length(ioks)
  
  if (noks < 5L) { # || noks / len < .2
    return(NULL)
  }
  
  # remove spikes outside of the half width
  im <- .Internal(which(ioks == imax))

  if (im == 1L || im == noks) {
    return(NULL)
  }

  ds1   <- diff(ioks[1:im])
  ds2   <- diff(ioks[im:noks])
  ioks1 <- .Internal(which(ds1 < 5L))
  ioks2 <- .Internal(which(ds2 < 5L)) + im
  ioks  <- ioks[c(ioks1, im, ioks2)]
  noks  <- length(ioks)
  
  if (!noks) {
    return(NULL)
  }

  iend <- ioks[[noks]]
  ista <- ioks[[1]]
  fwhm <- ts[[iend]] - ts[[ista]]
  
  if (fwhm >= 25) {
    return(NULL)
  }
  
  if (fwhm / (ts[[len]] - ts[[1]]) > .667) { # blunt; was .8
    return(NULL)
  }

  # stopifnot(len == len0)
  
  list (fwhm = fwhm, ista = ista, iend = iend)
}


#' Find the baseline values of Y
#'
#' @param vals A vector of Y values.
#' @param vmax The maximum of vals.
#' @param perc The percentage to the base-peak intensity for cut-offs.
#' @param lwr A lower bound of intensities for defining baseline.
#' @param max_perc The maximum percentage of baseline levels.
#' @return A conservative baseline value for finding signal edges.
find_baseline <- function (vals, vmax, perc = .02, max_perc = .05, lwr = 2.0)
{
  # exclude plateau and paddings
  vs <- vals[vals > lwr]
  vmin <- min(vs)
  upr <- vmax * .25 + vmin
  vs  <- vs[(!is.na(vs)) & vs <= upr]
  nvs <- length(vs)
  
  if (nvs < 5L) { # no good baseline
    return(NULL)
    # return(vmax * perc)
    # return(bx + vmin)
  }
  
  # may check that vs are populated mostly on either ends of vals...

  bx   <- vmax * max_perc

  # some vs remains part of peaks but `diff(vs)` are expected to be small 
  min(mean(abs(diff(vs))) * 3 + vmin, bx)
}


#' Find peak edges by intensity threshold relative to the base peak
#'
#' @param ys A vector of intensity values.
#' @param sta The starting point (off-set) of the current ys (ys[ri]) in the
#'   whole ys.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param y_rpl A replacement value of intensities.
#' @param perc The percentage to the base-peak intensity for cut-offs.
#' @param min_n The minimum number of \code{ys} to define an LC gate.
find_lc_gates2 <- function (perc = .02, ys, sta, n_dia_scans = 6L, y_rpl = 2.0, 
                            min_n = 10L) # was 6L
{
  ymin <- max(ys, na.rm = TRUE) * perc
  ys[ys <= ymin] <- NA_real_
  ys <- fill_lc_gaps(ys, n_dia_scans = n_dia_scans, y_rpl = 2.0)
  
  ioks  <- .Internal(which(ys > 0))
  edges <- find_gate_edges(ioks)
  ups   <- edges[["ups"]] # in relative to ys
  dns   <- edges[["dns"]]
  nps   <- length(ups)
  
  # discard narrow gates
  if (FALSE) {
    if (nps) {
      oks <- dns - ups + 1L >= min_n
      ups <- ups[oks]
      dns <- dns[oks]
    }
  }
  
  if (!nps) {
    return(list(stas = NULL, ends = NULL, dn1 = NULL, upn = NULL, 
                first_bad = FALSE, last_bad = FALSE))
  }
  
  ###
  first_bad <- FALSE
  last_bad  <- FALSE
  stx <- sta - 1L
  upn <- ioks[[ups[[nps]]]] + stx
  dn1 <- ioks[[dns[[1]]]] + stx
  ###
  
  if (nps == 1L) {
    lenx <- sum(ys[ioks[ups:dns]] > y_rpl)
  }
  else {
    # rngs <- mapply(`:`, ups, dns, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    # lenx <- lapply(rngs, function (rng) sum(ys[ioks[rng]] > y_rpl))
    
    ###
    lenx <- vector("integer", nps)
    yhs  <- vector("numeric", nps) # peak heights
    rngs <- vector("list", nps)
    
    for (i in 1:nps) {
      rngs[[i]] <- rgi <- ups[[i]]:dns[[i]]
      ysi <- ys[ioks[rgi]]
      lenx[[i]] <- sum(ysi > y_rpl)
      yhs[[i]]  <- max(ysi, na.rm = TRUE)
    }
    
    yhs <- yhs - ymin # subtract the baseline
    yhm <- max(yhs)
    
    if (first_bad <- yhm / yhs[[1]] > 5) { # up-bending partial peak or small hump
      lenx[[1]] <- 0L
    }
    
    if (last_bad <- yhm / yhs[[nps]] > 5) {
      lenx[[nps]] <- 0L
    }
    ###
  }
  
  if (length(bads <- .Internal(which(lenx < min_n)))) {
    ups <- ups[-bads]
    dns <- dns[-bads]
  }
  
  # sta <- ioks[[upi]] # the starting index of gate-i in relative to ys
  # stas[[i]] <- sta + ioks[ups] - 1L # in relative to ys
  # ends[[i]] <- sta + ioks[dns] - 1L
  
  if (length(ups)) {
    list(stas = stx + ioks[ups], ends = stx + ioks[dns], 
         dn1 = dn1, upn = upn, first_bad = first_bad, last_bad = last_bad)
  }
  else {
    list(stas = NULL, ends = NULL, dn1 = NULL, upn = NULL, 
         first_bad = FALSE, last_bad = FALSE)
  }
}


#' Calculate area under a curve
#'
#' @param ys Y values around an LC peak. Should contain no NA values in the
#'   sequence.
#' @param ts Time values corresponding to \code{ys}.
#' @param rng The range indexes of \code{ys}.
#' @param yco An intensity cut-off.
#' @param ytot_co A cut-off in peak area.
#' @param err An error tolerance of retention times between max and median (in
#'   seconds).
#' @param min_n The minimum points across the peak profile of \code{ys} for
#'   consideration.
calcAUC <- function (ys, ts, rng, yco = 100, ytot_co = 2E5, min_n = 15L, 
                     err = 2.0)
{
  ## (1) find the apex index
  len  <- length(ys)
  
  if (len < min_n) {
    return(NULL) 
  }
  
  imax <- .Internal(which.max(ys))
  dts  <- diff(ts)
  csum <- cumsum(ys * c(dts, .5))
  tot  <- csum[[len]]
  mval <- tot / 2
  imed <- .Internal(which(csum >= mval))[[1]]
  
  ###
  # not good: need to first handle spikes...
  # rmed <- imed / len # may be spike
  # if (rmed <= .25 || rmed >= .75 || imed <= 3L || imed + 3L >= len || 
  #     imax <= 3L || imax + 3L >= len) { return(NULL) }
  ###
  
  tmed <- ts[[imed]]
  tmax <- ts[[imax]]
  dt1  <- tmed - tmax
  bigt <- abs(dt1) > err
  
  ##  handle spikes
  if (bigt) {
    count <- 3L
    vmed  <- mean(ys[max(1L, imed - 2L):min(len, imed + 2L)])
    vmax  <- mean(ys[max(1L, imax - 2L):min(len, imax + 2L)])
    
    while((vmax < vmed) && count) {
      ys[[imax]] <- vmax
      imax <- .Internal(which.max(ys))
      tmax <- ts[[imax]]
      dt1  <- tmed - tmax
      bigt <- abs(dt1) > err
      vmax <- mean(ys[max(1L, imax - 2L):min(len, imax + 2L)])
      
      count <- count - 1L
    }
  }
  
  ## (3) refine the peak range
  sk  <- if (bigt) if (dt1 > 0) 1L else -1L else 0L
  
  if (tot < ytot_co) {
    return(NULL)
  }
  
  ans_fw <- calcFWHM(ys, ts, yco = 100)
  fwhm   <- ans_fw[["fwhm"]]
  if (is.null(fwhm)) { return(NULL) } # flat peak
  ista   <- ans_fw[["ista"]]
  iend   <- ans_fw[["iend"]]

  if (sk == 1L) {
    w1 <- fwhm * 1.2
    w2 <- fwhm * 1.6
  }
  else if (sk == -1L) {
    w1 <- fwhm * 1.6
    w2 <- fwhm * 1.2
  }
  else {
    w1 <- fwhm * 1.2
    w2 <- fwhm * 1.2
  }
  
  ioks <- .Internal(which(ts >= (tmax - w1) & ts <= (tmax + w2)))
  iok1 <- ioks[[1]]
  idx  <- imax - iok1 + 1L
  len  <- length(ioks)
  iokn <- ioks[[len]]
  
  # to replace diff(ts) by a subset of dts or subtract csum ...
  tot <- cumsum(ys[ioks] * c(diff(ts[ioks]), .5))[[len]]

  if (tot < ytot_co) {
    return(NULL)
  }

  ok_convex <- check_peak_convex(
    ys = ys, i_sta = iok1, i_midsta = ista, i_midend = iend, i_end = iokn)

  if (!ok_convex) {
    return(NULL) 
  }
  
  return(list(area = tot, idx = idx, rng = rng[ioks], fwhm = fwhm))

  
  
  # calculate peak area (may or may not use `tot`)
  oks <- .Internal(which(ys >= ys[[imax]] * .03))
  ysub <- ys[oks]
  tsub <- ts[oks]
  area <- sum(ysub * c(diff(tsub), .5), na.rm = TRUE)
  
  list(area = area, idx = idx)
}


#' Fills LC gaps.
#'
#' @param ys A vector of intensity values.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param y_rpl A replacement value to fill the gaps of \code{ys}.
#' @examples
#' ys <- c(10,0,0,0,11,0,15,0,0,12,0,10,0,0,10)
#' mzion:::fill_lc_gaps(ys, 3)
#' mzion:::fill_lc_gaps(ys, 2)
#'
#' ys <- c(10,0,0,0,11,15,15,0,0,12,0,10,0,0,10)
#' mzion:::fill_lc_gaps(ys, 3)
#' mzion:::fill_lc_gaps(ys, 2)
fill_lc_gaps <- function (ys, n_dia_scans = 6L, y_rpl = 2.0)
{
  if (n_dia_scans <= 1L) {
    return(ys)
  }

  xs <- .Internal(which(ys > 0))
  ds <- diff(xs)
  ps <- .Internal(which(ds <= n_dia_scans & ds > 1L))
  ng <- length(ps)
  
  if (!ng) {
    return(ys)
  }

  for (i in seq_along(ps)) {
    p <- ps[[i]]
    rng <- (xs[[p]] + 1L):(xs[[p + 1L]] - 1L)
    ys[rng] <- y_rpl
  }
  
  ys
}


#' Check partial peaks
#'
#' @param ys A vector of intensity values
#' @param width The width of a peak (number of points).
#' @param min_p The minimum number of points on one side of a peak for excluding
#'   the call for a partial peak.
is_partial_peak <- function (ys, width = 200L, min_p = 10L)
{
  # "no-left/right" partial peaks; probably due to space-charge effects
  # does it handle abreast peaks? may skip this and handled by space charge...
  
  ys <- ys[!is.na(ys)]
  ys <- ys[ys >= max(ys) * .05] # in case of skewing tails
  # sum(diff(ys) > 0) / length(ys) < .3 # not working: local wobbling of ups and downs
  imax <- which.max(ys)
  nyi  <- length(ys)
  
  nyi > width && (imax < min_p || nyi < imax + min_p)
}


#' Make MS1 X and Y matrices across adjacent scans
#'
#' Allow adjacent values in \code{unv} and later collapse adjacent
#' columns/values.
#'
#' @param xs Vectors of m-over-z values in an ASCENDING order.
#' @param ys Vectors of intensity values corresponding to \code{xs}.
#' @param lwr The lower mass limit.
#' @param step The bin size in converting numeric m-over-z values to integers.
#' @param reord Logical; re-order data or not.
#' @param cleanup Logical; to perform data row clean-ups or not. Rows may drop
#'   at \code{cleanup = TRUE}.
#' @param sum_y Logical; to sum Y values or not. Either TRUE or FALSE for
#'   Thermo's and TRUE for collapsing PASEF MS2 slices.
#' @param coll Logical; collapse adjacent columns or not.
#' @param add_colnames Logical; if TRUE, add the indexes of mass bins to the
#'   column names of \code{matx}.
#' @param look_back Logical; look up the preceding MS bin or not.
#' @importFrom fastmatch %fin%
#' @examples
#' # Twos adjacent bins of xs: 392.1796, 392.1845
#' # xs longer than the universe
#' xs <- list(c(391.1883,391.2848,391.6627,391.6995,392.1646,392.1796,392.1845,
#'            392.2030,392.2887,392.6641,392.7812,393.0833,393.2975))
#' ys <- list(c(12827.41,337002.19,617819.69,18045.10,205851.53,15194.98,11318.61,
#'              12970.02,118604.48,75726.89,11676.51,23723.18,55749.93))
#' mzion:::collapse_mms1ints(xs, ys, lwr = 389.6529)
#'
#' xs <- list(c(400.6596,401.7076,402.1813,402.1944,402.1969,402.2094,402.5438,402.7112,403.1812,404.1777),
#'            c(400.6599,401.7075,402.1954,402.1975,402.7112,403.1822,404.2777))
#' ys <- list(c(24003.98,53431.96,110619.26,10988.55,12291.00,140045.06,67601.16,11413.04,21651.61,16686.06),
#'            c(10000.1,40000.1,20000.1,50000.1,2500.2,5000.1,30000.1))
#' mzion:::collapse_mms1ints(xs, ys, lwr = 400.1994)
#'
#' xs <- ys <- vector("list", 13L)
#' xs[[7]] <- 954.607849; xs[[8]] <- 954.630249; xs[[10]] <- 954.622925
#' ys[[7]] <- 15706.2627; ys[[8]] <- 19803.5879; ys[[10]] <- 31178.9648
#' mzion:::collapse_mms1ints(xs, ys, lwr = 951.089731)
collapse_mms1ints <- function (xs = NULL, ys = NULL, lwr = 115L, step = 1e-5, 
                               reord = FALSE, cleanup = FALSE, sum_y = FALSE, 
                               coll = TRUE, add_colnames = FALSE, 
                               look_back = FALSE)
{
  ### 
  # the utility is often called heavily;
  # DO NOT gc() that will slow things down
  ### 
  
  # mostly FALSE; otherwise cause drops in the number of data entries
  if (cleanup) {
    null_out <- list(x = NULL, y = NULL, u = NULL)
    
    # 1. all xs are NULL
    if (!any(oks <- lengths(xs) > 0L)) {
      return(null_out)
    }
    
    oks <- .Internal(which(oks))
    xs <- xs[oks]
    ys <- ys[oks]
    
    # 2. remove zero intensities
    oky <- lapply(ys, `>`, 0)
    
    for (i in seq_along(xs)) {
      oki <- .Internal(which(oky[[i]]))
      xs[[i]] <- xs[[i]][oki]
      ys[[i]] <- ys[[i]][oki]
    }
    
    # 2.1 checks again
    if (!any(oks <- lengths(xs) > 0L)) {
      return(null_out)
    }
    
    oks <- .Internal(which(oks))
    xs <- xs[oks]
    ys <- ys[oks]
  }
  
  if (reord) {
    lens <- lengths(xs)
    
    for (i in seq_along(xs)) {
      xi <- xs[[i]]
      
      if (lens[[i]] > 1L) {
        ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, 
                                   TRUE, xi))
        xs[[i]] <- xi[ord]
        ys[[i]] <- ys[[i]][ord]
      }
    }
  }
  
  ixs <- lapply(xs, index_mz, lwr, step)
  
  # 1. remove duplicated ixs and collapse the Y values under the same ixs
  #    Y values are only for de-isotoping, not for precursor intensities
  # 2. `sum_y` seems have little effect with Thermo's but PASEF
  if (sum_y) {
    for (i in seq_along(xs)) {
      ix <- ixs[[i]]
      x  <- xs[[i]]
      y  <- ys[[i]]
      ps <- .Internal(which(duplicated(ix)))
      
      if (l <- length(ps)) {
        for (j in 1:l) {
          okp <- .Internal(which(ix == ix[ps[[j]]]))
          y[okp[[1]]] <- sum(y[okp], na.rm = TRUE)
        }
        
        ixs[[i]] <- ix[-ps]
        xs[[i]]  <- x[-ps]
        ys[[i]]  <- y[-ps]
      }
    }
  }
  else {
    for (i in seq_along(xs)) {
      ix <- ixs[[i]]
      x  <- xs[[i]]
      y  <- ys[[i]]
      ps <- .Internal(which(duplicated(ix)))
      
      if (length(ps)) {
        ixs[[i]] <- ix[-ps]
        xs[[i]]  <- x[-ps]
        ys[[i]]  <- y[-ps]
      }

      # if (any(dups <- duplicated(ix))) {
      #   oks      <- .Internal(which(!dups))
      #   ixs[[i]] <- ix[oks]
      #   xs[[i]]  <- x[oks]
      #   ys[[i]]  <- y[oks]
      # }
    }
  }
  # rm(list = c("x", "y", "ix", "ps"))
  
  ## maps ixs vectors to unv (presence or absence)
  unv  <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv  <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups  <- lapply(ixs, function (x) unv %fin% x)
  
  # note one-to-one correspondence between ixs and xs
  xmat <- mapcoll_xyz(vals = xs, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)
  ymat <- mapcoll_xyz(vals = ys, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)
  
  # if (!coll) { return(list(x = xmat, y = ymat, u = unv)) }

  # collapses adjacent entries
  ps   <- find_gates(unv)
  lenp <- length(ps)
  
  # which(sapply(ps, `[[`, 1L) == which(unv == 167269)) # 6325
  # which(sapply(ps, `[[`, 2L) == which(unv == 167269))
  
  # all discrete values
  if (is.null(ps)) {
    return(list(x = xmat, y = ymat, u = unv))
  }
  
  # for recording matrix columns to be dropped
  ps1 <- vector("integer", lenp)
  
  # collapses matrix columns with +/-1 in bin indexes
  for (i in 1:lenp) {
    # test from 6323:6325
    # i= 6324; i = 6325
    c12 <- ps[[i]]
    c1  <- c12[[1]]
    c2  <- c12[[2]]
    
    # with values in both columns, X: 2 <- 1; Y: 1 + 2; c2 can be 0
    if (c2 == 0L) {
      next
    }
    
    if (FALSE) {
      xm1 <- xmat[, c1]
      xm2 <- xmat[, c2]
      rows1 <- is.na(xm1)
      rows2 <- is.na(xm2)
      rows  <- .Internal(which(!rows1 & rows2))
      xmat[rows, c2] <- xm1[rows]
      # Y values on both rows1 and rows2 will not be collapsed
      ymat[rows, c2] <- rowSums(ymat[rows, c1:c2, drop = FALSE], na.rm = TRUE)
      xmat[rows, c1] <- NA_real_ # to prevent value traversing in fwd and bwd looking
      ymat[rows, c1] <- NA_real_
    }
    
    rows <- .Internal(which(!is.na(xmat[, c1])))
    xmat[rows, c2] <- xmat[rows, c1]
    ymat[rows, c2] <- ymat[rows, c1] # rowSums(ymat[rows, c1:c2, drop = FALSE], na.rm = TRUE)
    # no need since columns traverse from left to right and c1 will be dropped
    xmat[rows, c1] <- NA_real_ # toggle on; to prevent value traversing in fwd and bwd looking
    ymat[rows, c1] <- NA_real_ # toggle on
    
    if (look_back && i < lenp && (af <- ps[[i+1]][[1]]) == (c2 + 1L)) {
      if (FALSE) { # disadvantageous on Y values
        xc2   <- xmat[, c2]
        xaf   <- xmat[, af]
        rows1 <- is.na(xc2)
        rows2 <- is.na(xaf)
        rows  <- .Internal(which(rows1 & !rows2))
        xmat[rows, c2] <- xaf[rows]
        # Y values on both rows1 and rows2 will not be collapsed
        ymat[rows, c2] <- rowSums(ymat[rows, c2:af, drop = FALSE], na.rm = TRUE)
      }

      ###
      # transfer only if tangent...
      if (FALSE) {
        rows2 <- .Internal(which(is.na(xmat[, c2]) & !is.na(xmat[, af])))
        # ys1 <- ymat[rows2, c2]
        ys2 <- ymat[rows2, af]
        gs2 <- find_lc_gates(ys = ymat[, af], n_dia_scans = 0L, y_rpl = 0)
        
        local({
          ymin <- max(ys, na.rm = TRUE) * perc
          ys[ys <= ymin] <- NA_real_
          ys <- fill_lc_gaps(ys, n_dia_scans = n_dia_scans, y_rpl = 2.0)
          
          ioks  <- .Internal(which(ys > 0))
          edges <- find_gate_edges(ioks)
          ups   <- edges[["ups"]] # in relative to ys
          dns   <- edges[["dns"]]
          nps   <- length(ups)
        })

        # get the ranges
        
        # find tangent ranges -> c2
        
        ys2[4:43]
        rows2[4:43] # 399:439
      }
      
      ###
      
      rows2 <- .Internal(which(is.na(xmat[, c2])))
      xmat[rows2, c2] <- xmat[rows2, af]
      ymat[rows2, c2] <- ymat[rows2, af] # rowSums(ymat[rows2, c2:af, drop = FALSE], na.rm = TRUE)
      # about ok to enable, but may be unnecessary
      xmat[rows2, af] <- NA_real_ # toggle on
      ymat[rows2, af] <- NA_real_ # toggle on
    }
    
    # colnames(xmat) <- colnames(ymat) <- unv; zy <- ymat[, c12]; plot(zy[, 2])
    ps1[[i]] <- c1
  }
  
  # note `-ps1` ok in that at least one ps1 is not 0
  # identical(xmat[, -c(0, 2:3), drop = FALSE], xmat[, -c(2:3), drop = FALSE])
  # !identical(xmat[, 0, drop = FALSE], xmat)
  
  # zy2 <- ymat[, c(11415, 11416)]
  
  xmat <- xmat[, -ps1, drop = FALSE]
  ymat <- ymat[, -ps1, drop = FALSE]
  unv  <- unv[-ps1]
  
  # possible additional all-NA columns by look_back
  if (look_back && 
      length(nas <- .Internal(which(colSums(is.na(xmat)) == lenx)))) {
    xmat <- xmat[, -nas, drop = FALSE]
    ymat <- ymat[, -nas, drop = FALSE]
    unv  <- unv[-nas]
  }
  
  list(x = xmat, y = ymat, u = unv)
}


#' Collapse MS1 matrices of X and Y
#'
#' @param matx A matrix of m-over-z values.
#' @param maty A matrix of intensity values.
#' @return A list of x: the weighted means of m-over-z values; y: the sum of
#'   intensities over flanking MS1s; n: the numbers of observed precursors.
calc_ms1xys <- function (matx, maty)
{
  nr <- nrow(matx)
  
  if (nr == 1L) {
    ysums <- maty[1, ]
  }
  else {
    ysums <- colSums(maty, na.rm = TRUE)
  }

  xmeans <- colSums(matx * maty, na.rm = TRUE) / ysums
  ns <- colSums(!is.na(maty), na.rm = TRUE)
  
  # use tallied ysums, not ysums/ns, for deisotoping
  list(x = xmeans, y = ysums, n = ns)
}


#' Finds mDDA precursors
#'
#' Averages of multiple MS1 scans.
#'
#' @param msx_moverzs Vectors of bracketing (e.g. +/-6 scans) full-spectrum MS1
#'   m-over-z values. Each \code{msx_moverzs} needs to be in an ascending order.
#' @param msx_ints Vectors of bracketing full-spectrum MS1 intensity values
#'   corresponding to \code{msx_moverzs}.
#' @param iso_ctr Vectors of isolation centers (e.g. top-12 scans) of MS2.
#' @param iso_lwr Vectors of isolation lowers of MS2.
#' @param iso_upr Vectors of isolation uppers of MS2. Both \code{iso_lwr} and
#'   \code{iso_upr} are used, in case of asymmetrical isolation to
#'   \code{iso_ctr}.
#' @param ppm Mass error tolerance.
#' @param step A step size.
#' @param maxn_precurs Maximum number of precursors for consideration.
#' @param max_ms1_charge Maximum charge state of precursors for consideration.
#' @param width_left The left-width of an extended MS1 isolation window to
#'   contain a full isotope envelops.
#' @param width_right The right-width of an extended MS1 isolation window to
#'   contain a full isotope envelops.
#' @param margin The margin of an m-over-z extended to an isolation window (to
#'   account for the imperfection of the voltage gates of isolation).
#' @param use_defpeaks Use default peak info or not.
#' @param look_back Logical; look up the preceding MS bin or not.
#' @param is_pasef Logical; is TIMS TOF data or not.
#' @param min_y Not yet used. The minimum intensity for consideration of
#'   deisotoping (guard against PASEF data). The setting should be MS platform
#'   dependent, e.g., a much smaller value with PASEF. Also note different
#'   settings between MS1 and MS2.
#' @inheritParams matchMS
#' @inheritParams find_ms1stat
#' @return A list. x: monoisotopic moverzs (weighted mean statistics); y:
#'   intensities (mean); z: charge states.
find_mdda_mms1s <- function (msx_moverzs = NULL, msx_ints = NULL, 
                             iso_ctr = NULL, iso_lwr = NULL, iso_upr = NULL, 
                             ppm = 10L, maxn_precurs = 5L, max_ms1_charge = 4L, 
                             n_fwd = 15L, n_bwd = 20L, offset_upr = 20L, 
                             offset_lwr = 30L, grad_isotope = 1.6, 
                             fct_iso2 = 3.0, use_defpeaks = FALSE, 
                             width_left = 2.01, width_right = 1.25, margin = .5, 
                             min_mass = 200L, step = ppm/1e6, 
                             look_back = FALSE, is_pasef = FALSE, min_y = 10, 
                             ms1_min_ratio = 0)
{
  if (!(len1 <- length(msx_moverzs))) {
    return(NULL)
  }
  
  if (!(len2 <- length(iso_ctr))) {
    return(NULL)
  }
  
  # ansx1 etc.: (6+1+6) MS1 frames subset by one MS2 isolation window
  # ansx2 etc.: collapsed MS1 within an isolation window for each MS2
  # ymats: matrices of MS1 intensity. Columns: masses; rows: MS1 frames
  # ymat0: delayed MS1 intensities for each MS1 frame; 
  #   the information is only available after the determination of mono m/z
  
  ansx1 <- ansy1 <- vector("list", len1)
  xs <- ys <- zs <- vector("list", len2)
  # ymats <- vector("list", len2)
  # ymat0 <- replicate(len2, rep_len(NA_real_, len1), simplify = FALSE)
  # ymat0 <- matrix(rep_len(NA_real_, len1 * len2), ncol = len2)
  
  ## by MS2 entries
  for (i in 1:len2) {
    m2  <- iso_ctr[[i]]
    lwi <- iso_lwr[[i]]
    upi <- iso_upr[[i]]
    lwr <- lwi - width_left
    upr <- upi + width_right
    
    # 1. gather e.g. +/-6 MS1s
    for (j in 1:len1) {
      x1s <- msx_moverzs[[j]]
      y1s <- msx_ints[[j]]
      oks <- .Internal(which(x1s > lwr & x1s < upr))
      ansx1[[j]] <- x1s[oks]
      ansy1[[j]] <- y1s[oks]
    }

    # 2. generate flanking matrices of MS1 X and Y
    # xs must be in an ascending order
    # rows: by flanking MS1 scans; columns: by moverzs bin indexes
    xys <- collapse_mms1ints(
      xs = ansx1, ys = ansy1, lwr = min_mass, step = step, reord = FALSE, 
      cleanup = FALSE, add_colnames = FALSE, look_back = FALSE)
    
    # 3. collapse MS1s for de-isotoping
    ans <- calc_ms1xys(xys[["x"]], xys[["y"]])
    ax  <- ans[["x"]]
    if (is.null(ax)) { next }
    ay <- ans[["y"]]
    an <- ans[["n"]]
    
    # 4. de-isotoping
    mic <- find_ms1stat(
      moverzs = ax, msxints = ay, n_ms1s = an, 
      center = m2, exclude_reporter_region = FALSE, ppm = ppm, ms_lev = 1L, 
      maxn_feats = maxn_precurs, max_charge = max_ms1_charge, n_fwd = n_fwd, 
      n_bwd = n_bwd, offset_upr = offset_upr, offset_lwr = offset_lwr, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks, is_pasef = is_pasef, 
      min_y = min_y, ms1_min_ratio = ms1_min_ratio)
    
    masses <- mic[["masses"]]
    if (!length(masses)) next
    intensities <- mic[["intensities"]]
    charges <- mic[["charges"]]

    # 5. subset by isolation window
    if (TRUE) {
      # mono-isotopic peaks are left-biased
      # masses left to the iso-window may still be offseted from the mono-isotopic 
      oks <- masses <= upi + margin 
      oks <- .Internal(which(oks))
      len <- length(oks)

      # may removes masses[[i]] left to the iso_lwr[[i]] that has little 
      # evidence within the isolation window
      
      # length(xs) drops by 1 if is.null(masses)
      if (len && len < length(masses)) {
        xs[[i]] <- masses[oks]
        ys[[i]] <- intensities[oks]
        zs[[i]] <- charges[oks]
      }
      else {
        # accepts all, e.g., all masses are greater than the upper bound
        # may consider the supermum/infimum to the iso_upr.
        xs[[i]] <- masses
        ys[[i]] <- intensities
        zs[[i]] <- charges
      }
    }
    else {
      xs[[i]] <- masses
      ys[[i]] <- intensities
      zs[[i]] <- charges
    }
    
    # PASEF reorder: 
    # xs: first within the isoWindow: from high to low intensity
    #  then outsie the isoWindow from high to low intensity
    if (is_pasef) {
      oksa <- masses >= lwi & masses <= upi
      oksb <- .Internal(which(!oksa))
      oksa <- .Internal(which(oksa))
      
      if (lena <- length(oksa)) {
        xsa <- masses[oksa]
        ysa <- intensities[oksa]
        zsa <- charges[oksa]
        
        if (lena > 1L) {
          orda <- .Internal(radixsort(na.last = TRUE, decreasing = TRUE, FALSE, 
                                      TRUE, ysa))
          xsa <- xsa[orda]
          ysa <- ysa[orda]
          zsa <- zsa[orda]
        }
      }
      else {
        xsa <- ysa <- zsa <- NULL
      }
      
      if (lenb <- length(oksb)) {
        xsb <- masses[oksb]
        ysb <- intensities[oksb]
        zsb <- charges[oksb]
        
        if (lenb > 1L) {
          ordb <- .Internal(radixsort(na.last = TRUE, decreasing = TRUE, FALSE, 
                                     TRUE, ysb))
          xsb <- xsb[ordb]
          ysb <- ysb[ordb]
          zsb <- zsb[ordb]
        }
      }
      else {
        xsb <- ysb <- zsb <- NULL
      }
      
      xs[[i]] <- c(xsa, xsb)
      ys[[i]] <- c(ysa, ysb)
      zs[[i]] <- c(zsa, zsb)
    }
  }
  
  list(x = xs, y = ys, z = zs)
}


#' Finds precursors using MS2 data.
#' 
#' For absolutely no precursor identification based on MS1 data.
#' 
#' @param charges MS2 charges.
#' @param iso_lwr The lower width of an isolation window.
#' @param iso_upr The upper width of an isolation window.
#' @inheritParams find_ms1stat
find_ms1byms2 <- function (moverzs = NULL, msxints = NULL, charges = NULL, 
                           center, iso_lwr, iso_upr)
{
  na_out <- list(ms1_moverzs = NA_real_, ms1_masses = NA_real_, 
                 ms1_charges = NA_integer_, ms1_ints = NA_real_)
  
  if (!(len_ms <- length(moverzs))) {
    return(na_out)
  }

  oks <- (moverzs >= iso_lwr) & (moverzs <= iso_upr)
  oks <- .Internal(which(oks))
  moverzs <- moverzs[oks]
  
  if (!(len_ms <- length(moverzs))) {
    return(na_out)
  }

  msxints <- msxints[oks]
  charges <- charges[oks]
  
  charges[is.null(charges)] <- NA_integer_
  charges[is.na(charges)]   <- 2L # arbitrary
  masses <- (moverzs - 1.00727647) * charges
  
  if (len_ms == 1L) {
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  }

  okc <- !is.na(charges)
  okc <- .Internal(which(okc))
  charges <- charges[okc]
  
  if (!(len_ms <- length(charges))) {
    return(na_out)
  }

  moverzs <- moverzs[okc]
  msxints <- msxints[okc]
  masses  <- masses[okc]
  
  if (len_ms == 1L) {
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  }

  idx <- .Internal(which.max(msxints))
  moverzs <- moverzs[idx]
  charges <- charges[idx]
  msxints <- msxints[idx]
  masses  <- masses[idx]
  
  list(ms1_moverzs = moverzs, ms1_masses = masses, 
       ms1_charges = charges, ms1_ints = msxints)
}


