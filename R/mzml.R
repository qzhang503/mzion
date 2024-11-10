#' Helper in preparing peaklists
#' 
#' @param filelist A list of peaklist files.
#' @param data_type A data type of either RAW or mzML.
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
                      ppm_ms1 = 10L, ppm_ms2 = 10L, 
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
  qs::qsave(data_type, file.path(mgf_path, "data_type.rds"), preset = "fast")
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
    is_dia <- attr(peakfile1, "is_dia", exact = TRUE)
    mzml_type <- attr(peakfile1, "mzml_type", exact = TRUE)
    peakfiles <- unlist(peakfiles)
  }
  else {
    peakfiles <- list.files(file.path(mgf_path, "temp_dir"), pattern = "\\.d\\.rds$")
    # peakfiles <- list.files(file.path(mgf_path, "temp_dir"), pattern = "\\.raw\\.rds$")
    peakfiles <- peakfiles[!grepl("^predeisoDDA_", peakfiles)]
    peakfiles <- peakfiles[!grepl("^pasefms1_", peakfiles)]
    # peakfiles <- peakfiles[grepl("^LFQ_", peakfiles)]
    is_dia <- FALSE
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
  n_cores <- max(min(n_pcs, ceiling(rams / if (is_pasef) 30 else 5), lenf), 1L) # 20 -> 30
  n_para  <- max(min(n_pcs, round(n_pcs / n_cores)), 1L)
  
  if (is_pasef) {
    if (n_mdda_flanks) {
      n_mdda_flanks <- 0L
      warning("Coerce to `n_mdda_flanks = 0` for PASEF data.")
    }
    
    multi_pasef <- lenf > n_cores
    
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
  else {
    file_indexes = seq_along(peakfiles)
    multi_pasef <- FALSE
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
    yco <- if (is_pasef) 10 else 100
    
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
            ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
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
            ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
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
        # if (is_pasef) { n_para <- max(floor(min(n_para, 16L) / n_cores), 1L) }
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
            ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
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
  raws <- unlist(raws, recursive = FALSE, use.names = TRUE)
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
#' @inheritParams deisoDDA
#' @inheritParams matchMS
hdeisoDDA <- function (filename, raw_id = 1L, 
                       out_path = NULL, mgf_path = NULL, temp_dir = NULL, 
                       mzml_type = "raw", ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       maxn_mdda_precurs = 5L, 
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
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
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
#' @inheritParams matchMS
deisoDDA <- function (filename = NULL, 
                      out_path = NULL, mgf_path = NULL, temp_dir = NULL, 
                      mzml_type = "raw", ppm_ms1 = 10L, ppm_ms2 = 10L, 
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
                      use_lfq_intensity = TRUE, ppm_ms1trace = 5L, 
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
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2, 
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
    
    if (FALSE) {
      dfx <- df |> dplyr::filter(orig_scan == 18208.14)
      data <- data.frame(x = dfx$msx_moverzs[[1]], y = dfx$msx_ints[[1]])
      data |>
        # dplyr::filter(x >= 690, x <= 699) |>
        ggplot2::ggplot() + 
        ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                              color = "gray", linewidth = .1)
    }
  }
  
  ## LFQ: replaces intensities with apex values
  # No MS1 info at maxn_mdda_precurs == 0L
  if (use_lfq_intensity && maxn_mdda_precurs) {
    # for MBR in proteoQ
    path_ms1 <- create_dir(file.path(out_path, "ms1data"))
    
    qs::qsave(df[with(df, ms_level == 1L), 
                 c("ret_time", "scan_num", "orig_scan", "msx_moverzs", "msx_ints")], 
              file.path(path_ms1, paste0("ms1full_", filename)), preset = "fast")

    df <- get_ms1xs_space(df)
    
    if (FALSE) {
      qs::qsave(df[with(df, ms_level == 1L), 
                   c("ret_time", "scan_num", "orig_scan", "ms1_moverz", "ms1_int", 
                     "ms1_charge")], 
                file.path(path_ms1, paste0("ms1subspace_", filename)), 
                preset = "fast")
    }

    step <- ppm_ms1trace * 1e-6
    rows1 <- which(df$ms_level == 1L)
    rt_gap <- estimate_rtgap(df$ret_time[rows1])

    # Set aside MS2 scans before the first MS1 (more likely with timsTOF)
    dfx <- if (row1  <- rows1[[1]] - 1L) df[1:row1, ] else NULL
    rm(list = c("rows1", "row1"))
    
    ans_prep <- pretraceXY(
      df[, c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", "ms_level", 
             "msx_moverzs", "msx_ints", "msx_charges", "orig_scan", "ret_time")], 
      from = min_mass, step = step, 
      # 1024L: more RAM, same speed
      n_chunks = ceiling(sum(df$ms_level == 1L) / 512L), 
      gap = rt_gap)
    dfs <- ans_prep$dfs
    df1s <- ans_prep$df1s
    gaps <- unlist(ans_prep$gaps, use.names = FALSE, recursive = FALSE)
    lenv <- length(df1s)
    gaps_bf <- ans_prep$gaps_bf <- c(0L, gaps[1:(lenv - 1L)])
    gaps_af <- ans_prep$gaps_af <- c(gaps[2:lenv], 0L)
    
    if (FALSE) {
      qs::qsave(ans_prep, file.path(path_ms1, paste0("msxspace_", filename)), 
                preset = "fast")
    }
    rm(list = c("ans_prep", "rt_gap"))

    cols <- c("ms_level", "ms1_moverz", "ms1_int", "orig_scan") # orig_scan for troubleshooting
    min_y <- if (is_pasef) 500 else 2e6

    if (TRUE) {
      cl <- parallel::makeCluster(getOption("cl.cores", 2L))
      out <- parallel::clusterMap(
        cl, htraceXY, 
        lapply(df1s, `[[`, "msx_moverzs"), 
        lapply(df1s, `[[`, "msx_ints"), 
        lapply(df1s, `[[`, "orig_scan"), 
        lapply(df1s, `[[`, "ret_time"), 
        lapply(dfs, `[`, cols), 
        gaps_bf, 
        gaps_af, 
        MoreArgs = list(
          n_dia_scans = n_dia_scans, from = min_mass, step = step, 
          y_perc = y_perc, yco = yco, look_back = TRUE, min_y = min_y
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
      parallel::stopCluster(cl)
    }
    else {
      vxs <- lapply(df1s, `[[`, "msx_moverzs")
      vys <- lapply(df1s, `[[`, "msx_ints")
      vss <- lapply(df1s, `[[`, "orig_scan") # for troubleshooting
      vts <- lapply(df1s, `[[`, "ret_time")
      vdf <- lapply(dfs, `[`, cols)
      out <- vector("list", lenv)

      for (i in 1:lenv) {
        out[[i]] <- htraceXY(
          xs = vxs[[i]], ys = vys[[i]], ss = vss[[i]], ts = vts[[i]], 
          df = vdf[[i]], gap_bf = gaps_bf[[i]], gap_af = gaps_af[[i]], 
          n_dia_scans = n_dia_scans, from = min_mass, step = step, 
          y_perc = y_perc, yco = yco, look_back = TRUE, min_y = min_y)
        
        if (FALSE) {
          if (all(lengths(vxs[[i]]) == 0L)) {
            li <- length(vxs[[i]]) - gaps_bf[[i]] - gaps_af[[i]] + 1L
            null <- rep_len(list(NULL), li)
            
            out[[i]] <- tibble::tibble(
              ms_level = vdf[[i]]$ms_level, 
              ms1_moverz = null, 
              ms1_int = null, 
              apex_scan_num = null)
          } else {
            out[[i]] <- htraceXY(
              xs = vxs[[i]], ys = vys[[i]], ss = vss[[i]], df = vdf[[i]], 
              gap_bf <- gaps_bf[[i]], gap_af = gaps_af[[i]], 
              # may be > n_mdda_flanks
              n_dia_scans = n_dia_scans, from = min_mass, step = step, 
              y_perc = y_perc, yco = yco)
          }
        }
      }
      rm(list = c("vxs", "vys", "vdf", "gaps", "lenv"))
    }

    out <- dplyr::bind_rows(out)
    
    if (!is.null(dfx)) {
      out <- dplyr::bind_rows(dfx[, cols], out)
    }
    
    if (nrow(df) == nrow(out)) {
      df[, cols] <- out[, cols]
    }
    else {
      stop("Developer: checks for row dropping.")
    }

    # is_list <- class(out$apex_scan_num) == "list" # must be
    df$apex_scan_num <- out$apex_scan_num
    rm(list = ("out"))
    
    df <- local({
      rows2 <- df$ms_level != 1L
      df2 <- df[rows2, ]
      empties <- lengths(df2$apex_scan_num) == 0L
      
      if (use_ms2scan_for_ms1apex <- TRUE) {
        df2$apex_scan_num[empties] <- as.list(df2$scan_num[empties])
      }
      else {
        df2$apex_scan_num[empties] <- as.list(NA_integer_)
      }

      # chimeric MS1 entries
      lens2m <-lengths(df2$ms1_moverz)
      lens2a <- lengths(df2$apex_scan_num)
      bads <- lens2m & (lens2m != lens2a)
      dfx <- df2[bads, ]
      dfx$apex_scan_num <- mapply(function (x, n) rep_len(x, n), 
                                  dfx$apex_scan_num, lens2m[bads],
                                  USE.NAMES = FALSE, SIMPLIFY = FALSE)
      df2[bads, ] <- dfx
      
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
  # maybe useful later
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
  
  if (maxn_mdda_precurs >= 0L) { # was 1
    bads <- lapply(df$ms1_mass, function (x) length(x) == 1L && is.na(x))
    bads <- unlist(bads)
    df <- df[!bads, ]
    
    oks <- lapply(df$ms1_mass, function (x) 
      .Internal(which(x >= min_mass & x <= max_mass)))
    # oks <- lapply(df$ms1_mass, function (x) .Internal(which(x >= 230 & x <= max_mass)))

    # better check all MS1 columns with "list" properties...
    # also check all equal: lengths(df$ms1_mass), lengths(apex_scan_num)...
    df$ms1_mass <- mapply(function (x, y) x[y], df$ms1_mass, oks, 
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_moverz <- mapply(function (x, y) x[y], df$ms1_moverz, oks, 
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_charge <- mapply(function (x, y) x[y], df$ms1_charge, oks, 
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_int <- mapply(function (x, y) x[y], df$ms1_int, oks, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$apex_scan_num <- mapply(function (x, y) x[y], df$apex_scan_num, oks, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    # ms1_mass may again contain numeric(0) following the above filtration
    df <- df[lengths(df$ms1_mass) > 0L, ]
  }
  # may be deleted; ms1_charge, ms1_mass are list not scalar anymore
  else {
    # debugging; expecting vectors of ms1_charge etc. but getting lists.
    if (is_pasef) {
      qs::qsave(df, file.path("~", paste0(filename, "_.rds")), preset = "fast")
    }
    
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


#' Obtains the MS1-X space
#'
#' Followed by de-isotoping and temporarily puts MS1-XYZ of MS2 scans to the
#' preceding MS1 scans.
#'
#' Note that the Y-values were averaged from franking MS1 scans, not the real
#' precursor intensities for for ascribing X-values.
#'
#' @param df A data frame.
get_ms1xs_space <- function (df)
{
  pos_levs <- getMSrowIndexes(df$ms_level)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  
  ms1_moverzs <- df[["ms1_moverz"]]
  ms1_ints <- df[["ms1_int"]]
  ms1_charges <- df[["ms1_charge"]]
  ms1_masses <- df[["ms1_mass"]]
  
  outx <- outy <- outz <- outm <- vector("list", length(ms1_stas))
  
  # i <- which(ms1_stas == 13187) # 2146
  for (i in seq_along(ms1_stas)) {
    rng1 <- ms1_stas[[i]]
    rng2 <- ms2_stas[[i]]:ms2_ends[[i]]
    
    # isolation windows can have overlaps -> 
    #   the same precursor at multiple windows -> duplicated MS1 entries
    xs <- .Internal(unlist(ms1_moverzs[rng2], recursive = FALSE, use.names = FALSE))
    ys <- .Internal(unlist(ms1_ints[rng2], recursive = FALSE, use.names = FALSE))
    zs <- .Internal(unlist(ms1_charges[rng2], recursive = FALSE, use.names = FALSE))
    ms <- .Internal(unlist(ms1_masses[rng2], recursive = FALSE, use.names = FALSE))
    nx <- length(xs)
    
    if (!nx) {
      next
    }

    if (nx > 1L) {
      ord <- order(xs)
      xs <- xs[ord]
      ys <- ys[ord]
      zs <- zs[ord]
      ms <- ms[ord]
      
      oks <- !duplicated(xs)
      xs <- xs[oks]
      ys <- ys[oks]
      zs <- zs[oks]
      ms <- ms[oks]
    }
    
    outx[[i]] <- xs
    outy[[i]] <- ys
    outz[[i]] <- zs
    outm[[i]] <- ms
  }
  
  df[["ms1_moverz"]][ms1_stas] <- outx
  df[["ms1_int"]][ms1_stas] <- outy
  df[["ms1_charge"]][ms1_stas] <- outz
  df[["ms1_mass"]][ms1_stas] <- outm
  
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
  mid_len <- max(as.integer(len/2), 1L)
  
  max_rt <- rts[len]
  mid_rt <- rts[mid_len]
  
  grs <- which(rts > mid_rt + d)
  upr <- if (length(grs)) grs[[1]] else len
  
  max(upr - mid_len, 1L)
}


#' Helper of \link{deisoDDA}.
#' 
#' @param debug Logical; in a debug mode or not.
#' @inheritParams deisoDDA
predeisoDDA <- function (filename = NULL, temp_dir = NULL, 
                         mzml_type = "raw", ppm_ms1 = 10L, ppm_ms2 = 10L, 
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
  # 0 if undertermined by Bruker metadata
  ms1_moverzs <- ans$ms1_moverzs
  ms1_ints <- ans$ms1_ints
  ms1_charges <- ans$ms1_charges
  
  scan_title <- ans$scan_title
  ms_level <- ans$ms_level
  ret_time <- ans$ret_time
  scan_num <- ans$scan_num
  orig_scan <- ans$orig_scan
  iso_ctr <- ans$iso_ctr
  iso_lwr <- ans$iso_lwr
  iso_upr <- ans$iso_upr
  
  raw_file <- ans$raw_file # scalar
  
  ## NULL with Thermo's data
  mobility <- ans$mobility
  slice_start = ans$slice_start
  slice_end = ans$slice_end

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
  
  msx_moverzs <- restmt[["xvals"]]
  msx_ints <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  rm(list = "restmt")
  
  if (deisotope_ms2) {
    message("Deisotoping MS2.")
    oks2 <- ms_level == 2L
    ans2 <- deisoDDAMS2(msx_moverzs[oks2], msx_ints[oks2], 
                        topn_ms2ions = topn_ms2ions, 
                        ppm_ms2_deisotope = ppm_ms2_deisotope, 
                        max_ms2_charge = max_ms2_charge, 
                        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
                        quant = quant, tmt_reporter_lower = tmt_reporter_lower, 
                        tmt_reporter_upper = tmt_reporter_upper, 
                        exclude_reporter_region = exclude_reporter_region, 
                        n_fwd = 10L, n_bwd = 10L, offset_upr = 30L, 
                        offset_lwr = 30L, n_para = n_para)
    msx_moverzs[oks2] <- ans2[["msx_moverzs"]]
    msx_ints[oks2] <- ans2[["msx_ints"]]
    msx_charges[oks2] <- ans2[["msx_charges"]]
    rm(list = c("ans2", "oks2"))
    
    message("Completed MS2 deisotoping.")
  }
  
  # MS1 X, Y and Z are NA from Mzion::readRAW and require deisotoping
  if (maxn_mdda_precurs) {
    message("Deisotoping MS1.")

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
      
      # `ms1_ints` are the sum over flankings, 
      #    not the corresponding MS1-Y of ms1_moverzs
      ms1_moverzs <- ans$ms1_moverzs
      ms1_masses <- ans$ms1_masses
      ms1_ints <- ans$ms1_ints
      ms1_charges <- ans$ms1_charges
      rm(list = "ans")
      message("Completed MS1 deisotoping at: ", Sys.time())
    }

    # look up MS2 for undetermined precursor charge states
    if (deisotope_ms2) {
      rows1 <- lapply(ms1_moverzs, is.null) # all MS1 and some MS2
      rows1 <- .Internal(unlist(rows1, recursive = FALSE, use.names = FALSE))
      rows2 <- ms_level == 2L
      rows <- .Internal(which(rows2 & rows1)) # MS2 without precursor info
      nrows <- length(rows)
      
      if (nrows) {
        # msx_moverzs[rows] and msx_ints[rows] not NULL; msx_charges[rows]: NULL
        # msx_charges[rows] <- NA_real_
        ans2 <- mapply(find_ms1byms2, 
                       moverzs = msx_moverzs[rows], msxints = msx_ints[rows], 
                       charges = msx_charges[rows], center = iso_ctr[rows], 
                       iso_lwr = iso_lwr[rows], iso_upr = iso_upr[rows], 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
        
        if (length(ans2) != length(rows)) {
          stop("Develeoper: checks for row drops.")
        }

        ans2 <- dplyr::bind_rows(ans2)
        ans2$ms1_moverzs <- as.list(ans2$ms1_moverzs)
        ans2$ms1_masses <- as.list(ans2$ms1_masses)
        ans2$ms1_ints <- as.list(ans2$ms1_ints)
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
        ms1_masses[rows] <- ans2$ms1_masses
        ms1_ints[rows] <- ans2$ms1_ints
        ms1_charges[rows] <- ans2$ms1_charges
        rm(list = c("ans2", "nas"))
      }
      
      message("Completed auxillary MS1 deisotoping.")
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
      split(ms2_moverzs, grps), split(ms2_ints, grps), 
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
    
    out <- list(
      msx_moverzs = outx,
      msx_ints = outy,
      msx_charges = outz)
    
    # msx_moverzs[rows2] <- outx
    # msx_ints[rows2] <- outy
    # msx_charges[rows2] <- outz
    # rm(list = c("out", "outx", "outy", "outz", "n_chunks", "cl", "grps", 
    #             "oks2", "ms2_moverzs", "ms2_ints"))
  }
  else {
    out <- getMS2xyz(
      ms2_moverzs, ms2_ints, 
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
    
    # no need to check rptr_moverzs, rptr_ints alignment since equal length
    # msx_moverzs[rows2] <- out[[1]]
    # msx_ints[rows2] <- out[[2]]
    # msx_charges[rows2] <- out[[3]]
    
    # rm(list = c("out", "oks2"))
  }
  
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


#' De-isotoping PASEF DDA-MS1
#' 
#' Depreciated; for PASEF with MS1 combining of all slices.
#' 
#' @param msx_moverzs Full-spectrum MS1 or MS2 m-over-z values.
#' @param msx_ints Full-spectrum MS1 or MS2 intensities.
#' @param ms1_fr MS1 frame numbers.
#' @param ms_level Vectors of MS levels.
#' @param orig_scan Original scan numbers.
#' @param iso_ctr A vector of isolation centers.
#' @param iso_lwr A vector of isolation lowers.
#' @param iso_upr A vector of isolation uppers.
#' @param slice_start A vector of indexes of slice starts.
#' @param slice_end A vector of indexes of slice ends.
#' @param mobility A vector of mobilities.
#' @param ms1_moverzs MS1 moverz values.
#' @param ms1_ints MS1 intensity values.
#' @param ms1_charges MS1 charge states.
#' @param filename A filename for logging.
#' @inheritParams matchMS
pasefMS1xyz0 <- function (msx_moverzs, msx_ints, ms1_fr, ms_level, orig_scan, 
                         iso_ctr = NULL, iso_lwr = NULL, iso_upr = NULL, 
                         slice_start = NULL, slice_end = NULL, mobility = NULL, 
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
      ms1_fr = ms1_fr, ms_level = ms_level, orig_scan = orig_scan,
      iso_ctr = iso_ctr, iso_lwr = iso_lwr, iso_upr = iso_upr, 
    )
  }

  if (!use_defpeaks) {
    ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", length(msx_moverzs))
  }
  
  # separate into MS1 and MS2 data
  oks1 <- ms_level == 1L
  oks2 <- .Internal(which(!oks1))
  oks1 <- .Internal(which(oks1))
  
  xs1 <- msx_moverzs[oks1]
  ys1 <- msx_ints[oks1]
  scans1 <- as.integer(ms1_fr[oks1])
  len1 <- length(ys1)
  
  if (debug) {
    df0 <- df0[oks2, ]
    df0 <- split(df0, df0$ms1_fr)
  }

  frs  <- ms1_fr[oks2]
  ctrs <- iso_ctr[oks2]
  lwrs <- iso_lwr[oks2]
  uprs <- iso_upr[oks2]
  
  ctrs <- split(ctrs, frs)
  lwrs <- split(lwrs, frs)
  uprs <- split(uprs, frs)
  rngs <- split(oks2, frs)
  frs  <- split(frs,  frs)
  
  # find the corresponding MS1 scan number /frame for an MS2 frame
  idxes <- match(as.integer(names(frs)), as.integer(scans1))

  for (i in seq_along(idxes)) {
    # i <- which(names(frs) == "9915"); i <- which(names(frs) == "18203")
    # can it be all NA since some MS1 scans were removed?
    if (is.na(row <- idxes[[i]])) {
      next
    }

    rng1 <- max(1L, row - n_mdda_flanks):min(len1, row + n_mdda_flanks)
    ctri <- ctrs[[i]]
    lwri <- lwrs[[i]]
    upri <- uprs[[i]]
    rng2 <- rngs[[i]]
    
    ###
    # may set a minimum Y of 500 if found multiple precursors and some >= 500
    ###
    
    ans <- find_mdda_mms1s(
      msx_moverzs = xs1[rng1], 
      msx_ints = ys1[rng1], 
      iso_ctr = ctri, iso_lwr = lwri, iso_upr = upri, 
      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
      max_ms1_charge = max_ms1_charge, n_fwd = 15L, n_bwd = 20L, 
      offset_upr = 20L, offset_lwr = 30L, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks, min_mass = min_mass, margin = 0, 
      width_left = 0.26, width_right = 0.25, 
      is_pasef = TRUE, min_y = 10, ms1_min_ratio = 0.05)
    
    # Precursor x, y and z values for each MS2
    xs <- ans[["x"]]
    ys <- ans[["y"]]
    zs <- ans[["z"]]
    
    if (FALSE) {
      for (j in seq_along(xs)) {
        xsj <- xs[[j]]
        p <- which.min(abs(xsj - ctri[[j]]))
        xs[[j]] <- xsj[[p]]
        ys[[j]] <- ys[[p]]
        zs[[j]] <- zs[[p]]
      }
    }

    if (length(xs) != length(rng2)) { stop("Check for entry-dropping.") }
    
    if (debug) {
      dfx <- dplyr::bind_cols(tibble::tibble(x = xs, y = ys, z = zs), df0[[i]]) |>
        dplyr::select(-c("msx_moverzs", "msx_ints", "ms1_fr", "ms_level"))
    }
    
    # empties <- lengths(xs) == 0L
    # xs[empties] <- NA_real_
    # ys[empties] <- NA_integer_
    # zs[empties] <- NA_integer_
    
    # updates corresponding MS1 X, Y and Z for each MS2
    oks <- .Internal(which(lengths(xs) > 0L))
    ms1_moverzs[rng2][oks] <- xs[oks]
    ms1_ints[rng2][oks] <- ys[oks]
    ms1_charges[rng2][oks] <- zs[oks]
  }
  
  ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
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
  rm(list = "pos_levs")
  
  # go from z = min_ms1_charge:max_ms1_charge first,  
  # then if (max_ms1_charge < 6) max_ms1_charge:6
  
  len1  <- length(ms1_stas)

  # get MS2 precursor XYZ values from multiple adjacent MS1 scans
  for (i in 1:len1) {
    rng1 <- ms1_stas[max(1L, i - n_mdda_flanks):min(len1, i + n_mdda_flanks)]
    rng2 <- ms2_stas[i]:ms2_ends[i]

    ans <- find_mdda_mms1s(
      msx_moverzs = msx_moverzs[rng1], msx_ints = msx_ints[rng1], 
      iso_ctr = iso_ctr[rng2], iso_lwr = iso_lwr[rng2], iso_upr = iso_upr[rng2],
      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
      max_ms1_charge = max_ms1_charge, n_fwd = 15L, n_bwd = 20L, 
      offset_upr = 20L, offset_lwr = 30L, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks, min_mass = min_mass)
    
    # Precursor X, Y and Z values for each MS2
    # ys: summed over flanking MS1 scans, not the true precursor intensity of xs
    xs <- ans[["x"]]
    ys <- ans[["y"]]
    zs <- ans[["z"]]

    if (length(xs) != length(rng2)) {
      stop("Check for entries drop in ", fun, " for ", filename) 
    }
    
    # updates corresponding MS1 x, y and z for each MS2
    if (length(oks <- .Internal(which(lengths(xs) > 0L)))) {
      ms1_moverzs[rng2][oks] <- xs[oks]
      ms1_charges[rng2][oks] <- zs[oks]
      ms1_ints[rng2][oks] <- ys[oks]
    }
  }

  ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  if (length(ms1_moverzs) != len0) {
    stop("Developer: check for entries dropping.")
  }

  message("Completed ", fun, " for ", filename)
  
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


#' Spreads MS data across adjacent scans.
#' 
#' Not used.
#' 
#' @param matx The data matrix of moverzs.
#' @param maty The data matrix of intensities.
#' @param matz The data matrix of charge states.
#' @param ns The number of observing spectra that have contributed to an MS
#'   feature.
#' @param ps The row positions (along LC scans) of features.
spreadMS_v1 <- function (matx, maty, matz, ns, ps)
{
  n_scans <- nrow(matx)
  n_masses <- ncol(matx)
  
  if (n_scans <= 2L)
    return(list(x = matx, y = maty, z = matz))
  
  gap1 <- 2L
  gap2 <- 2L
  gap0 <- 1L
  
  for (i in 1:n_masses) {
    ni <- ns[[i]]
    pi <- ps[[i]]
    
    rows0 <- pi[ni >= 3L]
    rows1 <- pi[ni == 1L]
    rows2 <- pi[ni == 2L]
    
    for (j in seq_along(rows1)) {
      rj1 <- rows1[[j]]
      rng1 <- max((rj1 - gap1), 1L):min((rj1 + gap1), n_scans)
      
      # what if !is.na(matx[rng1, i])??? do not overwrite...

      matx[rng1, i] <- matx[rj1, i]
      maty[rng1, i] <- maty[rj1, i]
      matz[rng1, i] <- matz[rj1, i]
    }
    
    for (j in seq_along(rows2)) {
      rj2 <- rows2[[j]]
      rng2 <- max((rj2 - gap2), 1L):min((rj2 + gap2), n_scans)
      matx[rng2, i] <- matx[rj2, i]
      maty[rng2, i] <- maty[rj2, i]
      matz[rng2, i] <- matz[rj2, i]
    }
    
    for (j in seq_along(rows0)) {
      rj0 <- rows0[[j]]
      rng0 <- max((rj0 - gap0), 1L):min((rj0 + gap0), n_scans)
      matx[rng0, i] <- matx[rj0, i]
      maty[rng0, i] <- maty[rj0, i]
      matz[rng0, i] <- matz[rj0, i]
    }
  }
  
  list(x = matx, y = maty, z = matz)
}


#' Combines adjacent MS2 traces.
#' 
#' Not used.
#' 
#' @param xs Vectors of MS2 moverzs.
#' @param ys Vectors of MS2 intensities.
#' @param zs Vectors of MS2 charge states.
#' @param ws Weights.
#' @inheritParams matchMS
#' @importFrom fastmatch %fin% 
comb_mstraces <- function (xs, ys, zs, ws, n_dia_ms2bins = 1L)
{
  if (!n_dia_ms2bins)
    return(list(x = xs, y = ys, z = zs))
  
  ansz <- ansy <- ansx <- vector("list", len <- length(xs))
  
  for (i in 1:n_dia_ms2bins) {
    ansx[[i]] <- xs[[i]]
    ansy[[i]] <- ys[[i]]
    ansz[[i]] <- zs[[i]]
  }
  
  sta <- n_dia_ms2bins + 1L
  end <- len - n_dia_ms2bins
  
  for (i in sta:end) {
    bf <- i - n_dia_ms2bins
    af <- i + n_dia_ms2bins
    rng <- bf:af
    # xbf <- .Internal(unlist(xs[bf:(i-1)], recursive = FALSE, use.names = FALSE))
    # xaf <- .Internal(unlist(xs[(i+1):af], recursive = FALSE, use.names = FALSE))
    # xcr <- xs[[i]]
    # ixbf <- index_mz(xbf, from, step)
    # ixaf <- index_mz(xaf, from, step)
    # ixcr <- index_mz(xcr, from, step)
    # oks_bf <- ixbf %fin% ixcr | (ixbf - 1L) %fin% ixcr | (ixbf + 1L) %fin% ixcr

    ansx[[i]] <- .Internal(unlist(xs[rng], recursive = FALSE, use.names = FALSE))
    ansy[[i]] <- .Internal(unlist(mapply(function (y, w) y * w, ys[rng], ws), 
                                  recursive = FALSE, use.names = FALSE))
    ansz[[i]] <- .Internal(unlist(zs[rng], recursive = FALSE, use.names = FALSE))
  }
  
  for (i in (end+1):len) {
    ansx[[i]] <- xs[[i]]
    ansy[[i]] <- ys[[i]]
    ansz[[i]] <- zs[[i]]
  }
  
  list(x = ansx, y = ansy, z = ansz)
}


#' Finds the positions of logical gates.
#' 
#' @param vals A logical vector.
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

  ups <- edges[["ups"]]
  dns <- edges[["dns"]]
  
  ps   <- mapply(function (x, y) x:y, ups, dns, SIMPLIFY = FALSE, 
                 USE.NAMES = FALSE)
  lens <- lengths(ps)
  oks  <- lens > 2L
  
  if (any(oks)) {
    oks <- .Internal(which(oks))
    ps0 <- ps[-oks]
    ps1 <- ps[oks]
    lens1 <- lens[oks]
    
    for (i in seq_along(ps1)) {
      psi <- ps1[[i]]
      
      if ((li <- lens1[[i]]) %% 2L) {
        psi <- c(ps1[[i]], 0L)
        li <- li + 1L
      }
      
      fcts <- rep(1:(li/2L), each = 2L)
      ps1[[i]] <- split(psi, fcts)
    }
    
    ps1 <- .Internal(unlist(ps1, recursive = FALSE, use.names = FALSE))
    ps <- c(ps0, ps1)
    
    ord <- lapply(ps, `[[`, 1)
    ord <- .Internal(unlist(ord, recursive = FALSE, use.names = FALSE))
    ord <- order(ord)
    ps <- ps[ord]
  }
  
  ps
}


#' Helper to find the positions of logical gates.
#' 
#' @param vals A logical vector.
find_gate_edges <- function (vals)
{
  if (length(vals) <= 1L) {
    return(NULL)
  }

  vec <- diff(vals, 1L) == 1L
  lenv <- length(vec)
  
  if (vec[lenv]) {
    vec <- c(vec, FALSE)
    lenv <- lenv + 1L
  }
  
  # all discrete
  if (!any(vec)) {
    return(NULL)
  }

  ds <- diff(vec)
  ups <- which(ds == 1L) + 1L
  dns <- which(ds == -1L) + 1L
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
  ps <- find_gates(unv)
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
#' @param y_rpl A replacement value of intensities.
#' @param max_perc The maximum percentage of baseline levels.
#' @examples
#' mzion:::find_lc_gates(xs = rep_len(100, 15), ys = c(10,0,0,0,11,15,15,0,0,12,0,10,0,0,10), ts = seq_len(15), n_dia_scans = 2)
#' mzion:::find_lc_gates(xs = rep_len(100, 18), ys = c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 1), 20, 50), ts = seq_len(18))
#' mzion:::find_lc_gates(xs = rep_len(100, 21), ys = c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 4), 20, 50), ts = seq_len(21))
#'
#' # all discrete
#' mzion:::find_lc_gates(xs = rep_len(100, 18), ys = c(rep(0, 5), 100, rep(0, 6), 200, rep(0, 4), 50), ts = seq_len(18))
#' 
#' # two passes
#' ys <- c(rep_len(0, 10), rep_len(100, 50), rep_len(0, 10), rep_len(100, 50), rep_len(3, 5), rep_len(50, 80), rep_len(0, 10))
#' xs <- rep(5, length(ys))
#' ts <- seq_along(ys)
#' mzion:::find_lc_gates(xs, ys, ts)
#' @importFrom fastmatch %fin%
#' @return Scan indexes of LC peaks.
find_lc_gates <- function (xs = NULL, ys = NULL, ts = NULL, n_dia_scans = 6L, 
                           y_rpl = 2.0, step = 5e-6, max_perc = .05)
{
  # if (n_dia_scans <= 0L) return(.Internal(which(ys > 0))) # should not occur
  ymax <- max(ys, na.rm = TRUE)
  ymin <- find_baseline(ys, ymax)
  ys[ys < ymin] <- NA_real_
  ys <- fill_lc_gaps(ys = ys, n_dia_scans = n_dia_scans, y_rpl = y_rpl)

  ioks <- .Internal(which(ys > 0))
  nx   <- length(ioks)
  i1hs <- ioks # for one-hit-wonders
  
  # one `one-hit wonder` across LC
  if (nx == 1L) {
    return(list(apex = ioks, yints = ys[ioks], ns = 1L, ranges = ioks, 
                xstas = ioks))
  }
  
  edges <- find_gate_edges(ioks)
  
  # all discrete `one-hit wonders`
  if (is.null(edges)) {
    return(list(apex = ioks, yints = ys[ioks], ns = rep_len(1L, nx), 
                ranges = ioks, xstas = ioks))
  }
  
  ups <- edges[["ups"]]
  dns <- edges[["dns"]]
  nps <- length(ups)
  
  stas <- ends <- vector("list", nps)
  
  for (i in seq_len(nps)) {
    # i = 6
    ri   <- ioks[ups[[i]]:dns[[i]]] # indexes in relative to ys
    i1hs <- i1hs[!i1hs %fin% ri]
    lenr <- length(ri)
    
    if (lenr <= 60L) {
      stas[[i]] <- ri[[1]] # the starting index of gate-i in relative to ys
      ends[[i]] <- ri[lenr]
    }
    else {
      # same X but different Ys will be treated as one peak
      if (FALSE) {
        ans <- find_lc_xedges(xs = xs[ri], sta = ri[[1]], step = step)
        stas[[i]] <- ans$stas
        ends[[i]] <- ans$ends
      }
      
      r_first <- ri[[1]]
      perc <- min(ymin/ymax, max_perc)
      ans <- lapply(c(perc, .06 * 1:10 + perc), find_lc_gates2, 
                    ys = ys[ri], sta = r_first, n_dia_scans = n_dia_scans, 
                    y_rpl = y_rpl)
      ans_stas <- lapply(ans, `[[`, "stas")
      n_stas   <- lengths(ans_stas)
      n_max    <- .Internal(which.max(n_stas))
      ans_stas <- ans_stas[[n_max]]
      ans_ends <- ans[[n_max]][["ends"]]
      n_stas   <- n_stas[[n_max]]
      
      if (n_stas <= 1) { # can be zero
        stas[[i]] <- r_first
        ends[[i]] <- ri[lenr]
      }
      else {
        stax <- endx <- vector("list", n_stas)
        stax[[1]] <-  ri[[1]]
        endx[[n_stas]] <- ri[lenr]
        
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
  }
  
  # flattening
  stas <- unlist(stas, recursive = FALSE, use.names = FALSE)
  ends <- unlist(ends, recursive = FALSE, use.names = FALSE)
  nps2 <- length(stas)
  
  yints  <- xapex <- xstas <- ranges <- vector("list", nps2)

  for (i in seq_along(stas)) {
    si <- stas[[i]]
    ranges[[i]] <- ixi <- si:ends[[i]]
    yts <- calcAUC(ys[ixi], ts[ixi])
    mi  <- yts[["idx"]]
    
    yints[[i]] <- yts[["area"]]
    xapex[[i]] <- si + mi - 1L
    xstas[[i]] <- si
  }
  
  ## outputs
  xapex <- .Internal(unlist(xapex, use.names = FALSE, recursive = FALSE))
  yints <- .Internal(unlist(yints, use.names = FALSE, recursive = FALSE))
  xstas <- .Internal(unlist(xstas, use.names = FALSE, recursive = FALSE))
  
  apexs  <- c(i1hs, xapex)
  yout   <- c(ys[i1hs], yints)
  ns     <- c(rep_len(1L, length(i1hs)), lengths(ranges))
  rout   <- c(i1hs, ranges)
  xstas  <- c(i1hs, xstas)

  if (length(apexs) > 1L) {
    ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, 
                               TRUE, apexs))
    apexs <- apexs[ord]
    yout  <- yout[ord]
    ns    <- ns[ord]
    rout  <- rout[ord]
    xstas <- xstas[ord]
  }

  list(apex = apexs, yints = yout, ns = ns, ranges = rout, xstas = xstas)
}


#' Find the baseline values of Y
#'
#' @param vals A vector of Y values.
#' @param vmax The maximum of vals.
#' @param perc The percentage to the base-peak intensity for cut-offs.
#' @param y_rpl A replacement value of intensities.
#' @param max_perc The maximum percentage of baseline levels.
find_baseline <- function (vals, vmax, perc = .02, max_perc = .05, y_rpl = 2.0)
{
  len <- length(vals)
  vco <- vmax * perc

  # exclude plateau and paddings
  vs  <- vals[(!is.na(vals)) & vals <= vmax * .10 & vals > y_rpl]
  
  if (length(vs) < 5L) {
    return(vco)
  }
  
  # some vs remains part of peaks but `diff(vs)` are expected to be small 
  base <- mean(abs(diff(vs))) * 3 + min(vs)
  min(base, vmax * max_perc)
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
                            min_n = 5L)
{
  ys[ys <= max(ys, na.rm = TRUE) * perc] <- NA_real_
  ys <- fill_lc_gaps(ys, n_dia_scans = n_dia_scans, y_rpl = 2.0)
  
  ioks  <- .Internal(which(ys > 0))
  edges <- find_gate_edges(ioks)
  ups <- edges[["ups"]] # in relative to ys
  dns <- edges[["dns"]]
  
  # discard narrow gates
  if (length(ups)) {
    oks <- dns - ups >= min_n
    ups <- ups[oks]
    dns <- dns[oks]
  }

  # sta <- ioks[[upi]] # the starting index of gate-i in relative to ys
  # stas[[i]] <- sta + ioks[ups] - 1L # in relative to ys
  # ends[[i]] <- sta + ioks[dns] - 1L
  
  list(stas = sta + ioks[ups] - 1L, ends = sta + ioks[dns] - 1L)
}


#' Find peak edges by intensity threshold relative to the base peak
#' 
#' Not yet used.
#' 
#' @param ys A vector of intensity values.
#' @param ps A vector of Y positions.
#' @param range A vector Y ranges.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @param y_rpl A replacement value of intensities.
find_lc_edges_bp <- function (ys, ps, range, n_dia_scans = 6L, y_rpl = 2.0)
{
  ymax <- max(ys, na.rm = TRUE)
  offp <- ps[[1]] - 1L
  offr <- range[[1]] - 1L
  
  ys[ys <= ymax * .03] <- NA_real_
  ys <- fill_lc_gaps(ys = ys, n_dia_scans = n_dia_scans, y_rpl = y_rpl)
  ioks <- .Internal(which(ys > 0))
  
  edges <- find_gate_edges(ioks)
}


#' Calculate area under a curve
#' 
#' @param ys Y values around an LC peak. Should be no NA values in the sequence.
#' @param ts Time values corresponding to \code{ys}.
calcAUC <- function (ys, ts)
{
  # find the apex index
  imax <- .Internal(which.max(ys))
  len  <- length(ys)
  csum <- cumsum(ys * c(diff(ts), .5))
  tot  <- csum[[len]]
  mval <- tot / 2
  imed <- .Internal(which(csum >= mval))[[1]]
  idx  <- if (abs(imax - imed) > 3) imed else imax

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
                               add_colnames = FALSE, look_back = FALSE)
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
  # 2. `sum_y` seems have little effect with Thermo's but necessary with PASEF
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
      x <- xs[[i]]
      y <- ys[[i]]
      oks <- .Internal(which(!duplicated(ix)))
      
      if (length(oks)) {
        ixs[[i]] <- ix[oks]
        xs[[i]]  <- x[oks]
        ys[[i]]  <- y[oks]
      }
    }
  }
  # rm(list = c("x", "y", "ix", "oks"))
  
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
  
  # collapses adjacent entries
  ps   <- find_gates(unv)
  lenp <- length(ps)
  
  # all discrete values
  if (is.null(ps)) {
    return(list(x = xmat, y = ymat, u = unv))
  }
  
  # for recording matrix columns to be dropped
  ps1 <- vector("list", lenp)
  
  # collapses matrix columns with +/-1 in bin indexes
  for (i in 1:lenp) {
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
    ymat[rows, c2] <- rowSums(ymat[rows, c1:c2, drop = FALSE], na.rm = TRUE)
    # no need since columns traverse from left to right and c1 will be dropped
    # xmat[rows, c1] <- NA_real_ # to prevent value traversing in fwd and bwd looking
    # ymat[rows, c] <- NA_real_
    
    if (look_back && i < lenp) {
      af <- ps[[i+1]][[1]]
      
      if ((af == c2 + 1L)) {
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
        
        rows2 <- .Internal(which(is.na(xmat[, c2])))
        xmat[rows2, c2] <- xmat[rows2, af]
        ymat[rows2, c2] <- rowSums(ymat[rows2, c2:af, drop = FALSE], na.rm = TRUE)
        # about ok to enable, but may be unnecessary
        # xmat[rows2, af] <- NA_real_
        # ymat[rows2, af] <- NA_real_
      }
    }
    
    ps1[[i]] <- c1
  }
  
  # note `-ps1` ok in that at least one ps1 is not 0
  # identical(xmat[, -c(0, 2:3), drop = FALSE], xmat[, -c(2:3), drop = FALSE])
  # !identical(xmat[, 0, drop = FALSE], xmat)
  
  ps1  <- unlist(ps1, recursive = FALSE, use.names = FALSE)
  xmat <- xmat[, -ps1, drop = FALSE]
  ymat <- ymat[, -ps1, drop = FALSE]
  unv  <- unv[-ps1]
  
  if (FALSE && look_back) {
    oks  <- colSums(!is.na(xmat)) > 0L
    xmat <- xmat[, oks]
    ymat <- ymat[, oks]
    unv  <- unv[oks]
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
    if (is.null(ax)) next
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
  charges[is.na(charges)] <- 2L # arbitrary
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
  masses <- masses[okc]
  
  if (len_ms == 1L) {
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  }

  idx <- .Internal(which.max(msxints))
  moverzs <- moverzs[idx]
  charges <- charges[idx]
  msxints <- msxints[idx]
  masses <- masses[idx]
  
  list(ms1_moverzs = moverzs, ms1_masses = masses, 
       ms1_charges = charges, ms1_ints = msxints)
}


