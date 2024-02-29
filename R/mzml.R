#' Helper in preparing mzML inputs.
#' 
#' @param filelist A list of mzML files.
#' @inheritParams load_mgfs
readmzML <- function (filelist = NULL, mgf_path = NULL, 
                      topn_ms2ions = 150L, topn_dia_ms2ions = 2400L, 
                      delayed_diams2_tracing = FALSE, 
                      maxn_dia_precurs = 1000L, 
                      n_dia_ms2bins = 1L, n_dia_scans = 4L, 
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
                      quant = "none", use_lfq_intensity = TRUE, digits = 4L)
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

  
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))

  if (TRUE) {
    peakfiles <- hloadMZML(filelist, mgf_path, temp_dir)
    is_dia <- attr(peakfiles[[1]], "is_dia", exact = TRUE)
    iso_width <- attr(peakfiles[[1]], "iso_width", exact = TRUE)
    peakfiles <- unlist(peakfiles)
    gc() # free up xml pointers
  }
  else {
    is_dia <- TRUE
    peakfiles <- "23aug2017_hela_serum_timecourse_wide_1a.raw.rds"
    iso_width <- 12.0054016
    
    is_dia <- FALSE
    peakfiles <- "20230108_AST_Neo1_DDA_UHG_HeLa_200ng_2th2p5ms_Cycle_01_20230808143253.raw.rds"
    iso_width <- 2.0
    
    is_dia <- FALSE
    peakfiles <- "01CPTAC3_Benchmarking_W_BI_20170508_BL_f02.raw.rds"
    iso_width <- 0.699999988
    
    # peakfiles <- qs::qread("~/peakfiles_bi_g1.rds")
  }

  lenf <- length(peakfiles)
  rams <- find_free_mem()/1024
  n_pcs <- detect_cores(64L) - 1L
  
  ### 
  # rams <- 24
  # n_pcs <- 15
  ###
  
  n_cores <- max(min(n_pcs, ceiling(rams/5L), lenf), 1L)
  r_cores <- round(n_pcs/n_cores)
  n_para <- max(min(n_pcs, r_cores), 1L)

  if (isTRUE(is_dia)) {
    message("Deisotoping DIA-MS.")

    # `else` is only slightly faster than `if` 
    if (n_cores <= 1L) {
      
      # need to average flanking spectra prior to deisotoping...
      
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
    message("Deisotoping DDA-MS.")

    if (n_cores <= 1L) {
      raws <- mapply(
        hdeisoDDA, 
        peakfiles, seq_along(peakfiles), 
        MoreArgs = list(
          mgf_path = mgf_path, 
          temp_dir = temp_dir, 
          ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
          maxn_mdda_precurs = maxn_mdda_precurs, 
          topn_ms2ions = topn_ms2ions, 
          n_mdda_flanks = n_mdda_flanks, 
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
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          use_defpeaks = use_defpeaks, 
          n_para = n_para
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      raws <- parallel::clusterMap(
        cl, hdeisoDDA, 
        peakfiles, seq_along(peakfiles), 
        MoreArgs = list(
          mgf_path = mgf_path, 
          temp_dir = temp_dir, 
          ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
          maxn_mdda_precurs = maxn_mdda_precurs, 
          topn_ms2ions = topn_ms2ions, 
          n_mdda_flanks = n_mdda_flanks, 
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
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          use_defpeaks = use_defpeaks, 
          n_para = n_para
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      parallel::stopCluster(cl)
    }
  }

  message("Completed deisotoping at: ", Sys.time())
  raws <- unlist(raws, recursive = FALSE, use.names = TRUE)
  qs::qsave(raws, file.path(mgf_path, "raw_indexes.rds"), preset = "fast")
  
  type_acqu <- if (is_dia) "dia" else "dda"
}


#' Helper of \link{loadMZML}.
#'
#' @param temp_dir A temporary file folder.
#' @inheritParams readMGF
#' @inheritParams matchMS
hloadMZML <- function (filelist = NULL, mgf_path = NULL, temp_dir = NULL)
{
  # if (maxn_mdda_precurs > 1L && use_defpeaks) {
  #   warning("Default peaks not used at maxn_mdda_precurs > 1;", 
  #           "\nCoerce to use_defpeaks = FALSE.")
  #   use_defpeaks <- FALSE
  # }
  
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
  idx_pwiz <- which(unlist(lapply(softwaresC, xml2::xml_attr, "id")) == "pwiz")
  pwiz <- softwaresC[[idx_pwiz]]
  pwiz_ver <- xml2::xml_attr(pwiz, "version")
  pwiz_ver <- strsplit(pwiz_ver, ".", fixed = TRUE)[[1]]
  
  if ((len_pwiz <- length(pwiz_ver)) >= 2L)
    pwiz_ver_major <- as.numeric(paste0(pwiz_ver[[1]], ".", pwiz_ver[[2]]))
  else if (len_pwiz == 1L)
    pwiz_ver_major <- as.numeric(pwiz_ver)
  else {
    warning("Unknown MSConvert version.")
    pwiz_ver_major <- 3.0
  }
  
  if (pwiz_ver_major < 3)
    stop("Use MSConvert version >= 3.0.")
  
  raw_file <- local({
    idx_file <- which(xml2::xml_name(mzC) == "fileDescription")
    file_des <- mzC[[idx_file]]
    idx_srcl <- which(xml2::xml_name(xml2::xml_children(file_des)) == "sourceFileList")
    info_raw <- xml2::xml_children(file_des)[[idx_srcl]]
    idx_srcf <- which(xml2::xml_name(xml2::xml_children(info_raw)) == "sourceFile")
    info_fi  <- xml2::xml_children(info_raw)[[idx_srcf]]
    raw_file <- xml2::xml_attr(info_fi, "name")
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
  idx_sc <- which(id_nms == "scan")
  if (!length(idx_osc <- which(id_nms == "originalScan"))) idx_osc <- idx_sc
  
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
  
  if (!len)
    stop("No spectrum data found.")
  
  for (i in rng) {
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (xml2::xml_attr(xc[[idx_mslev]], "value") == "2") {
      idx_scan_lwr_2 <- grep("lowest observed m/z", xc)
      if (!length(idx_scan_lwr_2)) {
        warning("Fields of `lowest observed m/z` not found.")
        idx_scan_lwr_2 <- 8L
      }
      
      idx_scan_upr_2 <- grep("highest observed m/z", xc)
      if (!length(idx_scan_upr_2)) {
        warning("Fields of `highest observed m/z` not found.")
        idx_scan_upr_2 <- 9L
      }
      
      idx_precursor_2 <- grep("precursorList", xc)
      if (!length(idx_precursor_2)) {
        warning("Fields of `precursorList` not found.")
        idx_precursor_2 <- 12L
      }
      
      idx_scanList_2 <- grep("scanList", xc)
      if (!length(idx_scanList_2)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_2 <- 11L
      }
      
      idx_bin_2 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_2)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_2 <- 13L
      }
      
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      
      idx_rt_2 <- which(xml2::xml_name(scanList) == "scan")
      if (!length(idx_rt_2)) {
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
      
      idx_ms2_reso <- which(scanList_ret_attrs == "mass resolving power")
      if (length(idx_ms2_reso)) {
        ms2_reso <- scanList_ret[[idx_ms2_reso]]
        ms2_reso <- as.integer(xml2::xml_attr(ms2_reso, "value"))
      } else {
        warning("Fields of `mass resolving power` not found.")
        ms2_reso <- 120000L
      }
      
      # entire MS2 is empty
      if (length(xc) < idx_precursor_2)
        next
      
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
      
      idx_ms1int <- which(selion_nms == "peak intensity") # DIA: zero intensity
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
      idx_scan_lwr_1 <- grep("lowest observed m/z", xc)
      if (!length(idx_scan_lwr_1)) {
        warning("Fields of `lowest observed m/z` not found.")
        idx_scan_lwr_1 <- 8L
      }
      
      idx_scan_upr_1 <- grep("highest observed m/z", xc)
      if (!length(idx_scan_upr_1)) {
        warning("Fields of `highest observed m/z` not found.")
        idx_scan_upr_1 <- 9L
      }
      
      idx_scanList_1 <- grep("scanList", xc)
      if (!length(idx_scanList_1)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_1 <- 11L
      }
      
      idx_bin_1 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_1)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_1 <- 12L
      }
      
      local({
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
                     idx_scan_lwr_2 = 8L, idx_scan_upr_2 = 9L)
{
  len <- length(spec)
  
  ret_times <- orig_scans <- scan_nums <- scan_titles <- 
    iso_ctr <- iso_lwr <- iso_upr <- character(len)
  
  # msx_: both MS1 and MS2, differentiated by ms_lev
  # ms0_: precursor info by other peak-pickings, e.g., MSConvert
  
  ms0_moverzs <- ms0_ints <- ms0_charges <- 
    msx_moverzs <- msx_ints <- msx_charges <- vector("list", len)
  ms_levs <- msx_ns <- integer(len)
  
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
      if (length(xc) < idx_precursor_2)
        next
      
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
    
    if (length(msData) == 2L) {
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_moverzs[[i]] <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ints[[i]] <- readBin(r2, "double", n = msx_n, size = 8L)
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
  
  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Helper of \link{deisoDDA}.
#' 
#' @param raw_id A raw file id.
#' @param n_para The allowance of parallel processing. 
#' @inheritParams deisoDDA
hdeisoDDA <- function (filename, raw_id = 1L, mgf_path = NULL, temp_dir = NULL, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       maxn_mdda_precurs = 5L, 
                       topn_ms2ions = 150L, n_mdda_flanks = 6L, 
                       min_mass = 200L, max_mass = 4500L,
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       quant = "none", use_lfq_intensity = TRUE, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, use_defpeaks = FALSE, 
                       n_para = 1L)
{
  df <- deisoDDA(
    filename, 
    temp_dir = temp_dir, 
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
    maxn_mdda_precurs = maxn_mdda_precurs, 
    topn_ms2ions = topn_ms2ions, 
    n_mdda_flanks = n_mdda_flanks, 
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
    n_para = n_para)
  
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
#' @param n_para The allowance of parallel processing. 
#' @inheritParams matchMS
deisoDDA <- function (filename = NULL, temp_dir = NULL, 
                      ppm_ms1 = 10L, ppm_ms2 = 10L, 
                      maxn_mdda_precurs = 5L, topn_ms2ions = 150L, 
                      n_mdda_flanks = 6L, 
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
                      n_para = 1L)
{
  ###
  # msx_: full spectra of ms1 and ms2, differentiated by ms_lev
  ###
  
  # reads parsed peak lists
  ans <- qs::qread(file.path(temp_dir, filename))
  msx_moverzs <- ans$msx_moverzs
  msx_ints <- ans$msx_ints
  msx_ns <- ans$msx_ns
  ms1_moverzs <- ans$ms1_moverzs # by MSConvert
  ms1_ints <- ans$ms1_ints # by MSConvert
  ms1_charges <- ans$ms1_charges # by MSConvert
  scan_title <- ans$scan_title
  raw_file <- ans$raw_file # scalar
  ms_level <- ans$ms_level
  ret_time <- ans$ret_time 
  scan_num <- ans$scan_num
  orig_scan <- ans$orig_scan
  iso_ctr <- ans$iso_ctr
  iso_lwr <- ans$iso_lwr
  iso_upr <- ans$iso_upr
  rm(list = "ans")
  
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
  
  msx_charges <- vector("list", length(msx_moverzs))
  
  ##
  # try to subset early by min_ and max_ masses...
  # even at unknown z, not very likely to have z >= 2 at small m/z...
  ##
  
  if (deisotope_ms2) {
    if (n_para > 1L) {
      # may subset data at ms_level == 2L
      n_cores <- n_para
      n_chunks <- n_cores * 4L
      grps <- sep_vec(msx_moverzs, n_chunks)
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      out <- parallel::clusterMap(
        cl, getMS2xyz, 
        split(msx_moverzs, grps), split(msx_ints, grps), split(ms_level, grps), 
        MoreArgs = list(
          topn_ms2ions = topn_ms2ions, 
          max_ms2_charge = max_ms2_charge, 
          ppm_ms2_deisotope = ppm_ms2_deisotope, 
          grad_isotope = grad_isotope, 
          fct_iso2 = fct_iso2, 
          quant = quant, 
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region 
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
      parallel::stopCluster(cl)
      
      msx_moverzs <- unlist(lapply(out, function (x) x[[1]]), 
                            recursive = FALSE, use.names = FALSE)
      msx_ints <- unlist(lapply(out, function (x) x[[2]]), 
                         recursive = FALSE, use.names = FALSE)
      msx_charges <- unlist(lapply(out, function (x) x[[3]]), 
                            recursive = FALSE, use.names = FALSE)
      rm(list = c("out", "n_chunks", "cl", "grps"))
    }
    else {
      # may subset data at ms_level == 2L
      out <- getMS2xyz(
        msx_moverzs, msx_ints, ms_level, 
        topn_ms2ions = topn_ms2ions, 
        max_ms2_charge = max_ms2_charge, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, 
        quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region)
      
      msx_moverzs <- out[[1]]
      msx_ints <- out[[2]]
      msx_charges <- out[[3]]
      rm(list = "out")
    }
  }
  
  if (maxn_mdda_precurs) {
    if (FALSE) { # check RAM issue
      n_cores <- n_para
      n_chunks <- n_cores
      
      grps <- local({
        idxes_ms1 <- which(ms_level == 1L)
        diff_ms1 <- c(0L, diff(idxes_ms1))
        oks <- which(diff_ms1 > 1L) # non-consecutive MS1s
        ms1_stas <- idxes_ms1[oks - 1L]
        ms2_stas <- ms1_stas + 1L
        ms2_ends <- idxes_ms1[oks] - 1L
        
        brs <- floor(length(ms2_ends)/n_chunks) * 1:(n_chunks - 1L)
        grps <- findInterval(seq_along(ms_level), ms2_ends[brs])
      })
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      ans <- parallel::clusterMap(
        cl, getMS1xyz, 
        split(msx_moverzs, grps), split(msx_ints, grps), split(ms_level, grps), 
        split(iso_ctr, grps), split(iso_lwr, grps), split(ms1_moverzs, grps), 
        split(ms1_charges, grps), split(ms1_ints, grps), 
        MoreArgs = list(
          maxn_mdda_precurs = maxn_mdda_precurs, n_mdda_flanks = n_mdda_flanks,  
          topn_ms2ions = topn_ms2ions, quant = quant, 
          tmt_reporter_lower = tmt_reporter_lower, 
          tmt_reporter_upper = tmt_reporter_upper, 
          exclude_reporter_region = exclude_reporter_region, 
          max_ms1_charge = max_ms1_charge, ppm_ms1_deisotope = ppm_ms1_deisotope, 
          grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
          use_defpeaks = use_defpeaks
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
      parallel::stopCluster(cl)
      
      ms1_moverzs <- unlist(lapply(ans, function (x) x[["ms1_moverzs"]]), 
                            recursive = FALSE, use.names = FALSE)
      ms1_masses <- unlist(lapply(ans, function (x) x[["ms1_masses"]]), 
                           recursive = FALSE, use.names = FALSE)
      ms1_charges <- unlist(lapply(ans, function (x) x[["ms1_charges"]]), 
                            recursive = FALSE, use.names = FALSE)
      ms1_ints <- unlist(lapply(ans, function (x) x[["ms1_ints"]]), 
                         recursive = FALSE, use.names = FALSE)
      ms1_stas <- unlist(lapply(ans, function (x) x[["ms1_stas"]]), 
                         recursive = FALSE, use.names = FALSE)
      ms2_stas <- unlist(lapply(ans, function (x) x[["ms2_stas"]]), 
                         recursive = FALSE, use.names = FALSE)
      ms2_ends <- unlist(lapply(ans, function (x) x[["ms2_ends"]]), 
                         recursive = FALSE, use.names = FALSE)
      rm(list = c("ans", "grps", "cl"))
      gc()
    }
    else {
      ans <- getMS1xyz(
        msx_moverzs = msx_moverzs, msx_ints = msx_ints, 
        ms_level = ms_level, iso_ctr = iso_ctr, iso_lwr = iso_lwr, 
        ms1_moverzs = ms1_moverzs, ms1_charges = ms1_charges, ms1_ints = ms1_ints, 
        maxn_mdda_precurs = maxn_mdda_precurs, n_mdda_flanks = n_mdda_flanks,  
        topn_ms2ions = topn_ms2ions, quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region, 
        max_ms1_charge = max_ms1_charge, ppm_ms1_deisotope = ppm_ms1_deisotope, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
        use_defpeaks = use_defpeaks)
      
      # list(NULL) at ms_level == 1L
      ms1_moverzs <- ans$ms1_moverzs
      ms1_masses <- ans$ms1_masses
      ms1_charges <- ans$ms1_charges
      ms1_ints <- ans$ms1_ints
      # for LFQ MS1
      # compiles ms1_ends later: multiple MS1s followed by multipe MS2s... 
      ms1_stas <- ans$ms1_stas
      ms2_stas <- ans$ms2_stas
      ms2_ends <- ans$ms2_ends
      rm(list = "ans")
    }
    
    # look up ms2 for undetermined precursor charge states
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
        
        ans2 <- dplyr::bind_rows(ans2)
        ans2$ms1_moverzs <- as.list(ans2$ms1_moverzs)
        ans2$ms1_masses <- as.list(ans2$ms1_masses)
        ans2$ms1_charges <- as.list(ans2$ms1_charges)
        ans2$ms1_ints <- as.list(ans2$ms1_ints)
        
        # outputs contain NA and revent them back to list(NULL)
        nas <- which(is.na(ans2$ms1_moverzs))
        
        if (length(nas)) {
          ans2$ms1_ints[nas] <- ans2$ms1_charges[nas] <- 
            ans2$ms1_masses[nas] <- ans2$ms1_moverzs[nas] <- list(NULL)
        }
        
        ms1_moverzs[rows] <- ans2$ms1_moverzs
        ms1_masses[rows] <- ans2$ms1_masses
        ms1_charges[rows] <- ans2$ms1_charges
        ms1_ints[rows] <- ans2$ms1_ints
        rm(list = c("ans2", "nas"))
      }
      
      rm(list = c("rows1", "rows2", "rows"))
    }
    
    ###
    # Up to this point, ms1_moverzs are list(NULL) at ms_level == 1L;
    # Not to trace all MS1 features by only the monoisotopic 
    #  that have been assigned to MS2 scans.
    ###
    
    # Pools precursors at ms_level == 2 to ms_level == 1 for LFQ MS1;
    for (i in seq_along(ms1_stas)) {
      rng1 <- ms1_stas[[i]]
      rng2 <- ms2_stas[i]:ms2_ends[i]
      
      # isolation windows can have overlaps -> 
      #   the same precursor at multiple windows -> duplicated MS1 entries
      xs <- .Internal(unlist(ms1_moverzs[rng2], recursive = FALSE, use.names = FALSE))
      ys <- .Internal(unlist(ms1_ints[rng2], recursive = FALSE, use.names = FALSE))
      zs <- .Internal(unlist(ms1_charges[rng2], recursive = FALSE, use.names = FALSE))
      ms <- .Internal(unlist(ms1_masses[rng2], recursive = FALSE, use.names = FALSE))
      
      if (length(xs)) {
        ord <- order(xs)
        ms1_moverzs[[rng1]] <- xs[ord]
        ms1_ints[[rng1]] <- ys[ord]
        ms1_charges[[rng1]] <- zs[ord]
        ms1_masses[[rng1]] <- ms[ord]
        rm(list = "ord")
      }
    }
    rm(list = c("rng1", "rng2", "xs", "ys", "zs", "ms"))
    # rm(list = c("ms1_stas", "ms2_stas", "ms2_ends"))
  }
  else {
    ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                         ms1_moverzs, ms1_charges, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  # msx_moverzs at ms_level == 1L correspond to full-spectrum ms1_moverzs
  # msx_charges at ms_level == 1L are list(NULL)
  df <- tibble::tibble(
    scan_title = scan_title,
    raw_file = raw_file,
    ms_level = ms_level, 
    ms1_moverz = ms1_moverzs, 
    ms1_mass = ms1_masses,
    ms1_charge = ms1_charges, 
    ms1_int = ms1_ints, 
    # mzML: ret_times in minutes; MGF: ret_times in seconds
    ret_time = ret_time, 
    scan_num = scan_num, 
    orig_scan = orig_scan,
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_charges = msx_charges, 
    msx_n = msx_ns, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints)

  ## LFQ: replaces intensities with apex values
  if (use_lfq_intensity) {
    # df[["orig_ms1_ints"]] <- df[["ms1_int"]]
    step = ppm_ms1 * 1e-6
    
    # later subset by +/- 2 mins...
    ans_prep <- prep_traceXY(
      df[, c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", "ms_level", 
             "msx_moverzs", "msx_ints", "msx_charges", "orig_scan")], 
      from = min_mass, step = step, 
      # 128L MS1 scans corresponds to ~ 2 mins in LC; not yet tested with Astral
      n_chunks = ceiling(sum(df$ms_level == 1L)/512L), # 1024L, more RAM, same speed 
      # included n_dia_scans in parameters later...
      gap = 128L, n_dia_scans = 4L)
    
    dfs <- ans_prep$dfs
    df1s <- ans_prep$df1s
    gaps <- ans_prep$gaps
    types <- ans_prep$types
    rm(list = "ans_prep")
    gc()
    
    cols <- c("ms_level", "ms1_moverz", "ms1_int")

    if (TRUE) {
      cl <- parallel::makeCluster(getOption("cl.cores", 2L))
      out <- parallel::clusterMap(
        cl, htraceXY, 
        lapply(df1s, `[[`, "msx_moverzs"), lapply(df1s, `[[`, "msx_ints"), 
        lapply(dfs, `[`, cols), gaps, types, 
        MoreArgs = list(
          n_dia_scans = 4L, from = min_mass, step = step
        ), SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = "dynamic")
      parallel::stopCluster(cl)
    }
    else { # slow
      valxs <- lapply(df1s, `[[`, "msx_moverzs")
      valys <- lapply(df1s, `[[`, "msx_ints")
      valdf <- lapply(dfs, `[`, cols)
      out <- vector("list", lenvs <- length(df1s))
      
      for (i in seq_along(df1s)) {
        out[[i]] <- htraceXY(
          xs = valxs[[i]], ys = valys[[i]], df = valdf[[i]], gap = gaps[[i]], 
          type = types[[i]], n_dia_scans = 4L, from = min_mass, step = step
        )
      }
      
      rm(list = c("valxs", "valys", "valdf"))
    }

    out <- dplyr::bind_rows(out)
    df[, cols] <- out
    rm(list = "out")
    
    ## Obtains ms1_int from the corresponding full MS1 (msx_ints)
    # cols <- c("ms1_moverz", "ms1_int", "msx_moverzs", "msx_ints")
    # df[rows, cols] <- getMS1Int(df[rows, cols], from = min_mass, step = step)
    
    # later use apex values to update retention times...
  }

  ## cleans up
  rows <- df$ms_level == 1L
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
  
  if (maxn_mdda_precurs >= 1L) {
    bads <- lapply(df$ms1_mass, function (x) length(x) == 1L && is.na(x))
    bads <- unlist(bads)
    df <- df[!bads, ]
    
    oks <- lapply(df$ms1_mass, function (x) .Internal(which(x >= min_mass & x <= max_mass)))

    df$ms1_mass <- mapply(function (x, y) x[y], df$ms1_mass, oks, 
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_moverz <- mapply(function (x, y) x[y], df$ms1_moverz, oks, 
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_charge <- mapply(function (x, y) x[y], df$ms1_charge, oks, 
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
    df$ms1_int <- mapply(function (x, y) x[y], df$ms1_int, oks, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
    # ms1_mass may again contain numeric(0) following the above filtration
    df <- df[lengths(df$ms1_mass) > 0L, ]
  }
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



#' Obtains the indexes of MS1 and MS2 starts and ends.
#' 
#' @param ms_level A vector of MS levles.
#' @param pad_nas Logical; if TRUE, adds padding values to keep the same length.
getMSrowIndexes <- function (ms_level, pad_nas = FALSE)
{
  idxes_ms1 <- which(ms_level == 1L)
  diff_ms1 <- c(0L, diff(idxes_ms1))
  oks <- which(diff_ms1 > 1L) # non-consecutive MS1s
  ms1_stas <- idxes_ms1[oks - 1L]
  ms2_stas <- ms1_stas + 1L
  ms2_ends <- idxes_ms1[oks] - 1L
  
  if (pad_nas) {
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


#' Deisotoping DDA-MS1.
#' 
#' @param ms_level Vectors of MS levels
#' @param iso_ctr A vector of isolation centers.
#' @param iso_lwr A vecor of isolation lowers.
#' @inheritParams find_mdda_mms1s
#' @inheritParams matchMS
getMS1xyz <- function (msx_moverzs = NULL, msx_ints = NULL, 
                        ms_level = NULL, iso_ctr = NULL, iso_lwr = NULL, 
                        ms1_moverzs = NULL, ms1_charges = NULL, ms1_ints = NULL, 
                        maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                        topn_ms2ions = 150L, quant = "none", 
                        tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                        exclude_reporter_region = FALSE, 
                        max_ms1_charge = 4L, ppm_ms1_deisotope = 8L, 
                        grad_isotope = 1.6, fct_iso2 = 3.0, 
                        use_defpeaks = FALSE)
{
  ## Low priority: no data filtration by scan_nums; 
  #  plus, as.integer(scan_nums) may be invalid with Bruker's
  #  better filter data by retention times
  
  if (!use_defpeaks) {
    ms1_moverzs <- ms1_charges <- ms1_ints <- vector("list", length(msx_moverzs))
  }
  
  pos_levs <- getMSrowIndexes(ms_level)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")

  # go from z = min_ms1_charge:max_ms1_charge first,  
  # then if (max_ms1_charge < 6) max_ms1_charge:6
  
  len <- length(ms1_stas)
  
  # get MS2 precursor xyz values from multiple adjacent MS1 scans
  for (i in 1:len) {
    rng1 <- ms1_stas[max(1L, i - n_mdda_flanks):min(len, i + n_mdda_flanks)]
    rng2 <- ms2_stas[i]:ms2_ends[i]
    
    ans <- find_mdda_mms1s(
      msx_moverzs = msx_moverzs[rng1], 
      msx_ints = msx_ints[rng1], 
      iso_ctr = iso_ctr[rng2], iso_lwr = iso_lwr[rng2], 
      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
      max_ms1_charge = max_ms1_charge, n_fwd = 20L, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks)

    # Precursor x, y and z values for each MS2
    xs <- ans[["x"]]
    ys <- ans[["y"]]
    zs <- ans[["z"]]
    
    # updates corresponding MS1 x, y and z for each MS2
    oks <- .Internal(which(lengths(xs) > 0L))
    ms1_moverzs[rng2][oks] <- xs[oks]
    ms1_ints[rng2][oks] <- ys[oks]
    ms1_charges[rng2][oks] <- zs[oks]
  }
  rm(list = c("rng1", "rng2", "xs", "ys", "zs", "oks"))
  
  ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  list(ms1_moverzs = ms1_moverzs, ms1_masses = ms1_masses, 
       ms1_charges = ms1_charges, ms1_ints = ms1_ints, 
       ms1_stas = ms1_stas, ms2_stas = ms2_stas, ms2_ends = ms2_ends)
}


#' Deisotoping DDA-MS2.
#' 
#' @param ms_level Vectors of MS levels
#' @inheritParams find_mdda_mms1s
#' @inheritParams matchMS
getMS2xyz <- function (msx_moverzs = NULL, msx_ints = NULL, ms_level = NULL, 
                        topn_ms2ions = 150L, quant = "none", 
                        tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                        exclude_reporter_region = FALSE, 
                        max_ms2_charge = 3L, ppm_ms2_deisotope = 10L, 
                        grad_isotope = 1.6, fct_iso2 = 3.0)
{
  if (!(len <- length(msx_moverzs)))
    return(NULL)
  
  msx_charges <- vector("list", len)
  is_tmt <- if (isTRUE(grepl("^tmt.*\\d+", quant))) TRUE else FALSE
  
  for (i in 1:len) {
    if (ms_level[[i]] == 2L) {
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
        n_fwd = 10L, offset_upr = 30L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
        order_mz = FALSE)
      
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
        max_charge = max_ms2_charge, n_fwd = 10L, 
        offset_upr = 30L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
        order_mz = FALSE)
      
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
        n_fwd = 20L, offset_upr = 30L, offset_lwr = 30L, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, order_mz = FALSE)
      
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
#' @param gap The size of gap.
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
#' @param neigh_start The starting point of neighbors for the data spreading.
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
  if (is.null(edges <- find_gate_edges(vals)))
    return(NULL)
  
  ups <- edges[["ups"]]
  dns <- edges[["dns"]]
  
  ps <- mapply(function (x, y) x:y, ups, dns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  lens <- lengths(ps)
  oks <- lens > 2L
  
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
  if (length(vals) <= 1L)
    return(NULL)
  
  vec <- diff(vals, 1L) == 1L
  lenv <- length(vec)
  
  if (vec[lenv]) {
    vec <- c(vec, FALSE)
    lenv <- lenv + 1L
  }
  
  # all discrete
  if (!any(vec)) 
    return(NULL)
  
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
  else if (lenu < lend)
    ups <- c(1L, ups)
  else if (lenu > lend) # should not occur with ensured trailing FALSE
    dns <- c(dns, ups[lenu])
  
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
      gates <- find_lc_gates(ansy[, i], n_dia_scans = n_dia_scans)
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
      gates <- find_lc_gates(ansy[, i], n_dia_scans = n_dia_scans)
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


#' Collapse MS data.
#'
#' Similar to \link{collapse_mms1ints} with the addition of zs.
#'
#' Values of \code{xs} can be smaller than \code{lwr} since no ms2_mass cut-offs
#' (some z undetermined).
#'
#' @param xs Vectors of moverzs at a given isolation center (\code{icenter}).
#'   Each entry corresponds to a full vector of MS1 or MS2.
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
    if (!any(oksx <- lengths(xs) > 0L))
      return(null_out)
    
    oksx <- .Internal(which(oksx))
    xs <- xs[oksx]
    ys <- ys[oksx]
    zs <- zs[oksx]
    rm(list = "oksx")
    
    # remove zero intensities
    for (i in seq_along(oky <- lapply(ys, `>`, 0))) {
      oki <- .Internal(which(oky[[i]]))
      xs[[i]] <- xs[[i]][oki]
      ys[[i]] <- ys[[i]][oki]
      zs[[i]] <- zs[[i]][oki]
    }
    rm(list = "oky", "oki")
    
    # does this again after ys removals
    if (!any(oks <- lengths(xs) > 0L))
      return(null_out)
    
    oks <- .Internal(which(oks))
    xs <- xs[oks]
    ys <- ys[oks]
    zs <- zs[oks]
    rm(list = "oks")
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
  rm(list = c("oks", "ix"))
  
  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups <- lapply(ixs, function (x) unv %in% x)
  
  xout <- mapcoll_xyz(vals = xs, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "xs", direct_out = direct_out)
  rm(list = "xs")
  gc()
  
  yout <- mapcoll_xyz(vals = ys, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "ys", direct_out = direct_out)
  rm(list = "ys")
  gc()
  
  zout <- mapcoll_xyz(vals = zs, ups = ups, lenx = lenx, lenu = lenu, 
                      temp_dir = temp_dir, icenter = icenter, ms_lev = ms_lev, 
                      type = "zs", direct_out = direct_out)
  rm(list = c("zs", "ups"))
  gc()
  
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
  if (is.null(ps))
    return(list(x = xout, y = yout, z = zout))

  ps2 <- lapply(ps, `[[`, 2)
  ps2 <- .Internal(unlist(ps2, recursive = FALSE, use.names = FALSE))
  
  for (i in 1:lenp) {
    c12 <- ps[[i]]
    c2 <- c12[[2]]
    
    # with values in both columns, simply overwrite: 1 <- 2; 
    # c2 can be 0
    if (!c2)
      next
    
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
#' @param lenx The length of \code{vals}.
#' @param lenu The number of entries in the universe.
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
  if (!length(vals))
    return(NULL)
  
  out <- if (type == "zs")
    rep_len(list(rep_len(NA_integer_, lenu)), lenx)
  else if (type %in% c("xs", "ys"))
    rep_len(list(rep_len(NA_real_, lenu)), lenx)
  else
    rep_len(list(rep_len(NA_real_, lenu)), lenx)
  
  for (i in 1:lenx) {
    cols <- ups[[i]]
    out[[i]][cols] <- vals[[i]]
  }
  
  out <- do.call(rbind, out)
  
  if (direct_out)
    return(out)
  
  out_name <- paste0(type, ms_lev, "_diauniv_", icenter, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  
  invisible(NULL)
}


#' Finds the gates of retention times.
#'
#' @param ys A vector of intensity value for the same (approximate) mass along
#'   LC.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @examples
#' library(mzion)
#'
#' find_lc_gates(c(10,0,0,0,11,15,15,0,0,12,0,10,0,0,10), 2)
#'
#' # find_lc_gates(c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 1), 20, 50))
#' # find_lc_gates(c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 4), 20, 50))
#'
#' # all discrete
#' # find_lc_gates(c(rep(0, 5), 100, rep(0, 6), 200, rep(0, 4), 50))
#' @return Scan indexes of LC peaks.
find_lc_gates <- function (ys, n_dia_scans = 4L)
{
  # should not occur
  # if (n_dia_scans <= 0L) return(.Internal(which(ys > 0)))
  
  ys <- fill_lc_gaps(ys, n_dia_scans)
  xs <- .Internal(which(ys > 0))
  nx <- length(xs)
  
  # a case of one one-hit wonder across LC
  if (nx == 1L)
    return(list(apex = xs, ns = 1L, ranges = xs))
  
  ysub <- ys[xs]
  edges <- find_gate_edges(xs)
  
  # all discrete one-hit wonders
  if (is.null(edges))
    return(list(apex = xs, ns = rep_len(1L, nx), ranges = xs))
  
  ups <- edges[["ups"]]
  dns <- edges[["dns"]]
  rm(list = "edges")
  
  xus <- xs[ups]
  xds <- xs[dns]
  widths <- xds - xus + 1L
  
  len <- length(ups)
  ps <- ranges <- vector("list", len)
  
  for (i in 1:len) {
    ui <- ups[[i]]
    di <- dns[[i]]
    ps[[i]] <- ui:di
    ranges[[i]] <- xs[ui:di]
  }

  # ps <- mapply(function (x, y) x:y, ups, dns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  xs1 <- xs[-.Internal(unlist(ps, recursive = FALSE, use.names = TRUE))]
  
  nps <- length(ps)
  xps <- rep_len(NA_integer_, nps)
  # yps <- rep_len(NA_real_, nps)
  
  for (i in seq_len(nps)) {
    pi <- ps[[i]]
    xi <- xs[pi]
    yi <- ysub[pi]
    mi <- which.max(yi)
    xps[[i]] <- xi[[1]] + mi - 1L
    # yps[[i]] <- yi[[mi]]
  }
  
  if (FALSE) {
    if (nps > 1L) {
      for (i in 2:nps) {
        ix <- i - 1
        
        if ((xus[[i]] - xds[[ix]]) <= n_dia_scans) {
          if (yps[[i]] >= yps[[ix]])
            xus[[ix]] <- NA_integer_
          else
            xus[[i]] <- NA_integer_
        }
      }
    }
    
    oks <- !is.na(xus)
    xps <- xps[oks]
  }
  
  list(apex = c(xs1, xps), ns = c(rep_len(1L, length(xs1)), widths), 
       ranges = c(xs1, ranges))
}


#' Fills LC gaps.
#'
#' @param ys A vector of intensity values.
#' @param n_dia_scans The number of adjacent MS scans for constructing a peak
#'   profile and thus for determining the apex scan number of an moverz value
#'   along LC.
#' @examples
#' ys <- c(10,0,0,0,11,0,15,0,0,12,0,10,0,0,10)
#' fill_lc_gaps(ys, 3)
#' fill_lc_gaps(ys, 2)
#'
#' ys <- c(10,0,0,0,11,15,15,0,0,12,0,10,0,0,10)
#' fill_lc_gaps(ys, 3)
#' fill_lc_gaps(ys, 2)
fill_lc_gaps <- function (ys, n_dia_scans = 4L)
{
  if (n_dia_scans <= 1L)
    return(ys)
  
  xs <- .Internal(which(ys > 0))
  ds <- diff(xs)
  oks <- ds <= n_dia_scans & ds > 1L
  ps <- .Internal(which(oks))
  ng <- length(ps)
  
  if (!ng)
    return(ys)
  
  for (i in seq_along(ps)) {
    p <- ps[i]
    rng <- (xs[p]+1L):(xs[p+1]-1L)
    ys[rng] <- 1.0
  }
  
  ys
}


#' Collapse MS1 intensities.
#'
#' Allow adjacent values in unv and later collapse adjacent columns/values.
#'
#' @param xs Vectors of SORTED m-over-z values.
#' @param ys Vectors of intensity values corresponding to xs.
#' @param lwr The lower mass limit.
#' @param step The bin size in converting numeric m-over-z values to integers.
#' @param reord Logical; re-order data or not.
#' @param coll Logical; to further collapse results or not.
#' @param cleanup Logical; to perform data row clean-ups or not. Rows may drop
#'   at \code{cleanup = TRUE}.
#' @importFrom fastmatch %fin%
#' @examples
#' # Twos adjacent bins of xs: 392.1796, 392.1845
#' # xs longer than the universe
#' xs <- list(c(391.1883,391.2848,391.6627,391.6995,392.1646,392.1796,392.1845,
#'            392.2030,392.2887,392.6641,392.7812,393.0833,393.2975))
#' ys <- list(c(12827.41,337002.19,617819.69,18045.10,205851.53,15194.98,11318.61,
#'              12970.02,118604.48,75726.89,11676.51,23723.18,55749.93))
#' # collapse_mms1ints(xs, ys, lwr = 389.6529)
#'
#' xs <- list(c(400.6596,401.7076,402.1813,402.1944,402.1969,402.2094,402.5438,402.7112,403.1812,404.1777),
#'            c(400.6599,401.7075,402.1954,402.1975,402.7112,403.1822,404.2777))
#' ys <- list(c(24003.98,53431.96,110619.26,10988.55,12291.00,140045.06,67601.16,11413.04,21651.61,16686.06),
#'            c(10000.1,40000.1,20000.1,50000.1,2500.2,5000.1,30000.1))
#' # collapse_mms1ints(xs, ys, lwr = 400.1994)
#'
#' xs <- ys <- vector("list", 13L)
#' xs[[7]] <- 954.607849; xs[[8]] <- 954.630249; xs[[10]] <- 954.622925
#' ys[[7]] <- 15706.2627; ys[[8]] <- 19803.5879; ys[[10]] <- 31178.9648
#' # collapse_mms1ints(xs, ys, lwr = 951.089731)
collapse_mms1ints <- function (xs = NULL, ys = NULL, lwr = 115L, step = 1e-5, 
                               reord = FALSE, coll = TRUE, cleanup = TRUE)
{
  null_out <- if (coll)
    list(x = NULL, y = NULL, n = NULL)
  else
    list(x = NULL, y = NULL)
  
  if (cleanup) {
    # 1. all xs are NULL
    if (!any(oks <- lengths(xs) > 0L))
      return(null_out)
    
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
    if (!any(oks <- lengths(xs) > 0L))
      return(null_out)
    
    oks <- .Internal(which(oks))
    xs <- xs[oks]
    ys <- ys[oks]
  }
  
  if (reord) {
    lens <- lengths(xs)
    ords <- lapply(xs, order)
    
    for (i in seq_along(xs)) {
      if (lens[[i]]) {
        ordi <- ords[[i]]
        xs[[i]] <- xs[[i]][ordi]
        ys[[i]] <- ys[[i]][ordi]
      }
    }
    # rm(list = "ords", "ordi")
  }

  ixs <- lapply(xs, index_mz, lwr, step)
  
  # 3. remove duplicated ixs
  for (i in seq_along(ixs)) {
    ix <- ixs[[i]]
    x <- xs[[i]]
    y <- ys[[i]]
    oks <- .Internal(which(!duplicated(ix)))
    ixs[[i]] <- ix[oks]
    xs[[i]]  <- x[oks]
    ys[[i]]  <- y[oks]
  }
  # rm(list = c("x", "y", "ix", "oks"))
  
  ## maps ixs vectors to unv (presence or absence)
  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups <- lapply(ixs, function (x) unv %in% x)
  
  # note one-to-one correspondence between ixs and xs
  xmat <- mapcoll_xyz(vals = xs, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)
  # rm(list = "xs")
  # gc()
  
  ymat <- mapcoll_xyz(vals = ys, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)
  # rm(list = c("ys", "ups"))
  # gc()
  
  ## collapses adjacent entries
  ps <- find_gates(unv)
  lenp <- length(ps)
  
  # all discrete values
  if (is.null(ps)) {
    if (coll)
      return(calc_ms1xys(xmat, ymat))
    else
      return(list(x = xmat, y = ymat))
  }
  
  ps2 <- vector("integer", lenp) # columns to be removed

  # collapses matrix columns with +/-1 in bin indexes
  for (i in 1:lenp) {
    c12 <- ps[[i]]
    c2 <- c12[[2]]
    
    # with values in both columns, simply overwrite: 1 <- 2; 
    # c2 can be 0
    if (!c2)
      next
    
    # with values in both columns but simply overwrite: 1 <- 2
    ps2[[i]] <- c2
    c1 <- c12[1]
    rows <- .Internal(which(!is.na(xmat[, c2])))
    xmat[rows, c1] <- xmat[rows, c2]
    ymat[rows, c1] <- ymat[rows, c2]
  }
  # gc()
  
  # note that at least one ps2 are not 0
  # identical(xmat[, -c(0, 2:3), drop = FALSE], xmat[, -c(2:3), drop = FALSE])
  # !identical(xmat[, 0, drop = FALSE], xmat)
  xmat <- xmat[, -ps2, drop = FALSE]
  ymat <- ymat[, -ps2, drop = FALSE]
  
  if (coll)
    calc_ms1xys(xmat, ymat)
  else
    list(x = xmat, y = ymat)
}


#' Makes MS1 matrices of moverzs and intensities
#'
#' @param matx A matrix of moverzs onto the data universe.
#' @param maty Vectors of intensities onto the data universe.
#' @return A list of x: the weighted means of moverzs; y: the means of
#'   intensities; n: the numbers of observations.
calc_ms1xys <- function (matx, maty)
{
  ysums <- colSums(maty, na.rm = TRUE)
  xmeans <- colSums(matx * maty, na.rm = TRUE)/ysums
  ns <- colSums(!is.na(maty), na.rm = TRUE)
  # use ysums, not ysums/ns, for deisotoping
  list(x = xmeans, y = ysums, n = ns)
}


#' Finds mDDA precursors.
#'
#' Averages of multiple MS1 scans.
#'
#' @param msx_moverzs Vectors of bracketing (e.g. +/-6 scans) full-spectrum
#'   moverzs.
#' @param msx_ints Vectors of bracketing full-spectrum intensities.
#' @param iso_ctr Vectors of isolation centers (e.g. top-12 scans) for MS2.
#' @param iso_lwr Vectors of isolation lowers for MS2.
#' @param ppm Mass error tolerance.
#' @param step A step size.
#' @param maxn_precurs Maximum number of precursors for consideration.
#' @param max_ms1_charge Maximum charge state of precursors for consideration.
#' @param width The width of an MS1 window. A wide window is used for containing
#'   isotope envelops.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries.
#' @param use_defpeaks Use default peak info or not.
#' @inheritParams matchMS
#' @return A list. x: monoisotopic moverzs (weighted mean statistics); y:
#'   intensities (mean); z: charge states.
find_mdda_mms1s <- function (msx_moverzs = NULL, msx_ints = NULL, 
                             iso_ctr = NULL, iso_lwr = NULL, 
                             ppm = 10L, maxn_precurs = 5L, max_ms1_charge = 4L, 
                             n_fwd = 20L, grad_isotope = 1.6, fct_iso2 = 3.0, 
                             use_defpeaks = FALSE, width = 2.01, step = ppm/1e6)
{
  # for all (6+1+6) MS1 frames subset by one MS2 iso-window
  ansx1 <- ansy1 <- vector("list", len1 <- length(msx_moverzs))
  # for all MS2s from averaged (6+1+6 -> 1) MS1s 
  ansn2 <- ansx2 <- ansy2 <- vector("list", len2 <- length(iso_ctr))
  
  # go through MS2 entries
  for (i in 1:len2) {
    m2  <- iso_ctr[[i]]
    lwr <- m2 - width
    upr <- m2 + width
    
    # gather e.g. +/-6 MS1s
    for (j in 1:len1) {
      x1s <- msx_moverzs[[j]]
      y1s <- msx_ints[[j]]
      oks <- .Internal(which(x1s > lwr & x1s < upr))
      ansx1[[j]] <- x1s[oks]
      ansy1[[j]] <- y1s[oks]
    }
    
    # collapsed MS1s
    # assume xs are already ordered from low to high
    ans <- collapse_mms1ints(xs = ansx1, ys = ansy1, lwr = lwr, step = step, 
                             reord = FALSE, coll = TRUE, cleanup = TRUE)
    ansx <- ans[["x"]] # the weighted-mean of precursor moverzs
    ansy <- ans[["y"]] # the mean of precursor intensities
    ansn <- ans[["n"]] # the numbers of precursor observations

    # assign precursor info to the corresponding MS2 entry
    if (is.null(ansx))
      next
    
    ansx2[[i]] <- ansx
    ansy2[[i]] <- ansy
    ansn2[[i]] <- ansn
  }
  rm(list = c("ansx1", "ansy1", "ansx", "ansy", "ansn", "ans", "len1", 
              "m2", "lwr", "upr", "x1s", "y1s", "oks"))

  ## Deisotopes precursors for each MS2 scan
  # ansx2[[i]] can be NULL (no precursor found in the isolation window)
  mics <- mapply(
    find_ms1stat, 
    moverzs = ansx2, msxints = ansy2, n_ms1s = ansn2, center = iso_ctr, 
    MoreArgs = list(
      exclude_reporter_region = FALSE, 
      ppm = ppm, ms_lev = 1L, maxn_feats = maxn_precurs, 
      max_charge = max_ms1_charge, n_fwd = n_fwd, offset_upr = 30L, 
      offset_lwr = 30L, order_mz = FALSE, grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2, use_defpeaks = use_defpeaks
    ), SIMPLIFY = FALSE, USE.NAMES = FALSE)

  # (2) subset by isolation window
  # ( `width = 2.01` contains isotope envelope and now need subsetting)
  xs <- ys <- zs <- vector("list", len2)
  
  for (i in seq_len(len2)) {
    mic <- mics[[i]]
    masses <- mic[["masses"]]
    intensities <- mic[["intensities"]]
    charges <- mic[["charges"]]
    
    m <- iso_ctr[[i]]
    w <- iso_lwr[[i]]
    oks <- masses > m - w & masses < m + w
    oks <- .Internal(which(oks))
    
    # accepts all
    if (!length(oks))
      oks <- seq_along(masses)
    
    xs[[i]] <- masses[oks]
    ys[[i]] <- intensities[oks]
    zs[[i]] <- charges[oks]
  }
  # impurities <- lapply(ys, function (x) x/sum(x))
  
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
  
  if (!(len_ms <- length(moverzs)))
    return(na_out)
  
  oks <- (moverzs >= center - iso_lwr) & (moverzs <= center + iso_upr)
  oks <- .Internal(which(oks))
  moverzs <- moverzs[oks]
  
  if (!(len_ms <- length(moverzs)))
    return(na_out)
  
  msxints <- msxints[oks]
  charges <- charges[oks]
  
  charges[is.null(charges)] <- NA_integer_
  charges[is.na(charges)] <- 2L # arbitrary
  masses <- (moverzs - 1.00727647) * charges
  
  if (len_ms == 1L)
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  
  okc <- !is.na(charges)
  okc <- .Internal(which(okc))
  charges <- charges[okc]
  
  if (!(len_ms <- length(charges)))
    return(na_out)
  
  moverzs <- moverzs[okc]
  msxints <- msxints[okc]
  masses <- masses[okc]
  
  if (len_ms == 1L)
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  
  idx <- .Internal(which.max(msxints))
  moverzs <- moverzs[idx]
  charges <- charges[idx]
  msxints <- msxints[idx]
  masses <- masses[idx]
  
  list(ms1_moverzs = moverzs, ms1_masses = masses, 
       ms1_charges = charges, ms1_ints = msxints)
}


