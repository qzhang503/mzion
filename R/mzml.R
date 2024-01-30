#' Helper in preparing mzML inputs.
#' 
#' @param filelist A list of mzML files.
#' @inheritParams load_mgfs
readmzML <- function (filelist = NULL, mgf_path = NULL, 
                      topn_ms2ions = 150L, topn_dia_ms2ions = 500L, 
                      maxn_dia_precurs = 1000L, n_dia_ms2s = 0L, 
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
                      ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                      quant = "none", digits = 4L)
{
  # - hloadMZML (helper)
  #   - loadMZML
  #     - extrDIA
  #     - extrDDA
  # 
  # [Y] is_dia?
  # - deisoDIA (helper)
  #   - deconvDIA
  # - subDIA_MS1
  # - htraceDIA (helper)
  #   - traceDIA
  # [Y] is_dda?
  # - hdeisoDDA
  #   - deisoDDA
  #     - deconvDDA1
  #     - deconvDDA2
  
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))
  
  peakfiles <- hloadMZML(
    filelist = filelist, 
    mgf_path = mgf_path, 
    temp_dir = temp_dir, 
    topn_ms2ions = topn_ms2ions, 
    min_mass = min_mass, 
    max_mass = max_mass, 
    min_ms2mass = min_ms2mass, 
    max_ms2mass = max_ms2mass, 
    min_ms1_charge = min_ms1_charge, 
    max_ms1_charge = max_ms1_charge, 
    min_scan_num = min_scan_num, 
    max_scan_num = max_scan_num, 
    min_ret_time = min_ret_time, 
    max_ret_time = max_ret_time, 
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
    tmt_reporter_lower = tmt_reporter_lower, 
    tmt_reporter_upper = tmt_reporter_upper, 
    exclude_reporter_region = exclude_reporter_region, 
    mgf_cutmzs = mgf_cutmzs, 
    mgf_cutpercs = mgf_cutpercs, 
    deisotope_ms2 = deisotope_ms2, 
    max_ms2_charge = max_ms2_charge, 
    use_defpeaks = use_defpeaks, 
    maxn_mdda_precurs = maxn_mdda_precurs, 
    n_mdda_flanks = n_mdda_flanks, 
    ppm_ms1_deisotope = ppm_ms1_deisotope, 
    ppm_ms2_deisotope = ppm_ms2_deisotope, 
    grad_isotope = grad_isotope, 
    fct_iso2 = fct_iso2, 
    quant = quant, 
    digits = digits)
  is_dia <- attr(peakfiles[[1]], "is_dia", exact = TRUE)
  iso_width <- attr(peakfiles[[1]], "iso_width", exact = TRUE)
  peakfiles <- unlist(peakfiles)
  # free up xml pointers
  gc()
  
  if (isTRUE(is_dia)) {
    message("Deisotoping DIA-MS.")

    if (iso_width > 2.5 && topn_ms2ions <= 500L) {
      topn_ms2ions <- as.integer(iso_width * 40)
      warning("Increase the maximum number of MS2 features to ", topn_ms2ions, 
              " for wide-window DIA.")
    }
    
    n_pcs <- detect_cores(64L) - 1L
    n_cores <- min(n_pcs, ceiling(find_free_mem()/1024/5L), length(peakfiles))
    n_para <- min(n_pcs, round(n_pcs/n_cores))
    
    if (n_cores == 1L) {
      dia_files <- lapply(
        peakfiles, deisoDIA,
        temp_dir = temp_dir, 
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
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      dia_files <- parallel::clusterApply(
        cl, peakfiles, deisoDIA,
        temp_dir = temp_dir, 
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

    len_files <- length(dia_files)
    dia_paths <- file.path(temp_dir, dia_files)
    
    if (len_files == 1L) {
      dfs <- lapply(dia_paths, subDIA_MS1)
    }
    else {
      n_cores <- max(min(parallel::detectCores() - 1L, len_files), 1L)
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      dfs <- parallel::clusterApply(cl, dia_paths, subDIA_MS1)
      parallel::stopCluster(cl)
    }
    
    raws <- mapply(
      htraceDIA, dfs, seq_along(dfs), 
      MoreArgs = list(
        mgf_path = mgf_path, 
        min_ret_time = min_ret_time, 
        max_ret_time = max_ret_time, 
        min_mass = min_mass, 
        max_mass = max_mass, 
        min_ms2mass = min_ms2mass, 
        max_ms2mass = max_ms2mass, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        n_dia_ms2s = n_dia_ms2s, 
        topn_ms2ions = topn_ms2ions, 
        mgf_cutmzs = mgf_cutmzs, 
        mgf_cutpercs = mgf_cutpercs, 
        quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region
      ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  else {
    message("Deisotoping DDA-MS.")
    
    n_pcs <- detect_cores(64L) - 1L
    n_cores <- min(n_pcs, ceiling(find_free_mem()/1024/5L), length(peakfiles))
    # n_para <- min(n_pcs, round(n_pcs/n_cores) * 2L)
    n_para <- min(n_pcs, round(n_pcs/n_cores))

    if (n_cores == 1L) {
      raws <- mapply(
        hdeisoDDA, 
        peakfiles, seq_along(peakfiles), 
        MoreArgs = list(
          mgf_path = mgf_path, 
          temp_dir = temp_dir, 
          ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
          maxn_mdda_precurs = maxn_mdda_precurs, 
          topn_ms2ions = topn_ms2ions, n_mdda_flanks = n_mdda_flanks, 
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
        cl, hdeisoDDA, peakfiles, seq_along(peakfiles), 
        MoreArgs = list(
          mgf_path = mgf_path, 
          temp_dir = temp_dir, 
          ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
          maxn_mdda_precurs = maxn_mdda_precurs, 
          topn_ms2ions = topn_ms2ions, n_mdda_flanks = n_mdda_flanks, 
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
  
  invisible(NULL)
}


#' Helper of \link{loadMZML}.
#'
#' @param temp_dir A temporary file folder.
#' @inheritParams readMGF
#' @inheritParams matchMS
hloadMZML <- function (filelist = NULL, mgf_path = NULL, temp_dir = NULL, 
                       topn_ms2ions = 150L, min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, quant = "none", 
                       digits = 4L)
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

  if (n_cores == 1L) {
    out_names <- vector("list", len)
    
    for (i in 1:len) {
      out_names[[i]] <- loadMZML(
        files[[i]], 
        temp_dir = temp_dir, 
        topn_ms2ions = topn_ms2ions, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_mass = min_mass, 
        max_mass = max_mass, 
        min_ms2mass = min_ms2mass, 
        max_ms2mass = max_ms2mass, 
        min_ms1_charge = min_ms1_charge, 
        max_ms1_charge = max_ms1_charge, 
        min_scan_num = min_scan_num, 
        max_scan_num = max_scan_num, 
        min_ret_time = min_ret_time, 
        max_ret_time = max_ret_time, 
        deisotope_ms2 = deisotope_ms2, 
        max_ms2_charge = max_ms2_charge, 
        use_defpeaks = use_defpeaks, 
        maxn_mdda_precurs = maxn_mdda_precurs, 
        n_mdda_flanks = n_mdda_flanks, 
        ppm_ms1_deisotope = ppm_ms1_deisotope, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        grad_isotope = grad_isotope,
        fct_iso2 = fct_iso2,
        quant = quant, digits = digits)
    }
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out_names <- parallel::clusterApply(
      cl, files, loadMZML, 
      temp_dir = temp_dir, 
      topn_ms2ions = topn_ms2ions, 
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region, 
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass, 
      max_ms2mass = max_ms2mass, 
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge, 
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      deisotope_ms2 = deisotope_ms2, 
      max_ms2_charge = max_ms2_charge, 
      use_defpeaks = use_defpeaks, 
      maxn_mdda_precurs = maxn_mdda_precurs, 
      n_mdda_flanks = n_mdda_flanks, 
      ppm_ms1_deisotope = ppm_ms1_deisotope, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      grad_isotope = grad_isotope,
      fct_iso2 = fct_iso2,
      quant = quant, digits = digits)
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
#' @inheritParams matchMS
#' @inheritParams load_mgfs
loadMZML <- function (xml_file = NULL, temp_dir = NULL, topn_ms2ions = 150L, 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, 
                      ppm_ms1 = 10L, ppm_ms2 = 10L, 
                      min_mass = 200L, max_mass = 4500L, 
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      min_ms1_charge = 2L, max_ms1_charge = 4L, 
                      min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                      min_ret_time = 0, max_ret_time = Inf, 
                      deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                      use_defpeaks = FALSE, 
                      maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                      ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                      grad_isotope = 1.6, fct_iso2 = 3.0, quant = "none", 
                      digits = 4L)
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
                       quant = "none", 
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
                      quant = "none", 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, use_defpeaks = FALSE, 
                      n_para = 1L)
{
  # msx_: full spectra of ms1 and ms2, differentiated by ms_lev
  
  ans <- qs::qread(file.path(temp_dir, filename))
  msx_moverzs <- ans$msx_moverzs
  msx_ints <- ans$msx_ints
  msx_ns <- ans$msx_ns
  ms1_moverzs <- ans$ms1_moverzs
  ms1_ints <- ans$ms1_ints
  ms1_charges <- ans$ms1_charges
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
  
  if (deisotope_ms2) {
    if (n_para > 1L) {
      # may subset data at ms_level == 2L
      
      n_cores <- n_para
      n_chunks <- n_cores * 4L
      grps <- sep_vec(msx_moverzs, n_chunks)
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      out <- parallel::clusterMap(
        cl, deconvDDA2, 
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
      rm(list = "out")
    }
    else {
      # may subset data at ms_level == 2L
      
      out <- deconvDDA2(
        msx_moverzs, msx_ints, ms_level, 
        topn_ms2ions = topn_ms2ions, 
        max_ms2_charge = max_ms2_charge, 
        ppm_ms2_deisotope = ppm_ms2_deisotope, 
        grad_isotope = grad_isotope, 
        fct_iso2 = fct_iso2, 
        quant = quant, 
        tmt_reporter_lower = tmt_reporter_lower, 
        tmt_reporter_upper = tmt_reporter_upper, 
        exclude_reporter_region = exclude_reporter_region 
      )
      msx_moverzs <- out[[1]]
      msx_ints <- out[[2]]
      msx_charges <- out[[3]]
      rm(list = "out")
    }
  }

  # not yet parallel...
  
  if (maxn_mdda_precurs) {
    ans <- deconvDDA1(
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
    
    ms1_moverzs <- ans$ms1_moverzs
    ms1_masses <- ans$ms1_masses
    ms1_charges <- ans$ms1_charges
    ms1_ints <- ans$ms1_ints
    
    # look up ms2 for undetermined precursor charge states
    if (deisotope_ms2) {
      rows1 <- lapply(ms1_moverzs, is.null)
      rows1 <- unlist(rows1, recursive = FALSE, use.names = FALSE)
      rows2 <- ms_level == 2L
      rows <- rows1 & rows2
      
      if (any(rows)) {
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
        
        ms1_moverzs[rows] <- ans2$ms1_moverzs
        ms1_masses[rows] <- ans2$ms1_masses
        ms1_charges[rows] <- ans2$ms1_charges
        ms1_ints[rows] <- ans2$ms1_ints
      }
      
      rm(list = c("rows1", "rows2", "rows", "ans2"))
    }
  }
  else {
    ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                         ms1_moverzs, ms1_charges, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }

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
    ms2_moverzs = msx_moverzs, 
    ms2_ints = msx_ints, 
    ms2_charges = msx_charges, 
    ms2_n = msx_ns, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints, 
  )
  
  df <- df[with(df, ms_level != 1L), ]
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
    
    oks <- lapply(df$ms1_mass, function (x) x >= min_mass & x <= max_mass)
    
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
  
  df
}


#' Deisotoping MS1.
#' 
#' @param ms_level Vectors of MS levels
#' @param iso_ctr A vector of isolation centers.
#' @param iso_lwr A vecor of isolation lowers.
#' @inheritParams find_mdda_mms1s
#' @inheritParams matchMS
deconvDDA1 <- function (msx_moverzs = NULL, msx_ints = NULL, 
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
  
  idxes_ms1 <- which(ms_level == 1L)
  diff_ms1 <- c(0L, diff(idxes_ms1))
  oks <- which(diff_ms1 > 1L) # non-consecutive MS1s
  ms1_stas <- idxes_ms1[oks - 1L]
  ms2_stas <- ms1_stas + 1L
  ms2_ends <- idxes_ms1[oks] - 1L
  rm(list = c("oks", "diff_ms1", "idxes_ms1"))

  # go from z = min_ms1_charge:max_ms1_charge first 
  # if (max_ms1_charge < 6) max_ms1_charge:6
  
  len <- length(ms1_stas)
  
  for (i in 1:len) {
    stas1 <- ms1_stas[max(1L, i - n_mdda_flanks):min(len, i + n_mdda_flanks)]
    stas2 <- ms2_stas[i]
    ends2 <- ms2_ends[i]
    rng <- stas2:ends2
    
    ans <- find_mdda_mms1s(
      msx_moverzs = msx_moverzs[stas1], msx_ints = msx_ints[stas1], 
      ms1_moverzs = ms1_moverzs[rng], ms1_charges = ms1_charges[rng], 
      ms1_ints = ms1_ints[rng], 
      iso_ctr = iso_ctr[rng], iso_lwr = iso_lwr[rng], 
      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_mdda_precurs, 
      max_ms1_charge = max_ms1_charge, n_fwd = 20L, 
      grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
      use_defpeaks = use_defpeaks)
    
    ms1_moverzs[rng] <- ans[[1]]
    ms1_charges[rng] <- ans[[2]]
    ms1_ints[rng] <- ans[[3]]
  }
  # rm(list = c("stas1", "stas2", "ends2", "ms1_stas", "ms2_stas", "ms2_ends"))
  
  ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  list(ms1_moverzs = ms1_moverzs, ms1_masses = ms1_masses, 
       ms1_charges = ms1_charges, ms1_ints = ms1_ints)
}


#' Deisotopes DDA-MS2
#' 
#' @param ms_level Vectors of MS levels
#' @inheritParams find_mdda_mms1s
#' @inheritParams matchMS
deconvDDA2 <- function (msx_moverzs = NULL, msx_ints = NULL, ms_level = NULL, 
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
      mic <- find_ms1stat(
        moverzs = msx_moverzs[[i]], msxints = msx_ints[[i]], center = 0, 
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
#' @inheritParams deisoDIA
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
  
  out <- list(
    # either MS1 or MS2
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    
    scan_title = scan_titles,
    raw_file = raw_file, # single
    ms_level = ms_levs, 
    # mzML: ret_times in minutes; MGF: in seconds
    ret_time = as.numeric(ret_times) * 60, 
    scan_num = as.integer(scan_nums), 
    orig_scan = orig_scans,
    iso_ctr = as.numeric(iso_ctr), 
    iso_lwr = as.numeric(iso_lwr), 
    iso_upr = as.numeric(iso_upr), 
    demux = as.integer(demux)
  )
  
  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  invisible(out_name)
}


#' Helper of \link{deconvDIA}.
#' 
#' @param filename A peaklist filename.
#' @param temp_dir A temp_dir to the filename.
#' @param n_para The allowance of parallel processing.
#' @inheritParams extrDDA
#' @inheritParams matchMS
deisoDIA <- function (filename = NULL, temp_dir = NULL, 
                      ppm_ms1 = 10L, ppm_ms2 = 10L, 
                      min_mass = 200L, # max_mass = 4500L,
                      min_ms2mass = 115L, # max_ms2mass = 4500L, 
                      maxn_dia_precurs = 500L, topn_dia_ms2ions = 500L, 
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
  
  if (n_para > 1L) {
    # n_cores <- max(min(parallel::detectCores() - 1L, 128L), 1L)
    n_cores <- n_para
    n_chunks <- n_cores * 4L
    grps <- sep_vec(msx_moverzs, n_chunks)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out <- parallel::clusterMap(
      cl, deconvDIA, 
      split(msx_moverzs, grps), split(msx_ints, grps), split(ms_level, grps), 
      MoreArgs = list(
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
    rm(list = "out")
    
    msx_moverzs <- unlist(msx_moverzs, recursive = FALSE, use.names = FALSE)
    msx_ints <- unlist(msx_ints, recursive = FALSE, use.names = FALSE)
    msx_charges <- unlist(msx_charges, recursive = FALSE, use.names = FALSE)
    
    ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", length(msx_moverzs))
  }
  else {
    out <- deconvDIA(
      msx_moverzs, msx_ints, ms_level, 
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
    )
    
    msx_moverzs <- out[[1]]
    msx_ints <- out[[2]]
    msx_charges <- out[[3]]
    rm(list = "out")
    
    ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", length(msx_moverzs))
  }

  df <- tibble::tibble(
    scan_title = scan_title,
    raw_file = raw_file,
    ms_level = ms_level, 
    ret_time = ret_time, 
    scan_num = scan_num, 
    orig_scan = orig_scan,
    # either MS1 or MS2
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_charges = msx_charges, 
    msx_ns = msx_ns, 
    
    # MS1 placeholder
    ms1_moverzs = ms1_moverzs, 
    ms1_ints = ms1_ints, 
    # MS1: vector or empty vector for MS1; MS2: NULL
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
#' @param topn_dia_ms2ions Top-N DIA-MS2 features for deisotoping.
#' @inheritParams deisoDIA
#' @inheritParams matchMS
deconvDIA <- function (msx_moverzs = NULL, msx_ints = NULL, ms_level = NULL, 
                       maxn_dia_precurs = 500L, topn_dia_ms2ions = 500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       grad_isotope = 1.6, fct_iso2 = 3.0, quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE)
{
  if (!(len <- length(msx_moverzs)))
    return(NULL)
  
  is_tmt <- if (isTRUE(grepl("^tmt.*\\d+", quant))) TRUE else FALSE
  # if (is_tmt) stop("TMT not yet supported with DIA workflows.")
  
  msx_charges <- vector("list", len)
  
  for (i in 1:len) {
    ms_lev <- ms_level[[i]]
    
    if (ms_lev == 2L && deisotope_ms2) {
      mic <- find_ms1stat(moverzs = msx_moverzs[[i]], msxints = msx_ints[[i]], 
                          center = 0, exclude_reporter_region = is_tmt, 
                          tmt_reporter_lower = tmt_reporter_lower, 
                          tmt_reporter_upper = tmt_reporter_upper, 
                          is_dda = FALSE, ppm = ppm_ms2_deisotope, 
                          ms_lev = 2L, maxn_feats = topn_dia_ms2ions, 
                          max_charge = max_ms2_charge, n_fwd = 10L, 
                          offset_upr = 30L, offset_lwr = 30L, 
                          grad_isotope = grad_isotope, fct_iso2 = fct_iso2, 
                          order_mz = FALSE)
    }
    else if (ms_lev == 1L) {
      mic <- find_ms1stat(moverzs = msx_moverzs[[i]], msxints = msx_ints[[i]], 
                          center = 0, exclude_reporter_region = FALSE, 
                          is_dda = FALSE, ppm = ppm_ms1_deisotope, 
                          ms_lev = 1L, 
                          # by the width of MS1: 395 - 1005???
                          maxn_feats = maxn_dia_precurs, 
                          max_charge = max_ms1_charge, 
                          n_fwd = 20L, offset_upr = 30L, offset_lwr = 30L, 
                          grad_isotope = grad_isotope, 
                          fct_iso2 = fct_iso2, order_mz = FALSE)
    }
    
    msx_moverzs[[i]] <- mic[["masses"]]
    msx_ints[[i]] <- mic[["intensities"]]
    msx_charges[[i]] <- mic[["charges"]] # MS2: NULL; MS1: integer vectors
  }
  
  list(msx_moverzs = msx_moverzs, msx_ints = msx_ints, msx_charges = msx_charges)
}


#' Get precursor mass candidates for each DIA window.
#' 
#' Subset by MS2 isolation windows.
#' 
#' @param filename A filename with prepending path.
subDIA_MS1 <- function (filename = NULL)
{
  df <- qs::qread(filename)
  
  # correspondence of row indexes between MS1 and MS2 scans
  idxes_ms1 <- which(df$ms_level == 1L)
  
  if (length(idxes_ms1) == 1L) {
    ms1_stas <- idxes_ms1
    ms2_stas <- ms1_stas + 1L
    ms2_ends <- nrow(df)
  }
  else {
    diff_ms1 <- c(0L, diff(idxes_ms1))
    oks <- which(diff_ms1 > 1L)
    ms1_stas <- idxes_ms1[oks - 1L]
    ms2_stas <- ms1_stas + 1L
    ms2_ends <- idxes_ms1[oks] - 1L
    rm(list = c("oks", "diff_ms1", "idxes_ms1"))
  }
  
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
    
    # the first center not the smallest: up then down
    if ((imin <- which.min(cents)) > 1L) {
      edg_b <- c(cents[[1]] - half, cents[1:(imin-1L)] + half)
      cuts_b <- findInterval(xs, edg_b)
      xs_b <- split(xs, cuts_b)
      ys_b <- split(ys, cuts_b)
      css_b <- split(css, cuts_b)
      xs_b <- xs_b[-1]
      ys_b <- ys_b[-1]
      css_b <- css_b[-1]
      rows_b <- as.integer(names(xs_b))
      
      edg_a <- c(cents[[imin]] - half, cents[imin:length(cents)] + half)
      cuts_a <- findInterval(xs, edg_a)
      xs_a <- split(xs, cuts_a)
      ys_a <- split(ys, cuts_a)
      css_a <- split(css, cuts_a)
      xs_a <- xs_a[-(nitv <- length(xs_a))]
      ys_a <- ys_a[-nitv]
      css_a <- css_a[-nitv]
      rows_a <- as.integer(names(xs_a))
      
      rows <- c(as.integer(names(xs_b)), as.integer(names(xs_a)) + imin - 1L)
      df$ms1_moverzs[rows2][rows] <- c(xs_b, xs_a)
      df$ms1_ints[rows2][rows] <- c(ys_b, ys_a)
      df$ms1_charges[rows2][rows] <- c(css_b, css_a)
    }
    else {
      edges <- c(cents[[1]] - half, cents + half)
      cuts <- findInterval(xs, edges)
      xs_cuts <- split(xs, cuts)
      ys_cuts <- split(ys, cuts)
      css_cuts <- split(css, cuts)
      
      # excludes precursors outside the isolation window
      rows <- as.integer(names(xs_cuts))
      oks <- rows & rows <= length(rows2)
      rows <- rows[oks]
      df$ms1_moverzs[rows2][rows] <- xs_cuts[oks]
      df$ms1_ints[rows2][rows] <- ys_cuts[oks]
      df$ms1_charges[rows2][rows] <- css_cuts[oks]
    }
  }
  
  df
}


#' Traces DIA-MS2 against MS1.
#'
#' For a single MS1 range, e.g., 592.5 to 604.5.
#'
#' @param df A data frame.
#' @param step1 A step size for MS1.
#' @param step2 A step size for MS2.
#' @inheritParams matchMS
traceDIA <- function (df, min_mass = 200L, min_ms2mass = 115L, step1 = 1E-5, 
                        step2 = 1E-5, n_dia_ms2s = 0L)
{
  # n_dia_ms2s > 1L not yet optimized, longer sequences ms2s than topn_ms2_ions...
  
  mat1 <- trace_ms(xs = df[["ms1_moverz"]], ys = df[["ms1_int"]], 
                   zs = df[["ms1_charge"]], from = min_mass, step = step1)
  mat1x <- mat1[["x"]]
  mat1y <- mat1[["y"]]
  mat1z <- mat1[["z"]]
  rm(list = "mat1")
  
  if (!n_dia_ms2s) {
    for (j in 1:nrow(df)) {
      m1xj <- mat1x[j, ]
      m1yj <- mat1y[j, ]
      m1zj <- mat1z[j, ]
      oks1 <- !is.na(m1xj)
      m1xj <- m1xj[oks1]
      m1yj <- m1yj[oks1]
      m1zj <- m1zj[oks1]
      
      df$ms1_moverz[[j]] <- m1xj
      df$ms1_int[[j]] <- m1yj
      df$ms1_charge[[j]] <- m1zj
    }
    
    return(df)
  }
  
  mat2 <- trace_ms(xs = df[["ms2_moverzs"]], ys = df[["ms2_ints"]], 
                   zs = df[["ms2_charges"]], from = min_ms2mass, step = step2)
  mat2x <- mat2[["x"]]
  mat2y <- mat2[["y"]]
  mat2z <- mat2[["z"]]
  rm(list = "mat2")
  
  len <- nrow(df)
  len2 <- len - n_dia_ms2s
  
  for (j in 1:len) {
    m1xj <- mat1x[j, ]
    m1yj <- mat1y[j, ]
    m1zj <- mat1z[j, ]
    oks1 <- !is.na(m1xj)
    m1xj <- m1xj[oks1]
    m1yj <- m1yj[oks1]
    m1zj <- m1zj[oks1]
    
    sta <- if (j > n_dia_ms2s) j - n_dia_ms2s + 1L else 1L
    end <- if (j < len2) j + n_dia_ms2s - 1L else len
    js <- sta:end
    m2xj <- mat2x[js, ]
    m2yj <- mat2y[js, ]
    m2zj <- mat2z[js, ]
    
    if (n_dia_ms2s == 1L) {
      oks2 <- !is.na(m2xj)
      m2xj <- m2xj[oks2]
      m2yj <- m2yj[oks2]
      m2zj <- m2zj[oks2]
    }
    else {
      ysum <- colSums(m2yj, na.rm = TRUE)
      oks2 <- ysum > 0
      ysum <- ysum[oks2]
      m2xj <- colSums(m2xj[, oks2] * m2yj[, oks2], na.rm = TRUE)/ysum
      m2yj <- ysum
      m2zj <- colMeans(m2zj[, oks2], na.rm = TRUE)
      m2zj <- as.integer(m2zj)
    }
    
    df$ms1_moverz[[j]] <- m1xj
    df$ms1_int[[j]] <- m1yj
    df$ms1_charge[[j]] <- m1zj
    df$ms2_moverzs[[j]] <- m2xj
    df$ms2_ints[[j]] <- m2yj
    df$ms2_charges[[j]] <- m2zj
  }
  
  df
}


#' Helper of tracing DIA.
#' 
#' @param df A data frame.
#' @param raw_id A raw file id. 
#' @inheritParams load_mgfs
htraceDIA <- function (df, raw_id, mgf_path, 
                       min_ret_time = 0L, max_ret_time = Inf, 
                       min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       ppm_ms1 = 20L, ppm_ms2 = 20L, n_dia_ms2s = 0L, 
                       topn_ms2ions = 150L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE)
{
  ### cleanDIAMS
  df <- df[with(df, ms_level != 1L), ]
  # NULL: no MS1 in the bins
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
  
  ## Trace MS1 - MS2
  cols <- c("ms1_moverz", "ms1_int", "ms1_charge", 
            "ms2_moverzs", "ms2_ints", "ms2_charges")
  
  if (!all(cols %in% names(df)))
    stop("Not all required columns found for tracing MS features.")
  
  dfs <- split(df, df$iso_ctr)
  rm(list = "df")
  n_cores <- max(min(parallel::detectCores() - 1L, ceiling(length(dfs)/2L)), 1L)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  dfs <- parallel::clusterApply(cl, dfs, traceDIA, 
                                min_mass = min_mass, 
                                min_ms2mass = min_ms2mass, 
                                step1 = ppm_ms1/1e6, step2 = ppm_ms2/1e6, 
                                n_dia_ms2s = n_dia_ms2s)
  parallel::stopCluster(cl)
  df <- dplyr::bind_rows(dfs)
  rm(list = "dfs")
  
  ### 
  lens <- lengths(df$ms1_moverz)
  df <- df[lens > 0L, ] # 69701
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
  
  if (any(oks <- lens > 2L)) {
    ps0 <- ps[!oks]
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


#' Finds the gates of retention times.
#' 
#' @param ys A vector of corresponding signal intensities.
#' @param gap The size of gap for LC peaks collapsion.
#' 
#' @examples
#' library(mzion)
#' 
#' # find_lc_gates(c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 1), 20, 50))
#' # find_lc_gates(c(rep(0, 7), 100, 101, rep(0, 2), seq(200, 500, 100), rep(0, 4), 20, 50))
#' 
#' # all discrete
#' # find_lc_gates(c(rep(0, 5), 100, rep(0, 6), 200, rep(0, 4), 50))
find_lc_gates <- function (ys, gap = 4L)
{
  xs <- which(ys > 0)
  yxs <- ys[xs]
  
  # all discrete
  if (is.null(edges <- find_gate_edges(xs)))
    return(.Internal(which(ys > 0)))
  
  ups <- edges[["ups"]]
  dns <- edges[["dns"]]
  rm(list = "edges")
  
  # all discrete (already handled above)
  if (FALSE && is.null(ups)) {
    xus <- xs
    len <- length(xs)
    
    if (len > 1L) {
      for (i in 2:len) {
        ix <- i - 1
        gi <- xus[[i]] - xus[[ix]]
        
        if ((!is.na(gi)) && gi <= gap) {
          if (yxs[[i]] >= yxs[[ix]])
            xus[[ix]] <- NA_integer_
          else
            xus[[i]] <- NA_integer_
        }
      }
    }
    
    return(xus[!is.na(xus)])
  }
  
  xus <- xs[ups]
  xds <- xs[dns]
  ps <- mapply(function (x, y) x:y, ups, dns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  len <- length(ps)
  xps <- rep_len(NA_integer_, len)
  ymaxs <- rep_len(NA_real_, len)
  
  for (i in seq_len(len)) {
    pi <- ps[[i]]
    xi <- xs[pi]
    yi <- yxs[pi]
    mi <- which.max(yi)
    xps[[i]] <- xi[[1]] + mi - 1L
    ymaxs[[i]] <- yi[[mi]]
  }
  
  if (len > 1L) {
    for (i in 2:len) {
      ix <- i - 1
      gi <- xus[[i]] - xds[[ix]]
      
      if (gi <= gap) {
        if (ymaxs[[i]] >= ymaxs[[ix]])
          xus[[ix]] <- NA_integer_
        else
          xus[[i]] <- NA_integer_
      }
    }
  }
  
  oks <- !is.na(xus)
  # list(xvals = xps[oks], yvals = ymaxs[oks])
  xps[oks]
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


#' Traces MS peaks.
#' 
#' For both moverzs and intensities.
#' 
#' @param xs Vectors of moverzs.
#' @param ys Vectors of intensities.
#' @param zs Vectors of charge states.
#' @param from The starting value for mass binning.
#' @param step A step size for mass binning.
#' @param gap The gap size for collapsing moverzs.
#' @param reord Logical; re-order data or not.
trace_ms <- function (xs, ys, zs, from = 115L, step = 1E-5, gap = 4L, 
                      reord = TRUE)
{
  if (reord) {
    ords <- lapply(xs, order)
    xs <- mapply(function (x, ord) x[ord], xs, ords, 
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ys <- mapply(function (x, ord) x[ord], ys, ords, 
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    zs <- mapply(function (x, ord) x[ord], zs, ords, 
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    rm(list = "ords")
  }
  
  ans <- collapse_xyz(xs = xs, ys = ys, zs = zs, lwr = from, step = step)
  ansx <- ans[["xvals"]]
  ansy <- ans[["yvals"]]
  ansz <- ans[["zvals"]]
  rm(list = c("ans"))
  
  nrc <- dim(ansy)
  nr <- nrc[[1]]
  nc <- nrc[[2]]
  null_charge <- null_mass <- rep_len(NA_integer_, nr)
  null_inty <- rep_len(NA_real_, nr)
  
  for (i in 1:nc) {
    rows <- find_lc_gates(ansy[, i], gap = gap)
    x1 <- ansx[rows, i]
    y1 <- ansy[rows, i]
    z1 <- ansz[rows, i]
    
    ansx[, i] <- null_mass
    ansy[, i] <- null_inty
    ansz[, i] <- null_charge
    ansx[rows, i] <- x1
    ansy[rows, i] <- y1
    ansz[rows, i] <- z1
  }
  
  list(x = ansx, y = ansy, z = ansz)
}


#' Collapse MS1 intensities.
#' 
#' Allow adjacent values in unv and later collapse adjacent columns/values.
#' 
#' @param xs Vectors of SORTED m-over-z values.
#' @param ys Vectors of intensity values corresponding to xs.
#' @param zs Vectors of charge-state vaues corresponding to xs.
#' @param lwr The lower mass limit.
#' @param step The bin size in converting numeric m-over-z values to integers.
#' @param coll Logical; to further collapse results or not.
#' @importFrom fastmatch %fin%
#' @examples
#' # Twos adjacent bins of xs: 392.1796, 392.1845
#' # xs longer than the universe
#' xs <- list(c(391.1883,391.2848,391.6627,391.6995,392.1646,392.1796,392.1845,
#'            392.2030,392.2887,392.6641,392.7812,393.0833,393.2975))
#' ys <- list(c(12827.41,337002.19,617819.69,18045.10,205851.53,15194.98,11318.61,
#'              12970.02,118604.48,75726.89,11676.51,23723.18,55749.93))
#' # collapse_mms1ints(xs, ys, 389.6529)
#' 
#' xs <- list(c(400.6596,401.7076,402.1813,402.1944,402.1969,402.2094,402.5438,402.7112,403.1812,404.1777), 
#'            c(400.6599,401.7075,402.1954,402.1975,402.7112,403.1822,404.2777))
#' ys <- list(c(24003.98,53431.96,110619.26,10988.55,12291.00,140045.06,67601.16,11413.04,21651.61,16686.06), 
#'            c(10000.1,40000.1,20000.1,50000.1,2500.2,5000.1,30000.1))
#' # collapse_mms1ints(xs, ys, 400.1994)
#' 
#' xs <- ys <- vector("list", 13L)
#' xs[[7]] <- 954.607849; xs[[8]] <- 954.630249; xs[[10]] <- 954.622925
#' ys[[7]] <- 15706.2627; ys[[8]] <- 19803.5879; ys[[10]] <- 31178.9648
#' # collapse_mms1ints(xs, ys, 951.089731)
collapse_mms1ints <- function (xs = NULL, ys = NULL, zs = NULL, lwr = 115L, 
                               step = 1e-5, coll = TRUE)
{
  null_out <- if (coll)
    list(xvals = NULL, yvals = NULL, n = NULL)
  else
    list(xvals = NULL, yvals = NULL)
  
  # all xs are NULL
  if (!any(oks <- lengths(xs) > 0L))
    return(null_out)
  
  xs <- xs[oks]
  ys <- ys[oks]
  
  # remove zero intensities
  oky <- lapply(ys, `>`, 0)
  xs <- mapply(function (x, i) x[i], xs, oky, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ys <- mapply(function (x, i) x[i], ys, oky, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  if (!any(oks <- lengths(xs) > 0L))
    return(null_out)
  
  xs <- xs[oks]
  ys <- ys[oks]
  ixs <- lapply(xs, index_mz, lwr, step)
  
  # remove duplicated ixs
  for (i in seq_along(ixs)) {
    ix <- ixs[[i]]
    x <- xs[[i]]
    y <- ys[[i]]
    oks <- !duplicated(ix)
    ixs[[i]] <- ix[oks]
    xs[[i]]  <- x[oks]
    ys[[i]]  <- y[oks]
  }
  
  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  
  # map ixs onto unv (presence or absence)
  # note one-to-one correspondence between ixs and xs
  xps <- lapply(ixs, function (x) unv %fin% x)
  x0 <- y0 <- rep(NA_real_, length(unv))
  
  xout <- mapply(function (i, v) {
    x0[i] <- v
    x0
  }, xps, xs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  yout <- mapply(function (i, v) {
    y0[i] <- v
    y0
  }, xps, ys, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  xout <- do.call(rbind, xout)
  yout <- do.call(rbind, yout)
  
  # all discrete values
  if (is.null(ps <- find_gates(unv))) {
    if (coll)
      return(calc_ms1xys(xout, yout))
    else
      return(list(xvals = xout, yvals = yout))
  }
  
  ps1 <- lapply(ps, `[[`, 1)
  ps2 <- lapply(ps, `[[`, 2)
  ps1 <- .Internal(unlist(ps1, recursive = FALSE, use.names = FALSE))
  ps2 <- .Internal(unlist(ps2, recursive = FALSE, use.names = FALSE))
  
  for (i in seq_along(ps)) {
    cols <- ps[[i]]
    
    # possible to have values in both columns but simply overwrite: 1 <- 2
    if (c2 <- cols[2]) {
      c1 <- cols[1]
      oks <- !is.na(xout[, c2])
      xout[oks, c1] <- xout[oks, c2]
      yout[oks, c1] <- yout[oks, c2]
    }
  }
  
  xout <- xout[, -ps2, drop = FALSE]
  yout <- yout[, -ps2, drop = FALSE]
  
  if (coll)
    calc_ms1xys(xout, yout)
  else
    list(xvals = xout, yvals = yout)
}


#' Collapse MS data.
#' 
#' Almost identical to \link{collapse_mms1ints} with the addition of zs.
#' 
#' @param zs Vectors of charge states.
#' @inheritParams collapse_mms1ints
#' @importFrom fastmatch %fin%
collapse_xyz <- function (xs = NULL, ys = NULL, zs = NULL, lwr = 115L, 
                          step = 1e-5)
{
  null_out <- list(xvals = NULL, yvals = NULL, zvals = NULL)
  
  # all xs are NULL
  if (!any(oks <- lengths(xs) > 0L))
    return(null_out)
  
  xs <- xs[oks]
  ys <- ys[oks]
  zs <- zs[oks]
  
  # remove zero intensities
  oky <- lapply(ys, `>`, 0)
  xs <- mapply(function (x, i) x[i], xs, oky, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ys <- mapply(function (x, i) x[i], ys, oky, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  zs <- mapply(function (x, i) x[i], zs, oky, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  if (!any(oks <- lengths(xs) > 0L))
    return(null_out)
  
  xs <- xs[oks]
  ys <- ys[oks]
  zs <- zs[oks]
  ixs <- lapply(xs, index_mz, lwr, step)
  
  # remove duplicated ixs
  for (i in seq_along(ixs)) {
    ix <- ixs[[i]]
    x <- xs[[i]]
    y <- ys[[i]]
    z <- zs[[i]]
    oks <- !duplicated(ix)
    ixs[[i]] <- ix[oks]
    xs[[i]]  <- x[oks]
    ys[[i]]  <- y[oks]
    zs[[i]]  <- z[oks]
  }
  
  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  
  xps <- lapply(ixs, function (x) unv %fin% x)
  z0 <- x0 <- y0 <- rep(NA_real_, length(unv))
  
  xout <- mapply(function (i, v) {
    x0[i] <- v
    x0
  }, xps, xs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  yout <- mapply(function (i, v) {
    y0[i] <- v
    y0
  }, xps, ys, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  zout <- mapply(function (i, v) {
    z0[i] <- v
    z0
  }, xps, zs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  xout <- do.call(rbind, xout)
  yout <- do.call(rbind, yout)
  zout <- do.call(rbind, zout)
  
  # all discrete values
  if (is.null(ps <- find_gates(unv)))
    return(list(xvals = xout, yvals = yout, zvals = zout))
  
  ps1 <- lapply(ps, `[[`, 1)
  ps2 <- lapply(ps, `[[`, 2)
  ps1 <- .Internal(unlist(ps1, recursive = FALSE, use.names = FALSE))
  ps2 <- .Internal(unlist(ps2, recursive = FALSE, use.names = FALSE))
  
  for (i in seq_along(ps)) {
    cols <- ps[[i]]
    
    # possible to have values in both columns but simply overwrite: 1 <- 2
    if (c2 <- cols[2]) {
      c1 <- cols[1]
      oks <- !is.na(xout[, c2])
      xout[oks, c1] <- xout[oks, c2]
      yout[oks, c1] <- yout[oks, c2]
      zout[oks, c1] <- zout[oks, c2]
    }
  }
  
  xout <- xout[, -ps2, drop = FALSE]
  yout <- yout[, -ps2, drop = FALSE]
  zout <- zout[, -ps2, drop = FALSE]
  
  list(xvals = xout, yvals = yout, zvals = zout)
}


#' Makes MS1 matrixes of moverzs and intensities
#' 
#' @param xs Vectors of moverzs
#' @param ys Vectors of intensities
calc_ms1xys <- function (xs, ys)
{
  ysum <- colSums(ys, na.rm = TRUE)
  xmeans <- colSums(xs * ys, na.rm = TRUE)/ysum
  n <- colSums(!is.na(ys), na.rm = TRUE)
  list(xvals = xmeans, yvals = ysum, n = n)
}


#' Finds mDDA precursors.
#'
#' Averages of multiple MS1 scans.
#'
#' @param msx_moverzs Vectors of moverzs (full spectrum).
#' @param msx_ints Vectors of intensities (full spectrum.
#' @param ms1_moverzs Vectors of MS1 moverzs (associated with MS2).
#' @param ms1_charges Vectors of MS1 charge states (associated with MS2).
#' @param ms1_ints Vectors of MS1 intensities (associated with MS2).
#' @param iso_ctr Vectors of isolation centers.
#' @param iso_lwr Vectors of isolation lowers.
#' @param ppm Mass error tolerance.
#' @param step A step size. 
#' @param maxn_precurs Maximum number of precursors for consideration.
#' @param max_ms1_charge Maximum charge state of precursors for consideration.
#' @param width The width of an MS1 window. A wide window is used for containing
#'   isotope envelops.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries.
#' @param use_defpeaks Use default peak info or not.
#' @inheritParams matchMS
find_mdda_mms1s <- function (msx_moverzs = NULL, msx_ints = NULL, 
                             ms1_moverzs = NULL, ms1_charges = NULL, ms1_ints = NULL, 
                             iso_ctr = NULL, iso_lwr = NULL, 
                             ppm = 10L, maxn_precurs = 5L, max_ms1_charge = 4L, 
                             n_fwd = 20L, grad_isotope = 1.6, fct_iso2 = 3.0, 
                             use_defpeaks = FALSE, width = 2.01, step = ppm/1e6)
{
  # for all (6+1+6) MS1 frames subset by one MS2 iso-window
  ansx1 <- ansy1 <- vector("list", len1 <- length(msx_moverzs))
  # for all MS2s from averaged (6+1+6 -> 1) MS1s 
  ansn2 <- ansx2 <- ansy2 <- vector("list", len2 <- length(iso_ctr))
  
  # go through MS2
  for (i in 1:len2) {
    m2  <- iso_ctr[[i]]
    lwr <- m2 - width
    upr <- m2 + width
    
    # go through MS1s
    for (j in 1:len1) {
      x1s <- msx_moverzs[[j]]
      y1s <- msx_ints[[j]]
      oks <- x1s > lwr & x1s < upr
      ansx1[[j]] <- x1s[oks]
      ansy1[[j]] <- y1s[oks]
    }
    
    # collapse MS1s
    ans <- collapse_mms1ints(xs = ansx1, ys = ansy1, lwr =  lwr, step = step)
    ansx <- ans[["xvals"]]
    ansy <- ans[["yvals"]]
    ansn <- ans[["n"]]
    
    if (!is.null(ansx)) {
      ansx2[[i]] <- ansx
      ansy2[[i]] <- ansy
      ansn2[[i]] <- ansn
    }
  }
  
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
    ), 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  masses <- lapply(mics, `[[`, "masses")
  charges <- lapply(mics, `[[`, "charges")
  intensities <- lapply(mics, `[[`, "intensities")
  
  # (2) subset by isolation window
  # ( `width = 2.01` contains isotope envelope and now need subsetting)
  moks <- mapply(function (x, m, w) {
    if (any(oks <- x > m - w & x < m + w)) oks else rep_len(TRUE, length(x))
  }, x = masses, m = iso_ctr, w = iso_lwr, 
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  # masses may be NULL and the corresponding moks is logical(0)
  masses <- mapply(function (x, i) x[i], x = masses, i = moks, 
                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
  charges <- mapply(function (x, i) x[i], x = charges, i = moks, 
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  intensities <- mapply(function (x, i) x[i], x = intensities, i = moks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # impurities <- lapply(intensities, function (x) x/sum(x))
  rm(list = "moks")
  
  # (3) updates
  rows <- lengths(masses) > 0L
  ms1_moverzs[rows] <- masses[rows]
  ms1_charges[rows] <- charges[rows]
  ms1_ints[rows] <- intensities[rows]
  
  list (ms1_moverzs = ms1_moverzs, 
        ms1_charges = ms1_charges, 
        ms1_ints = ms1_ints)
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
  
  if (!(len_ms <- length(moverzs))) return(na_out)
  
  oks <- (moverzs >= center - iso_lwr) & (moverzs <= center + iso_upr)
  moverzs <- moverzs[oks]
  if (!(len_ms <- length(moverzs))) return(na_out)
  msxints <- msxints[oks]
  charges <- charges[oks]
  masses <- (moverzs - 1.00727647) * charges
  
  if (len_ms == 1L)
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  
  okc <- !is.na(charges)
  charges <- charges[okc]
  if (!(len_ms <- length(charges))) return(na_out)
  moverzs <- moverzs[okc]
  msxints <- msxints[okc]
  masses <- masses[okc]
  
  if (len_ms == 1L)
    return(list(ms1_moverzs = moverzs, ms1_masses = masses, 
                ms1_charges = charges, ms1_ints = msxints))
  
  idx <- which.max(msxints)
  moverzs <- moverzs[idx]
  charges <- charges[idx]
  msxints <- msxints[idx]
  masses <- masses[idx]
  
  list(ms1_moverzs = moverzs, ms1_masses = masses, 
       ms1_charges = charges, ms1_ints = msxints)
}


#' Separates a vector into groups.
#' 
#' @param vec A vector.
#' @param n_chunks The number of chunks.
sep_vec <- function (vec, n_chunks)
{
  seqs <- seq_along(vec)
  labs <- levels(cut(seqs, n_chunks))
  lwrs <- floor(as.numeric( sub("\\((.+),.*", "\\1", labs)))
  grps <- findInterval(seqs, lwrs)
}

