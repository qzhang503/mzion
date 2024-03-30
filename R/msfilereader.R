readRAW <- function (filepath = NULL, filelist = NULL, out_path = NULL, 
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
                     use_defpeaks = FALSE, maxn_dia_precurs = 1000L, 
                     maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                     ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                     quant = "none") 
{
  # filepath = mgf_path
  local({
    fn_suffix <- tolower(gsub("^.*\\.([^.]*)$", "\\1", filelist[[1]]))
    if (fn_suffix == "raw") {
      data_format <- "Thermo-RAW"
      mgf_format <- "Mzion"
    }
    
    qs::qsave(list(data_format = data_format, mgf_format = mgf_format), 
              file.path(filepath, "info_format.rds"), preset = "fast")
  })
  
  len <- length(filelist)
  n_cores <- min(len, detect_cores(32L))
  
  if (n_cores == 1L) {
    raw_files <- vector("list", len)
    # filepeaks <- lapply(filelist[[1]], exeReadRAW, filepath)
    filepeaks <- list("01CPTAC3_Benchmarking_W_BI_20170508_BL_f02.peaks")
    
    i = 1
    file <- filepeaks[[i]]

    for (i in 1:len) {
      
    }
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    raw_files <- parallel::clusterMap(cl, proc_raws, 
                                      1:len, filelist, 
                                      MoreArgs = list(filepath = filepath, 
                                                      raw_file = raw_file), 
                                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    parallel::stopCluster(cl)
  }
  
}


#' Helper in processing MGF entries in chunks.
#'
#' @param lines Lines of MGF.
#' @inheritParams proc_mgf_chunks
proc_raws <- function (lines, topn_ms2ions = 150L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       
                       
                       type_mgf = "mzion_thermo", n_bf_begin = 0L, 
                       n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L,
                       n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L,
                       n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                       sep_pepmass = " ", nfields_pepmass = 2L, raw_file = NULL, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 1000L, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       quant = "none", digits = 4L) 
{
  options(digits = 9L)
  
  ###
  if (FALSE) {
    message("Parsing '", file, "'.")
    lines <- stringi::stri_read_lines(file)
    basename <- gsub("\\.[^.]*$", "", file)
    
    idx_scans <- .Internal(which(stringi::stri_cmp_eq(lines, "SCAN")))
    scan_nums <- lines[idx_scans + 1L]
    ret_times <- as.numeric(lines[idx_scans + 3L]) * 60
    low_masses <- as.numeric(lines[idx_scans + 5L])
    high_masses <- as.numeric(lines[idx_scans + 7L])
    iso_ctr <- as.numeric(lines[idx_scans + 9L])
    ms_levels <- as.integer(lines[idx_scans + 11L]) # 0 at exception
    iso_width <- as.numeric(lines[idx_scans + 13L])
    msx_ns <- as.integer(lines[idx_scans + 15L])
    
    msx_moverzs <- mapply(
      function (s, n) stringi::stri_split_fixed(s, pattern = ",", n = n, 
                                                simplify = TRUE), 
      lines[idx_scans + 17L], msx_ns, 
      SIMPLIFY = TRUE, USE.NAMES = FALSE)
    
    msx_ints <- mapply(
      function (s, n) stringi::stri_split_fixed(s, pattern = ",", n = n, 
                                                simplify = TRUE), 
      lines[idx_scans + 19L], msx_ns, 
      SIMPLIFY = TRUE, USE.NAMES = FALSE)
    
    scan_titles <- lines[idx_scans + 21L]
    
    
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
  }
  ###
  
  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends   <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))
  
  ## MS1 
  # (1) m-over-z and intensity
  ms1s <- stringi::stri_replace_first_fixed(lines[begins + n_to_pepmass], 
                                            "PEPMASS=", "")
  ms1s <- lapply(ms1s, stringi::stri_split_fixed, pattern = sep_pepmass, 
                 n = nfields_pepmass, simplify = TRUE)
  ms1_moverzs <- lapply(ms1s, function (x) as.numeric(x[, 1]))
  ms1_moverzs <- .Internal(unlist(ms1_moverzs, recursive = FALSE, use.names = FALSE))
  # not as.integer; intensity may be > .Machine$integer.max (2147483647)
  ms1_ints <- lapply(ms1s, function (x) as.numeric(x[, 2]))
  ms1_ints <- .Internal(unlist(ms1_ints, recursive = FALSE, use.names = FALSE))
  rm(list = c("ms1s"))
  
  # (2) retention time
  ret_times <- stringi::stri_replace_first_fixed(lines[begins + n_to_rt], "RTINSECONDS=", "")
  ret_times <- as.numeric(ret_times)
  
  # (3) MS1 charges and masses
  ms1_charges <- stringi::stri_replace_first_fixed(lines[begins + n_to_charge], "CHARGE=", "")
  ms1_charges <- gsub("[\\+]$", "", ms1_charges)
  if (FALSE) {
    if ((polarity <- gsub(".*([\\+-])$", "\\1", ms1_charges[[1]])) == "+")
      ms1_charges <- gsub("[\\+]$", "", ms1_charges)
    else if (polarity == "-")
      ms1_charges <- gsub("-$", "", ms1_charges)
  }
  ms1_charges <- as.integer(ms1_charges)
  
  # timsTOF no CHARGE line -> NAs introduced by coercion
  ms1_masses <- mapply(function (x, y) x * y - y * 1.00727647, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)
  
  ans_rts <- reset_rettimes(ret_times = ret_times, min_ret_time = min_ret_time, 
                            max_ret_time = max_ret_time)
  min_ret_time <- ans_rts$min_ret_time
  max_ret_time <- ans_rts$max_ret_time
  
  rows <- (ms1_charges >= min_ms1_charge & ms1_charges <= max_ms1_charge & 
             ret_times >= min_ret_time & ret_times <= max_ret_time & 
             ms1_masses >= min_mass & ms1_masses <= max_mass & 
             # timsTOF: no MS1 masses
             !is.na(ms1_masses))
  
  # timsTOF may have undetermined charge states
  if (length(na_rows <- .Internal(which(is.na(rows))))) 
    rows[na_rows] <- FALSE
  
  begins <- begins[rows]
  ends <- ends[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- round(ret_times[rows], digits = 2L)
  ms1_masses <- ms1_masses[rows]
  rm(list = "rows")
  
  # Others
  scan_titles <- 
    stringi::stri_replace_first_fixed(lines[begins + n_to_title], "TITLE=", "")
  
  if (type_mgf %in% c("msconv_thermo", "msconv_pasef")) {
    raw_files <- 
      stringi::stri_replace_first_regex(scan_titles, "^.* File:\"([^\"]+)\".*", "$1")
    scan_nums <- 
      stringi::stri_replace_first_regex(scan_titles, 
                                        "^.*\\.(\\d+)\\.\\d+\\.\\d+ File:\".*", 
                                        "$1")
  } 
  else if (type_mgf == "pd") {
    raw_files <- gsub("^.*File: \"([^\"]+)\".*", "\\1", scan_titles)
    raw_files <- gsub("\\\\", "/", raw_files)
    raw_files <- gsub("^.*/(.*)", "\\1", raw_files)
    scan_nums <- gsub("^.* scans: \"([0-9]+)\"$", "\\1", scan_titles)
  } 
  else if (type_mgf == "default_pasef") {
    # one raw_file one .d file guaranteed
    raw_files <- rep_len(raw_file, length.out = length(begins))
    scan_nums <- stringi::stri_replace_first_fixed(lines[begins + n_to_scan], "RAWSCANS=", "")
  } 
  else {
    stop("Unknown MGF format.")
  }
  
  ## Low priority: no data subsetting by scan_nums; use retention time instead
  
  ## MS2
  # (-1L: one line above "END IONS")
  ms2s <- mapply(function (x, y) lines[(x + n_hdr) : (y - 1L)], 
                 begins, ends, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ms2s <- lapply(ms2s, stringi::stri_split_fixed, pattern = sep_ms2s, 
                 n = nfields_ms2s, simplify = TRUE)
  ms2_moverzs <- lapply(ms2s, function (x) as.numeric(x[, 1]))
  # not as.integer; intensity may be > .Machine$integer.max
  ms2_ints <- lapply(ms2s, function (x) as.numeric(x[, 2]))
  rm(list = c("ms2s"))
  
  ms2_charges <- vector("list", length(ms2_moverzs))
  is_tmt <- if (grepl("^tmt.*\\d+", quant)) TRUE else FALSE
  
  if (deisotope_ms2) {
    mics <- mapply(find_ms1stat, moverzs = ms2_moverzs, msxints = ms2_ints, 
                   MoreArgs = list(
                     center = 0, 
                     exclude_reporter_region = is_tmt, 
                     tmt_reporter_lower = tmt_reporter_lower, 
                     tmt_reporter_upper = tmt_reporter_upper, 
                     ppm = ppm_ms2_deisotope, 
                     ms_lev = 2L, maxn_feats = topn_ms2ions, 
                     max_charge = max_ms2_charge, n_fwd = 10L, 
                     offset_upr = 30L, offset_lwr = 30L, order_mz = FALSE
                   ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ms2_moverzs <- lapply(mics, `[[`, "masses")
    ms2_ints <- lapply(mics, `[[`, "intensities")
    ms2_charges <- lapply(mics, `[[`, "charges")
  }
  
  # extract the TMT region of MS2 moverz and intensity
  # (also convert reporter-ion intensities to integers)
  restmt <- extract_mgf_rptrs(ms2_moverzs, 
                              ms2_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  
  ms2_moverzs <- restmt[["xvals"]]
  ms2_ints <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  rm(list = "restmt")
  
  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = ms2_moverzs, 
                          ms2_ints = ms2_ints, 
                          ms2_charges = ms2_charges, 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  
  ms2_moverzs <- mz_n_int[["ms2_moverzs"]]
  ms2_ints <- mz_n_int[["ms2_ints"]]
  ms2_charges <- mz_n_int[["ms2_charges"]]
  lens <- mz_n_int[["lens"]]
  rm(list = "mz_n_int")
  
  df <- tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_files,
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_int = ms1_ints,
    ms1_charge = ms1_charges,
    ret_time = ret_times,
    scan_num = scan_nums,
    ms2_moverzs= ms2_moverzs,
    ms2_ints = ms2_ints,
    ms2_charges = ms2_charges, 
    ms2_n = lens, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints, )
}





#' Helper in the execution of CSharp-compiled ReadRAW.exe
#' 
#' @param raw_file A RAW file name.
#' @param type The type of MSFileReader request.
#' @param mgf_path A file path to peak lists or RAW MS data.
exeReadRAW <- function(raw_file, mgf_path, type = "scans")
{
  if (FALSE) {
    type <- "scans"
    mgf_path <- "~/mzion/BI_G1_test/mgf"
    mgf_path <- find_dir(mgf_path)
    raw_file <- file.path(mgf_path, "01CPTAC3_Benchmarking_W_BI_20170508_BL_f02.raw")
  }

  raw_full <- file.path(mgf_path, raw_file)
  sys_path <- system.file("extdata", package = "mzion")
  acceptMSFileReaderLicense(sys_path)
  exe <- file.path(sys_path, name_exe = "ReadRAW.exe")
  mono <- if (Sys.info()['sysname'] %in% c("Darwin", "Linux")) TRUE else FALSE
  
  peaks <- tempfile(tmpdir = mgf_path, fileext = ".peaks")
  stdout <- tempfile(tmpdir = mgf_path, fileext = ".stdout")
  stderr <- tempfile(tmpdir = mgf_path, fileext = ".stderr" )

  if (mono) {
    rvs <- system2(Sys.which("mono"), 
                   args = c(shQuote(exe),
                            shQuote(raw_full),
                            type, 
                            shQuote(peaks)),
                   stdout = stdout,
                   stderr = stderr)
  } else{
    rvs <- system2(exe, 
                   args = c(shQuote(raw_full), 
                            type, 
                            shQuote(peaks)), 
                   stdout = stdout,
                   stderr = stderr)
  }
  
  if (rvs || !file.exists(peaks))
    stop("Fail to process ", raw_file, ".")

  raw_name <- readLines(peaks, 2L)[[2]]
  new_name <- paste0(gsub("\\.[^.]*$", "", raw_name), ".peaks")
  raw_temp <- fs::path_file(peaks)
  file.rename(file.path(mgf_path, raw_temp), file.path(mgf_path, new_name))
  unlink(c(stdout, stderr))
  
  new_name
}


#' License agreement.
#' 
#' @param sys_path The system path of mzion.
acceptMSFileReaderLicense <- function(sys_path) 
{
  license  <- file.path(sys_path, name_exe = "RawFileReaderLicense.txt")
  
  if (!dir.exists(libdir <- tools::R_user_dir("mzion", which = 'cache')))
    dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(file.path(libdir, "rspn.txt")))
    return(TRUE)
  
  msg <- "By selecting YES, you are accepting the Thermo's License agreement."
  file.show(license)
  response <- readline(prompt = sprintf("Accept '%s'? [y/n]: ", license))
  
  if (tolower(response) != "y")
    stop("You have to accept the License agreement!")
  
  rspnFile <- file.path(libdir, "rspn.txt")
  fileConn <- file(rspnFile)
  writeLines(paste(msg, paste0(date()), "Response = TRUE", sep = "\n"), fileConn)
  close(fileConn)
  
  return(TRUE)
}

