#' Reads Thermo RAW files.
#' 
#' @param mgf_path A file path RAW MS files.
#' @param filelist A list of RAW MS files.
readRAW <- function (mgf_path = NULL, filelist = NULL) 
{
  sys_path <- system.file("extdata", package = "mzion")
  acceptMSFileReaderLicense(sys_path)
  
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))
  
  local({
    fn_suffix <- tolower(gsub("^.*\\.([^.]*)$", "\\1", filelist[[1]]))
    
    if (fn_suffix == "raw") {
      data_format <- "Thermo-RAW"
      mgf_format <- "Mzion"
    }
    
    qs::qsave(list(data_format = data_format, mgf_format = mgf_format), 
              file.path(mgf_path, "info_format.rds"), preset = "fast")
  })
  
  len <- length(filelist)
  n_cores <- min(len, detect_cores(32L))
  
  if (n_cores == 1L) {
    filenames <- lapply(filelist, proc_raws, mgf_path, temp_dir)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    filenames <- parallel::clusterApply(cl, filelist, proc_raws, mgf_path, temp_dir)
    parallel::stopCluster(cl)
  }
  
  filenames
}


#' Helper in processing MGF entries in chunks.
#'
#' @param raw_file The file name of RAW MS data.
#' @param mgf_path The file path of MGF files.
#' @param temp_dir A file path for temporary files. 
proc_raws <- function (raw_file, mgf_path, temp_dir) 
{
  options(digits = 9L)
  
  filepeak <- exeReadRAW(raw_file, mgf_path)
  pathfile <- file.path(mgf_path, filepeak)
  lines <- readLines(pathfile)
  # DON'T; see the help document; memory allocation or access error
  # lines <- stringi::stri_read_lines(pathfile)

  idx_scans <- .Internal(which(stringi::stri_cmp_eq(lines, "SCAN")))
  scan_nums <- lines[idx_scans + 1L]
  ret_times <- as.numeric(lines[idx_scans + 3L]) * 60
  ms_levels <- as.integer(lines[idx_scans + 5L]) # 0 at exception
  iso_ctrs <- as.numeric(lines[idx_scans + 7L]) # 0 if MS1
  iso_widths <- as.numeric(lines[idx_scans + 9L]) # large values if MS1
  msx_ns <- as.integer(lines[idx_scans + 11L])
  
  msx_moverzs <- mapply(
    function (s, n) 
      as.numeric(
        stringi::stri_split_fixed(s, pattern = ",", n = n, simplify = TRUE)), 
    lines[idx_scans + 13L], msx_ns, 
    SIMPLIFY = TRUE, USE.NAMES = FALSE)
  
  msx_ints <- mapply(
    function (s, n) 
      as.numeric(
        stringi::stri_split_fixed(s, pattern = ",", n = n, simplify = TRUE)), 
    lines[idx_scans + 15L], msx_ns, 
    SIMPLIFY = TRUE, USE.NAMES = FALSE)
  
  scan_titles <- lines[idx_scans + 17L]
  half_widths <- iso_widths / 2
  iso_lwrs <- iso_ctrs - half_widths
  iso_uprs <- iso_ctrs + half_widths
  
  len <- length(msx_moverzs)
  na_ints <- rep_len(NA_integer_, len)
  na_reals <- rep_len(NA_real_, len)
  
  out <- list(
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    ms1_moverzs = na_reals, 
    ms1_charges = na_ints, 
    ms1_ints = na_reals, 
    
    scan_title = scan_titles,
    raw_file = raw_file, # scalar
    ms_level = ms_levels, 
    ret_time = ret_times, 
    scan_num = as.integer(scan_nums), 
    orig_scan = scan_nums,
    iso_ctr = as.numeric(iso_ctrs), 
    iso_lwr = as.numeric(iso_lwrs), 
    iso_upr = as.numeric(iso_uprs)
  )
  
  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  
  idx_first_ms1 <- which(ms_levels == 1L)[[1]]
  iso_ctr1 <- iso_ctrs[[idx_first_ms1]]
  attr(out_name, "is_dia") <- if (iso_ctr1 < 0) TRUE else FALSE
  attr(out_name, "mzml_type") <- "raw"
  
  invisible(out_name)
}


#' Helper in executing CSharp-compiled ReadRAW.exe
#' 
#' Requires mono for Linux or Mac OS.
#' 
#' @param raw_file A RAW file name.
#' @param mgf_path A file path RAW MS files.
exeReadRAW <- function(raw_file, mgf_path)
{
  sys_path <- system.file("extdata", package = "mzion")
  exe <- file.path(sys_path, name_exe = "ReadRAW.exe")
  mono <- if (Sys.info()['sysname'] %in% c("Darwin", "Linux")) TRUE else FALSE
  
  peaks  <- tempfile(tmpdir = mgf_path, fileext = ".peaks")
  stdout <- tempfile(tmpdir = mgf_path, fileext = ".stdout")
  stderr <- tempfile(tmpdir = mgf_path, fileext = ".stderr" )
  raw_full <- file.path(mgf_path, raw_file)
  
  if (mono) {
    rvs <- system2(Sys.which("mono"), 
                   args = c(shQuote(exe),
                            shQuote(raw_full),
                            shQuote(peaks)),
                   stdout = stdout,
                   stderr = stderr)
  }
  else{
    rvs <- system2(exe, 
                   args = c(shQuote(raw_full), 
                            shQuote(peaks)), 
                   stdout = stdout,
                   stderr = stderr)
  }
  
  if (rvs || !file.exists(peaks)) {
    stop("Fail to process ", raw_file, ".")
  }
  
  out_name <- paste0(gsub("\\.[^.]*$", "", raw_file), ".peaks")
  file.rename(peaks, file.path(mgf_path, out_name))
  unlink(c(stdout, stderr))
  
  out_name
}


#' License agreement.
#' 
#' @param sys_path The system path of Mzion.
acceptMSFileReaderLicense <- function(sys_path) 
{
  license  <- file.path(sys_path, name_exe = "RawFileReaderLicense.txt")
  libdir <- create_dir(tools::R_user_dir("mzion", which = 'cache'))

  if (file.exists(file.path(libdir, "rspn.txt")))
    return(TRUE)
  
  msg <- "By selecting YES, you are accepting the Thermo's License agreement."
  file.show(license)
  response <- readline(prompt = sprintf("Accept '%s'? [y/n]: ", license))
  
  if (tolower(response) != "y")
    stop("Need to accept the license agreement for processing Thermo's RAW.")
  
  rspnFile <- file.path(libdir, "rspn.txt")
  fileConn <- file(rspnFile)
  writeLines(paste(msg, paste0(date()), "Response = TRUE", sep = "\n"), fileConn)
  close(fileConn)
  
  NULL
}


