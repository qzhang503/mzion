#' Reads Bruker's PASEF files.
#' 
#' @param mgf_path A file path RAW MS files.
#' @param filelist A list of RAW MS files.
#' @inheritParams matchMS
readPASEF <- function (mgf_path = NULL, filelist = NULL, topn_ms2ions = 150L) 
{
  sys_path <- system.file("extdata", package = "mzion")
  acceptBrukerLicense(sys_path)
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))
  
  qs::qsave(list(data_format = "Bruker-RAW", mgf_format = "Mzion"), 
            file.path(mgf_path, "info_format.rds"), preset = "fast")

  len <- length(filelist)
  # n_cores <- min(len, detect_cores(4L))
  n_cores <- 1L
  
  if (n_cores == 1L) {
    filenames <- lapply(filelist, proc_pasefs, mgf_path, temp_dir, topn_ms2ions)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    filenames <- parallel::clusterApply(cl, filelist, proc_pasefs, mgf_path, 
                                        temp_dir, topn_ms2ions)
    parallel::stopCluster(cl)
  }
  
  filenames
}


#' Helper in processing MGF entries in chunks.
#'
#' @param raw_file The file name of RAW MS data.
#' @param mgf_path The file path of MGF files.
#' @param temp_dir A file path for temporary files. 
#' @inheritParams matchMS
proc_pasefs <- function (raw_file, mgf_path, temp_dir, topn_ms2ions = 150L) 
{
  options(digits = 9L)
  
  if (FALSE) {
    ans <- exeReadPASEF(raw_file, mgf_path, topn_ms2ions)
    spectra <- readLines(ans[[1]])
    # writeLines(spectra[1:12339321], file.path(mgf_path, raw_file, "spectra_sub.txt"))
    precursors <- readLines(ans[[2]])
  }
  
  lines <- readLines(file.path(mgf_path, raw_file, "spectra_sub.txt")) # 4035
  hdr <- lines[1:3]
  lines <- lines[-c(1:3)]
  idx_empts <- .Internal(which(stringi::stri_cmp_eq(lines, "TITLE"))) + 2L
  lines <- split(lines, cut(seq_along(lines), c(0, idx_empts)))
  
  n_frames <- length(lines)
  last_entry <- lines[[n_frames]]
  n_last <- length(last_entry)
  
  # in case all "" lines
  if (n_last <= 4L) {
    lines <- lines[-n_frames]
  }
  
  # adds a trailing line of ""
  if (last_entry[length(last_entry)] != "") {
    last_entry <- c(last_entry, "")
    lines[[n_frames]] <- last_entry
  }

  if (FALSE) {
    i <- 4023
    z <- extract_pasef_frame(lines[[i]])
    
    # ans <- vector("list", length(lines))
    ans <- vector("list", length(lines)) # 4023
    for (i in seq_along(ans)) {
      print(i)
      ans[[i]] <- extract_pasef_frame(lines[[i]])
      # print(i)
    }
    
    i
  }

  n_cores <- detect_cores(16L)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  ans <- 
    parallel::clusterApply(cl, chunksplit(lines, n_cores), hextract_pasef_frame)
  parallel::stopCluster(cl)
  
  if (FALSE) {
    # idx_xs <- .Internal(which(stringi::stri_cmp_eq(lines, "X")))
    # idx_ys <- .Internal(which(stringi::stri_cmp_eq(lines, "Y")))
    idx_npeaks <- .Internal(which(stringi::stri_cmp_eq(lines, "NPEAKS")))
    idx_mobs <- .Internal(which(stringi::stri_cmp_eq(lines, "MOBILITY")))
    
    # shared across scans
    idx_scans <- .Internal(which(stringi::stri_cmp_eq(lines, "SCAN"))) # 268
    ret_times <- as.numeric(lines[idx_scans + 3L]) * 60
    ms_levels <- as.integer(lines[idx_scans + 5L])
    scan_titles <- lines[idx_scans + 7L]
    idx_empties <- idx_scans + 8L
    
    ln_chunks <- split(lines, cut(seq_along(lines), c(0, idx_empties)))
  }
  
  # lines <- readLines(file.path(mgf_path, raw_file, "spectra.txt"))
  # precursors <- readr::read_tsv(file.path(mgf_path, raw_file, "ms2precursors.txt"))

  
  # collapse MS1 scans under the same frame...
  
  
  
  # out_name <- paste0(gsub("\\.[^.]*$", "", raw_file), ".peaks")
  # file.rename(spectra, file.path(mgf_path, out_name))
  # file.rename(precursors, file.path(mgf_path, out_name))
  # unlink(c(stdout, stderr))
  # out_name
  
  # pathfile <- file.path(mgf_path, filepeak)
  # lines <- readLines(pathfile)

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
        stringi::stri_split_fixed(s, pattern = " ", n = n, simplify = TRUE)), 
    lines[idx_scans + 13L], msx_ns, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  msx_ints <- mapply(
    function (s, n) 
      as.numeric(
        stringi::stri_split_fixed(s, pattern = " ", n = n, simplify = TRUE)), 
    lines[idx_scans + 15L], msx_ns, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
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


#' Extracts data from a PASEF frame.
#' 
#' @param data A PASEF frame (with multiple scans).
extract_pasef_frame <- function (data)
{
  len <- length(data)
  
  # fields common within a frame
  title <- data[len-1]
  ms_lev <- as.integer(data[len-3])
  ret_time <- as.numeric(data[len-5])
  scan_num <- as.integer(data[len-7])
  
  # fields specific to each scan in a frame
  idx_npeaks <- .Internal(which(stringi::stri_cmp_eq(data, "NPEAKS"))) + 1L
  idx_mobs <- idx_npeaks + 2L
  idx_ys <- idx_npeaks - 2L
  idx_xs <- idx_ys - 2L
  
  npeaks <- as.integer(data[idx_npeaks])
  mobils <- as.numeric(data[idx_mobs])
  mobils <- sum(mobils)/length(mobils)
  
  xs <- mapply(
    function (s, n) 
      as.numeric(
        stringi::stri_split_fixed(s, pattern = " ", n = n, simplify = TRUE)), 
    data[idx_xs], npeaks, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ys <- mapply(
    function (s, n) 
      as.numeric(
        stringi::stri_split_fixed(s, pattern = " ", n = n, simplify = TRUE)), 
    data[idx_ys], npeaks, 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  yco <- if (ms_lev == 1L) 75 else 10
  
  for (i in seq_along(xs)) {
    ysi <- ys[[i]]
    xsi <- xs[[i]]
    oki <- ysi >= yco
    ys[[i]] <- ysi[oki] # works at rhs numeric(0)
    xs[[i]] <- xsi[oki]
  }
  
  oks <- lengths(xs) > 0L
  xs <- xs[oks]
  ys <- ys[oks]
  
  if (!length(xs)) {
    return(
      list(title = title,
           ms_lev = ms_lev,
           ret_time = ret_time,
           scan_num = scan_num, 
           mobility = mobils, 
           xs = 0, ys = 0))
  }
  
  xs <- unlist(xs, recursive = FALSE, use.names = FALSE)
  ys <- unlist(ys, recursive = FALSE, use.names = FALSE)
  ord <- order(xs)
  xs <- xs[ord]
  ys <- ys[ord]
  
  if (anyDuplicated(xs)) {
    ys <- split(ys, xs)
    xs <- split(xs, xs)
    dups <- lengths(xs) > 1L
    xs[dups] <- lapply(xs[dups], `[[`, 1)
    ys[dups] <- lapply(ys[dups], sum, na.rm = TRUE)
    xs <- unlist(xs, recursive = FALSE, use.names = FALSE)
    ys <- unlist(ys, recursive = FALSE, use.names = FALSE)
  }
  
  list(title = title,
       ms_lev = ms_lev,
       ret_time = ret_time,
       scan_num = scan_num, 
       mobility = mobils, 
       xs = xs, ys = ys)
}


#' Helper of \link{extract_pasef_frame}
#' 
#' @param mdata A list of PASEF frame.
hextract_pasef_frame <- function (mdata)
{
  lapply(mdata, extract_pasef_frame)
}


#' Helper in executing Cpp-compiled timsdataSampleCpp.exe
#' 
#' Requires mono for Linux or Mac OS.
#' 
#' @param raw_file A RAW file name.
#' @param mgf_path A file path RAW MS files.
exeReadPASEF <- function(raw_file, mgf_path, topn_ms2ions = 150L)
{
  sys_path <- system.file("extdata", package = "mzion")
  exe <- file.path(sys_path, name_exe = "timsdataSampleCpp.exe")
  mono <- if (Sys.info()['sysname'] %in% c("Darwin", "Linux")) TRUE else FALSE
  
  spectra <- file.path(mgf_path, raw_file, "spectra.txt")
  precursors <- file.path(mgf_path, raw_file, "ms2precursors.txt")
  # stdout <- tempfile(tmpdir = mgf_path, fileext = ".stdout")
  # stderr <- tempfile(tmpdir = mgf_path, fileext = ".stderr" )
  raw_full <- file.path(mgf_path, raw_file)

  if (!file.exists(analysis <- file.path(mgf_path, raw_file, "analysis.tdf"))) {
    stop("File not found: ", analysis)
  }
  
  if (mono) {
    rvs <- system2(Sys.which("mono"), 
                   args = c(shQuote(exe),
                            shQuote(raw_full),
                            topn_ms2ions
                            ),
                   stdout = stdout,
                   stderr = stderr)
  }
  else{
    rvs <- system2(exe, 
                   args = c(shQuote(raw_full), 
                            topn_ms2ions
                            ), 
                   stdout = stdout,
                   stderr = stderr)
  }
  
  if (rvs) {
    stop("Fail to process ", raw_file, ".")
  }
  
  # stitch together spectra and precursors -> .peaks
  
  list(spectra, precursors)

  # out_name <- paste0(gsub("\\.[^.]*$", "", raw_file), ".peaks")
  # file.rename(spectra, file.path(mgf_path, out_name))
  # file.rename(precursors, file.path(mgf_path, out_name))
  # unlink(c(stdout, stderr))
  # out_name
}




#' License agreement.
#' 
#' Bruker.
#' 
#' @param sys_path The system path of Mzion.
acceptBrukerLicense <- function(sys_path) 
{
  license  <- file.path(sys_path, name_exe = "THIRD-PARTY-LICENSE-README.txt")
  libdir <- create_dir(tools::R_user_dir("mzion", which = 'cache'))
  
  if (file.exists(file.path(libdir, "rspn_bruker.txt")))
    return(TRUE)
  
  msg <- "By selecting YES, you are accepting the Bruker's License agreement."
  file.show(license)
  response <- readline(prompt = sprintf("Accept '%s'? [y/n]: ", license))
  
  if (tolower(response) != "y")
    stop("Need to accept the license agreement for processing Bruker's PASEF.")
  
  rspnFile <- file.path(libdir, "rspn_bruker.txt")
  fileConn <- file(rspnFile)
  writeLines(paste(msg, paste0(date()), "Response = TRUE", sep = "\n"), fileConn)
  close(fileConn)
  
  NULL
}

