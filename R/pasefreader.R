#' Reads Bruker's PASEF files.
#'
#' @param mgf_path A file path to .d folders.
#' @param filelist A list of .d folders.
#' @param bypass_rawexe Logical; to bypass \link{exeReadPASEF} processing of
#'   RAW .d files or not.
#' @inheritParams matchMS
readPASEF <- function (mgf_path = NULL, filelist = NULL, topn_ms2ions = 150L, 
                       bypass_rawexe = FALSE) 
{
  sys_path <- system.file("extdata", package = "mzion")
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))
  logs     <- file.path(temp_dir, "log.txt")
  acceptBrukerLicense(sys_path)
  
  qs::qsave(list(data_format = "Bruker-RAW", mgf_format = "Mzion"), 
            file.path(mgf_path, "info_format.rds"), preset = "fast")

  len   <- length(filelist)
  n_pcs <- detect_cores(64L)
  ram_units <- max(floor(find_free_mem() / 1024 / 10), 1L)
  
  if (ram_units <= 1L) {
    n_cores <- 1L
  }
  else {
    n_cores <- find_min_ncores(len, min(len, detect_cores(16L), ram_units))
  }
  
  if (!bypass_rawexe) {
    if (n_cores <= 1L) {
      lapply(filelist, exeReadPASEF, mgf_path)
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores), 
                                  outfile = logs)
      parallel::clusterApply(cl, filelist, exeReadPASEF, mgf_path)
      parallel::stopCluster(cl)
    }
  }
  
  # n_cores2 <- min(len, detect_cores(8L), ram_units)
  # n_cores2 <- find_min_ncores(len, n_cores2)
  n_cores2 <- 1L
  # 90 min per file at * 4; was * 2
  n_para   <- min(max(floor(ram_units * 4 / n_cores2), 1L), n_pcs - 1L)

  message("Compiling RAW PASEF peak lists at: ", Sys.time())
  
  if (n_cores2 <= 1L) {
    filenames <- lapply(filelist, hproc_pasefs, mgf_path = mgf_path, 
                        temp_dir = temp_dir, n_para = n_para)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores2), outfile = logs)
    filenames <- parallel::clusterApply(cl, filelist, hproc_pasefs, 
      mgf_path = mgf_path, temp_dir = temp_dir, n_para = n_para)
    parallel::stopCluster(cl)
  }
  
  filenames
}


#' Processes single .d file
#'
#' By each .d file.
#' 
#' @param raw_file The name of a .d folder.
#' @param mgf_path A file path to .d folders.
#' @param temp_dir A file path for temporary files.
#' @param n_para The number of maximum allowed parallel processes.
#' @param debug Debugging mode or not.
hproc_pasefs <- function (raw_file, mgf_path, temp_dir, n_para = 1L, 
                          debug = FALSE)
{
  options(digits = 9L)
  out_name <- paste0(raw_file, ".rds")
  logs <- file.path(temp_dir, paste0("log_", raw_file, ".txt"))
  path <- file.path(mgf_path, raw_file)
  files <- list.files(path, pattern = "^spectra_\\d+\\.txt$")
  n_files <- length(files)

  if (!n_files) {
    stop("No spectra data found under ", path, ".")
  }
  
  files  <- files[order(file.size(file.path(path, files)), decreasing = TRUE)]
  
  ## (1) combine information of MS2 isolation windows
  comb_pasef_ms1ms2iso(path)
  
  ## (2) process spectra_[...].txt
  if (n_para <= 1L) {
    if (n_files == 1L) {
      out <- hextract_pasef(files, path = path)
    }
    else {
      out <- lapply(files, hextract_pasef, path = path) |>
        dplyr::bind_rows()
    }
  }
  else {
    cl  <- parallel::makeCluster(getOption("cl.cores", n_para), outfile = logs)
    out <- parallel::clusterApplyLB(cl, files, hextract_pasef, path = path)
    parallel::stopCluster(cl)
    out <- dplyr::bind_rows(out)
  }

  ## (3) outputs
  #  stagger between MS1 and MS2 frames
  out <- out[with(out, order(ms1_fr, slice_start, ms_level)), ]
  out$scan_num <- seq_along(out$scan_num)
  
  # PASEF full MS1
  rows <- out$slice_start == 0L
  qs::qsave(out[rows, ], file.path(temp_dir, paste0("pasefms1_", out_name)), 
            preset = "fast")
  out <- out[!rows, ]

  lens <- lengths(out$msx_moverzs)
  oks  <- out[["ms_level"]] == 1L | (lens >= 25L & lens <= 2000L)
  out  <- out[oks, ]
  out <- as.list(out)
  out$raw_file <- raw_file # scalar
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  
  message("Completed PASEF peak list processing: ", out_name)

  attr(out_name, "is_dia") <- FALSE
  attr(out_name, "mzml_type") <- "raw" # a token for Mzion de-isotoping

  invisible(out_name)
}


#' Splits single \code{spec_[...].txt} by frames
#' 
#' @param file A file name of \code{spec_[...].txt}.
#' @param path A path to a \code{.d} folder.
split_paseflines_by_frames <- function (file, path)
{
  if (!grepl("^spectra_\\d+\\.txt$", file)) {
    stop("The input file needs to be \"spectra_[0-9].txt\" etc.")
  }
  
  lines   <- readr::read_lines(file.path(path, file))
  lines   <- lines[-length(lines)] # two empty lines at the end
  empties <- .Internal(which(lines == ""))
  lines   <- split(lines, cut(seq_along(lines), c(0, empties)))
  lines   <- unname(lines)
}


#' Combines \code{precursors.txt} and \code{ms2precursors.txt}
#' 
#' @param path A file path to a \code{.d} folder.
comb_pasef_ms1ms2iso <- function (path)
{
  parents <- 
    readr::read_tsv(file.path(path, "precursors.txt"), 
                    col_types = cols(
                      Id = col_integer(), 
                      MonoisotpoicMz = col_number(),
                      Charge = col_integer(),
                      Intensity = col_integer(),
                      Parent = col_integer()
                    ), show_col_types = FALSE) |>
    dplyr::rename(Precursor = Id)
  
  iso_info <- 
    readr::read_tsv(file.path(path, "ms2precursors.txt"), 
                    col_types = cols(
                      Frame = col_integer(), 
                      ScanNumBegin = col_integer(),
                      ScanNumEnd = col_integer(), 
                      IsolationMz = col_number(), 
                      IsolationWidth = col_number(),
                      Precursor = col_integer()
                    ), show_col_types = FALSE) |>
    dplyr::left_join(parents, by = "Precursor") |>
    dplyr::rename(MS2Frame = Frame, MS1Frame = Parent) |>
    reloc_col_after("MS1Frame", "MS2Frame") # |> add_pasef_precursors()
  
  readr::write_tsv(iso_info, file.path(path, "ms1ms2iso.txt"))
}


#' Adds PASEF precursor IDs to \code{iso_info}
#'
#' Not used; only for cross-validation with the Bruker's grouping.
#'
#' @param iso_info A data frame contains MS2 isolation information.
#' @param keys The column keys in \code{iso_info} for combinations of unique
#'   identifiers.
add_pasef_precursors <- function (
    iso_info, keys = c("MS1Frame", "ScanNumBegin", "IsolationMz"))
{
  if (TRUE) {
    iso_info <- iso_info |>
      tidyr::unite(uid, all_of(keys), sep = ".", remove = FALSE) |>
      dplyr::group_by(uid) |>
      dplyr::mutate(PrecursorID = dplyr::cur_group_id()) |>
      dplyr::ungroup()
  }
  else {
    iso_info <- iso_info |>
      tidyr::unite(uid, all_of(c("MS1Frame", "ScanNumBegin")), sep = ".", 
                   remove = FALSE)
    iso_info  <- split(iso_info, iso_info$IsolationMz)
    
    for (i in seq_along(iso_info)) {
      iso_info[[i]]$PrecursorID <- as.numeric(i)
    }
    
    lens <- iso_info |>
      lapply(function (x) length(unique(x$uid))) |> 
      unlist(recursive = FALSE, use.names = FALSE)
    
    if (any(oks2 <- lens > 1L)) {
      df2s <- iso_info[oks2] |> 
        lapply(function (x) {
          # may check minmax(MS1Frame) distance
          subids <- match(x$uid, names(split(x$uid, x$uid)))
          
          # add padding '0' to '1' -> '01'; and 10 remains the same
          # better with seq_len(len) / 10^nchar(as.character(len))
          if (!all(ten_less <- subids < 10)) {
            subids[ten_less] <- paste0("0", subids[ten_less])
          }
          
          x$PrecursorID <- as.numeric(paste0(x$PrecursorID, ".", subids))
          x
        })
      
      iso_info <- dplyr::bind_rows(
        dplyr::bind_rows(iso_info[!oks2]), 
        dplyr::bind_rows(df2s))
    }
    else {
      iso_info <- dplyr::bind_rows(iso_info)
    }
    
    iso_info <- iso_info |>
      dplyr::select(-c("uid")) |>
      dplyr::arrange(MS2Frame, ScanNumBegin)
    
    pids <- iso_info$PrecursorID
    iso_info$PrecursorID <- match(pids, names(split(pids, pids)))
  }

  iso_info
}


#' Separates PASEF MS2Info by a range of frames
#'
#' @param iso_info A data frame containing ms2info.
#' @param mdata Line data of MS frames (containing both MS1 and MS2).
#' @param offset_fr The offset from an empty line to the line of a FRAME number.
bracket_pasef_ms2info <- function (iso_info, mdata, offset_fr = 3L)
{
  len <- length(mdata) # number of frames
  d1  <- mdata[[1]]    # the first FRAME in the current chunk of frames
  dn  <- mdata[[len]]  # the last  FRAME in the current chunk of frames
  
  sta <- as.integer(d1[length(d1) - offset_fr])
  end <- as.integer(dn[length(dn) - offset_fr])

  # always: values MS2Frame > MS1Frame
  iso_info[with(iso_info, MS2Frame >= sta & MS2Frame <= end), ]
}


#' Extracts data of PASEF frames
#'
#' (1) Pairs MS2 data and iso_info by frame IDs; (2) Collapses data by Precursor
#' groups.
#' 
#' @param file A file name of \code{spectra_[...].txt}.
#' @param path The full path to a .d folder.
#' @param keys The key names of output lists.
#' @param step A step size for mass binning.
hextract_pasef <- function (
    file = "spectra_1.txt", path = NULL, 
    keys = c("msx_moverzs", "msx_ints", "msx_ns", "ms1_fr", 
             "ms1_moverzs", "ms1_ints", "ms1_charges", "scan_title", "ms_level", 
             "ret_time", "scan_num", "orig_scan", "slice_start", "slice_end", 
             "iso_ctr", "iso_lwr", "iso_upr", "mobility"), 
    step = 1.6e-5)
{
  ## (1.a) split line data by frames
  mdata <- split_paseflines_by_frames(file, path)
  
  # Y
  # X
  # SLICE
  # NPEAKS
  # MOBILITY
  # ...
  # FRAME
  # TIME
  # MSORDER
  # \n
  
  ## (1.b) bracket `iso_info` to the current chunk of MS data
  if (!file.exists(iso_file <- file.path(path, "ms1ms2iso.txt"))) {
    stop("File not found ", iso_file, ".")
  }
  
  iso_info <- readr::read_tsv(
    iso_file, 
    col_types = cols(
      MS2Frame = col_integer(),
      MS1Frame = col_integer(), 
      ScanNumBegin = col_integer(), 
      ScanNumEnd = col_integer(), 
      Precursor = col_integer(), 
      Charge = col_integer(), 
      Intensity = col_integer())) |>
    bracket_pasef_ms2info(mdata = mdata, offset_fr = 3L)

  lens    <- lengths(mdata) # numbers of lines in each frame
  ms_levs <- mapply(function (x, n) x[n - 1L], mdata, lens, 
                    SIMPLIFY = TRUE, USE.NAMES = FALSE) |> as.integer()
  ms_frs  <- mapply(function (x, n) x[n - 3L], mdata, lens, 
                    SIMPLIFY = TRUE, USE.NAMES = FALSE) |> as.integer()

  oks1 <- ms_levs == 1L
  oks2 <- .Internal(which(!oks1))
  oks1 <- .Internal(which(oks1))
  
  ## (2) collapse scans into slices by ScanNumBegin:ScanNumEnd (MS1 and MS2)
  #  MS1 slices may not have corresponding MS2 data
  #  MS2 at the same ScanNumBegin:ScanNumEnd and MS1 frame can correspond to 
  #   different species (m/z)

  if (length(oks2)) {
    out2 <- extract_pasefms2(
      mdata = mdata[oks2], ms_frs = ms_frs[oks2], iso_info = iso_info, 
      ms_col = "MS2Frame", keys = keys, title = path, step = step)
    empty2 <- FALSE
  }
  else {
    message("No MS2 spectra for ", file.path(path, file))
    out2 <- vector("list", length(keys))
    names(out2) <- keys
    empty2 <- TRUE
  }
  message("Completed PASEF MS2 extraction: ", file)
  
  if (length(oks1)) {
    iso_info_1 <- 
      unique(iso_info[, c("MS1Frame", "ScanNumBegin", "ScanNumEnd")]) |>
      dplyr::arrange(MS1Frame, ScanNumBegin)
    
    temp_1 <- extract_pasefms1(
      mdata = mdata[oks1], ms_frs = ms_frs[oks1], iso_info = iso_info_1, 
      ms_col = "MS1Frame", keys = keys, title = path, step = step)
    out1 <- temp_1$ms1_slices
    ms1_full <- temp_1$ms1_full
    empty1 <- FALSE
  }
  else {
    message("No MS1 spectra for ", file.path(path, file))
    out1 <- ms1_full <- vector("list", length(keys))
    names(ms1_full) <- names(out1) <- keys
    empty1 <- TRUE
  }
  message("Completed PASEF MS1 extraction: ", file)

  if (identical(names(out1), names(out2))) {
    out <- mapply(`c`, ms1_full, out1, out2, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  }
  else {
    stop("Developer: uneven number of columns between MS1 and MS2 data.")
  }
  
  if (empty1 && empty2) {
    return(NULL)
  }
  
  # when out2 is empty; out is out1 and orig_scan is integer (may already fixed)
  if (!is.character(out$orig_scan)) {
    out$orig_scan <- as.character(out$orig_scan)
  }
  
  df <- tibble::tibble(
    msx_moverzs = out$msx_moverzs, 
    msx_ints = out$msx_ints, 
    msx_ns = out$msx_ns, 
    ms1_fr = out$ms1_fr, 
    ms1_moverzs = out$ms1_moverzs, 
    ms1_ints = out$ms1_ints, 
    ms1_charges = out$ms1_charges, 
    scan_title = out$scan_title, 
    ms_level = out$ms_level, 
    ret_time = out$ret_time, 
    scan_num = out$scan_num, 
    orig_scan = out$orig_scan, 
    slice_start = out$slice_start,
    slice_end = out$slice_end,
    iso_ctr = out$iso_ctr, 
    iso_lwr = out$iso_lwr, 
    iso_upr = out$iso_upr, 
    mobility = out$mobility, )
  
  if (ncol(df) != length(out)) {
    stop("Developer: checks for column drops.")
  }
  
  message("Done ", file.path(path, file))
  
  df
}


#' Extracts PASEF MS1 data.
#'
#' @param mdata A list of \emph{MS1} PASEF frames. Each list entry corresponding
#'   to one MS1 frame. Note that each frame contains multiple mobility slices.
#' @param lens The number of lines in each MS2 frame.
#' @param iso_info The information of MS1 slices.
#' @param ms_col The criteria column in \code{iso_info} for subsetting.
#' @param keys The key names of output lists.
#' @param title The title (.d file  name).
#' @param step A step size for mass binning.
extract_pasefms1 <- function (mdata, ms_frs, iso_info, ms_col = "MS1Frame", 
                              keys, title = "", step = 1.6e-5)
{
  ## (1) mutually subset `mdata` and `iso_info` by MS1 frames
  res <- subset_pasefms(
    mdata = mdata, ms_frs = ms_frs, iso_info = iso_info, ms_col = ms_col)
  mdata <- res[["mdata"]]
  ms_frs <- res[["ms_frs"]]
  iso_info <- res[["iso_info"]]
  rm(list = "res")
  
  ## (2) convert line data to lists for each frame
  ans <- lapply(mdata,  extract_pasef_frame)
  rm(list = "mdata")
  
  # assume no empty frames (e.g. at no additional intensity filtration)
  if (FALSE) {
    oks <- lengths(ans) > 0L
    if (!all(oks)) {
      oks <- which(oks)
      ans <- ans[oks]
      iso_info <- iso_info[oks]
      ms_frs <- ms_frs[oks]
    }
    rm(list = c("oks"))
  }

  if (FALSE) {
    scans <- sapply(ans, `[[`, "scan_num")
    # i <- which(scans == 9915) # i <- 114
    i <- which(scans == 18203L) # i <- 82; i <- 1L
    ai <- ans[[i]]
    df <- data.frame(x = ai$msx_moverzs, y = ai$msx_ints)
    dfx <- df |> dplyr::filter(x >= 807.5, x <= 808.5)
    ggplot2::ggplot() + 
      ggplot2::geom_segment(dfx, mapping = aes(x = x, y = y, xend = x, yend = 0), 
                            color = "gray", linewidth = .1) # + 
    # scale_y_continuous(limits = c(0, 1000))
    tempdata1 <- centroid_pasefms1(
      xs = ai$msx_moverzs, ys = ai$msx_ints, reso = 60000, maxn = 2000L, 
      grad = .667, fct_reso = 1.5, ms_lev = 1L, tol = .10)
    df <- data.frame(x = tempdata1$x, y = tempdata1$y)
    dfx2 <- df |> dplyr::filter(x >= 807.5, x <= 808.5)
  }
  
  ## (3) collapse scans -> slices frame-wisely by ScanNumBegin:ScanNumEnd ranges
  lena <- length(ans)
  if (lena != length(iso_info)) {
    stop("Developer: mismatches in the lengths of MS data and metadata.")
  }

  # (3.1) collapse by full MS1
  X <- Y <- vector("list", lena)
  ret_times <- vector("numeric", lena)
  
  for (i in seq_along(ans)) {
    ansi <- ans[[i]]
    ret_times[[i]] <- ansi$ret_time
    xys  <- collapse_pasef_xys(ansi$msx_moverzs, ansi$msx_ints)

    tempd <- centroid_pasefms1(
      xs = xys[["x"]], ys = xys[["y"]], reso = 60000L, maxn = 2000L, 
      grad = .667, fct_reso = 1.5, ms_lev = 1L, tol = .10)
    X[[i]] <- tempd[["x"]]
    Y[[i]] <- tempd[["y"]]
  }
  rm(list = c("ansi", "xys", "tempd"))
  
  lenx <- length(X)
  na_intsx <- rep_len(NA_integer_, lenx)
  na_realsx <- rep_len(NA_real_, lenx)
  
  ms1_full <- list(
    msx_moverzs = X, 
    msx_ints = Y, 
    slice_start = rep_len(0L, lenx),
    slice_end = rep_len(10000L, lenx), # an arbitrary large integer
    msx_ns = lengths(Y), 
    scan_num = as.numeric(ms_frs), 
    mobility = rep_len(0.0, lenx),
    ms1_fr = ms_frs,
    ret_time = ret_times,
    orig_scan = as.character(ms_frs), 
    scan_title = paste0(title, "; MS1: ", ms_frs), 
    ms_level = rep_len(1L, lenx), 
    ms1_moverzs = na_realsx, 
    ms1_charges = na_intsx, 
    ms1_ints = na_intsx, 
    iso_ctr = na_realsx, 
    iso_lwr = na_realsx, 
    iso_upr = na_realsx)
  
  rm(list = c("lenx", "na_intsx", "na_realsx", "X", "Y"))

  # (3.2) collapse by mobility slices (in accordance with MS2)
  # `msx_moverzs = 0` for slices in metadata but not in MS data
  ans <- mapply(collapse_pasef_ms1scans, ans, iso_info)
  
  # assume no empty entries (e.g. at no additional intensity filtration)
  if (FALSE) {
    oks <- lengths(ans) > 0L
    ans <- ans[oks]
    iso_info <- iso_info[oks]
  }

  ## (4) flatten slices in frames (need to handle that all ans are flat)
  #  nrow(iso_info[[i]]) should == the number of slices in each frame
  #  over-killing to maintain 1-to-1 correspondence, but good tidiness
  nrows <- unlist(lapply(iso_info, nrow))
  
  # check for slices in metadata but not in MS data
  if (length(bads <- which(lengths(ans) < nrows))) {
    iso_info[bads] <- purge_iso_info(ans[bads], iso_info[bads])
    nrows[bads] <- lapply(iso_info[bads], nrow) |>
      unlist(use.names = FALSE, recursive = FALSE)
  }
  
  if (any(nrows > 1L)) {
    ans <- unlist(ans, recursive = FALSE, use.names = FALSE)
  }
  
  if (sum(nrows) != length(ans)) {
    stop("Develooper: mismatched between MS data and metadata.")
  }
  
  # remove slices with `msx_moverzs = 0`
  bads <- lapply(ans, function (x) {
    ms <- x$msx_moverzs
    length(ms) == 1L && ms == 0
  }) |>
    unlist(use.names = FALSE, recursive = FALSE)
  ans <- ans[!bads]

  ## (5) centroid MS1 data
  for (i in seq_along(ans)) {
    ai <- ans[[i]]
    tempdata1 <- centroid_pasefms1(
      xs = ai$msx_moverzs, ys = ai$msx_ints, reso = 60000, maxn = 2000L, 
      grad = .667, fct_reso = 1.5, ms_lev = 1L, tol = .10)
    ans[[i]]$msx_moverzs <- tempdata1[["x"]]
    ans[[i]]$msx_ints <- tempdata1[["y"]]
  }
  rm(list = c("ai", "tempdata1", "bads"))
  
  nms1 <- names(ans[[1]])
  out1 <- vector("list", length(nms1))
  for (i in seq_along(out1)) {
    out1[[i]] <- lapply(ans, `[[`, i)
  }
  names(out1) <- nms1
  
  out1[["ret_time"]] <- 
    unlist(out1[["ret_time"]], recursive = FALSE, use.names = FALSE)
  out1[["scan_num"]] <- 
    unlist(out1[["scan_num"]], recursive = FALSE, use.names = FALSE)
  out1[["ms1_fr"]] <- 
    unlist(out1[["ms1_fr"]], recursive = FALSE, use.names = FALSE)
  out1[["orig_scan"]] <- as.character(out1[["ms1_fr"]])
  out1[["mobility"]] <- 
    unlist(out1[["mobility"]], recursive = FALSE, use.names = FALSE)
  out1[["msx_ns"]] <- lengths(out1[["msx_moverzs"]])
  out1[["slice_rng"]] <- 
    unlist(out1[["slice_rng"]], recursive = FALSE, use.names = FALSE)
  out1[["slice_start"]] <- 
    unlist(out1[["slice_start"]], recursive = FALSE, use.names = FALSE)
  out1[["slice_end"]] <- 
    unlist(out1[["slice_end"]], recursive = FALSE, use.names = FALSE)
  
  len1 <- length(out1[["msx_ns"]])
  na_ints1 <- rep_len(NA_integer_, len1)
  na_reals1 <- rep_len(NA_real_, len1)
  
  out1[["scan_title"]]  <- 
    paste0(title, "; MS1: ", out1[["ms1_fr"]], "[", out1[["slice_rng"]], "]")
  out1[["slice_rng"]] <- NULL
  
  out1[["ms_level"]]    <- rep_len(1L, len1)
  out1[["ms1_moverzs"]] <- na_reals1
  out1[["ms1_charges"]] <- na_ints1
  out1[["ms1_ints"]]    <- na_ints1
  out1[["iso_ctr"]]     <- na_reals1
  out1[["iso_lwr"]]     <- na_reals1
  out1[["iso_upr"]]     <- na_reals1

  if (length(out1) != length(keys) || !all(keys %in% names(out1))) {
    stop("Developer: mismatched columns in PASEF-MS1 data.")
  }
  
  list(ms1_slices = out1[keys], ms1_full = ms1_full[keys])
}


#' Extracts PASEF MS1 data.
#'
#' @param mdata A list of \emph{MS2} PASEF frames. Each list entry
#'   corresponding to one MS2 frame.
#' @param ms_frs The frame numbers corresponding to msdata.
#' @param lens The number of lines in each MS2 frame.
#' @param iso_info The information of MS2 isolation windows etc.
#' @param ms_col The criteria column in \code{iso_info} for subsetting.
#' @param keys The key names of output lists.
#' @param title The title (.d file  name).
#' @param step A step size for mass binning.
extract_pasefms2 <- function (mdata, ms_frs, iso_info, ms_col = "MS2Frame", 
                              keys, title = "", step = 1.6e-5)
{
  ## (0) One `Precursor` ID can only occur in one `MS1Frame` -> 
  #  different `Precursor` IDs for the same precursor at different `MS1Frame`s
  if (FALSE) {
    isx <- unique(iso_info[, c("MS1Frame", "Precursor")])
    isx <- split(isx, isx$Precursor)
    nrs <- sapply(isx, nrow)
    all(nrs == 1L)
  }
  
  ## (1) mutually subset `mdata` and `iso_info` by MS2 frames
  # for each MS2 frame: 
  #  the slice numbers are in an ascending order for `mdata`;
  #  the ScanNumBegin values are also in an ascending order for each `iso_info`
  res <- subset_pasefms(
    mdata = mdata, ms_frs = ms_frs, iso_info = iso_info, ms_col = ms_col)
  mdata    <- res[["mdata"]]
  ms_frs   <- res[["ms_frs"]]
  iso_info <- res[["iso_info"]]
  rm(list = "res")

  # each entry corresponds to one MS2 frame (multiple scans in each frame)
  if (FALSE) {
    rng <- 1149:1158
    ans <- lapply(mdata[rng], extract_pasef_frame)
    ans <- mapply(collapse_pasef_ms2scans, ans, iso_info[rng])
    ans <- unlist(ans, recursive = FALSE, use.names = FALSE)
    precursors <- lapply(ans, `[[`, "precursor") |>
      unlist(recursive = FALSE, use.names = FALSE)
    
    ax2 <- split(ans, precursors)
    z2 <- ax2[[which(names(ax2) == 44101)]] # 44080
    z <- group_ms2pasef_by_precursors(z2)
  }
  
  ## (2) convert line data to lists for each MS2 frame
  ans <- lapply(mdata, extract_pasef_frame)
  rm(list = "mdata")
  
  # assume no empty frames (e.g. at no additional intensity filtration)
  if (FALSE) {
    oks <- lengths(ans) > 0L
    if (!all(oks)) {
      oks <- which(oks)
      ans <- ans[oks]
      iso_info <- iso_info[oks]
      ms_frs <- ms_frs[oks]
    }
    rm(list = c("oks"))
  }
  
  ## (3) collapse scans -> slices frame-wisely by ScanNumBegin:ScanNumEnd ranges
  # ms1_moverzs == 0: undertermined by Bruker metadata
  lena <- length(ans)
  if (lena != length(iso_info)) {
    stop("Developer: mismatches in the lengths of MS2 data and metadata.")
  }
  ans <- mapply(collapse_pasef_ms2scans, ans, iso_info)
  
  # assume no empty entries (e.g. at no additional intensity filtration)
  if (FALSE) {
    oks <- lengths(ans) > 0L
    ans <- ans[oks]
    iso_info <- iso_info[oks]
  }
  
  ## (4) flatten slices in frames (need to handle cases that all ans are flat)
  #  nrow(iso_info[[i]]) indicates the number of slices in each frame
  nrows <- unname(unlist(lapply(iso_info, nrow)))
  
  # with slices in metadata but missing from MS2 data
  if (length(bads <- which(lengths(ans) < nrows))) {
    iso_info[bads] <- purge_iso_info(ans[bads], iso_info[bads])
    nrows[bads] <- lapply(iso_info[bads], nrow) |>
      unlist(use.names = FALSE, recursive = FALSE)
  }
  
  if (any(nrows > 1L)) {
    ans <- unlist(ans, recursive = FALSE, use.names = FALSE)
  }
  
  if (sum(nrows) != length(ans)) {
    stop("Developer: mismatched between MS2 data and metadata.")
  }

  ## (5) group MS2Frames by the same `Precursor` (must have the same MS1Frame)
  precursors <- lapply(ans, `[[`, "precursor") |>
    unlist(recursive = FALSE, use.names = FALSE)
  
  if (FALSE) {
    ax2 <- split(ans, precursors)
    ai2 <- ax2[[which(names(ax2) == "89235")]] # 89235: precursor ID
    z <- group_ms2pasef_by_precursors(ai2, step = 1.6e-5)
    
    a <- centroid_pasefms2(
      xs = z$msx_moverzs, ys = z$msx_ints, reso = 30000, maxn = 500L, 
      grad = .25, fct_reso = 1.25, ms_lev = 2L, keep_ohw = TRUE, 
      tol = .1)
    dfa <- data.frame(x = a$x, y = a$y)
    dfa |>
      dplyr::filter(x >= 690, x <= 699) |>
      ggplot2::ggplot() + 
      ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                            color = "gray", linewidth = .1)
    
    dfx <- data.frame(x = z$msx_moverzs, y = z$msx_ints)
    dfx |>
      dplyr::filter(x >= 988.8, x <= 989.0) |>
      ggplot2::ggplot() + 
      ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                            color = "gray", linewidth = .1)
  }
  
  # Note: for multiple MS2Frames (73002.01, 73003.01, ...) under the same 
  #  `Precursor`, the first scan_num (73002.01) was used after the grouping.
  ans <- lapply(split(ans, precursors), group_ms2pasef_by_precursors, 
                step = step, title = title)
  rm(list = "precursors")

  ## (6) centroid MS2 data
  for (i in seq_along(ans)) {
    ai <- ans[[i]]
    ac <- centroid_pasefms2(
      xs = ai$msx_moverzs, ys = ai$msx_ints, reso = 30000L, maxn = 500L, 
      grad = .25, fct_reso = 1.25, ms_lev = 2L, keep_ohw = TRUE, tol = .1)
    mx <- ac[["x"]]
    my <- as.integer(ac[["y"]])

    if (!is.null(mx)) {
      ans[[i]]$msx_moverzs <- mx
      ans[[i]]$msx_ints <- my
      ans[[i]]$msx_ns <- length(my)
    }
  }
  rm(list = c("ai", "ac", "mx", "my"))
  
  ## (6) outputs
  nms2 <- names(ans[[1]])
  out2 <- vector("list", length(nms2))
  
  for (i in seq_along(out2)) {
    out2[[i]] <- lapply(ans, `[[`, i)
  }
  names(out2) <- nms2
  
  out2[["scan_title"]] <- 
    unlist(out2[["scan_title"]], recursive = FALSE, use.names = FALSE)
  out2[["ms_level"]] <- 
    rep_len(2L, length(out2[["scan_title"]]))
  out2[["msx_ns"]] <- 
    unlist(out2[["msx_ns"]], recursive = FALSE, use.names = FALSE)
  out2[["ms1_fr"]] <- 
    unlist(out2[["ms1_fr"]], recursive = FALSE, use.names = FALSE)
  out2[["slice_start"]] <- 
    unlist(out2[["slice_start"]], recursive = FALSE, use.names = FALSE)
  out2[["slice_end"]] <- 
    unlist(out2[["slice_end"]], recursive = FALSE, use.names = FALSE)
  out2[["ret_time"]] <- 
    unlist(out2[["ret_time"]], recursive = FALSE, use.names = FALSE)
  out2[["scan_num"]] <- 
    unlist(out2[["scan_num"]], recursive = FALSE, use.names = FALSE)
  out2[["orig_scan"]] <- 
    unlist(out2[["orig_scan"]], recursive = FALSE, use.names = FALSE)
  out2[["iso_ctr"]] <- 
    unlist(out2[["iso_ctr"]], recursive = FALSE, use.names = FALSE)
  out2[["iso_lwr"]] <- 
    unlist(out2[["iso_lwr"]], recursive = FALSE, use.names = FALSE)
  out2[["iso_upr"]] <- 
    unlist(out2[["iso_upr"]], recursive = FALSE, use.names = FALSE)
  out2[["mobility"]] <- 
    unlist(out2[["mobility"]], recursive = FALSE, use.names = FALSE)
  out2[["ms1_moverzs"]] <- 
    unlist(out2[["ms1_moverzs"]], recursive = FALSE, use.names = FALSE)
  out2[["ms1_charges"]] <- 
    unlist(out2[["ms1_charges"]], recursive = FALSE, use.names = FALSE)
  out2[["ms1_ints"]] <- 
    unlist(out2[["ms1_ints"]], recursive = FALSE, use.names = FALSE)
  
  if (length(out2) != length(keys) || !all(keys %in% names(out2))) {
    stop("Developer: mismatched columns in PASEF-MS2 data.")
  }
  
  out2 <- out2[keys]
}


#' Removes \code{iso_info} rows not found in MS data.
#' 
#' @param msdata MS data.
#' @param iso_info A subset of metadata paired with \code{msdata}.
purge_iso_info <- function (msdata, iso_info) 
{
  for (i in seq_along(msdata)) {
    data <- msdata[[i]]
    isoi <- iso_info[[i]]
    stai <- isoi$ScanNumBegin
    # endi <- isoi$ScanNumEnd
    
    stax <- lapply(data, `[[`, "slice_start") |>
      unlist(recursive = FALSE, use.names = FALSE)
    
    iso_info[[i]] <- isoi[stai %in% stax, ]
  }
  
  iso_info
}


#' Subset PASEF MS data by frame numbers
#'
#' @param mdata A list of PASEF frames. Each list entry corresponding to one MS
#'   frame.
#' @param ms_frs MS frame numbers corresponding to \code{mdata}.
#' @param iso_info The metadata os slice starts etc.
#' @param ms_col The criteria column for subsetting.
#' @seealso \link{bracket_pasef_ms2info} that brackets \code{iso_info} by the
#'   first and the last frames of \code{mdata}.
subset_pasefms <- function (mdata, ms_frs, iso_info, ms_col = "MS1Frame")
{
  if (!ms_col %in% c("MS1Frame", "MS2Frame")) {
    stop("Key column needs to be either `MS1Frame` or `MS2Frame`.")
  }
  
  if (!ms_col %in% names(iso_info)) {
    stop("Key column not found: ", ms_col)
  }
  
  iso_info <- iso_info[iso_info[[ms_col]] %in% ms_frs, ]
  iso_info <- split(iso_info, iso_info[[ms_col]])
  
  mts <- match(ms_frs, as.integer(names(iso_info)))
  oks_ms <- which(!is.na(mts))
  mdata  <- mdata[oks_ms]
  ms_frs <- ms_frs[oks_ms]
  
  ## already satisfied
  # iso_info <- iso_info[mts[oks_ms]]
  
  list(mdata = mdata, ms_frs = ms_frs, iso_info = iso_info)
}


#' Sum peak area
#'
#' @param ys A sub vector of intensity values around a peak.
#' @param imax The index of the peak position in \code{ys}.
#' @param grad A threshold of gradient between two adjacent peaks (more abundant
#'   over less abundant).
sum_pasef_ms1 <- function (ys, imax, grad = .667)
{
  len <- length(ys)

  if (len == 1L) {
    return(ys)
  }
  
  # not quite for MS2, but perhaps not matter much if an MS2 only has two peaks
  if (len == 2L) {
    return(sum(ys))
  }
  
  yval <- ys[[imax]]

  if (imax < len) {
    pr1  <- imax + 1L
    yr1  <- ys[pr1]
    yval <- yval + yr1
    
    if (pr1 < len) {
      for (j in (pr1 + 1L):len) {
        yr2 <- ys[[j]]
        
        if (yr1 / yr2 >= grad) { # || yr2 <= 100
          yval <- yval + yr2
          yr1  <- yr2
        }
        else {
          d    <- yr1 * grad
          yval <- yval + d
          yr1  <- yr2 - d
        }
      }
    }
  }

  if (imax > 1L) {
    pl1  <- imax - 1L
    yl1  <- ys[[pl1]]
    yval <- yval + yl1
    
    if (pl1 > 1L) {
      for (j in (pl1 - 1L):1) {
        yl2 <- ys[[j]]
        
        if (yl1 / yl2 >= grad) { # || yl2 <= 100
          yval <- yval + yl2
          yl1  <- yl2
        }
        else {
          d    <- yl1 * grad
          yval <- yval + d
          yl1  <- yl2 - d
        }
      }
    }
  }

  as.integer(yval)
}


#' Specialty summing of monotonically decreased or increased intensities.
#' 
#' @param xs A vector of moverz values
#' @param ys A vector of intensity values.
#' @param reso2 Adjusted resolution.
#' @param down Logical; are the \code{ys} monotonically decreased or not.
sum_mono_yints <- function (xs, ys, reso2 = 24000, down = TRUE)
{
  len <- length(ys)
  
  if (!len) {
    return(NULL)
  }
  
  xout <- yout <- NULL
  
  while (len) {
    if (len == 1L) {
      xout <- c(xout, xs)
      yout <- c(yout, ys)
      break
    }
    
    if (down) {
      x   <- xs[[1]]
      oks <- .Internal(which(xs <= x + x / reso2))
    }
    else {
      x   <- xs[[len]]
      oks <- .Internal(which(xs >= x - x / reso2))
    }
    
    xout <- c(xout, x)
    yout <- c(yout, sum(ys[oks]))
    xs   <- xs[-oks]
    ys   <- ys[-oks]
    len  <- length(ys)
  }
  
  if (down) {
    list(x = xout, y = yout)
  }
  else {
    rng <- length(yout):1
    list(x = xout[rng], y = yout[rng])
  }
}


#' Sum PASEF MS1 peak area
#' 
#' @param xs A vector of ascending m-over-z values.
#' @param ys A vector of intensity values.
#' @param reso The resolution of a peak.
#' @param maxn The maximum number of peaks.
#' @param ymin The minimum Y values for considering in peak centroiding.
#' @param tol The tolerance of Y for defining a peak profile.
#' @param grad A threshold of gradient between two adjacent peaks (more abundant
#'   over less abundant).
#' @param fct_reso A scaling factor for instrument resolution.
#' @param ms_lev The level of MS.
#' @importFrom fastmatch %fin%
#' @examples
#' # example code
#' mzion:::centroid_pasefms1(c(500), c(10))
#' mzion:::centroid_pasefms1(c(500, 500.01), c(1, 10))
#' mzion:::centroid_pasefms1(c(500, 500.01, 500.2), c(1, 10, 20))
#' mzion:::centroid_pasefms1(c(500, 500.01, 500.2), c(10, 5, 2))
#' 
#' # ignore the trailing half peak
#' mzion:::centroid_pasefms1(500 + .01 * 0:3, c(1, 10, 2, 5))
#' 
#' mzion:::centroid_pasefms1(500 + .01 * 0:4, c(1, 10, 2, 5, 3))
#' 
#' # trailing max at the last position
#' mzion:::centroid_pasefms1(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,33,23,22,41))
#' 
#' mzion:::centroid_pasefms1(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,15,23,30,41))
centroid_pasefms1 <- function (xs, ys, reso = 60000, maxn = 2000L, ymin = 50L, 
                               grad = .667, fct_reso = 1.5, ms_lev = 1L, 
                               tol = .10)
{
  len <- length(ys)
  
  if (!len) {
    return(NULL)
  }
  
  if (len == 1L) {
    return(list(x = xs, y = ys))
  }
  
  reso2 <- reso / fct_reso
  
  # rising edge
  ds1 <- diff(ys) > 0
  
  # only two peaks
  if (length(ds1) == 1L) {
    return(list(x = if (ds1) xs[[2]] else xs[[1]], y = sum(ys)))
  }
  
  # all falling
  if (!any(ds1)) {
    return(list(x = xs[[1]], y = sum(ys)))
  }
  
  # falling edge
  ps <- .Internal(which(diff(ds1) == -1L)) + 1L
  ps <- ps[ys[ps] >= ymin] # optional
  ns <- length(ps)
  
  # no falling
  if (!ns) {
    return(list(x = xs[[len]], y = sum(ys)))
  }
  
  if (ns == 1L) {
    return(list(x = xs[[ps]], y = ys[[ps - 1L]] + ys[[ps]] + ys[[ps + 1L]]))
  }
  
  ns    <- min(ns, maxn)
  xvals <- vector("numeric", ns)
  yvals <- vector("integer", ns)
  ord   <- 
    .Internal(radixsort(na.last = TRUE, decreasing = TRUE, FALSE, TRUE, ys[ps]))
  
  ct <- 0L
  for (i in 1:ns) {
    if (ct == maxn)
      break
    
    oi <- ord[[i]]
    p  <- ps[oi]
    
    if (is.na(p))
      next
    
    x <- xs[[p]]
    w <- x / reso2
    idxes <- .Internal(which(xs >= x - w & xs <= x + w))
    ysubs <- ys[idxes]
    imax  <- .Internal(which(idxes == p)) # relative to ysubs
    
    yvals[[i]] <- sum_pasef_ms1(ys = ysubs, imax = imax, grad = grad)
    xvals[[i]] <- xs[[p]]
    
    # ps[ps %fin% idxes] <- NA_integer_ # slower?
    rng  <- max(1L, oi - 5L):min(len, oi + 5L)
    nbrs <- .Internal(which(ps[rng] %in% idxes))
    ps[rng[nbrs]] <- NA_integer_
    
    ct <- ct + 1L
  }
  
  oks <- yvals > 0L
  xvals <- xvals[oks]
  yvals <- yvals[oks]
  
  if (length(yvals)) {
    ord <- 
      .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, xvals))
    xvals <- xvals[ord]
    yvals <- yvals[ord]
  }
  
  list (x = xvals, y = yvals)
}


#' Sum PASEF MS2 peak area
#' 
#' @param xs A vector of ascending m-over-z values.
#' @param ys A vector of intensity values.
#' @param reso The resolution of a peak.
#' @param maxn The maximum number of peaks.
#' @param ymin The minimum Y values for considering in peak centroiding.
#' @param tol The tolerance of Y for defining a peak profile.
#' @param grad A threshold of gradient between two adjacent peaks (more abundant
#'   over less abundant).
#' @param fct_reso A scaling factor for instrument resolution.
#' @param ms_lev The level of MS.
#' @param keep_ohw Logical; Keep one-hit-wonders or not.
#' @importFrom fastmatch %fin%
centroid_pasefms2 <- function (xs, ys, reso = 30000, maxn = 500L, ymin = 10L, 
                              grad = .25, fct_reso = 1.25, ms_lev = 2L, 
                              keep_ohw = TRUE, tol = .10)
{
  len <- length(ys)
  
  if (!len) {
    return(NULL)
  }
  
  if (len == 1L) {
    return(list(x = xs, y = ys))
  }
  
  ## (1) find peaks (local maximal)
  reso2 <- reso / fct_reso
  
  # rising edges
  ds1 <- diff(ys) > 0
  
  # only two peaks
  if (length(ds1) == 1L) {
    if (ms_lev == 1L) {
      return(list(x = if (ds1) xs[[2]] else xs[[1]], y = sum(ys)))
    }
    else {
      if (xs[[2]] <= xs[[1]] * (1 + reso2) /reso2) {
        return(list(x = if (ys[[2]] > ys[[1]]) xs[[2]] else xs[[1]], y = sum(ys)))
      }
      else {
        return(list(x = xs, y = ys))
      }
    }
  }
  
  # all falling
  if (!any(ds1)) {
    if (ms_lev == 1L) {
      return(list(x = xs[[1]], y = sum(ys)))
    }
    else {
      return(sum_mono_yints(xs = xs, ys = ys, reso2 = reso2, down = TRUE))
    }
  }
  
  # falling edge
  ps <- .Internal(which(diff(ds1) == -1L)) + 1L
  ps <- ps[ys[ps] >= ymin] # optional
  ns <- length(ps)
  
  # no falling
  if (!ns) {
    if (ms_lev == 1L) {
      return(list(x = xs[[len]], y = sum(ys)))
    }
    else {
      return(sum_mono_yints(xs = xs, ys = ys, reso2 = reso2, down = FALSE))
    }
  }
  
  if (ns == 1L) {
    if (ms_lev == 1L) {
      return(list(x = xs[[ps]], y = ys[[ps - 1L]] + ys[[ps]] + ys[[ps + 1L]]))
    }
    else {
      # may later collapse signals for adjacent peaks...
      # return(list(x = xs, y = ys))
    }
  }
  
  ns    <- min(ns, maxn)
  xvals <- vector("numeric", ns)
  yvals <- vector("integer", ns)
  ord   <- 
    .Internal(radixsort(na.last = TRUE, decreasing = TRUE, FALSE, TRUE, ys[ps]))
  i1hs  <- NULL # tracking one-hit-wonders
  
  ## (2) integrate peak profiles
  ct <- 0L
  for (i in 1:ns) {
    if (ct == maxn)
      break
    
    oi <- ord[[i]]
    p  <- ps[oi] # relative to ys
    
    if (is.na(p))
      next
    
    x <- xs[[p]]
    w <- x / reso2
    idxes <- .Internal(which(xs >= x - w & xs <= x + w))
    
    if (keep_ohw) {
      i1hs <- c(i1hs, idxes)
    }
    
    ysubs <- ys[idxes]
    imax  <- .Internal(which(idxes == p)) # relative to ysubs
    
    yvals[[i]] <- sum_pasef_ms1(ys = ysubs, imax = imax, grad = grad)
    xvals[[i]] <- xs[[p]]
    
    ps[ps %fin% idxes] <- NA_integer_ # slower?
    # rng  <- max(1L, oi - 5L):min(len, oi + 5L)
    # nbrs <- .Internal(which(ps[rng] %in% idxes))
    # ps[rng[nbrs]] <- NA_integer_
    
    ct <- ct + 1L
  }
  
  if (keep_ohw) {
    ohw <- .Internal(which(!(1:len) %fin% i1hs))
    xvals <- c(xvals, xs[ohw])
    yvals <- c(yvals, ys[ohw])
  }
  
  
  oks <- yvals > 0L
  xvals <- xvals[oks]
  yvals <- yvals[oks]
  
  if (length(yvals)) {
    ord <- 
      .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, xvals))
    xvals <- xvals[ord]
    yvals <- yvals[ord]
  }
  
  list (x = xvals, y = yvals)
}


#' Collapses PASEF MS2 data (belonging to the same \code{MS1Frame}) by
#' \code{Precursor} IDs.
#'
#' The same \code{Precursor} ID must have the same \code{MS1Frame}. The same
#' precursor at different \code{MS1Frame}s has different \code{Precursor} IDs.
#'
#' @param dat MS2 data under the same precursor ID.
#' @param lwr A lower bound as the starting point in mass binning.
#' @param step A step size for mass binning.
#' @param title The title (.d file  name).
group_ms2pasef_by_precursors <- function (dat, lwr = 115L, step = 1.6e-5, 
                                          title = "")
{
  oks <- .Internal(which(lapply(dat, `[[`, "msx_ns") > 0L))
  len <- length(oks)
  
  if (!len) {
    return(NULL)
  }

  if (len < length(dat)) {
    dat <- dat[oks]
  }

  xys <- collapse_mms1ints(
    lapply(dat, `[[`, "msx_moverzs"), lapply(dat, `[[`, "msx_ints"), 
    lwr = lwr, step = step, reord = FALSE, cleanup = FALSE, 
    sum_y = TRUE, add_colnames = FALSE, look_back = TRUE)
  xys <- calc_ms1xys(xys[["x"]], xys[["y"]])

  msx_moverzs <- xys[["x"]]
  msx_ints <- as.integer(xys[["y"]])
  # ns <- xys$n
  msx_ns <- length(msx_ints)
  
  ms1_frs  <- lapply(dat, `[[`, "ms1_fr") # should be all the same
  ret_time <- lapply(dat, `[[`, "ret_time")
  ret_time <- .Internal(unlist(ret_time, recursive = FALSE, use.names = FALSE))
  ret_time <- sum(ret_time)/len
  
  scan_num_all <- lapply(dat, `[[`, "scan_num")
  scan_num_all <- 
    .Internal(unlist(scan_num_all, recursive = FALSE, use.names = FALSE))
  scan_num_all <- paste0(scan_num_all, collapse = ",")
  
  mobility <- lapply(dat, `[[`, "mobility")
  mobility <- .Internal(unlist(mobility, recursive = FALSE, use.names = FALSE))
  mobility <- sum(mobility)/len
  
  slice_stas <- lapply(dat, `[[`, "slice_start")
  slice_ends <- lapply(dat, `[[`, "slice_end")
  slice_stas <- 
    .Internal(unlist(slice_stas, recursive = FALSE, use.names = FALSE))
  slice_ends <- 
    .Internal(unlist(slice_ends, recursive = FALSE, use.names = FALSE))
  
  ###
  if (length(unique(slice_stas)) > 1L) {
    stop("Multiple starts for the same precursor ID.")
  }
  if (length(unique(slice_ends)) > 1L) {
    stop("Multiple ends for the same precursor ID.")
  }
  ###
  
  slice_rng <- lapply(dat, `[[`, "slice_rng")
  slice_rng <- 
    .Internal(unlist(slice_rng, recursive = FALSE, use.names = FALSE))
  slice_rng <- paste0(slice_rng, collapse = ",")

  # uses the fist one as an surrogate for staggering with MS1 scan numbers
  # sapply(dat, `[[`, "scan_num")
  # [1] 73002.01 73003.01 73004.01 ...
  dat_1 <- dat[[1]]
  scan_num <- dat_1[["scan_num"]]
  ms1_fr <- ms1_frs[[1]]
  scan_title <- paste0(title, "; MS1: ", ms1_fr, 
                       "; scans: ", scan_num_all, 
                       "[", slice_rng, "]", "; Cmpd: ", dat_1[["precursor"]])
  iso_ctr <- dat_1[["iso_ctr"]]
  iso_lwr <- dat_1[["iso_lwr"]]
  iso_upr <- dat_1[["iso_upr"]]

  list(
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    scan_title = scan_title,
    ms1_fr = ms1_fr,
    slice_start = min(slice_stas), 
    slice_end = max(slice_ends), 
    ret_time = ret_time,
    scan_num = scan_num, # use the first scan_num under the same Precursor ID 
    orig_scan = scan_num_all, 
    iso_ctr = dat_1[["iso_ctr"]],
    iso_lwr = dat_1[["iso_lwr"]],
    iso_upr = dat_1[["iso_upr"]], 
    mobility = mobility, 
    ms1_moverzs = dat_1[["ms1_moverzs"]], 
    ms1_ints = dat_1[["ms1_ints"]], 
    ms1_charges = dat_1[["ms1_charges"]])
}


#' Collapses slices in an MS1 frame.
#'
#' @param data One frame of MS1 data before combining slices by IM groups.
#' @param iso_info A subset of MS1 isolation information for the current
#'   \code{frame}.
#' @param min_ms2n The minimum number of MS2 in a slices for considerations.
#' @return A vector of lists. Each list corresponding to a slice of MS1 between
#'   ScanNumBegin and ScanNumEnd in a frame.
collapse_pasef_ms1scans <- function (data, iso_info, min_ms2n = 0L)
{
  if (!length(data)) {
    stop("MS1 line data cannot be empty.")
  }
  
  xs <- data$msx_moverzs
  ys <- data$msx_ints
  slices <- data$slices
  mobils <- data$mobility
  ret_times <- data$ret_time
  ms1_fr <- data$scan_num
  stas <- iso_info$ScanNumBegin
  ends <- iso_info$ScanNumEnd
  len  <- length(stas)
  out  <- vector("list", len)

  if (!len) {
    stop("MS metadata cannot be empty.")
  }
  
  for (i in 1:len) {
    stai <- stas[[i]]
    endi <- ends[[i]]
    oksi <- .Internal(which(slices >= stai & slices <= endi))
    
    # slices may be in metadata but not in MS data
    if (length(oksi)) {
      ansi <- collapse_pasef_xys(xs[oksi], ys[oksi])
      ansx <- ansi[["x"]]
      ansy <- ansi[["y"]]
      
      out[[i]] <- list(
        msx_moverzs = ansx, 
        msx_ints = ansy,
        slice_rng = paste0(stai, "-", endi), 
        slice_start = stai, 
        slice_end = endi, 
        msx_ns = length(ansy), 
        scan_num = ms1_fr + i / 100, # simply assume < 100 slices
        mobility = sum(mobils[oksi]) / length(oksi), 
        ms1_fr = ms1_fr,
        ret_time = ret_times)
    }
    else {
      ansi <- collapse_pasef_xys(xs[oksi], ys[oksi])
      ansx <- ansi[["x"]]
      ansy <- ansi[["y"]]
      
      out[[i]] <- list(
        msx_moverzs = 0, # prevent entry drops
        msx_ints = 0L,
        slice_rng = paste0(stai, "-", endi), 
        slice_start = stai, 
        slice_end = endi, 
        msx_ns = length(ansy), 
        scan_num = ms1_fr + i / 100,
        mobility = 0, 
        
        ms1_fr = ms1_fr,
        ret_time = ret_times)
    }
  }
  
  out
}


#' Collapses slices in an MS2 frame.
#'
#' @param data One frame of MS2 data before combining slices by IM groups.
#' @param iso_info A subset of MS2 isolation information at the current MS2
#'   \code{frame}.
#' @param min_ms2n The minimum number of MS2 in a slices for considerations.
#' @return A vector of lists. Each list corresponding to a slice of MS2 between
#'   ScanNumBegin and ScanNumEnd in a frame, with additional fields of
#'   IsolationMz, IsolationWidth, Precursor etc.
collapse_pasef_ms2scans <- function (data, iso_info, min_ms2n = 0L)
{
  if (!length(data)) {
    return(NULL)
  }
  
  ## (1) combine slices by ranges of ScanNumBegin and ScanNumEnd
  breaks <- findInterval(data$slices, iso_info$ScanNumEnd, left.open = TRUE)
  slices <- split(data$slices, breaks)
  len    <- length(slices)
  
  if (!len) {
    return(NULL)
  }
  
  # removes `iso_info` rows without matched scan ranges to data$slices
  if (len < nrow(iso_info)) {
    iso_info <- iso_info[as.integer(names(slices)) + 1L, ]
  }
  
  xys <- mapply(collapse_pasef_xys, 
                split(data$msx_moverzs, breaks), 
                split(data$msx_ints, breaks), 
                SIMPLIFY = FALSE, USE.NAMES = FALSE)
  msx_moverzs <- lapply(xys, `[[`, 1)
  msx_ints    <- lapply(xys, `[[`, 2)
  mobs        <- split(data$mobility, breaks)
  mobs        <- lapply(mobs, function (x) sum(x) / length(x))
  
  ## (2) clean ups
  if (FALSE) {
    msx_ns <- lengths(msx_moverzs)
    
    # if (any(msx_ns) == 0L) {
    #   stop("Not expecting `msx_ns == 0`.")
    # }
    
    oks <- .Internal(which(msx_ns >= min_ms2n))
    len <- length(oks)
    
    if (!len) {
      return(NULL)
    }
    
    msx_ns <- msx_ns[oks]
    iso_ctrs <- iso_info$IsolationMz[oks]
    iso_widths <- iso_info$IsolationWidth[oks]
    precursors <- iso_info$Precursor[oks]
    ms1_moverzs <- iso_info$MonoisotpoicMz[oks]
    ms1_ints <- iso_info$Intensity[oks]
    ms1_charges <- iso_info$Charge[oks]
    ms1_frs <- iso_info$MS1Frame[oks]
    msx_moverzs <- msx_moverzs[oks]
    msx_ints <- msx_ints[oks]
    slices <- slices[oks]
    mobs <- mobs[oks]
  }
  
  msx_ns      <- lengths(msx_moverzs)
  iso_ctrs    <- iso_info$IsolationMz
  iso_widths  <- iso_info$IsolationWidth
  precursors  <- iso_info$Precursor
  ms1_moverzs <- iso_info$MonoisotpoicMz
  ms1_ints    <- iso_info$Intensity
  ms1_charges <- iso_info$Charge
  ms1_frs     <- iso_info$MS1Frame
  
  ## MonoisotpoicMz is almost always -.8 and -0.1 lower than IsolationMz: 
  ##  as it is weighted-mean of isolated features & 13C dist are right-censored
  ## IsolationWidth is 2 at ~ IsolationMz < 720 and 3 at IsolationMz > 800.
  ## The wide IsolationWidth may be intended to include lower m/z values in 
  #   deisotoping, but Mzion has its own extension (+/-2 around IsolationMz).
  
  iso_uprs <- iso_ctrs
  iso_lwrs <- iso_uprs - iso_widths / 2
  iso_ctrs <- (iso_uprs + iso_lwrs) / 2

  ## (3) outputs
  out <- vector("list", len)
  slice_stas <- iso_info$ScanNumBegin
  slice_ends <- iso_info$ScanNumEnd
  
  if (length(slice_stas) != len) {
    stop("Developer: check for mismatched MS2 slices.")
  }

  ret_times  <- rep_len(data$ret_time, len)
  # 73009.01 if >= 10 slices, but 73009.1 if < 10 slices
  #  -> and 73009.1 will later become 73009.10
  # scan_nums  <- data$scan_num + seq_len(len) / 10^nchar(as.character(len))
  scan_nums  <- data$scan_num + seq_len(len) / 100 # simply assume < 100 slices
  
  rngs <- mapply(function (x, y) paste0(x, "-", y), 
                 slice_stas, slice_ends, 
                 SIMPLIFY = TRUE, USE.NAMES = FALSE)
  
  for (i in 1:len) {
    out[[i]] <- list(
      msx_moverzs = msx_moverzs[[i]], 
      msx_ints = msx_ints[[i]],
      slice_rng = rngs[[i]], 
      slice_start = slice_stas[[i]], 
      slice_end = slice_ends[[i]], 
      msx_ns = msx_ns[[i]],
      scan_num = scan_nums[[i]],
      mobility = mobs[[i]], 
      ms1_fr = ms1_frs[[i]],
      ret_time = ret_times[[i]],
      iso_ctr = iso_ctrs[[i]],
      iso_lwr = iso_lwrs[[i]],
      iso_upr = iso_uprs[[i]],
      precursor = precursors[[i]], 
      ms1_moverzs = ms1_moverzs[[i]], 
      ms1_ints = ms1_ints[[i]], 
      ms1_charges = ms1_charges[[i]])
  }
  
  out  
}


#' Extracts MS1 or MS2 data from a PASEF frame
#'
#' @param data Line data in a PASEF frame (with multiple slices).
#' @param coll Logical; collapses MS data or not.Should be FALSE. Only TRUE at
#'   debugging.
extract_pasef_frame <- function (data, coll = FALSE)
{
  len <- length(data) # the number of lines
  
  ## (1) fields common within a frame
  # ms_lev <- as.integer(data[len - 1L])
  ret_time <- as.numeric(data[len - 2L])
  scan_num <- as.integer(data[len - 3L])
  
  ## fields specific to each slice in a frame
  len2 <- len - 4L
  data <- data[1:len2]
  rngs <- 5 * 0:(len2 / 5L - 1L)
  
  slices <- as.integer(data[rngs + 3L])
  npeaks <- as.integer(data[rngs + 4L])
  mobils <- as.numeric(data[rngs + 5L])
  
  ys <- data[rngs + 1L]
  xs <- data[rngs + 2L]
  ys <- stringi::stri_split_fixed(ys, pattern = " ", npeaks, simplify = FALSE)
  xs <- stringi::stri_split_fixed(xs, pattern = " ", npeaks, simplify = FALSE)
  ys <- lapply(ys, as.integer)
  xs <- lapply(xs, as.numeric)

  ## (2) clean up by Y values (to be disabled to avoid potential mismatches)
  # 10L is the default minimum value by PASEF
  for (i in seq_along(xs)) {
    ysi <- ys[[i]]
    oki <- .Internal(which(ysi >= 10L))
    
    if (length(oki) < npeaks[[i]]) {
      ys[[i]] <- ysi[oki] # works at rhs numeric(0)
      xs[[i]] <- xs[[i]][oki]
    }
  }
  
  ## (3) clean up by X values
  oks  <- .Internal(which(lengths(xs) > 0L))
  noks <- length(oks)
  
  if (!noks) {
    return(NULL)
  }
  
  if (noks < length(xs)) {
    xs <- xs[oks]
    ys <- ys[oks]
    slices <- slices[oks]
    mobils <- mobils[oks]
  }
  
  ###
  # Not to collapse MS1: 
  # lose info on the mobility dimension and increase entropy in precursors
  ###
  
  # collapses all MS slices in a Frame; should only be used for debugging
  if (coll) {
    xys <- collapse_pasef_xys(xs, ys)
    xs <- xys$x
    ys <- xys$y
    
    list(ret_time = ret_time,
         scan_num = scan_num,
         mobility = sum(mobils)/length(mobils),
         # vectors below
         msx_moverzs = xs, 
         msx_ints = ys)
  }
  else {
    list(ret_time = ret_time,
         scan_num = scan_num,
         # vectors below
         msx_moverzs = xs, 
         msx_ints = ys, 
         mobility = mobils,
         slices = slices)
  }
}


#' Collapse PASEF X and Y values within an MS1 or MS2 frame.
#' 
#' @param xs Slices of m-over-z values.
#' @param ys Slices of intensity values.
collapse_pasef_xys <- function (xs, ys)
{
  # number of slices
  len <- length(xs)
  
  if (!len) {
    return(NULL)
  }
  
  xs <- .Internal(unlist(xs, recursive = FALSE, use.names = FALSE))
  ys <- .Internal(unlist(ys, recursive = FALSE, use.names = FALSE))

  if (length(xs) > 1L) {
    ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, 
                               TRUE, xs))
    xs <- xs[ord]
    ys <- ys[ord]
  }

  if (anyDuplicated(xs)) {
    ys <- split(ys, xs)
    xs <- split(xs, xs)
    dups <- lengths(xs) > 1L
    xs[dups] <- lapply(xs[dups], `[[`, 1)
    ys[dups] <- lapply(ys[dups], sum, na.rm = TRUE)
    xs <- .Internal(unlist(xs, recursive = FALSE, use.names = FALSE))
    ys <- .Internal(unlist(ys, recursive = FALSE, use.names = FALSE))
  }
  
  list(x = xs, y = ys)
}


#' Helper in executing Cpp-compiled timsdataSampleCpp.exe
#' 
#' Requires mono for Linux or Mac OS.
#' 
#' @param raw_file A RAW file name.
#' @param mgf_path A file path RAW MS files.
exeReadPASEF <- function(raw_file, mgf_path)
{
  analysis <- file.path(mgf_path, raw_file, "analysis.tdf")
  if (!file.exists(analysis)) {
    stop("File not found: ", analysis)
  }
  
  sys_path <- system.file("extdata", package = "mzion")
  exe <- file.path(sys_path, name_exe = "timsdataSampleCpp.exe")
  mono <- if (Sys.info()['sysname'] %in% c("Darwin", "Linux")) TRUE else FALSE
  
  raw_full <- file.path(mgf_path, raw_file)
  stdout   <- tempfile(tmpdir = mgf_path, fileext = ".stdout")
  stderr   <- tempfile(tmpdir = mgf_path, fileext = ".stderr" )

  if (mono) {
    rvs <- system2(Sys.which("mono"), 
                   args = c(shQuote(exe), shQuote(raw_full)),
                   stdout = stdout,
                   stderr = stderr)
  }
  else{
    rvs <- system2(exe, 
                   args = c(shQuote(raw_full)), 
                   stdout = stdout,
                   stderr = stderr)
  }
  
  if (rvs) {
    stop("Fail to process ", raw_file, ".")
  }
  
  message("\tComplete ReadRAW: ", raw_file, " at ", Sys.time())
  unlink(c(stdout, stderr))
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


