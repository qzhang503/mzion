#' Reads Bruker's PASEF files.
#'
#' @param mgf_path A file path RAW MS files.
#' @param filelist A list of RAW MS files.
#' @param bypass_rawexe Logical; to bypass \link{exeReadPASEF} processing of
#'   RAW .d files or not.
#' @inheritParams matchMS
readPASEF <- function (mgf_path = NULL, filelist = NULL, topn_ms2ions = 150L, 
                       bypass_rawexe = FALSE) 
{
  sys_path <- system.file("extdata", package = "mzion")
  acceptBrukerLicense(sys_path)
  temp_dir <- create_dir(file.path(mgf_path, "temp_dir"))
  logs <- file.path(temp_dir, "log.txt")
  
  qs::qsave(list(data_format = "Bruker-RAW", mgf_format = "Mzion"), 
            file.path(mgf_path, "info_format.rds"), preset = "fast")

  len <- length(filelist)
  n_pcs <- detect_cores(64L)
  ram_units <- max(floor(find_free_mem()/1024/20), 1L)

  if (ram_units <= 1L) {
    n_cores <- 1L
  }
  else {
    n_cores <- find_min_ncores(len, min(len, detect_cores(16L), ram_units, 8L))
  }
  
  if (!bypass_rawexe) {
    if (n_cores <= 1L) {
      lapply(filelist, exeReadPASEF, mgf_path)
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores), outfile = logs)
      parallel::clusterApply(cl, filelist, exeReadPASEF, mgf_path)
      parallel::stopCluster(cl)
    }
  }

  # can be increased if change the MS1 intensity cut-off in exeReadPASEF to 100
  n_cores2 <- min(len, detect_cores(2L), ram_units)
  n_cores2 <- find_min_ncores(len, n_cores2)
  n_para <- max(floor(ram_units/n_cores2), 1L)
  # n_para <- floor(16L/n_cores2)
  # n_para <- if (len > 1L) 1L else max(min(floor(n_pcs/n_cores2/2.8), 16L), 1L)
  
  message("Extracting RAW PASEF peak lists at: ", Sys.time())

  if (n_cores2 <= 1L) {
    filenames <- lapply(filelist, proc_pasefs, mgf_path = mgf_path, 
                        temp_dir = temp_dir, n_para = n_para)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores2), outfile = logs)
    filenames <- parallel::clusterApply(cl, filelist, proc_pasefs, 
      mgf_path = mgf_path, temp_dir = temp_dir, n_para = n_para)
    parallel::stopCluster(cl)
  }
  
  filenames
}


#' Helper in processing MGF entries in chunks.
#'
#' By each .tdf file.
#' 
#' @param raw_file The file name of RAW MS data.
#' @param temp_dir A file path for temporary files.
#' @param n_para The number of maximum allowed parallel processes.
#' @param debug Debugging mode or not.
#' @inheritParams matchMS
proc_pasefs <- function (raw_file, mgf_path, temp_dir, n_para = 1L, 
                         debug = FALSE)
{
  options(digits = 9, warn = 1)
  logs <- file.path(temp_dir, "log.txt")
  
  # MS1 and MS2 spectra
  lines <- readr::read_lines(file.path(mgf_path, raw_file, "spectra.txt"), 
                             skip = 3L, num_threads = n_para)
  idx_empties <- .Internal(which(stringi::stri_cmp_eq(lines, "TITLE"))) + 2L
  lines <- split(lines, cut(seq_along(lines), c(0, idx_empties)))
  lines <- unname(lines) # 'cut' ranges being the names
  
  # Information of MS2 isolation windows
  parents <- 
    readr::read_tsv(file.path(mgf_path, raw_file, "precursors.txt"), 
                    col_types = cols(
                      Id = col_integer(), 
                      MonoisotpoicMz = col_number(),
                      Charge = col_integer(),
                      Intensity = col_integer(),
                      Parent = col_integer()
                    ), show_col_types = FALSE) |>
    dplyr::rename(Precursor = Id)
  
  iso_info <- 
    readr::read_tsv(file.path(mgf_path, raw_file, "ms2precursors.txt"), 
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
  
  if (n_para <= 1L) {
    out <- hextract_pasef(lines, iso_info)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_para), outfile = logs)
    lines <- chunksplit(lines, n_para)
    iso_info <- lapply(lines, sep_pasef_ms2info, iso_info)

    ans <- parallel::clusterMap(
      cl, hextract_pasef, 
      lines, 
      iso_info, 
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    parallel::stopCluster(cl)

    nms <- names(ans[[1]])
    n_col <- length(ans[[1]])
    out <- vector("list", n_col)
    for (i in seq_len(n_col)) {
      out[[i]] <- unlist(lapply(ans, `[[`, i), recursive = FALSE, 
                         use.names = FALSE)
    }
    names(out) <- nms
  }
  
  # good for code safety
  if (debug) {
    if (!is.numeric(out$scan_num))
      stop("Anticipating numeric values for `scan_num`.")
    
    ord <- order(out$scan_num)
    
    if (!identical(out$scan_num, out$scan_num[ord]))
      stop("Anticipating ascending `scan_num`.")
  }

  ###
  ##  may trigger data filtraiton by `scan_num` here...
  ###

  out$scan_num <- seq_along(out$scan_num)
  out$raw_file <- raw_file # scalar
  out_name <- paste0(raw_file, ".rds")
  qs::qsave(out, file.path(temp_dir, out_name), preset = "fast")
  message("Completed PASEF processing: ", out_name)

  attr(out_name, "is_dia") <- FALSE
  attr(out_name, "mzml_type") <- "raw" # a token for Mzion de-isotoping
  
  invisible(out_name)
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


#' Separates PASEF MS2Info
#' 
#' @param data Line data of MS frames (containing both MS1 and MS2).
#' @param iso_info Data frame of ms2info.
sep_pasef_ms2info <- function (data, iso_info)
{
  len <- length(data) # number of FRAMEs
  d_1 <- data[[1]] # the first FRAME in the chunk
  d_n <- data[[len]] # the last FRAME in the chunk
  sta <- as.integer(d_1[match("FRAME", d_1) + 1L])
  end <- as.integer(d_n[match("FRAME", d_n) + 1L])
  
  # values: MS2Frame > MS1Frame
  iso_info[with(iso_info, MS2Frame >= sta & MS2Frame <= end), ]
}


#' Helper of \link{extract_pasef_frame}
#'
#' (1) Pairs MS2 data and iso_info by frame IDs; (2) Collapses data by Precursor
#' groups.
#'
#' @param mdata A list of PASEF frames. Each list entry corresponding to one
#'   frame (over the ion-mobility dimension).
#' @param iso_info The information of MS2 isolation windows etc.
#' @param keys The key names of output lists.
#' @param step A step size for mass binning.
hextract_pasef <- function (
    mdata, iso_info, 
    keys = c("msx_moverzs", "msx_ints", "msx_ns", "ms1_fr", 
             "ms1_moverzs", "ms1_ints", "ms1_charges", 
             "scan_title", "ms_level", "ret_time", "scan_num", 
             "orig_scan", "iso_ctr", "iso_lwr", "iso_upr", "mobility"), 
    step = 1.6e-5)
{
  lens <- lengths(mdata) # numbers of lines in each frame
  ms_levs <- mapply(function (x, n) x[n - 3], mdata, lens, 
                    SIMPLIFY = TRUE, USE.NAMES = FALSE) |>
    as.integer()
  oks1 <- ms_levs == 1L
  oks2 <- .Internal(which(!oks1))
  
  ## MS1 and MS2
  out2 <- extract_pasef_ms2(ms2data = mdata[oks2], lens2 = lens[oks2], 
                            iso_info = iso_info, keys = keys, step = step)
  out1 <- extract_pasef_ms1(mdata = mdata[oks1], keys = keys)
  
  ## put together
  if (FALSE) {
    df2 <- tibble::tibble(
      msx_moverzs = out2$msx_moverzs, 
      msx_ints = out2$msx_ints, 
      msx_ns = out2$msx_ns, 
      ms1_fr = out2$ms1_fr, 
      ms1_moverzs = out2$ms1_moverzs, 
      ms1_ints = out2$ms1_ints, 
      ms1_charges = out2$ms1_charges, 
      scan_title = out2$scan_title, 
      ms_level = out2$ms_level, 
      ret_time = out2$ret_time, 
      scan_num = out2$scan_num, 
      orig_scan = out2$orig_scan, 
      iso_ctr = out2$iso_ctr, 
      iso_lwr = out2$iso_lwr, 
      iso_upr = out2$iso_upr, 
      mobility = out2$mobility, )
    df2x <- df2 |> dplyr::filter(ms1_fr == 9920)
    z2 <- df2x[39, ]
  }
  
  out <- mapply(`c`, out1, out2, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  names(out) <- names(out1)

  ord <- order(out[["scan_num"]])
  for (i in seq_along(out)) {
    out[[i]] <- out[[i]][ord]
  }

  out
}


#' Extracts PASEF MS1 data.
#'
#' @param mdata A list of \emph{MS1} PASEF frames. Each list entry corresponding
#'   to one MS1 frame.
#' @param keys The key names of output lists.
extract_pasef_ms1 <- function (mdata, keys)
{
  # 10 -> 20 -> 50; 100 too much
  ans1 <- lapply(mdata,  extract_pasef_frame, ms_lev = 1L, ymin = 50)
  ans1 <- ans1[lengths(ans1) > 0L]
  
  if (FALSE) {
    # qs::qsave(ans1, "~/ans1_extract_pasef_ms1.rds", preset = "fast")
    ans1 <- qs::qread("~/ans1_extract_pasef_ms1.rds")
    scans <- sapply(ans1, `[[`, "scan_num")
    # i <- which(scans == 9920)
    i <- which(scans == 9915) # i <- 114
    ai <- ans1[[i]]
    df <- data.frame(x = ai$msx_moverzs, y = ai$msx_ints)
    dfx <- df |> dplyr::filter(x >= 926, x <= 928)
    ggplot2::ggplot() + 
      ggplot2::geom_segment(dfx, mapping = aes(x = x, y = y, xend = x, yend = 0), 
                            color = "gray", linewidth = .1)
    
    tempdata1 <- centroid_pasefms(ai$msx_moverzs, ai$msx_ints)
    df <- data.frame(x = tempdata1$x, y = tempdata1$y)
    dfx <- df |> dplyr::filter(x >= 460, x <= 463)
  }

  for (i in seq_along(ans1)) {
    ai <- ans1[[i]]
    # may use weighted-mean... 
    tempdata1 <- centroid_pasefms(ai$msx_moverzs, ai$msx_ints)
    ans1[[i]]$msx_moverzs <- tempdata1[["x"]]
    ans1[[i]]$msx_ints <- tempdata1[["y"]]
  }
  rm(list = c("ai", "tempdata1"))

  nms1 <- names(ans1[[1]])
  out1 <- vector("list", length(nms1))
  for (i in seq_along(out1)) {
    out1[[i]] <- lapply(ans1, `[[`, i)
  }
  names(out1) <- nms1

  out1[["scan_title"]] <- 
    unlist(out1[["scan_title"]], recursive = FALSE, use.names = FALSE)
  out1[["ms_level"]] <- 
    unlist(out1[["ms_level"]], recursive = FALSE, use.names = FALSE)
  out1[["ret_time"]] <- 
    unlist(out1[["ret_time"]], recursive = FALSE, use.names = FALSE)
  out1[["scan_num"]] <- 
    unlist(out1[["scan_num"]], recursive = FALSE, use.names = FALSE)
  out1[["orig_scan"]] <- out1[["scan_num"]]
  out1[["mobility"]] <- 
    unlist(out1[["mobility"]], recursive = FALSE, use.names = FALSE)
  out1[["msx_ns"]] <- lengths(out1[["msx_moverzs"]])
  
  len1 <- length(out1[["msx_ns"]])
  na_ints1 <- rep_len(NA_integer_, len1)
  na_reals1 <- rep_len(NA_real_, len1)
  out1[["ms1_fr"]] <- na_ints1
  out1[["ms1_moverzs"]] <- na_reals1
  out1[["ms1_charges"]] <- na_ints1
  out1[["ms1_ints"]] <- na_ints1
  out1[["iso_ctr"]] <- na_reals1
  out1[["iso_lwr"]] <- na_reals1
  out1[["iso_upr"]] <- na_reals1
  
  out1[["scan_num"]] <- as.numeric(out1[["scan_num"]])
  
  if (length(keys) != length(out1)) {
    stop("Developer: unequal numbers of columns in PASEF-MS1 data.")
  }
  
  message("Completed PASEF MS1 extraction.")
  
  out1[keys]
}


#' Extracts PASEF MS1 data.
#'
#' @param ms2data A list of \emph{MS2} PASEF frames. Each list entry
#'   corresponding to one MS2 frame.
#' @param lens2 The number of lines in each MS2 frame.
#' @param iso_info The information of MS2 isolation windows etc.
#' @param keys The key names of output lists.
#' @param step A step size for mass binning.
extract_pasef_ms2 <- function (ms2data, lens2, iso_info, keys, step = 1.6e-5)
{
  # (1) pair MS2 data with `iso_info` by FRAME
  ms2frs <- mapply(function (x, n) x[n - 7L], ms2data, lens2, 
                   SIMPLIFY = TRUE, USE.NAMES = FALSE) |>
    as.integer()
  iso_info <- iso_info[with(iso_info, MS2Frame %in% ms2frs), ]
  iso_info <- split(iso_info, iso_info$MS2Frame) 
  
  if (!identical(as.integer(names(iso_info)), ms2frs)) {
    stop("Developer: mismatches in PASEF frame IDs.")
  }
  
  if (length(ms2data) > length(iso_info)) {
    ok_frs <- ms2frs %in% names(iso_info)
    ms2data <- ms2data[ok_frs]
    ms2frs <- ms2frs[ok_frs]
    rm(list = "ok_frs")
  }
  
  # each entry corresponds to one frame (multiple slices in each frame)
  if (FALSE) {
    rng <- 1149:1158
    ans2 <- mapply(add_pasef_ms2iso, ms2data[rng], iso_info[rng])
    ans2 <- unlist(ans2, recursive = FALSE, use.names = FALSE)
    precursors <- lapply(ans2, `[[`, "precursor") |>
      unlist(recursive = FALSE, use.names = FALSE)

    ax2 <- split(ans2, precursors)
    z2 <- ax2[[which(names(ax2) == 44101)]] # 44080
    z <- group_ms2pasef_by_precursors(z2)
  }
  
  ans2 <- mapply(add_pasef_ms2iso, ms2data, iso_info)
  ans2 <- ans2[lengths(ans2) > 0L]
  # flattens the slices in each frame
  ans2 <- unlist(ans2, recursive = FALSE, use.names = FALSE)
  
  ## (2) group MS2 slices by precursors
  precursors <- lapply(ans2, `[[`, "precursor") |>
    unlist(recursive = FALSE, use.names = FALSE)
  ans2 <- lapply(split(ans2, precursors), group_ms2pasef_by_precursors, 
                 step = step)
  
  nms2 <- names(ans2[[1]])
  out2 <- vector("list", length(nms2))
  
  for (i in seq_along(out2)) {
    out2[[i]] <- lapply(ans2, `[[`, i)
  }
  names(out2) <- nms2

  out2[["scan_title"]] <- 
    unlist(out2[["scan_title"]], recursive = FALSE, use.names = FALSE)
  out2[["ms_level"]] <- 
    unlist(out2[["ms_level"]], recursive = FALSE, use.names = FALSE)
  out2[["msx_ns"]] <- 
    unlist(out2[["msx_ns"]], recursive = FALSE, use.names = FALSE)
  out2[["ms1_fr"]] <- 
    unlist(out2[["ms1_fr"]], recursive = FALSE, use.names = FALSE)
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
  
  # out2[["scan_num"]] <- as.numeric(out2[["scan_num"]])
  
  if (length(keys) != length(out2)) {
    stop("Developer: unequal numbers of columns in PASEF MS2 data.")
  }

  message("Completed PASEF MS2 extraction.")
  
  out2 <- out2[keys]
}


#' Sum PASEF MS1 peak area
#' 
#' @param xs A vector of ascending m-over-z values.
#' @param ys A vector of intensity values.
#' @param reso The resolution of a peak.
#' @param maxn The maximum number of peaks.
#' @param ymin The minimum Y values for considering in peak centroiding.
#' @param tol The tolerance of Y for defining a peak profile.
#' @importFrom fastmatch %fin%
#' @examples
#' # example code
#' mzion:::centroid_pasefms(c(500), c(10))
#' mzion:::centroid_pasefms(c(500, 500.01), c(1, 10))
#' mzion:::centroid_pasefms(c(500, 500.01, 500.2), c(1, 10, 20))
#' mzion:::centroid_pasefms(c(500, 500.01, 500.2), c(10, 5, 2))
#' 
#' # ignore the trailing half peak
#' mzion:::centroid_pasefms(500 + .01 * 0:3, c(1, 10, 2, 5))
#' 
#' mzion:::centroid_pasefms(500 + .01 * 0:4, c(1, 10, 2, 5, 3))
#' 
#' # trailing max at the last position
#' mzion:::centroid_pasefms(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,33,23,22,41))
#' 
#' mzion:::centroid_pasefms(c(231.0019,371.1024,519.1426,542.3826,599.9552), c(12,15,23,30,41))
centroid_pasefms <- function (xs, ys, reso = 60000, maxn = 2000L, ymin = 100L, 
                              tol = .10)
{
  len <- length(ys)
  
  if (!len) {
    return(NULL)
  }

  if (len == 1L) {
    return(list(x = xs, y = ys))
  }

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
  
  ns <- min(ns, maxn)
  xvals <- vector("numeric", ns)
  yvals <- vector("integer", ns)
  ord <- order(ys[ps], decreasing = TRUE)

  ct <- 0L
  for (i in 1:ns) {
    if (ct == maxn)
      break
    
    oi <- ord[[i]]
    p <- ps[oi]
    
    if (is.na(p))
      next
    
    x <- xs[[p]]
    w <- x / reso * 1.5
    idxes <- .Internal(which(xs >= x - w & xs <= x + w))
    ysubs <- ys[idxes]
    imax <- .Internal(which(idxes == p)) # relative to ysubs
    
    yvals[[i]] <- sum_pasef_ms1(ysubs, imax)
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
    ord <- order(xvals)
    xvals <- xvals[ord]
    yvals <- yvals[ord]
  }
  
  list (x = xvals, y = yvals)
}


#' Sum peak area
#' 
#' @param ys A sub vector of intensity values around a peak.
#' @param imax The index of the peak position in \code{ys}.
sum_pasef_ms1 <- function (ys, imax)
{
  len <- length(ys)

  if (len == 1L) {
    return(ys)
  }
  
  if (len == 2L) {
    return(sum(ys))
  }
  
  yval <- ys[[imax]]

  if (imax < len) {
    pr1 <- imax + 1L
    yr1 <- ys[pr1]
    yval <- yval + yr1
    
    if (pr1 < len) {
      for (j in (pr1 + 1L):len) {
        yr2 <- ys[[j]]
        
        if (yr1 / yr2 >= .667) { # || yr2 <= 100
          yval <- yval + yr2
          yr1 <- yr2
        }
        else {
          d <- yr1 * .667
          yval <- yval + d
          yr1 <- yr2 - d
        }
      }
    }
  }

  if (imax > 1L) {
    pl1 <- imax - 1L
    yl1 <- ys[[pl1]]
    yval <- yval + yl1
    
    if (pl1 > 1L) {
      for (j in (pl1 - 1L):1) {
        yl2 <- ys[[j]]
        
        if (yl1 / yl2 >= .667) { # || yl2 <= 100
          yval <- yval + yl2
          yl1 <- yl2
        }
        else {
          d <- yl1 * .667
          yval <- yval + d
          yl1 <- yl2 - d
        }
      }
    }
  }

  as.integer(yval)
}


#' Collapses timsTOF X and Y values.
#' 
#' Not used.
#' 
#' @param xs A vector of m-over-z values.
#' @param ys A vector of intensity values.
#' @param maxn_peaks The maximum number of peaks for consideration.
#' @param yco The noise levels of intensity values.
#' @param ymin The minimum value of intensity for consideration as a peak.
#' @param tol The tolerance for masses binning.
collapse_rawtims_xys <- function (xs, ys, maxn_peaks = 1000L, yco = 20L, 
                                  ymin = 30L, tol = 2e-5)
{
  oks <- .Internal(which(ys > yco))
  xs <- xs[oks]
  ys <- ys[oks] - yco
  
  nx <- length(xs)
  maxn_peaks <- min(nx, maxn_peaks)
  peaks <- rep_len(NA_real_, maxn_peaks)
  intens <- rep_len(NA_integer_, maxn_peaks)
  
  p <- 1L
  while(p <= maxn_peaks) {
    imax <- .Internal(which.max(ys))
    ymax <- ys[imax]
    
    if (ymax < ymin) {
      break
    }
    
    is_peak <- if (imax > 1L && imax < nx) {
      ys[imax - 1L] < ymax && ys[imax + 1L] < ymax
    }
    else {
      FALSE
    }
    
    if (!is_peak) {
      ys[imax] <- 0L
      p <- p + 1L
      next
    }
    
    xmax <- xs[imax]
    upr <- min(nx, imax + 5L)
    gap <- min(diff(xs[imax:upr]))
    
    # next gather the peaks within +/- 20ppm & >= 3% or 5% of base peak
    xwin <- xmax * tol
    nbin <- as.integer(xwin / gap)
    
    rng <- max(1L, imax - nbin):min(nx, imax + nbin)
    
    peaks[p] <- xmax
    intens[p] <- sum(ys[rng])
    ys[rng] <- 0
    
    # also look into extended window of +/- 20 ppm based on intens[p]
    # if <= 3% -> 0

    # xmax - xwin
    # xmax + xwin
    
    p <- p + 1L
  }
  
  # intensity only approximately decreasing because of sum(ys[rng])
  oks2 <- .Internal(which(!is.na(peaks)))
  peaks <- peaks[oks2]
  intens <- intens[oks2]
  
  ord <- order(peaks)
  out <- list(xs = peaks[ord], ys = intens[ord])
}


#' Collapses PASEF MS2 data by precursor IDs.
#' 
#' @param dat MS2 data under the same precursor ID.
#' @param lwr A lower bound as the starting point in mass binning.
#' @param step A step size for mass binning.
group_ms2pasef_by_precursors <- function (dat, lwr = 115L, step = 1.6e-5)
{
  oks <- .Internal(which(lapply(dat, `[[`, "msx_ns") > 0L))
  len <- length(oks)
  
  if (!len)
    return(NULL)

  dat <- dat[oks]
  
  xys <- collapse_mms1ints(
    lapply(dat, `[[`, "msx_moverzs"), lapply(dat, `[[`, "msx_ints"), 
    lwr = lwr, step = step, reord = FALSE, cleanup = FALSE, 
    sum_y = TRUE, add_colnames = FALSE)
  xys <- calc_ms1xys(xys[["x"]], xys[["y"]])

  msx_moverzs <- xys$x
  msx_ints <- xys$y
  # ns <- xys$n
  msx_ns <- length(msx_ints)
  
  ms1_frs  <- lapply(dat, `[[`, "ms1_fr") # should be all the same
  ret_time <- lapply(dat, `[[`, "ret_time")
  ret_time <- unlist(ret_time, recursive = FALSE, use.names = FALSE)
  ret_time <- sum(ret_time)/len
  
  scan_num_all <- lapply(dat, `[[`, "scan_num")
  scan_num_all <- unlist(scan_num_all, recursive = FALSE, use.names = FALSE)
  scan_num_all <- paste0(scan_num_all, collapse = ",")
  
  mobility <- lapply(dat, `[[`, "mobility")
  mobility <- unlist(mobility, recursive = FALSE, use.names = FALSE)
  mobility <- sum(mobility)/len
  
  rng_slice <- lapply(dat, `[[`, "rng_slice")
  rng_slice <- unlist(rng_slice, recursive = FALSE, use.names = FALSE)
  rng_slice <- paste0(rng_slice, collapse = ",")
  
  dat_1 <- dat[[1]]
  
  # uses the fist one for staggering with MS1 scan numbers
  scan_num <- dat_1[["scan_num"]]
  scan_title <- dat_1[["scan_title"]]
  ms1_fr <- ms1_frs[[1]]
  scan_title <- paste0(scan_title, "; MS1: ", ms1_fr, 
                       "; scans: ", scan_num_all, 
                       "[", rng_slice, "]", "; Cmpd: ", dat_1[["precursor"]])
  ms_level <- dat_1[["ms_level"]]
  iso_ctr <- dat_1[["iso_ctr"]]
  iso_lwr <- dat_1[["iso_lwr"]]
  iso_upr <- dat_1[["iso_upr"]]

  list(
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    msx_ns = msx_ns,
    scan_title = scan_title,
    ms_level = dat_1[["ms_level"]], 
    ms1_fr = ms1_fr,
    ret_time = ret_time,
    scan_num = scan_num,
    orig_scan = scan_num_all, 
    iso_ctr = dat_1[["iso_ctr"]],
    iso_lwr = dat_1[["iso_lwr"]],
    iso_upr = dat_1[["iso_upr"]], 
    mobility = mobility, 
    ms1_moverzs = dat_1[["ms1_moverzs"]], 
    ms1_ints = dat_1[["ms1_ints"]], 
    ms1_charges = dat_1[["ms1_charges"]])
}


#' Add isolation information to an MS2 frame.
#'
#' @param data Lines of MS2 data in a \code{frame}.
#' @param iso_info A subset of MS2 isolation information at the current
#'   \code{frame}.
#' @param min_ms2n The minimum number of MS2 in a slices for considerations.
#' @return A vector of lists. Each list corresponding to a slice of MS2 between
#'   ScanNumBegin and ScanNumEnd in a frame, with additional fields of
#'   IsolationMz, IsolationWidth, Precursor etc.
add_pasef_ms2iso <- function (data, iso_info, min_ms2n = 0L)
{
  options(warn = 1)
  
  ## results from one MS2 frame
  ans <- extract_pasef_frame(data, ms_lev = 2L, ymin = 10)
  
  if (!length(ans)) {
    return(NULL)
  }
  
  ## Combine data slices by ranges of ScanNumBegin and ScanNumEnd
  breaks <- findInterval(ans$slices, iso_info$ScanNumEnd)
  slices <- split(ans$slices, breaks)
  len <- length(slices)
  
  if (!len) {
    return(NULL)
  }
  
  # removes iso_info rows without matched scan ranges to ans$slices
  if (len < nrow(iso_info)) {
    iso_info <- iso_info[as.integer(names(slices)) + 1L, ]
  }
  
  xys <- mapply(collapse_pasef_xys, 
                split(ans$msx_moverzs, breaks), 
                split(ans$msx_ints, breaks), 
                SIMPLIFY = FALSE, USE.NAMES = FALSE)
  msx_moverzs <- lapply(xys, `[[`, 1)
  msx_ints <- lapply(xys, `[[`, 2)
  mobs <- split(ans$mobility, breaks)
  mobs <- lapply(mobs, function (x) sum(x) / length(x))
  
  ## clean ups
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
  else {
    msx_ns <- lengths(msx_moverzs)
    iso_ctrs <- iso_info$IsolationMz
    iso_widths <- iso_info$IsolationWidth
    precursors <- iso_info$Precursor
    ms1_moverzs <- iso_info$MonoisotpoicMz
    ms1_ints <- iso_info$Intensity
    ms1_charges <- iso_info$Charge
    ms1_frs <- iso_info$MS1Frame
  }

  ## outputs
  out <- vector("list", len)
  rngs <- mapply(function (x, y) paste0(x, "-", y), 
                 lapply(slices, `[[`, 1), # starts
                 lapply(slices, function (x) x[length(x)]), # ends
                 SIMPLIFY = TRUE, USE.NAMES = FALSE)
  scan_titles <- rep_len(ans$scan_title, len)
  ms_levels <- rep_len(ans$ms_level, len)
  ret_times <- rep_len(ans$ret_time, len)
  scan_nums <- ans$scan_num + seq_len(len) / 10^nchar(as.character(len))
  
  ## MonoisotpoicMz is almost always -.8 and -0.1 lower than IsolationMz?
  ## IsolationWidth is 2 at ~ IsolationMz < 720 and 3 at IsolationMz > 800.
  ## The wide IsolationWidth may be intended to include lower m/z values in 
  #   deisotoping, but Mzion has its own extension (+/-2 around IsolationMz).

  # iso_lwrs <- iso_ctrs - iso_widths
  # iso_uprs <- iso_ctrs
  half_widths <- iso_widths / 2
  iso_lwrs <- iso_ctrs - half_widths
  iso_uprs <- iso_ctrs + half_widths

  for (i in 1:len) {
    out[[i]] <- list(
      msx_moverzs = msx_moverzs[[i]], 
      msx_ints = msx_ints[[i]],
      rng_slice = rngs[[i]], 
      msx_ns = msx_ns[[i]],
      scan_title = scan_titles[[i]],
      ms_level = ms_levels[[i]],
      ret_time = ret_times[[i]],
      scan_num = scan_nums[[i]],
      ms1_fr = ms1_frs[[i]],
      iso_ctr = iso_ctrs[[i]],
      iso_lwr = iso_lwrs[[i]],
      iso_upr = iso_uprs[[i]],
      mobility = mobs[[i]], 
      precursor = precursors[[i]], 
      ms1_moverzs = ms1_moverzs[[i]], 
      ms1_ints = ms1_ints[[i]], 
      ms1_charges = ms1_charges[[i]]
    )
  }
  
  out  
}


#' Extracts MS1 or MS2 data from a PASEF frame.
#' 
#' @param data A PASEF frame (with multiple slices).
#' @param ms_lev The level of MS data.
#' @param ymin The cut-off of intensity.
#' @param ymax The maximum intensity.
extract_pasef_frame <- function (data, ms_lev = 1L, ymin = 10, ymax = 1E7)
{
  len <- length(data) # the number of lines
  
  # fields common within a frame
  title <- data[len - 1L]
  # ms_lev <- as.integer(data[len - 3L])
  ret_time <- as.numeric(data[len - 5L])
  scan_num <- as.numeric(data[len - 7L])
  
  # fields specific to each slice in a frame
  idx_slices <- .Internal(which(stringi::stri_cmp_eq(data, "SLICE"))) + 1L
  slices <- as.integer(data[idx_slices])
  npeaks <- as.integer(data[idx_slices + 2L])
  mobils <- as.numeric(data[idx_slices + 4L])
  idx_ys <- idx_slices - 2L
  idx_xs <- idx_ys - 1L

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
  
  # clean up by Y values
  for (i in seq_along(xs)) {
    ysi <- ys[[i]]
    oki <- .Internal(which(ysi >= ymin))
    
    if (length(oki)) {
      ys[[i]] <- ysi[oki] # works at rhs numeric(0)
      xs[[i]] <- xs[[i]][oki]
    }
  }
  
  # clean up by X values
  oks <- .Internal(which(lengths(xs) > 0L))
  
  if (!length(oks)) {
    return(NULL)
  }
  
  xs <- xs[oks]
  ys <- ys[oks]
  slices <- slices[oks]
  mobils <- mobils[oks]
  # npeaks <- npeaks[oks]
  # idx_xs <- idx_xs[oks]
  # idx_ys <- idx_ys[oks]

  # need to collapse MS2 by IM slices
  if (ms_lev == 2L) {
    return(list(scan_title = title,
                ms_level = ms_lev,
                ret_time = ret_time,
                scan_num = scan_num,
                
                # vectors below
                msx_moverzs = xs, 
                msx_ints = ys, 
                mobility = mobils,
                slices = slices))
  }
  
  # collapses all MS1 scans in a Frame
  xys <- collapse_pasef_xys(xs, ys)
  xs <- xys$x
  ys <- xys$y

  list(scan_title = title,
       ms_level = ms_lev,
       ret_time = ret_time,
       scan_num = scan_num,
       mobility = sum(mobils)/length(mobils),
       
       # vectors below
       msx_moverzs = xs, 
       msx_ints = ys)
}


#' Collapse PASEF X and Y values within an MS1 or MS2 frame.
#' 
#' @param xs Slices of m-over-z values.
#' @param ys Slices of intensity values.
collapse_pasef_xys <- function (xs, ys)
{
  if (!length(xs))
    return(NULL)
  
  xs <- unlist(xs, recursive = FALSE, use.names = FALSE)
  ys <- unlist(ys, recursive = FALSE, use.names = FALSE)
  
  if (length(xs) > 1L) {
    ord <- order(xs)
    xs <- xs[ord]
    ys <- ys[ord]
  }

  if (anyDuplicated(xs)) {
    ys <- split(ys, xs)
    xs <- split(xs, xs)
    dups <- lengths(xs) > 1L
    xs[dups] <- lapply(xs[dups], `[[`, 1)
    ys[dups] <- lapply(ys[dups], sum, na.rm = TRUE)
    xs <- unlist(xs, recursive = FALSE, use.names = FALSE)
    ys <- unlist(ys, recursive = FALSE, use.names = FALSE)
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
  
  spectra <- file.path(mgf_path, raw_file, "spectra.txt")
  precursors <- file.path(mgf_path, raw_file, "ms2precursors.txt")
  raw_full <- file.path(mgf_path, raw_file)
  stdout <- tempfile(tmpdir = mgf_path, fileext = ".stdout")
  stderr <- tempfile(tmpdir = mgf_path, fileext = ".stderr" )

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
  
  message("Complete ReadRAW: ", raw_file)
  unlink(c(stdout, stderr))
  
  list(spectra, precursors)
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


