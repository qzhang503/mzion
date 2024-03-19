#' Subsets full MS by the universe of monoisotopic moverzs.
#'
#' Some non-monoisotopic moverzs can be removed.
#'
#' @param xs Vectors of full-spectrum m/z values.
#' @param ys Vectors of full-spectrum intensity values.
#' @param ms1s Vectors of MS2-specific MS1 m/z values (corresponding to
#'   \code{ms_level == 2L}).
#' @param from The starting m/z value for calculating bin indexes.
#' @param step The size of bins.
#' @param gap The gap of MS1 scans.
#' @importFrom fastmatch %fin%
subMSfull <- function (xs, ys, ms1s, from = 200L, step = 1E-5, gap = 256L)
{
  lens <- lengths(xs)
  len <- length(xs)
  yout <- xout <- vector("list", len)

  # for each xi, keeps entries found in ms1s.
  if (gap >= len) {
    ms1s <- .Internal(unlist(ms1s, recursive = FALSE, use.names = FALSE))
    ms1s <- index_mz(ms1s, from, step)
    ms1s <- sort(unique(ms1s))

    for (i in 1:len) {
      xi <- xs[[i]]
      li <- lens[[i]]
      
      if (!li)
        next
      
      ix <- as.integer(ceiling(log(xi/from)/log(1+step)))
      ps0 <- fastmatch::fmatch(ix, ms1s)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s)
      i0 <- .Internal(which(!is.na(ps0)))
      i1 <- .Internal(which(!is.na(ps1)))
      i2 <- .Internal(which(!is.na(ps2)))
      
      i1 <- i1[!i1 %fin% i0]
      i2 <- i2[!(i2 %fin% i0 | i2 %fin% i1)]
      i012 <- c(i0, i1, i2)
      i012 <- sort(i012)
      xout[[i]] <- xi[i012]
      yout[[i]] <- ys[[i]][i012]
    }
  }
  else {
    for (i in 1:len) {
      xi <- xs[[i]]
      li <- lens[[i]]
      
      if (!li)
        next
      
      sta <- max(1, i - gap)
      end <- min(i + gap, len)
      ms1s_sub <- .Internal(unlist(ms1s[sta:end], recursive = FALSE, use.names = FALSE))
      ms1s_sub <- index_mz(ms1s_sub, from, step)
      ms1s_sub <- sort(unique(ms1s_sub))
      
      ix <- as.integer(ceiling(log(xi/from)/log(1+step)))
      ps0 <- fastmatch::fmatch(ix, ms1s_sub)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s_sub)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s_sub)
      i0 <- .Internal(which(!is.na(ps0)))
      i1 <- .Internal(which(!is.na(ps1)))
      i2 <- .Internal(which(!is.na(ps2)))
      
      i1 <- i1[!i1 %fin% i0]
      i2 <- i2[!(i2 %fin% i0 | i2 %fin% i1)]
      i012 <- c(i0, i1, i2)
      i012 <- sort(i012) # logical(0) if all NA
      xout[[i]] <- xi[i012]
      yout[[i]] <- ys[[i]][i012]
    }
  }

  list(x = xout, y = yout)
}


#' Splits \code{df} for parallel tracing.
#' 
#' @param df A data frame of staggering MS1 and MS2 entries.
#' @param from The starting point for mass binning.
#' @param step The step size for mass binning.
#' @param n_chunks The number of chunks.
#' @param gap A gap size for forward and backward looking of precursors.
pretraceXY <- function (df, from = 200L, step = 1e-5, n_chunks = 4L, gap = 256L)
{
  # msx_moverzs at ms_level == 1L: full-spectrum ms1_moverzs
  # msx_charges at ms_level == 1L are list(NULL)
  
  cols1 <- c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", 
             "msx_moverzs", "msx_ints", "msx_charges", "orig_scan")
  rows1 <- df[["ms_level"]] == 1L
  df1 <- df[rows1, cols1]
  len1 <- nrow(df1)

  # Remove non-essential MS1 x and y values
  ans <- subMSfull(
    xs = df1$msx_moverzs, ys = df1$msx_ints, ms1s = df1$ms1_moverz, 
    from = from, step = step, gap = gap)
  
  df1$msx_moverzs <- ans$x # moverzs in ans$x are sorted
  df1$msx_ints <- ans$y
  rm(list = "ans")
  gc()
  
  # at least two chunks
  if (n_chunks <= 1L)
    n_chunks <- 2L
  
  df1s <- chunksplit(df1, n_chunks, type = "row")
  end1s <- cumsum(lapply(df1s, nrow))
  sta1s <- c(1L, end1s[1:(n_chunks-1L)] + 1L)
  
  # Adds 2-min gaps before and after
  gaps <- lapply(df1s, function (x) ceiling(min(gap, nrow(x)/2L)))
  df1s_bf <- df1s_af <- vector("list", n_chunks)
  
  for (i in 2:n_chunks) {
    df1s_bf[[i]] <- head(df1s[[i]], gaps[[i]])
  }
  
  for (i in 1:(n_chunks - 1L)) {
    df1s_af[[i]] <- tail(df1s[[i]], gaps[[i]])
  }
  
  df1s[[1]] <- dplyr::bind_rows(df1s[[1]], df1s_bf[[2]])
  df1s[[n_chunks]] <- dplyr::bind_rows(df1s_af[[n_chunks-1]], df1s[[n_chunks]])
  
  if (n_chunks > 2L) {
    for (i in 2:(n_chunks-1L))
      df1s[[i]] <- dplyr::bind_rows(df1s_af[[i-1]], df1s[[i]], df1s_bf[[i+1]])
  }
  rm(list = c("df1s_bf", "df1s_af", "df1"))
  
  ##  Splits `df` with bracketing entries
  # values: row indexes in `df`, length == nrow(df1) at pad_nas = TRUE
  ms1_stas <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)$ms1_stas
  cols <- c("ms_level", "ms1_moverz", "ms1_int")
  ms1_stax <- dfs <- vector("list", n_chunks)
  
  for (i in 1:n_chunks) {
    stai <- sta1s[[i]] # stai-th MS1 scan
    ms1_stax[[i]] <- ms1_stas[stai] # the corresponding row index in df
  }

  for (i in 1:(n_chunks - 1L)) {
    rowx <- ms1_stax[[i]]:(ms1_stax[[i+1]] - 1L) # ms2_endx can be NA
    dfs[[i]] <- df[rowx, cols] # both df1 and df2 data
  }
  
  dfs[[n_chunks]] <- df[ms1_stax[[n_chunks]]:nrow(df), cols]

  list(dfs = dfs, df1s = df1s, gaps = gaps)
}


#' Helper of \link{traceXY}.
#'
#' @param xs Vectors of full-spectrum MS1 m/z values.
#' @param ys Vectors of full-spectrum MS1 intensities.
#' @param ss Vectors of MS1 scan numbers.
#' @param df A data frame of staggering MS1 and MS2 in the same range of
#'   \code{xs}.
#' @param gap_bf A preceding gap size.
#' @param gap_af A following gap size.
#' @param from The starting point for mass binning.
#' @param step The step size for mass binning.
#' @inheritParams matchMS
htraceXY <- function (xs, ys, ss, df, gap_bf = 256L, gap_af = 256L, 
                      n_mdda_flanks = 6L, from = 200L, step = 1E5)
{
  if (all(lengths(xs) == 0L)) {
    # length(xs) - gap_bf - gap_af - nrow(df) == 0L
    li <- nrow(df)
    null <- rep_len(list(NULL), li)
    
    return(
      tibble::tibble(
        ms_level = df$ms_level, 
        ms1_moverz = null, 
        ms1_int = null, 
        apex_scan_num = null)
    )
  }
  
  mat <- traceXY(
    xs = xs, ys = ys, ss = ss, n_mdda_flanks = n_mdda_flanks, from = from, 
    step = step, reord = FALSE, cleanup = FALSE, # otherwise rows drop
    replace_ms1_by_apex = TRUE)
  
  matx <- mat[["x"]]
  maty <- mat[["y"]]
  nr <- nrow(matx)
  nc <- ncol(matx)
  
  ## apes, rngs and scans: each vector corresponds to a column of mass
  #  apes not ordered by orig_scan, one-hit-wonders goes first
  apes <- mat[["p"]]
  # rngs <- mat[["range"]]
  
  # apes correspond to the row numbers of matx and ss to orig_scan
  scan_apexs <- vector("list", nc)
  
  for (i in 1:nc) {
    rows <- apes[[i]]
    scan_apexs[[i]] <- as.integer(ss[rows])
  }

  rm(list = "mat")
  gc()
  
  if (gap_bf) {
    if (gap_af) { # middle
      sta <- gap_bf + 1L
      end <- nr - gap_af
    }
    else { # last
      sta <- gap_bf + 1L
      end <- nr
    }
  }
  else { # first
    sta <- 1L
    end <- nr - gap_af
  }
  
  if (FALSE) {
    matx <- matx[sta:end, ]
    maty <- maty[sta:end, ]
    ss <- ss[sta:end]
    # look current
    df <- updateMS1Int(df = df, matx = matx, maty = maty, from = from, 
                       step = step)
  }
  else {
    # look both before and after scans
    df <- updateMS1Int2(
      df = df, matx = matx, maty = maty, row_sta = sta, row_end = end, 
      scan_apexs = scan_apexs, from = from, step = step)
  }
}


#' Helper of MS1 tracing.
#'
#' @param xs Vectors of full-spectrum MS1 moverzs.
#' @param ys Vectors of full-spectrum MS1 intensities.
#' @param ss Vectors of MS1 scan numbers.
#' @param step Step size.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @param reord Logical; re-order data or not.
#' @param cleanup Logical; to clean up xs, ys and zs or not. Set the value to
#'   FALSE to maintain one-to-one correspondence between input (data frame) and
#'   the outputs. This will help, e.g., keep track of scan numbers in the input.
#' @param replace_ms1_by_apex Logical; if TRUE, fill all entries within a gate
#'   by its apex values.
#' @inheritParams matchMS
traceXY <- function (xs, ys, ss, n_mdda_flanks = 6L, from = 115L, step = 1E-5, 
                     reord = TRUE, cleanup = FALSE, replace_ms1_by_apex = FALSE)
{
  lens <- lengths(xs)
  
  if (all(lens == 0L)) {
    # return(list(x = NULL, y = NULL, n = NULL, p = NULL, range = NULL))
  }

  if (reord) {
    for (i in seq_along(xs)) {
      xi <- xs[[i]]
      
      if (lens[[i]]) {
        ord <- order(xi)
        xs[[i]] <- xi[ord]
        ys[[i]] <- ys[[i]][ord]
      }
    }
    rm(list = c("xi", "ord"))
  }
  
  ## collapses MS data by the indexes of mass bins; 
  # two matrix outputs; rows: scans; columns: masses or intensities
  
  # xs can be numeric(0)?
  # cleanup = FALSE; otherwise rows drop
  # often coll == cleanup
  ans <- collapse_mms1ints(
    xs = xs, ys = ys, lwr = from, step = step, reord = FALSE, coll = FALSE, 
    cleanup = FALSE, add_colnames = TRUE)
  
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  rm(list = c("ans"))
  gc()
  
  ## traces MS data matrices across LC scans; rows: scans; columns: masses
  nrc <- dim(ansy)
  nr <- nrc[[1]]
  nc <- nrc[[2]]
  rm(list = "nrc")
  
  if (nr != length(xs)) {
    stop("Developer: rows drop during MS1 tracing.")
  }
  
  xmat <- ymat <- matrix(rep_len(NA_real_, nc * nr), ncol = nc)
  colnames(xmat) <- colnames(ymat) <- colnames(ansx)
  rownames(xmat) <- rownames(ymat) <- ss
  ranges <- apexes <- ns <- vector("list", nc)
  
  if (replace_ms1_by_apex) {
    for (i in 1:nc) {
      # removes peaks at intensity < 2% of base peak
      yi <- ansy[, i]
      oks <- .Internal(which(!is.na(yi)))
      yoks <- yi[oks]
      yoks[yoks < max(yoks) * .02] <- NA_real_
      yi[oks] <- yoks
      gates <- find_lc_gates(yi, n_dia_scans = n_mdda_flanks)
      # ansy[, i] <- yi

      # one-hit-wonders go first, not ordered scans
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]

      for (j in seq_along(rows)) {
        rgj <- rngs[[j]]
        rwj <- rows[[j]]
        xmat[rgj, i] <- ansx[rwj, i]
        ymat[rgj, i] <- ansy[rwj, i]
      }
    }
  }
  else {
    for (i in 1:nc) {
      gates <- find_lc_gates(ansy[, i], n_dia_scans = n_mdda_flanks)
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]
      
      xmat[rows, i] <- ansx[rows, i]
      ymat[rows, i] <- ansy[rows, i]
    }
  }
  
  rm(list = c("ansx", "ansy"))
  gc()
  
  list(x = xmat, y = ymat, n = ns, p = apexes, range = ranges)
}


#' Updates MS1 intensity with apex values.
#'
#' @param df A data frame.
#' @param matx The matrix of moverzs Y: by masses; X: by LC scans.
#' @param maty The matrix of intensities. Y: by masses; X: by LC scans.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
updateMS1Int <- function (df, matx, maty, from = 200L, step = 1E-5)
                          
{
  nrow <- nrow(matx)
  
  pos_levs <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")
  
  for (i in seq_along(ms2_stas)) {
    ms2sta <- ms2_stas[[i]]
    
    if (is.na(ms2sta))
      next
    
    ms2end <- ms2_ends[[i]]
    df2 <- df[ms2sta:ms2end, ]
    
    xs <- matx[i, ]
    ys <- maty[i, ]
    oks <- .Internal(which(!is.na(xs)))
    xs <- xs[oks]
    ys <- ys[oks]
    ixs <- as.integer(ceiling(log(xs/from)/log(1+step)))

    for (j in 1:nrow(df2)) {
      xsj <- df2[["ms1_moverz"]][[j]]
      ixsj <- as.integer(ceiling(log(xsj/from)/log(1+step)))
      ps0 <- match(ixsj, ixs)
      
      if (any(nas0 <- is.na(ps0))) {
        ps1 <- match(ixsj + 1L, ixs)
        ps2 <- match(ixsj - 1L, ixs)
        
        if (length(i0 <- which(!nas0))) {
          df2[["ms1_int"]][[j]][i0] <- ys[ps0[!is.na(ps0)]]
        }
        
        if (length(i1 <- which(!is.na(ps1)))) {
          df2[["ms1_int"]][[j]][i1] <- ys[ps1[!is.na(ps1)]]
        }
        
        if (length(i2 <- which(!is.na(ps2)))) {
          df2[["ms1_int"]][[j]][i2] <- ys[ps2[!is.na(ps2)]]
        }
      }
      else {
        df2[["ms1_int"]][[j]] <- ys[ps0]
      }
    }
    
    df[["ms1_int"]][ms2sta:ms2end] <- df2[["ms1_int"]]
  }
  
  df
}


#' Updates MS1 intensity with apex values.
#'
#' Including forward and backward looking.
#'
#' @param df A data frame.
#' @param matx The matrix of moverzs Y: by masses; X: by LC scans.
#' @param maty The matrix of intensities. Y: by masses; X: by LC scans.
#' @param row_sta The starting row of \code{matx}.
#' @param row_end The ending row of \code{matx}.
#' @param scan_apexs The vectors of apex scan number. Each vector corresponds to
#'   a mass in \code{matx}.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @importFrom fastmatch %fin%
updateMS1Int2 <- function (df, matx, maty, row_sta, row_end, scan_apexs, 
                           from = 200L, step = 1E-5)
  
{
  # to update the apexs, vector since can be chimeric precursors
  df$apex_scan_num <- vector("list", nrow(df))
  
  nrow <- nrow(matx)
  unv <- as.integer(colnames(matx))
  ss <- as.integer(rownames(matx))
  
  pos_levs <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")
  gap <- row_sta - 1L
  gap2 <- 12L

  ###
  # rowi: the current row index along matx
  # scan: the current scan number (of rowi); 13588
  ###
  
  for (i in seq_along(ms2_stas)) { # the same as by ms1_stas
    # i = 527 - gap; which(rownames(matx) == 13588)
    # i = 100;
    ms2sta <- ms2_stas[[i]]
    if (is.na(ms2sta)) next
    ms2end <- ms2_ends[[i]]
    df2 <- df[ms2sta:ms2end, ]
    
    rowi <- gap + i
    scan <- ss[[rowi]]
    # rows <- max(1L, rowi - gap):min(rowi + gap2, nrow)

    for (j in 1:nrow(df2)) {
      # j = 8
      x2s <- df2[["ms1_moverz"]][[j]] # precursor masses associated with an MS2 scan
      nx <- length(x2s) # nx > 1 with a chimeric spectrum
      if (!nx) next
      ix2s <- as.integer(ceiling(log(x2s/from)/log(1+step)))
      ks <- lapply(ix2s, function (x) which(abs(x - unv) <= 1L))
      
      for (m in 1:nx) {
        k <- ks[[m]] # the k-th column
        
        if (!length(k)) {
          next
        }

        k <- k[[1]] # can have two adjacent matches
        
        # ix <- unv[[k]]; # a <- maty[, which(colnames(matx) == unv[[k]])]; a <- maty[, 3218:3219]
        # apex scan numbers; note: scan_apexs not ordered: one-hit-wonders first
        apexs <- scan_apexs[[k]] # all apexs (scan numbers) under the k-th mass column
        ds <- abs(apexs - scan)
        p1 <- .Internal(which.min(ds))
        ap1 <- apexs[[p1]] # the nearest apex scan; 13525
        ok1 <- .Internal(which(ss == ap1))
        y1 <- maty[ok1, k]

        # checks neighbors
        if (length(apexs) > 1L) {
          ps <- which_topx2(-ds, 2L)
          p2 <- ps[ps != p1]
          ap2 <- apexs[[p2]] # the second nearest apex scan; 13686
          ok2 <- .Internal(which(ss == ap2))
          y2 <- maty[ok2, k]

          # later checks peak spacing and peak width (for small satellite peaks)
          if ((y1 < y2) && (abs(ap2 - ap1) <= 200L)) { # 200L somewhat arbitrary
            y1 <- y2
            df2$apex_scan_num[[j]][[m]] <- ap2 # ss[ok2]
          }
          else {
            df2$apex_scan_num[[j]][[m]] <- ap1 # ss[ok1]
          }
        }

        df2[["ms1_int"]][[j]][m] <- max(df2[["ms1_int"]][[j]][m], y1)
        # df2[["ms1_int"]][[j]][m] <- y1
      }
    }

    ms2rng <- ms2sta:ms2end
    df[["apex_scan_num"]][ms2rng] <- df2[["apex_scan_num"]]
    df[["ms1_int"]][ms2rng] <- df2[["ms1_int"]]
  }
  
  df
}


#' Helper of MS1 tracing.
#' 
#' Not yet used.
#' 
#' @param df MS1 data.
#' @param step Step size.
#' @param n_ms1peakpicking_flanks The number of flanking MS1 scans for precursor
#'   peak picking.
#' @param replace_ms1_by_apex Logical; if TRUE, fill all entries within a gate
#'   by its apex values.
#' @param filename A peaklist filename.
#' @param temp_dir A temp_dir to the filename.
traceMS1 <- function (df, min_mass = 200L, step = 8E-6, 
                      n_ms1peakpicking_flanks = 4L, replace_ms1_by_apex = TRUE, 
                      filename = NULL, temp_dir = NULL)
  
{
  rows <- which(lengths(df$ms1_moverz) == 0L)
  
  if (length(rows)) {
    df$ms1_int[rows] <- df$ms1_mass[rows] <- df$ms1_moverz[rows] <- list(NA_real_)
    df$ms1_charge[rows] <- list(NA_integer_)
  }
  
  for (i in nrow(df)) {
    xi <- df[["ms1_moverz"]][[i]]
    yi <- df[["ms1_int"]][[i]]
    zi <- df[["ms1_charge"]][[i]]
    mi <- df[["ms1_mass"]][[i]]
    oki <- .Internal(which(yi > 0))
    
    if (length(oki)) {
      df[["ms1_moverz"]][[i]] <- xi[oki]
      df[["ms1_int"]][[i]] <- yi[oki]
      df[["ms1_charge"]][[i]] <- zi[oki]
      df[["ms1_mass"]][[i]] <- mi[oki]
    }
  }
  rm(list = c("xi", "yi", "zi", "mi", "oki"))
  
  mat <- traceLCMS(
    xs = df[["ms1_moverz"]], 
    ys = df[["ms1_int"]], 
    zs = df[["ms1_charge"]], 
    n_dia_scans = n_ms1peakpicking_flanks, 
    from = min_mass, 
    step = step, 
    reord = FALSE, # already ordered
    cleanup = FALSE, # already cleaned
    replace_ms1_by_apex = replace_ms1_by_apex, 
    direct_out = TRUE, 
    temp_dir = temp_dir)
}


#' Gets MS1 intensity values from full MS1 data.
#' 
#' Not yet used.
#' 
#' @param df1 Data frame corresponding to \code{ms_level == 1L}.
#' @param from The starting point for mass binning
#' @param step A step size.
#' @param set_missing_zero Logical; if TRUE, set 0-intensity for peaks not
#'   found.
getMS1Int <- function (df1, from = 200L, step = 8E-6, set_missing_zero = FALSE)
{
  for (i in 1:nrow(df1)) {
    xs1 <- df1$ms1_moverz[[i]]
    len <- length(xs1)
    
    if (!len)
      next
    
    ys <- df1$msx_ints[[i]]
    ixxs <- as.integer(ceiling(log(df1$msx_moverzs[[i]]/from)/log(1+step)))
    ixs1 <- as.integer(ceiling(log(xs1/from)/log(1+step)))
    ps0 <- match(ixs1, ixxs)
    oks <- !is.na(ps0)
    
    if (all(oks)) {
      df1[["ms1_int"]][[i]] <- ys[ps0]
    }
    else {
      ps1 <- match(ixs1 + 1L, ixxs)
      ps2 <- match(ixs1 - 1L, ixxs)
      i1 <- .Internal(which(!is.na(ps1)))
      i2 <- .Internal(which(!is.na(ps2)))
      i0 <- .Internal(which(oks))
      
      if (length(i0)) {
        df1[["ms1_int"]][[i]][i0] <- ys[ps0[i0]]
      }
      
      if (length(i1)) {
        df1[["ms1_int"]][[i]][i1] <- ys[ps1[i1]]
      }
      
      if (length(i2)) {
        df1[["ms1_int"]][[i]][i2] <- ys[ps2[i2]]
      }
      
      # not matched at all
      if (set_missing_zero) {
        i012 <- c(i0, i1, i2)
        
        if (length(i012) < len) {
          bads <- !1:len %in% i012
          df1[["ms1_int"]][[i]][bads] <- 0
        }
      }
    }
  }
  
  df1
}


