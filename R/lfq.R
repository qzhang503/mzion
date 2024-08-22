#' Subsets full MS1 by the universe of monoisotopic MS1 moverzs.
#'
#' Non-monoisotopic (13C_n) MS1-X are also removed.
#'
#' @param xs Vectors of MS1-X values.
#' @param ys Vectors of MS1-Y values.
#' @param ms1s Vectors of MS1 m/z values that are linked to the mono-isotopic
#'   DDA-MS2 space.
#' @param from The starting m/z value for calculating bin indexes.
#' @param step The size of bins.
#' @param gap The gap of MS1 scans.
#' @importFrom fastmatch %fin%
subMSfull <- function (xs, ys, ms1s, from = 200L, step = 1E-5, gap = 256L)
{
  lens <- lengths(xs)
  len  <- length(xs)
  yout <- xout <- vector("list", len)

  # for each xi, keeps entries found in ms1s.
  if (gap >= len) {
    ms1s <- .Internal(unlist(ms1s, recursive = FALSE, use.names = FALSE))
    ms1s <- index_mz(ms1s, from, step)
    ms1s <- sort(unique(ms1s))

    for (i in 1:len) {
      xi <- xs[[i]]
      li <- lens[[i]]
      if (!li) next
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
      if (!li) next
      sta <- max(1, i - gap)
      end <- min(i + gap, len)
      # the possible universe of mono-isotopic MS1 m/z values over a range
      ms1s_sub <- .Internal(unlist(ms1s[sta:end], recursive = FALSE, 
                                   use.names = FALSE))
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
#' @return dfs[[i]]:a chunk of MS1 and MS2 data (no bracketing scans); 
#' df1s[[i]]: the reference MS1 data with leading and trailing scans; 
#' gaps[[i]]: the number of bracketing MS1 and MS2 scans between chunks.
pretraceXY <- function (df, from = 200L, step = 1e-5, n_chunks = 4L, gap = 256L)
{
  # msx_moverzs at ms_level == 1L: full-spectrum MS1
  # msx_charges at ms_level == 1L are list(NULL)
  cols1 <- c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", 
             "msx_moverzs", "msx_ints", "msx_charges", 
             "orig_scan", "ret_time")
  df1  <- df[with(df, ms_level == 1L), cols1]
  len1 <- nrow(df1)

  # Remove non-essential (e.g., non mono-isotopic) MS1 x and y values
  ans <- subMSfull(
    xs = df1$msx_moverzs, ys = df1$msx_ints, ms1s = df1$ms1_moverz, 
    from = from, step = step, gap = gap)
  # moverzs in ans$x are in an ascending order
  df1$msx_moverzs <- ans$x
  df1$msx_ints <- ans$y
  rm(list = "ans")

  # at least two chunks
  if (n_chunks <= 1L) {
    n_chunks <- 2L
  }

  df1$ms1_mass <- df1$ms1_moverz <- df1$ms1_int <- df1$ms1_charge <- NULL
  df1s  <- chunksplit(df1, n_chunks, type = "row")
  min_rts <- unlist(lapply(df1s, function (x) x$ret_time[[1]]))
  max_rts <- unlist(lapply(df1s, function (x) x$ret_time[[nrow(x)]]))
  end1s <- cumsum(lapply(df1s, nrow))
  sta1s <- c(1L, end1s[1:(n_chunks - 1L)] + 1L)
  
  # Adds 2-min gaps before and after
  gaps <- lapply(df1s, function (x) ceiling(min(gap, nrow(x) / 2L)))
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
    for (i in 2:(n_chunks-1L)) {
      df1s[[i]] <- dplyr::bind_rows(df1s_af[[i-1]], df1s[[i]], df1s_bf[[i+1]])
    }
  }
  rm(list = c("df1s_bf", "df1s_af", "df1"))
  
  ##  Splits `df` with bracketing entries
  # values: row indexes in `df`, length == nrow(df1) at pad_nas = TRUE
  ms1_stas <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)$ms1_stas
  cols <- c("ms_level", "ms1_moverz", "ms1_int", "orig_scan") # orig_scan for troubleshooting
  ms1_stax <- vector("integer", n_chunks)
  dfs <- vector("list", n_chunks)
  
  for (i in 1:n_chunks) {
    stai <- sta1s[[i]] # stai-th MS1 scan
    ms1_stax[[i]] <- ms1_stas[stai] # the corresponding row index in df
  }

  for (i in 1:(n_chunks - 1L)) {
    rowx <- ms1_stax[[i]]:(ms1_stax[[i+1]] - 1L) # ms2_endx can be NA
    dfs[[i]] <- df[rowx, cols] # both df1 and df2 data
  }
  dfs[[n_chunks]] <- df[ms1_stax[[n_chunks]]:nrow(df), cols]

  # dfs[[i]]:  a chunk MS1 and MS2 data (no bracketing scans)
  # df1s[[i]]: the reference MS1 data + leading and trailing MS1 scans; 
  # gaps[[i]]: the number of bracketing MS1 and MS2 scans
  list(dfs = dfs, df1s = df1s, gaps = gaps, min_rts = min_rts, max_rts = max_rts)
}


#' Helper of \link{traceXY}.
#'
#' @param xs Vectors of full-spectrum MS1 m/z values.
#' @param ys Vectors of full-spectrum MS1 intensities.
#' @param ss Vectors of MS1 (original) scan numbers.
#' @param ts Vectors of MS1 retention times (for calculating area-under-a-peak).
#' @param df A data frame of MS1 and MS2 corresponding to \code{xs}.
#' @param gap_bf A preceding gap size.
#' @param gap_af A following gap size.
#' @param from The starting point for mass binning.
#' @param step The step size for mass binning.
#' @param yco The cut-off in y values.
#' @param y_perc The cut-off in intensity values in relative to the base peak.
#' @param look_back Logical; look up the preceding MS bin or not.
#' @param min_y The cut-off of intensity values.
#' @inheritParams matchMS
htraceXY <- function (xs, ys, ss, ts, df, gap_bf = 256L, gap_af = 256L, 
                      n_mdda_flanks = 6L, from = 200L, step = 5E-6, 
                      y_perc = .01, yco = 500, look_back = TRUE, min_y = 0)
{
  if (all(lengths(xs) == 0L)) {
    # length(xs) - gap_bf - gap_af - nrow(df) == 0L
    null <- rep_len(list(NULL), nrow(df))
    
    return(
      tibble::tibble(
        ms_level = df$ms_level, 
        ms1_moverz = null, 
        ms1_int = null, 
        apex_scan_num = null))
  }
  
  mat <- traceXY(
    xs = xs, ys = ys, ss = ss, ts = ts, n_mdda_flanks = n_mdda_flanks, 
    from = from, step = step, reord = FALSE, cleanup = FALSE, # otherwise rows drop
    replace_ms1_by_apex = TRUE, y_perc = y_perc, yco = yco, 
    look_back = look_back)
  
  ## apes, rngs and scans: each vector corresponds to a column of mass
  matx <- mat[["x"]]
  maty <- mat[["y"]]
  apes <- mat[["p"]]
  rngs <- mat[["range"]] # for debugging
  nr <- nrow(matx)
  nc <- ncol(matx)

  # apes correspond to the row numbers of matx and ss to orig_scan
  rt_apexs <- scan_apexs <- vector("list", nc)
  for (i in 1:nc) {
    ai <- apes[[i]]
    scan_apexs[[i]] <- as.integer(ss[ai])
    rt_apexs[[i]]   <- as.integer(ts[ai])
  }

  # i = 6921
  # apes[[i]][[24]]
  # ss[apes[[i]]][[24]] # scan 114404
  # matx[rngs[[i]][[24]], i]
  # maty[rngs[[i]][[24]], i]
  
  # i = 5053; scan_apexs[[i]]
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
  
  # look both before and after scans
  df <- updateMS1Int2(
    df = df, matx = matx, maty = maty, row_sta = sta, row_end = end, 
    scan_apexs = scan_apexs, rt_apexs = rt_apexs, rngs = rngs, 
    from = from, step = step, min_y = min_y)
}


#' Helper of MS1 tracing.
#'
#' @param xs Lists of full-spectrum MS1 moverzs vectors.
#' @param ys Lists of full-spectrum MS1 intensities vectors.
#' @param ss Vectors of MS1 scan numbers.
#' @param ts Vectors of MS1 retention times (for calculating area-under-a-peak).
#' @param step Step size.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @param reord Logical; re-order data or not.
#' @param cleanup Logical; to clean up xs, ys and zs or not. Set the value to
#'   FALSE to maintain one-to-one correspondence between input (data frame) and
#'   the outputs. This will help, e.g., keep track of scan numbers in the input.
#' @param replace_ms1_by_apex Logical; if TRUE, fill all entries within a gate
#'   by its apex values.
#' @param yco The cut-off in y values. An arbitrary small number. The MS1 peaks
#'   with high-resolution Thermo's instruments are not a continuum and the
#'   maximum number of MS1 peaks seem \eqn{\le 1500}.
#' @param y_perc The cut-off in intensity values in relative to the base peak.
#' @param look_back Logical; look up the preceding MS bin or not.
#' @inheritParams matchMS
traceXY <- function (xs, ys, ss, ts, n_mdda_flanks = 6L, from = 200L, 
                     step = 8E-6, reord = TRUE, cleanup = FALSE, 
                     replace_ms1_by_apex = FALSE,y_perc = .01, yco = 500, 
                     look_back = TRUE)
{
  lens <- lengths(xs)
  
  if (all(lens == 0L)) {
    return(list(x = NULL, y = NULL, n = NULL, p = NULL, range = NULL))
  }

  if (reord) {
    for (i in seq_along(xs)) {
      xi <- xs[[i]]
      
      if (lens[[i]] > 1L) {
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
  
  # may change to step = 7E-6
  # max(mass_delta) under a column can be up to 3 * step with look_back = TRUE
  
  ans <- collapse_mms1ints(
    xs = xs, ys = ys, lwr = from, step = step, reord = FALSE, cleanup = FALSE, 
    add_colnames = TRUE, look_back = look_back)
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  nr <- nrow(ansy)
  nc <- ncol(ansy)
  rm(list = c("ans"))
  
  if (FALSE) {
    rng <- 250:800
    i <- which(colnames(ansx) == index_mz(500.908, from, step) +  1) # 2761
    i <- which(colnames(ansx) == index_mz(500.908, from, step) - 1) # 2760
    ss[rng]
    
    plot(ansx[rng, i])
    plot(ansy[rng, i])
    
    plot(ansx[, i])
    plot(ansy[, i])
    (max(ansx[, i]) - min(ansx[, i])) / min(ansx[, i])
    summary(ansx[, i])
  }
  
  if (nr != length(xs)) {
    stop("Developer: check for row dropping during MS1 tracing.")
  }
  
  ## traces MS1 data matrices across LC scans; rows: scans; columns: masses
  xmat <- ymat <- matrix(rep_len(NA_real_, nc * nr), ncol = nc)
  colnames(xmat) <- colnames(ymat) <- colnames(ansx) # bin indexes of masses
  rownames(xmat) <- rownames(ymat) <- ss # scan numbers
  ranges <- apexes <- ns <- vector("list", nc)
  
  if (replace_ms1_by_apex) {
    for (i in 1:nc) {
      # i <- which(colnames(ansx) == index_mz(677.3303, from, step) + 1) # 7938
      # i <- 10083
      xi <- ansx[, i]
      yi <- ansy[, i]
      oks <- .Internal(which(!is.na(yi)))
      yoks <- yi[oks]
      # may be unnecessary, e.g., Thermo's MS1 peak distributions are discrete
      yoks[yoks < yco] <- NA_real_
      yi[oks] <- yoks
      # plot(yi[466:781])
      # plot(xi[466:781])
      
      # ss[466:781]
      gates <- 
        find_lc_gates(xs = xi, ys = yi, ts = ts, n_dia_scans = n_mdda_flanks)
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]
      yints <- gates[["yints"]]
      
      mx <- lapply(rngs, function (rng) mean(xi[rng], na.rm = TRUE))

      for (j in seq_along(rows)) {
        rgj <- rngs[[j]]
        # rwj <- rows[[j]]
        # xmat[rgj, i] <- ansx[rwj, i]
        xmat[rgj, i] <- mx[[j]]
        ymat[rgj, i] <- yints[[j]]
      }
    }
  }
  else {
    for (i in 1:nc) {
      gates <- find_lc_gates(ys = ansy[, i], ts = ts, n_dia_scans = n_mdda_flanks)
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]
      xmat[rows, i] <- ansx[rows, i]
      ymat[rows, i] <- ansy[rows, i]
    }
  }
  
  list(x = xmat, y = ymat, n = ns, p = apexes, range = ranges)
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
#' @param rt_apexs The vectors of apex retention time. Each vector corresponds
#'   to a mass in \code{matx}.
#' @param rngs The ranges of row indexes for each column in \code{matx}. The
#'   \code{rngs[[i]]} corresponds to \code{matx[, i]}. The argument is currently
#'   used for debugging.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @param min_y The cut-off of intensity values.
#' @importFrom fastmatch %fin%
updateMS1Int2 <- function (df, matx, maty, row_sta, row_end, scan_apexs, 
                           rt_apexs, rngs, from = 200L, step = 5E-6, min_y = 0)
  
{
  # to update the apexs, vector since can be chimeric precursors
  df$apex_scan_num <- vector("list", nrow(df))
  
  nrow <- nrow(matx)
  unv  <- as.integer(colnames(matx))
  ss   <- as.integer(rownames(matx))
  
  pos_levs <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")
  gap <- row_sta - 1L

  ###
  # rowi: the current row index along matx
  # scan: the current scan number (of rowi); 13588
  ###
  
  for (i in seq_along(ms2_stas)) { # the same as by ms1_stas
    if (FALSE) {
      df$orig_scan[ms1_stas] # look for orig_scan around 34009 ->
      i <- 190
      ms2sta <- ms2_stas[[i]]
      ms2end <- ms2_ends[[i]]
      df2 <- df[ms2sta:ms2end, ]
    }
    
    # i = 232; which(rownames(matx) == 18949) - gap -> i
    ms2sta <- ms2_stas[[i]]
    if (is.na(ms2sta)) next
    ms2end <- ms2_ends[[i]]
    df2 <- df[ms2sta:ms2end, ]
    
    rowi <- gap + i # row number in matx, maty
    scan <- ss[[rowi]] # the MS1 scan number at the current row
    
    for (j in 1:nrow(df2)) {
      # j <- 7
      x1s <- df2[["ms1_moverz"]][[j]] # MS1 masses associated with an MS2 scan
      nx <- length(x1s) # nx > 1 with a chimeric spectrum
      if (!nx) next
      ix1s <- as.integer(ceiling(log(x1s/from)/log(1+step)))
      ks <- lapply(ix1s, function (x) which(abs(x - unv) <= 1L))
      
      # go through chimeric precursors
      for (m in 1:nx) {
        k <- ks[[m]] # the k-th column
        if (!length(k)) next
        x1m <- x1s[[m]]

        # two adjacent matches: unv[k] <-> ix1s[[m]] - 1 and ix1s[[m]] + 1
        if (length(k) > 1L) {
          ka <- k[[1]]
          ans1 <- find_apex_scan(k = ka, xs_k = matx[, ka], ys_k = maty[, ka], 
                                 apexs_k = scan_apexs[[ka]], 
                                 rts_k = rt_apexs[[ka]], rngs_k = rngs[[ka]], 
                                 xm = x1m, scan = scan, ss = ss, step = step, 
                                 min_y = min_y)
          ap1a <- ans1$ap
          x1a  <- ans1$x
          y1a  <- ans1$y
          fra  <- ans1$from
          toa  <- ans1$to

          kb <- k[[2]]
          ans2 <- find_apex_scan(k = kb, xs_k = matx[, kb], ys_k = maty[, kb], 
                                 apexs_k = scan_apexs[[kb]], 
                                 rts_k = rt_apexs[[kb]], rngs_k = rngs[[kb]], 
                                 xm = x1m, scan = scan, ss = ss, step = step, 
                                 min_y = min_y)
          ap1b <- ans2$ap
          x1b  <- ans2$x
          y1b  <- ans2$y
          frb <- ans2$from
          tob <- ans2$to
          
          # (1) by mass errors
          erra <- x1a / x1m - 1
          errb <- x1b / x1m - 1
          # big_err <- if (abs(erra - errb) > 3e-6) TRUE else FALSE
          erra <- abs(erra)
          errb <- abs(errb)
          
          if (erra < errb) {
            is_a <- TRUE
            errx <- errb - erra
            err0 <- erra
          }
          else {
            is_a <- FALSE
            errx <- erra - errb
            err0 <- errb
          }
          
          # (2) by inclusiveness
          # can also often be discriminated by Y but not always as
          #  a large Y may be due to a irregular fat peak
          
          # `a` is mostly a subset of `b`
          if (fra + 2L >= frb && toa <= tob + 2L) {
            b_incl_a <- TRUE
          }
          else {
            b_incl_a <- FALSE
          }
          
          # `b` is mostly a subset of `a`
          if (frb + 2L >= fra && tob <= toa + 2L) {
            a_incl_b <- TRUE
          }
          else {
            a_incl_b <- FALSE
          }

          # (2.x) by scan number differences
          if (FALSE) {
            dsa <- abs(ap1a - scan)
            dsb <- abs(ap1b - scan)
            
            if (dsa < dsb) {
              ds1 <- dsb
              ds0 <- dsa
              is_da <- TRUE
            }
            else {
              ds1 <- dsa
              ds0 <- dsb
              is_da <- FALSE
            }
          }

          # (3) if `scan` is within any of the two peak ranges
          #     a scan cannot be both `oka` and `okb`
          # oka <- scan > fra + 1L && scan < toa - 1L
          # okb <- scan > frb + 1L && scan < tob - 1L
          
          if (errx > 3e-6) {
            if (is_a) {
              ap1 <- ap1a
              y1  <- y1a
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
            }
          }
          else if (a_incl_b) {
            ap1 <- ap1a
            y1  <- y1a
          }
          else if (b_incl_a) {
            ap1 <- ap1b
            y1  <- y1b
          }
          else if (FALSE && ds0 * 10 <= ds1) { # 10x closer ds0 <= 100 && 
            if (is_da) {
              ap1 <- ap1a
              y1  <- y1a
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
            }
          }
          else {
            if (y1a > y1b) {
              ap1 <- ap1a
              y1  <- y1a
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
            }
          }
        }
        else {
          ka <- k[[1]]
          ans1 <- find_apex_scan(k = ka, xs_k = matx[, ka], ys_k = maty[, ka], 
                                 apexs_k = scan_apexs[[ka]], 
                                 rts_k = rt_apexs[[ka]], rngs_k = rngs[[ka]], 
                                 xm = x1m, scan = scan, ss = ss, step = step, 
                                 min_y = min_y)
          ap1 <- ans1$ap
          y1  <- ans1$y
        }
        
        df2$apex_scan_num[[j]][[m]] <- ap1
        df2[["ms1_int"]][[j]][m] <- y1
      }
    }

    ms2rng <- ms2sta:ms2end
    df[["apex_scan_num"]][ms2rng] <- df2[["apex_scan_num"]]
    df[["ms1_int"]][ms2rng] <- df2[["ms1_int"]]
  }
  
  rows <- lengths(df$apex_scan_num) > 0L
  df$apex_scan_num[rows] <- 
    lapply(df$apex_scan_num[rows], unlist, recursive = FALSE, use.names = FALSE)
  
  df
}


#' Find the apex under a mass column
#'
#' @param k The k-th column in \code{matx} or \code{maty}.
#' @param xs_k The X values under \code{matx[, k]}.
#' @param ys_k The Y values under \code{maty[, k]}.
#' @param apexs_k All apexes scan numbers (\code{scan_apexs[[k]]}) under column
#'   k.
#' @param rts_k All apexes retention times (\code{scan_apexs[[k]]}) under column
#'   k.
#' @param rngs_k All apexes scan ranges (\code{scan_apexs[[k]]}) under column k.
#' @param xm The experimental X value for tracing.
#' @param scan The scan number of an MS2 event corresponding to the \code{xm}.
#' @param ss All of the scan numbers.
#' @param rngs The scan ranges of apexes under column k.
#' @param step A step size.
#' @param min_y The cut-off of intensity values.
#' @importFrom fastmatch %fin%
find_apex_scan <- function (k, xs_k, ys_k, apexs_k, rts_k, rngs_k, xm, scan, ss, 
                            step = 6E-6, min_y = 0)
{
  # (1) subset apexs by MS1 mass tolerance
  ok_xs  <- .Internal(which(abs(xs_k - xm) / xm <= step))
  ok_aps <- .Internal(which(apexs_k %fin% ss[ok_xs]))

  if (length(ok_aps)) { # && length(ok_xs)
    aps <- apexs_k[ok_aps]
    rtx <- rts_k[ok_aps]
    rgx <- rngs_k[ok_aps]
  }
  else {
    aps <- apexs_k
    rtx <- rts_k
    rgx <- rngs_k
  }
  # abs(xs_k[names(xs_k) %in% aps] - xm) / xm
  
  # (2) remove one-hit-wonders and spikes
  lens <- lengths(rgx)
  oks1 <- .Internal(which(lens > 10L))
  oks2 <- .Internal(which(lens > 5L))
  
  if (length(oks1)) {
    rtx  <- rtx[oks1]
    rgx  <- rgx[oks1]
    aps  <- aps[oks1]
    lens <- lens[oks1]
  }
  else if (length(oks2)) {
    rtx  <- rtx[oks2]
    rgx  <- rgx[oks2]
    aps  <- aps[oks2]
    lens <- lens[oks2]
  }
  
  # (3) subset apexes by distances between apex_scan and the triggering MS2 scan
  ds <- abs(aps - scan)
  p1 <- .Internal(which.min(ds))
  d1 <- ds[[p1]]
  px <- max(1L, p1 - 3L):min(length(ds), p1 + 3L)
  ds <- ds[px]

  # at least one (the most centered) peak is within 2500 scans
  if (d1 <= 2500) {
    oks_d <- ds <= 2500
    px <- px[oks_d]
    ds <- ds[oks_d]
  }
  
  aps  <- aps[px]
  rtx  <- rtx[px]
  rgx  <- rgx[px]
  lens <- lens[px]

  # by intensity cut-off
  oks <- .Internal(which(ss %fin% aps))
  xs  <- xs_k[oks]
  ys  <- ys_k[oks]
  
  ok_ys  <- .Internal(which(ys >= min_y))
  len_ys <- length(ok_ys)
  
  if (len_ys && len_ys < length(ys)) {
    xs <- xs[ok_ys]
    ys <- ys[ok_ys]
    aps  <- aps[ok_ys]
    ds   <- ds[ok_ys]
    rgx  <- rgx[ok_ys]
    lens <- lens[ok_ys]
  }
  
  if (len_ys == 1) {
    ssx <- ss[rgx[[1]]]
    return(list(x = xs[[1]], y = ys[[1]], ap = aps[[1]], from = ssx[[1]], to = ssx[[lens]]))
  }

  # the closest may be a spike... 
  # look for more "regular" peak nearby...
  
  # ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, ds))
  # topa <- ds[[ord[[1]]]]
  topa <- .Internal(which.min(ds))
  ya   <- ys[[topa]]
  xa   <- xs[[topa]]
  apa  <- aps[[topa]]
  da   <- ds[[topa]]
  rga  <- rgx[[topa]]
  ssa  <- ss[rga]
  lena <- lens[[topa]]
  # oka  <- apa %in% ssa # must be

  topb <- .Internal(which.max(ys))
  yb   <- ys[[topb]]
  xb   <- xs[[topb]]
  apb  <- aps[[topb]]
  db   <- ds[[topb]]
  rgb  <- rgx[[topb]]
  ssb  <- ss[rgb]
  lenb <- lens[[topb]]
  # okb  <- apb %in% ssb # must be

  if (topa == topb) {
    return(list(x = xa, y = ya, ap = apa, from = ssa[[1]], to = ssa[[lena]]))
  }
  
  # compare peak widths and distances
  lenr <- lena / lenb
  
  if (lenr > .67 || lenr < 1.5) {
    if (db > da + 200) {
      return(list(x = xa, y = ya, ap = apa, from = ssa[[1]], to = ssa[[lena]]))
    }
    else {
      return(list(x = xb, y = yb, ap = apb, from = ssb[[1]], to = ssb[[lenb]]))
    }
  }
  else {
    # use the fatter peak (assume the leaner one is a spike)
    return(list(x = xb, y = yb, ap = apb, from = ssb[[1]], to = ssb[[lenb]]))
  }
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
#' @param min_mass A minimum mass.
traceMS1 <- function (df, min_mass = 200L, step = 8E-6, 
                      n_ms1peakpicking_flanks = 4L, replace_ms1_by_apex = TRUE, 
                      filename = NULL, temp_dir = NULL)
  
{
  rows <- which(lengths(df$ms1_moverz) == 0L)
  
  if (length(rows)) {
    df$ms1_int[rows] <- df$ms1_mass[rows] <- df$ms1_moverz[rows] <- 
      list(NA_real_)
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
    
    if (!len) {
      next
    }

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


