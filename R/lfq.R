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
      
      ix  <- index_mz(xi, from, step)
      ps0 <- fastmatch::fmatch(ix, ms1s)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s)
      i0  <- .Internal(which(!is.na(ps0)))
      i1  <- .Internal(which(!is.na(ps1)))
      i2  <- .Internal(which(!is.na(ps2)))
      
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
      
      ix  <- index_mz(xi, from, step)
      ps0 <- fastmatch::fmatch(ix, ms1s_sub)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s_sub)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s_sub)
      i0  <- .Internal(which(!is.na(ps0)))
      i1  <- .Internal(which(!is.na(ps1)))
      i2  <- .Internal(which(!is.na(ps2)))
      
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
  df1   <- df[with(df, ms_level == 1L), cols1]
  len1  <- nrow(df1)

  # Remove non-essential (e.g., non mono-isotopic) MS1 x and y values
  ans <- subMSfull(
    xs = df1$msx_moverzs, ys = df1$msx_ints, ms1s = df1$ms1_moverz, 
    from = from, step = step, gap = gap)
  # moverzs in ans$x are in an ascending order
  df1$msx_moverzs <- ans$x
  df1$msx_ints    <- ans$y
  rm(list = "ans")

  # at least two chunks
  if (n_chunks <= 1L) {
    n_chunks <- 2L
  }

  df1$ms1_mass <- df1$ms1_moverz <- df1$ms1_int <- df1$ms1_charge <- NULL
  df1s    <- chunksplit(df1, n_chunks, type = "row")
  min_rts <- unlist(lapply(df1s, function (x) x$ret_time[[1]]))
  max_rts <- unlist(lapply(df1s, function (x) x$ret_time[[nrow(x)]]))
  end1s   <- cumsum(lapply(df1s, nrow))
  sta1s   <- c(1L, end1s[1:(n_chunks - 1L)] + 1L)
  
  # Adds 2-min gaps before and after
  gaps    <- lapply(df1s, function (x) ceiling(min(gap, nrow(x) / 2L)))
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
#' @param path_ms1 The path for saving the output of MS1 data.
#' @param out_name A file name of output.
#' @param out_cols The output column keys.
#' @inheritParams matchMS
htraceXY <- function (xs, ys, ss, ts, df, gap_bf = 256L, gap_af = 256L, 
                      out_name = NULL, n_dia_scans = 4L, from = 200L, 
                      step = 5E-6, y_perc = .01, yco = 500, look_back = TRUE, 
                      min_y = 0, path_ms1 = NULL, 
                      out_cols = c("ms_level", "ms1_moverz", "ms1_int", 
                                   "orig_scan", "apex_scan_num", 
                                   "apex_xs1", "apex_ys1", "apex_bin1", 
                                   "apex_ps1", "apex_ts1", 
                                   "apex_xs2", "apex_ys2", "apex_bin2", 
                                   "apex_ps2", "apex_ts2"))
{
  # length(xs) - gap_bf - gap_af - nrow(df) == 0L
  if (all(lengths(xs) == 0L)) {
    null <- rep_len(list(NULL), nrow(df))
    
    # not saving file `out_name`
    df <- tibble::tibble(
      ms_level = df$ms_level, 
      ms1_moverz = null, 
      ms1_int = null, 
      orig_scan = df$orig_scan, 
      # orig_scan = null, # don't should be `character` not `list`
      apex_scan_num = null, 
      apex_xs1 = null, 
      apex_ys1 = null, 
      apex_bin1 = null,
      apex_ps1 = null, 
      apex_ts1 = null,
      apex_xs2 = null, 
      apex_ys2 = null, 
      apex_bin2 = null,
      apex_ps2 = null, 
      apex_ts2 = null)
    
    return(df[, out_cols])
  }
  
  # may contain all NA columns
  mat <- traceXY(
    xs = xs, ys = ys, ss = ss, ts = ts, n_dia_scans = n_dia_scans, 
    from = from, step = step, reord = FALSE, cleanup = FALSE, # otherwise rows drop
    replace_ms1_by_apex = TRUE, y_perc = y_perc, yco = yco, 
    look_back = look_back)
  
  ## apes, rngs and scans: each vector corresponds to a column of mass; 
  #  their indexes corresponding to the row number in matx
  matx <- mat[["x"]]
  maty <- mat[["y"]]
  apes <- mat[["p"]]
  rngs <- mat[["range"]] # for debugging
  ns   <- mat[["n"]]

  if (length(bads <- which(lengths(apes) == 0L))) {
    matx <- matx[, -bads, drop = FALSE]
    maty <- maty[, -bads, drop = FALSE]
    apes <- apes[-bads]
    rngs <- rngs[-bads]
    ns   <- ns[-bads]
  }
  
  # not saving args
  if (!length(apes)) {
    df <- init_lfq_df(df)
    return(df[, out_cols])
  }
  
  nr <- nrow(matx)
  nc <- ncol(matx)

  # ss == rownames(matx) <-> orig_scan: ss[apes[[i]]]; how about with timsTOF?
  # apx = apes[[i]] <-> a row number of matx; rownames(matx)[apx] -> ss[apx]
  rt_apexs <- scan_apexs <- vector("list", nc)
  for (i in 1:nc) {
    apx <- apes[[i]]
    scan_apexs[[i]] <- as.integer(ss[apx])
    rt_apexs[[i]]   <- ts[apx] # as.integer()
  }

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
  args <- list(
    df = df, matx = matx, maty = maty, row_sta = sta, row_end = end, 
    scan_apexs = scan_apexs, rt_apexs = rt_apexs, rngs = rngs, 
    from = from, step = step, min_y = min_y)

  df <- do.call(updateMS1Int2, args)
  
  attr(df, "row_sta") <- sta
  attr(df, "row_end") <- end
  attr(df, "from")    <- from
  attr(df, "step")    <- step
  attr(df, "min_y")   <- min_y
  
  df <- df[, out_cols]
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
traceXY <- function (xs, ys, ss, ts, n_dia_scans = 4L, from = 200L, 
                     step = 8E-6, reord = TRUE, cleanup = FALSE, 
                     replace_ms1_by_apex = FALSE, y_perc = .01, yco = 500, 
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
        ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, xi))
        xs[[i]] <- xi[ord]
        ys[[i]] <- ys[[i]][ord]
      }
    }
    # rm(list = c("xi", "ord"))
  }
  
  ## collapses MS data by the indexes of mass bins; 
  #  two matrix outputs; rows: scans; columns: masses or intensities
  
  # xs can be numeric(0)?
  # cleanup = FALSE; otherwise rows drop
  
  # may change to step = 7E-6
  # max(mass_delta) under a column can be up to 3 * step with look_back = TRUE
  
  ans  <- collapse_mms1ints(
    xs = xs, ys = ys, lwr = from, step = step, reord = FALSE, cleanup = FALSE, 
    add_colnames = TRUE, look_back = look_back)
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  unv  <- ans[["u"]]
  # colnames(ansx) <- colnames(ansy) <- unv
  nr   <- nrow(ansy)
  nc   <- ncol(ansy)
  # rm(list = c("ans"))
  
  if (nr != length(xs)) {
    stop("Developer: check for row dropping during MS1 tracing.")
  }
  
  ## traces MS1 data matrices across LC scans; rows: scans; columns: masses
  xmat <- ymat <- matrix(rep_len(NA_real_, nc * nr), ncol = nc)
  ranges <- apexes <- ns <- vector("list", nc)
  
  if (replace_ms1_by_apex) {
    for (i in 1:nc) {
      # i <- which(unv == index_mz(627.8159, from, step) + 1)
      xi   <- ansx[, i]
      yi   <- ansy[, i]
      oks  <- .Internal(which(!is.na(yi)))
      yoks <- yi[oks]
      # may be unnecessary, e.g., Thermo's MS1 peak distributions are discrete
      yoks[yoks < yco] <- NA_real_ # does the same to xi?
      yi[oks]   <- yoks

      ## all NA or NaN if all yi < yco...
      # if (sum(is.na(yi)) + sum(is.nan(yi)) == nr) {
      #   ns[[i]] <- 0L; apexes[[i]] <- 0L; ranges[[i]] <- 0L
      # }
      
      if (FALSE) {
        data.frame(x = ts, y = yi) |>
          ggplot2::ggplot() + 
          ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                                color = "gray", linewidth = .1)
      }
      
      gates <- 
        find_lc_gates(xs = xi, ys = yi, ts = ts, n_dia_scans = n_dia_scans)
      
      if (is.null(gates)) {
        next
      }
      
      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]]     <-         gates[["ns"]]
      ranges[[i]] <- rngs <- gates[["ranges"]]
      yints       <-         gates[["yints"]]
      
      mx <- lapply(rngs, function (rng) mean(xi[rng], na.rm = TRUE)) # not median
      
      for (j in seq_along(rows)) {
        rgj <- rngs[[j]]
        xmat[rgj, i] <- mx[[j]]
        ymat[rgj, i] <- yints[[j]]
      }
    }
  }
  else {
    for (i in 1:nc) {
      gates <- find_lc_gates(ys = ansy[, i], ts = ts, n_dia_scans = n_dia_scans)
      apexes[[i]]   <- rows <- gates[["apex"]]
      ns[[i]]       <- gates[["ns"]] # number of observing scans
      ranges[[i]]   <- rngs <- gates[["ranges"]]
      xmat[rows, i] <- ansx[rows, i]
      ymat[rows, i] <- ansy[rows, i]
    }
  }
  
  colnames(xmat) <- colnames(ymat) <- unv # bin indexes of masses
  rownames(xmat) <- rownames(ymat) <- ss  # scan numbers
  
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
  df   <- init_lfq_df(df)
  nrow <- nrow(matx)
  unv  <- as.integer(colnames(matx))
  ss   <- as.integer(rownames(matx))
  
  rownames(matx) <- rownames(maty) <- NULL
  
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
  
  for (i in seq_along(ms1_stas)) {
    if (FALSE) {
      # look for penultimate scan number of psmQ.txt::pep_scan_num
      df$orig_scan[ms1_stas]
      i <- 344
      ms2sta <- ms2_stas[[i]]
      ms2end <- ms2_ends[[i]]
      ms2rng <- ms2sta:ms2end
      df2    <- df[ms2rng, ]
    }

    # i = 458; which(rownames(matx) == 78900) - gap -> i
    ms2sta <- ms2_stas[[i]]
    if (is.na(ms2sta)) { next }
    ms2end <- ms2_ends[[i]]
    ms2rng <- ms2sta:ms2end
    df2    <- df[ms2rng, ]
    rowi   <- gap + i    # row number in matx, maty
    scan   <- ss[[rowi]] # the MS1 scan number at the current row
    xs2    <- df2[["ms1_moverz"]]
    
    # by MS2 spectra
    for (j in 1:nrow(df2)) {
      # j <- 3L
      x1s  <- xs2[[j]] # MS1 masses associated with an MS2 scan
      nx   <- length(x1s) # nx > 1 at a chimeric spectrum
      if (!nx) { next }
      ix1s <- index_mz(x1s, from, step)
      ks   <- lapply(ix1s, function (x) which(abs(x - unv) <= 1L))
      
      # by chimeric precursors
      for (m in 1:nx) {
        k   <- ks[[m]] # the k-th column(s) in matx
        nk  <- length(k)
        if (!nk) { next }
        x1m <- x1s[[m]]

        # up to two with an adjacent match of ix1s[[m]] - 1 or ix1s[[m]] + 1
        if (nk > 1L) {
          ka <- k[[1]]
          kb <- k[[2]]
          
          if (FALSE) {
            data.frame(x = ss[1:nrow], y = maty[, ka]) |>
              ggplot2::ggplot() + 
              ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                                    color = "gray", linewidth = .1)
          }
          
          ans1 <- find_apex_scan(
            k = ka, xs_k = matx[, ka], ys_k = maty[, ka], 
            apexs_k = scan_apexs[[ka]], rts_k = rt_apexs[[ka]], 
            rngs_k = rngs[[ka]], xm = x1m, scan = scan, ss = ss, 
            step = step, min_y = min_y)
          # a tentative best
          ap1a <- ans1[["ap"]]
          x1a  <- ans1[["x"]]
          y1a  <- ans1[["y"]]
          fra  <- ans1[["from"]]
          toa  <- ans1[["to"]]
          
          # all candidates (low-quality entries in scan_apexs[[ka]] removed)
          xsa  <- ans1[["xs"]]
          ysa  <- ans1[["ys"]]
          apsa <- ans1[["aps"]]

          ans2 <- find_apex_scan(
            k = kb, xs_k = matx[, kb], ys_k = maty[, kb], 
            apexs_k = scan_apexs[[kb]], rts_k = rt_apexs[[kb]], 
            rngs_k = rngs[[kb]], xm = x1m, scan = scan, ss = ss, 
            step = step, min_y = min_y)
          ap1b <- ans2[["ap"]]
          x1b  <- ans2[["x"]]
          y1b  <- ans2[["y"]]
          frb  <- ans2[["from"]]
          tob  <- ans2[["to"]]
          
          xsb  <- ans2[["xs"]]
          ysb  <- ans2[["ys"]]
          apsb <- ans2[["aps"]]
          
          # (1) by mass errors
          erra <- abs(x1a / x1m - 1)
          errb <- abs(x1b / x1m - 1)
          
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
          # can also often be discriminated by Y but not always 
          #  as a large Y may be due to a irregular fat peak
          
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

          # (3) if a `scan` is within any of the two peak ranges
          #     the scan cannot be both `oka` and `okb`
          # oka <- scan > fra + 1L && scan < toa - 1L
          # okb <- scan > frb + 1L && scan < tob - 1L

          if (errx > 3e-6) {
            if (is_a) {
              ap1 <- ap1a
              y1  <- y1a
              k0  <- ka
              
              aps <- apsa
              xs  <- xsa
              ys  <- ysa
              
              kx  <- kb
              apx <- apsb
              xx  <- xsb
              yx  <- ysb
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
              k0  <- kb
              
              aps <- apsb
              xs  <- xsb
              ys  <- ysb
              
              kx  <- ka
              apx <- apsa
              xx  <- xsa
              yx  <- ysa
            }
          }
          else if (a_incl_b) {
            ap1 <- ap1a
            y1  <- y1a
            k0  <- ka
            
            aps <- apsa
            xs  <- xsa
            ys  <- ysa
            
            kx  <- kb
            apx <- apsb
            xx  <- xsb
            yx  <- ysb
          }
          else if (b_incl_a) {
            ap1 <- ap1b
            y1  <- y1b
            k0  <- kb
            
            aps <- apsb
            xs  <- xsb
            ys  <- ysb
            
            kx  <- ka
            apx <- apsa
            xx  <- xsa
            yx  <- ysa
          }
          else if (FALSE && ds0 * 10 <= ds1) { # 10x closer ds0 <= 100 && 
            if (is_da) {
              ap1 <- ap1a
              y1  <- y1a
              k0  <- ka
              
              aps <- apsa
              xs  <- xsa
              ys  <- ysa
              
              kx  <- kb
              apx <- apsb
              xx  <- xsb
              yx  <- ysb
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
              k0  <- kb
              
              aps <- apsb
              xs  <- xsb
              ys  <- ysb
              
              kx  <- ka
              apx <- apsa
              xx  <- xsa
              yx  <- ysa
            }
          }
          else {
            if (y1a > y1b) {
              ap1 <- ap1a
              y1  <- y1a
              k0  <- ka
              
              aps <- apsa
              xs  <- xsa
              ys  <- ysa
              
              kx  <- kb
              apx <- apsb
              xx  <- xsb
              yx  <- ysb
            }
            else {
              ap1 <- ap1b
              y1  <- y1b
              k0  <- kb
              
              aps <- apsb
              xs  <- xsb
              ys  <- ysb
              
              kx  <- ka
              apx <- apsa
              xx  <- xsa
              yx  <- ysa
            }
          }
          
          apexs_k0 <- scan_apexs[[k0]] # aps: after clean-ups
          rts_k0   <- rt_apexs[[k0]]
          rngs_k0  <- rngs[[k0]]
          xs_k0    <- matx[, k0]
          ys_k0    <- maty[, k0]
          
          apexs_kx <- scan_apexs[[kx]]
          rts_kx   <- rt_apexs[[kx]]
        }
        else {
          k0 <- k[[1]]
          apexs_k0 <- scan_apexs[[k0]]
          rts_k0   <- rt_apexs[[k0]]
          rngs_k0  <- rngs[[k0]]
          xs_k0    <- matx[, k0]
          ys_k0    <- maty[, k0]
          
          ans1 <- find_apex_scan(
            k = k0, xs_k = xs_k0, ys_k = ys_k0, 
            apexs_k = apexs_k0, rts_k = rts_k0, 
            rngs_k = rngs_k0, xm = x1m, scan = scan, ss = ss, 
            step = step, min_y = min_y)
          ap1 <- ans1[["ap"]]
          y1  <- ans1[["y"]]
          aps <- ans1[["aps"]]
          xs  <- ans1[["xs"]]
          ys  <- ans1[["ys"]]
          
          kx  <- NA_integer_ # an indicator
        }
        
        # scalar
        df2[["apex_scan_num"]][[j]][[m]] <- ap1
        df2[["ms1_int"]][[j]][[m]]       <- y1
        df2[["apex_bin1"]][[j]][[m]]     <- k0
        
        # vector
        oks <- .Internal(which(apexs_k0 %fin% aps))
        df2[["apex_ps1"]][[j]][[m]] <- apexs_k0[oks] # the same as `aps`
        df2[["apex_ts1"]][[j]][[m]] <- rts_k0[oks]
        
        # list
        df2[["apex_xs1"]][[j]][[m]]   <- xs
        df2[["apex_ys1"]][[j]][[m]]   <- ys
        # df2[["rngs_all"]][[j]][[m]]  <- rngs_k0[oks]
        
        if (!is.na(kx)) {
          df2[["apex_bin2"]][[j]][[m]] <- kx
          
          # vector
          okx <- .Internal(which(apexs_kx %fin% apx))
          df2[["apex_ps2"]][[j]][[m]] <- apexs_kx[okx]
          df2[["apex_ts2"]][[j]][[m]] <- rts_kx[okx]
          
          # list
          df2[["apex_xs2"]][[j]][[m]]   <- xx
          df2[["apex_ys2"]][[j]][[m]]   <- yx
        }
      }
    }
    
    df[["apex_scan_num"]][ms2rng] <- df2[["apex_scan_num"]]
    df[["ms1_int"]][ms2rng]       <- df2[["ms1_int"]]
    df[["apex_bin1"]][ms2rng]     <- df2[["apex_bin1"]]
    df[["apex_ps1"]][ms2rng]      <- df2[["apex_ps1"]]
    df[["apex_ts1"]][ms2rng]      <- df2[["apex_ts1"]]
    df[["apex_xs1"]][ms2rng]      <- df2[["apex_xs1"]]
    df[["apex_ys1"]][ms2rng]      <- df2[["apex_ys1"]]
    df[["apex_bin2"]][ms2rng]     <- df2[["apex_bin2"]]
    df[["apex_ps2"]][ms2rng]      <- df2[["apex_ps2"]]
    df[["apex_ts2"]][ms2rng]      <- df2[["apex_ts2"]]
    df[["apex_xs2"]][ms2rng]      <- df2[["apex_xs2"]]
    df[["apex_ys2"]][ms2rng]      <- df2[["apex_ys2"]]
    # df[["rngs_all"]][ms2rng] <- df2[["rngs_all"]]
  }
  
  df
}


#' Helper to initiate LFQ output columns
#' 
#' @param df A data frame.
init_lfq_df <- function (df)
{
  nr <- nrow(df)
  
  df[["apex_bin2"]] <- df[["apex_ys2"]] <- df[["apex_xs2"]] <- 
    df[["apex_ts2"]] <- df[["apex_ps2"]] <- 
    df[["apex_bin1"]] <- df[["apex_ys1"]] <- df[["apex_xs1"]] <- 
    df[["apex_ts1"]] <- df[["apex_ps1"]] <- 
    df[["apex_scan_num"]] <- # df[["rngs_all"]] <- 
    vector("list", nr)
  
  # replicate the lengths by the number of chimeric precursors
  rows2 <- which(df$ms_level == 2L)
  lens2 <- lengths(df$ms1_int[rows2])
  
  df[["apex_scan_num"]][rows2] <- df[["apex_bin1"]][rows2] <- 
    df[["apex_bin2"]][rows2] <- lapply(lens2, function (x) rep_len(0L, x))
  
  df[["apex_ps1"]][rows2] <- df[["apex_ts1"]][rows2] <- # df[["rngs_all"]][rows2] <- 
    df[["apex_xs1"]][rows2] <- df[["apex_ys1"]][rows2] <- 
    df[["apex_ps2"]][rows2] <- df[["apex_ts2"]][rows2] <- 
    df[["apex_xs2"]][rows2] <- df[["apex_ys2"]][rows2] <- 
    lapply(lens2, function (x) rep_len(list(NULL), x))
  
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
  aps  <- apexs_k
  rtx  <- rts_k
  rgx  <- rngs_k
  lens <- lengths(rgx)
  
  # (1) first-pass spike removals: in case that a major peak is just out of 
  # the mass tolerance and the spike is within the bound
  if (length(oks1 <- .Internal(which(lens > 5L)))) {
    aps  <- aps[oks1]
    rtx  <- rtx[oks1]
    rgx  <- rgx[oks1]
    lens <- lens[oks1]
  }
  
  # (2) subset apexs by MS1 mass tolerance
  if (length(ok_xs  <- .Internal(which(abs(xs_k - xm) / xm <= step)))) {
    if (length(ok_aps <- .Internal(which(aps %fin% ss[ok_xs])))) {
      aps  <- aps[ok_aps]
      rtx  <- rtx[ok_aps]
      rgx  <- rgx[ok_aps]
      lens <- lens[ok_aps]
    }
  }

  # abs(xs_k[names(xs_k) %in% aps] - xm) / xm
  
  # (3) remove one-hit-wonders and spikes
  if (length(oks3 <- .Internal(which(lens > 15L)))) { # may also include 20
    aps  <- aps[oks3]
    rtx  <- rtx[oks3]
    rgx  <- rgx[oks3]
    lens <- lens[oks3]
  }
  else if (length(oks4 <- .Internal(which(lens > 10L)))) {
    aps  <- aps[oks4]
    rtx  <- rtx[oks4]
    rgx  <- rgx[oks4]
    lens <- lens[oks4]
  }
  
  if (FALSE) {
    if (length(lens) > 3L) {
      oksn <- which_topx2(lens, 3L)
      aps  <- aps[oksn]
      rtx  <- rtx[oksn]
      rgx  <- rgx[oksn]
      lens <- lens[oksn]
    }
  }

  # (4) subset apexes by distances between apex_scan and the triggering MS2 scan
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
    x1  <- xs[[1]]
    y1  <- ys[[1]]
    ap1 <- aps[[1]]
    ssx <- ss[rgx[[1]]]
    return(list(x = x1, y = y1, ap = ap1, from = ssx[[1]], 
                to = ssx[[lens]], aps = ap1, xs = x1, ys = y1))
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
    return(list(x = xa, y = ya, ap = apa, from = ssa[[1]], to = ssa[[lena]], 
                aps = aps, xs = xs, ys = ys))
  }
  
  # The trigger MS2 scan approximately within the scans of a MS1 peak profile
  if (FALSE) {
    if (scan >= ssa[[1]] - 3L && scan <= ssa[[lena]] + 3L) {
      return(list(x = xa, y = ya, ap = apa, from = ssa[[1]], to = ssa[[lena]], 
                  aps = aps, xs = xs, ys = ys))
    }
    
    if (scan >= ssb[[1]] - 3L && scan <= ssb[[lenb]] + 3L) {
      return(list(x = xb, y = yb, ap = apb, from = ssb[[1]], to = ssb[[lenb]], 
                  aps = aps, xs = xs, ys = ys))
    }
  }

  # compare peak distances
  if (db > da + 200) {
    return(list(x = xa, y = ya, ap = apa, from = ssa[[1]], to = ssa[[lena]], 
                aps = aps, xs = xs, ys = ys))
  }
  else {
    return(list(x = xb, y = yb, ap = apb, from = ssb[[1]], to = ssb[[lenb]], 
                aps = aps, xs = xs, ys = ys))
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

    ys   <- df1$msx_ints[[i]]
    ixxs <- index_mz(df1$msx_moverzs[[i]], from, step)
    ixs1 <- index_mz(xs1, from, step)
    ps0  <- match(ixs1, ixxs)
    oks  <- !is.na(ps0)
    
    if (all(oks)) {
      df1[["ms1_int"]][[i]] <- ys[ps0]
    }
    else {
      ps1 <- match(ixs1 + 1L, ixxs)
      ps2 <- match(ixs1 - 1L, ixxs)
      i1  <- .Internal(which(!is.na(ps1)))
      i2  <- .Internal(which(!is.na(ps2)))
      i0  <- .Internal(which(oks))
      
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



