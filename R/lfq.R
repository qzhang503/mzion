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
  if (n_chunks <= 1L)
    n_chunks <- 2L
  
  df1$ms1_mass <- df1$ms1_moverz <- df1$ms1_int <- df1$ms1_charge <- NULL
  df1s  <- chunksplit(df1, n_chunks, type = "row")
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
    for (i in 2:(n_chunks-1L))
      df1s[[i]] <- dplyr::bind_rows(df1s_af[[i-1]], df1s[[i]], df1s_bf[[i+1]])
  }
  rm(list = c("df1s_bf", "df1s_af", "df1"))
  
  ##  Splits `df` with bracketing entries
  # values: row indexes in `df`, length == nrow(df1) at pad_nas = TRUE
  ms1_stas <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)$ms1_stas
  # cols <- c("ms_level", "ms1_moverz", "ms1_int")
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
  list(dfs = dfs, df1s = df1s, gaps = gaps)
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
#' @inheritParams matchMS
htraceXY <- function (xs, ys, ss, ts, df, gap_bf = 256L, gap_af = 256L, 
                      n_mdda_flanks = 6L, from = 200L, step = 1E5, 
                      y_perc = .01, yco = 500, look_back = TRUE)
{
  if (all(lengths(xs) == 0L)) {
    # length(xs) - gap_bf - gap_af - nrow(df) == 0L
    null <- rep_len(list(NULL), nrow(df))
    
    return(
      tibble::tibble(
        ms_level = df$ms_level, 
        ms1_moverz = null, 
        ms1_int = null, 
        apex_scan_num = null)
    )
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
  scan_apexs <- vector("list", nc)
  for (i in 1:nc) {
    scan_apexs[[i]] <- as.integer(ss[apes[[i]]])
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
    scan_apexs = scan_apexs, rngs = rngs, from = from, step = step)
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
traceXY <- function (xs, ys, ss, ts, n_mdda_flanks = 6L, from = 115L, 
                     step = 8E-6, reord = TRUE, cleanup = FALSE, 
                     replace_ms1_by_apex = FALSE,y_perc = .01, yco = 100, 
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
  ans <- collapse_mms1ints(
    xs = xs, ys = ys, lwr = from, step = step, reord = FALSE, cleanup = FALSE, 
    add_colnames = TRUE, look_back = look_back)
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  nr <- nrow(ansy)
  nc <- ncol(ansy)
  rm(list = c("ans"))
  
  if (FALSE) {
    rng <- 466:615
    i <- which(colnames(ansx) == index_mz(483.7299, from, step) +  1) # 3663
    i <- which(colnames(ansx) == index_mz(483.7299, from, step) - 1) # 3662
    ss[rng]
    
    plot(ansx[rng, i])
    plot(ansy[rng, i])
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
      # i <- which(colnames(xmat) == index_mz(935.0734, from, step)) # 10287
      yi <- ansy[, i]
      oks <- .Internal(which(!is.na(yi)))
      yoks <- yi[oks]
      # may be unnecessary, e.g., Thermo's MS1 peak distributions are discrete
      yoks[yoks < yco] <- NA_real_
      yi[oks] <- yoks
      gates <- find_lc_gates(ys = yi, ts = ts, n_dia_scans = n_mdda_flanks)

      apexes[[i]] <- rows <- gates[["apex"]]
      ns[[i]] <- gates[["ns"]] # number of observing scans
      ranges[[i]] <- rngs <- gates[["ranges"]]
      yints <- gates[["yints"]]
      
      for (j in seq_along(rows)) {
        rwj <- rows[[j]]
        rgj <- rngs[[j]]
        xmat[rgj, i] <- ansx[rwj, i]
        # ymat[rgj, i] <- ansy[rwj, i]
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
#' @param rngs The ranges of row indexes for each column in \code{matx}. The
#'   \code{rngs[[i]]} corresponds to \code{matx[, i]}. The argument is currently
#'   used for debugging.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @importFrom fastmatch %fin%
updateMS1Int2 <- function (df, matx, maty, row_sta, row_end, scan_apexs, rngs, 
                           from = 200L, step = 1E-5)
  
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
      df$orig_scan[ms1_stas] # look for orig_scan around 18950 ->
      i = 232
      ms2sta <- ms2_stas[[i]]
      ms2end <- ms2_ends[[i]]
      df2 <- df[ms2sta:ms2end, ]
    }
    
    # i = 232; which(rownames(matx) == 18949) # 482 - gap = 232
    ms2sta <- ms2_stas[[i]]
    if (is.na(ms2sta)) next
    ms2end <- ms2_ends[[i]]
    df2 <- df[ms2sta:ms2end, ]
    
    rowi <- gap + i # row number in matx, maty
    scan <- ss[[rowi]] # the MS1 scan number at the current row
    
    for (j in 1:nrow(df2)) {
      # j <- 3
      x2s <- df2[["ms1_moverz"]][[j]] # MS1 masses associated with an MS2 scan
      nx <- length(x2s) # nx > 1 with a chimeric spectrum
      if (!nx) next
      ix2s <- as.integer(ceiling(log(x2s/from)/log(1+step)))
      ks <- lapply(ix2s, function (x) which(abs(x - unv) <= 1L))
      
      # go through chimeric precursors
      for (m in 1:nx) {
        k <- ks[[m]] # the k-th column
        
        if (!length(k)) {
          next
        }

        # if (any(k %in% c(3662, 3664))) { print(k) }

        # two adjacent matches: unv[k] <-> ix2s[[m]] - 1 and ix2s[[m]] + 1
        if (length(k) > 1L) {
          # mk <- lapply(k, function (x) median(matx[, x], na.rm = TRUE))
          # mk <- .Internal(unlist(mk, use.names = FALSE, recursive = FALSE))
          # dx <- abs(x2s[[m]] - mk)
          
          # the first column
          ka <- k[[1]]
          apexa <- scan_apexs[[ka]] # all apexs (scan numbers) under the k-th mass column
          dsa <- abs(apexa - scan)

          p1a  <- .Internal(which.min(dsa))
          pxa  <- max(1L, p1a - 3L):min(length(dsa), p1a + 3L)
          
          # at least one (the most centered) peak is within 2500 scans
          if (dsa[p1a] <= 2500) {
            pxa <- pxa[dsa[pxa] <= 2500]
          }

          apa  <- apexa[pxa]
          oksa <- .Internal(which(ss %in% apa))
          ya   <- maty[oksa, ka]
          topa <- .Internal(which.max(ya))
          
          p1a  <- pxa[topa]
          y1a  <- ya[[topa]]
          rnga <- rngs[[ka]][[p1a]]
          ap1a <- apa[[topa]]
          ok1a <- oksa[topa]
          
          # the 2nd column
          kb <- k[[2]]
          apexb <- scan_apexs[[kb]]
          dsb <- abs(apexb - scan)
          
          p1b  <- .Internal(which.min(dsb))
          pxb  <- max(1L, p1b - 3L):min(length(dsb), p1b + 3L)
          
          if (dsb[p1b] <= 2500) {
            pxb <- pxb[dsb[pxb] <= 2500]
          }
          
          apb  <- apexb[pxb]
          oksb <- .Internal(which(ss %in% apb))
          yb   <- maty[oksb, kb]
          topb <- .Internal(which.max(yb))
          
          p1b  <- pxb[topb]
          y1b  <- yb[[topb]]
          rngb <- rngs[[kb]][[p1b]]
          ap1b <- apb[[topb]]
          ok1b <- oksb[topb]
          
          # the best of the two
          if (y1a > y1b) {
            k <- ka
            apexs <- apexa
            ds <- dsa
            p1 <- p1a
            ap1 <- ap1a
            ok1 <- ok1a
            y1 <- y1a
          }
          else {
            k <- kb
            apexs <- apexb
            ds <- dsb
            p1 <- p1b
            ap1 <- ap1b
            ok1 <- ok1b
            y1 <- y1b
          }
        }
        else {
          apexs <- scan_apexs[[k]] # all apexs (scan numbers) under the k-th mass column
          ds <- abs(apexs - scan)

          if (FALSE) {
            p1 <- .Internal(which.min(ds))
            ap1 <- apexs[[p1]] # the nearest apex scan number; 13525
            ok1 <- .Internal(which(ss == ap1))
            y1 <- maty[ok1, k]
          }
          
          p1  <- .Internal(which.min(ds))
          px  <- max(1L, p1 - 3L):min(length(ds), p1 + 3L)
          
          if (ds[p1] <= 2500) {
            px <- px[ds[px] <= 2500]
          }
          
          ap  <- apexs[px]
          oks <- .Internal(which(ss %in% ap))
          y   <- maty[oks, k]
          top <- .Internal(which.max(y))
          
          p1  <- px[top]
          y1  <- y[[top]]
          rng <- rngs[[k]][[p1]]
          ap1 <- ap[[top]]
          ok1 <- oks[top]
        }
        
        df2$apex_scan_num[[j]][[m]] <- ap1
        df2[["ms1_int"]][[j]][m] <- y1
        
        # k <- k[[1]] # can have two adjacent matches: ix2s - 1 and ix2s + 1
        # ix <- unv[[k]]; # a <- maty[, which(colnames(matx) == unv[[k]])]; a <- maty[, 3218:3219]
        
        if (FALSE) {
          # checks neighbors
          if (length(apexs) > 1L) {
            ps <- which_topx2(-ds, 2L)
            p2 <- ps[ps != p1]
            ap2 <- apexs[[p2]] # the second nearest apex scan; 13686
            ok2 <- .Internal(which(ss == ap2))
            y2 <- maty[ok2, k]
            
            # later checks peak spacing and peak width (for small satellite peaks)
            if ((y1 < y2) && (abs(ap2 - ap1) <= 200L)) { # 200L somewhat arbitrary
              df2$apex_scan_num[[j]][[m]] <- ap2 # ss[ok2]
              df2[["ms1_int"]][[j]][m] <- y2
            }
            else {
              df2$apex_scan_num[[j]][[m]] <- ap1 # ss[ok1]
              df2[["ms1_int"]][[j]][m] <- y1
            }
          }
          else {
            df2$apex_scan_num[[j]][[m]] <- ap1 # ss[ok1]
            df2[["ms1_int"]][[j]][m] <- y1
          }
        }

        
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


