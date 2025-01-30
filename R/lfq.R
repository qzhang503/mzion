#' Subsets full MS1 by the universe of monoisotopic MS1 moverzs.
#'
#' Non-monoisotopic (13C_n) MS1-X are also removed.
#'
#' @param xs Vectors of MS1-X values.
#' @param ys Vectors of MS1-Y values.
#' @param ts Vectors of MS1 retention times.
#' @param ms1s Vectors of MS1 m/z values that are linked to the mono-isotopic
#'   DDA-MS2 space.
#' @param from The starting m/z value for calculating bin indexes.
#' @param step The size of bins.
#' @param rt_tol The tolerance in retention times for m-over-zs looking up.
#' @param rt_gap An approximate gap size for forward and backward looking of
#'   precursors in a +/-3-min retention times. The argument is for fast
#'   subsetting without looking up the entire sequence. Note that at the early
#'   retention times of LC, there may be many more MS1 scans and thus a
#'   \code{gap = 500} may encompassing a narrower retention window (e.g. 2-min).
#'   This may not be critical since there are often fewer features in early RT.
#' @importFrom fastmatch %fin%
subMSfull <- function (xs, ys, ts, ms1s, from = 200L, step = 1E-5, 
                       rt_gap = 500L, rt_tol = 180)
{
  lens <- lengths(xs)
  len  <- length(xs)
  yout <- xout <- vector("list", len)
  
  # for each xi, keeps entries found in ms1s.
  if (rt_gap >= len) {
    ms1s <- .Internal(unlist(ms1s, recursive = FALSE, use.names = FALSE))
    ms1s <- index_mz(ms1s, from, step)
    ms1s <- sort(unique(ms1s))
    
    for (i in 1:len) {
      xi <- xs[[i]]
      li <- lens[[i]]
      if (!li) { next }
      
      ix  <- index_mz(xi, from, step)
      ps0 <- fastmatch::fmatch(ix, ms1s)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s)

      i012 <- .Internal(which((!is.na(ps0)) | (!is.na(ps1)) | (!is.na(ps2))))
      xout[[i]] <- xi[i012]
      yout[[i]] <- ys[[i]][i012]
    }
  }
  else {
    for (i in 1:len) {
      xi <- xs[[i]]
      li <- lens[[i]]
      if (!li) { next }
      
      ti <- ts[[i]]
      ### slow
      # rows <- which(ts >= ti - rt_tol & ts <= ti + rt_tol)
      # sta <- rows[[1]]
      # end <- rows[[length(rows)]]
      ### 
      tx   <- ts[max(1L, i - rt_gap):min(len, i + rt_gap)]
      rows <- which(tx >= ti - rt_tol & tx <= ti + rt_tol)
      off  <- max(i - rt_gap - 1L, 0L)
      sta  <- rows[[1]] + off
      end  <- rows[[length(rows)]] + off
      # sta <- max(1, i - rt_gap)
      # end <- min(i + rt_gap, len)

      # the possible universe of mono-isotopic MS1 m/z values over a range
      ms1s_sub <- 
        .Internal(unlist(ms1s[sta:end], recursive = FALSE, use.names = FALSE))
      ms1s_sub <- index_mz(ms1s_sub, from, step)
      ms1s_sub <- sort(unique(ms1s_sub))
      
      ix  <- index_mz(xi, from, step)
      ps0 <- fastmatch::fmatch(ix, ms1s_sub)
      ps1 <- fastmatch::fmatch(ix + 1L, ms1s_sub)
      ps2 <- fastmatch::fmatch(ix - 1L, ms1s_sub)

      i012 <- .Internal(which((!is.na(ps0)) | (!is.na(ps1)) | (!is.na(ps2))))
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
#' @return dfs[[i]]:a chunk of MS1 and MS2 data (no bracketing scans);
#'   df1s[[i]]: the reference MS1 data with leading and trailing scans;
#'   gaps[[i]]: the number of bracketing MS1 and MS2 scans between chunks.
#' @inheritParams subMSfull
pretraceXY <- function (df, from = 200L, step = 1e-5, rt_gap = 500L, 
                        rt_tol = 180)
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
    xs = df1$msx_moverzs, ys = df1$msx_ints, ts = df1$ret_time,
    ms1s = df1$ms1_moverz, from = from, step = step, 
    rt_gap = rt_gap, rt_tol = rt_tol)
  # moverzs in ans$x are in an ascending order
  df1$msx_moverzs <- ans$x
  df1$msx_ints    <- ans$y
  rm(list = "ans")

  df1$ms1_mass <- df1$ms1_moverz <- df1$ms1_int <- df1$ms1_charge <- NULL
  
  df1s <- sep_df1_byRTs(df1, col_rt = "ret_time")
  gaps_bf  <- attr(df1s, "gaps_bf",  exact = TRUE)
  gaps_af  <- attr(df1s, "gaps_af",  exact = TRUE)
  min_rts  <- attr(df1s, "min_rts",  exact = TRUE)
  max_rts  <- attr(df1s, "max_rts",  exact = TRUE)
  end1s    <- attr(df1s, "end1s",    exact = TRUE)
  sta1s    <- attr(df1s, "sta1s",    exact = TRUE)
  n_chunks <- attr(df1s, "n_chunks", exact = TRUE)
  
  ##  Splits `df` with bracketing entries
  # values: row indexes in `df`, length == nrow(df1) at pad_nas = TRUE
  if (n_chunks == 1L) {
    return(
      list(dfs = list(df), df1s = df1s, gaps_bf = gaps_bf, gaps_af = gaps_af, 
           min_rts = min_rts, max_rts = max_rts))
  }
  
  ms1_stas <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)$ms1_stas
  cols     <- c("ms_level", "ms1_moverz", "ms1_int", "orig_scan")
  ms1_stax <- vector("integer", n_chunks)
  dfs      <- vector("list", n_chunks)
  
  for (i in 1:n_chunks) {
    stai <- sta1s[[i]] # stai-th MS1 scan
    ms1_stax[[i]] <- ms1_stas[stai] # the corresponding row index in df
  }

  for (i in 1:(n_chunks - 1L)) {
    rowx <- ms1_stax[[i]]:(ms1_stax[[i + 1L]] - 1L) # ms2_endx can be NA
    dfs[[i]] <- df[rowx, cols] # both df1 and df2 data
  }
  dfs[[n_chunks]] <- df[ms1_stax[[n_chunks]]:nrow(df), cols]

  # dfs[[i]]:  a chunk MS1 and MS2 data (no bracketing scans)
  # df1s[[i]]: the reference MS1 data + leading and trailing MS1 scans; 
  # gaps_bf[[i]]: the number of preceding MS1 and MS2 scans...
  list(dfs = dfs, df1s = df1s, gaps_bf = gaps_bf, gaps_af = gaps_af, # gaps = gaps, 
       min_rts = min_rts, max_rts = max_rts)
}


#' Separate MS1 data
#'
#' Segments of 2+3+2 mins
#'
#' @param df1 A data frame of MS1.
#' @param col_rt The key of retention time column.
#' @param rt_size The width of each LC retention times in seconds.
#' @param rt_margin The bracketing margin before and after an LC retention time
#'   window.
sep_df1_byRTs <- function (df1, col_rt = "ret_time", rt_size = 180, 
                           rt_margin = 120)
{
  # if (n_chunks <= 1L) { stop("Developer: need at least two chunks.") }
  rts   <- df1[[col_rt]]
  rtmin <- min(rts, na.rm = TRUE)
  rtmax <- max(rts, na.rm = TRUE)
  delta <- rtmax - rtmin
  n_chunks <- ceiling(delta / rt_size)
                      
  if (n_chunks < 2L) {
    df1s <- list(df1)
    attr(df1s, "gaps_bf") <- 0L
    attr(df1s, "gaps_af") <- 0L
    attr(df1s, "min_rts") <- rtmin
    attr(df1s, "max_rts") <- rtmax
    attr(df1s, "end1s") <- nrow(df1)
    attr(df1s, "sta1s") <- 1L
    attr(df1s, "n_chunks") <- 1L
    
    return(df1s)
  }
  
  width <- delta / n_chunks
  brs   <- rtmin + width * seq_len(n_chunks - 1L)
  df1s  <- vector("list", n_chunks)

  for (i in 1:n_chunks) {
    rows <- if (i == 1L) {
      rts <= brs[[i]]
    }
    else if (i == n_chunks) {
      rts > brs[[n_chunks - 1L]]
    }
    else {
      rts > brs[[i-1]] & rts <= brs[[i]]
    }
    
    df1s[[i]] <- df1[rows, ]
  }
  
  min_rts <- unlist(lapply(df1s, function (x) x[[col_rt]][[1]]))
  max_rts <- unlist(lapply(df1s, function (x) x[[col_rt]][[nrow(x)]]))
  end1s   <- unname(cumsum(lapply(df1s, nrow)))
  sta1s   <- c(1L, end1s[1:(n_chunks - 1L)] + 1L)
  
  # Adds 2-min gaps before and after
  gaps_bf <- gaps_af <- vector("integer", n_chunks)
  
  for (i in seq_along(df1s)) {
    dfi  <- df1s[[i]]
    rti  <- dfi[[col_rt]]
    tmin <- min_rts[[i]]
    tmax <- max_rts[[i]]
    nri  <- nrow(dfi)
    
    if (i < n_chunks) {
      rows <- which(rti > tmax - rt_margin)
      gaps_af[[i]] <- nri - rows[[1]] + 1L
    }
    
    if (i > 1L) {
      rows <- which(rti <= tmin + rt_margin)
      gaps_bf[[i]] <- rows[[length(rows)]]
    }
  }
  
  df1s_bf <- df1s_af <- vector("list", n_chunks)
  
  for (i in 1:n_chunks) {
    df1s_bf[[i]] <- head(df1s[[i]], gaps_bf[[i]])
    df1s_af[[i]] <- tail(df1s[[i]], gaps_af[[i]])
  }
  
  df1s[[1]] <- dplyr::bind_rows(df1s[[1]], df1s_bf[[2]])
  df1s[[n_chunks]] <- dplyr::bind_rows(df1s_af[[n_chunks-1]], df1s[[n_chunks]])
  
  if (n_chunks > 2L) {
    for (i in 2:(n_chunks-1L)) {
      df1s[[i]] <- dplyr::bind_rows(df1s_af[[i-1]], df1s[[i]], df1s_bf[[i+1]])
    }
  }
  
  gaps_bf2 <- gaps_af2 <- vector("integer", n_chunks)
  gaps_bf2[2:n_chunks] <- gaps_af[1:(n_chunks-1L)]
  gaps_af2[1:(n_chunks-1L)] <- gaps_bf[2:n_chunks]
  
  attr(df1s, "gaps_bf")  <- gaps_bf2
  attr(df1s, "gaps_af")  <- gaps_af2
  attr(df1s, "min_rts")  <- min_rts
  attr(df1s, "max_rts")  <- max_rts
  attr(df1s, "end1s")    <- end1s
  attr(df1s, "sta1s")    <- sta1s
  attr(df1s, "n_chunks") <- n_chunks
  
  df1s
}


#' Helper of MS1 tracing.
#'
#' @param xs Vectors of full-spectrum MS1 m/z values.
#' @param ys Vectors of full-spectrum MS1 intensities.
#' @param ss Vectors of MS1 (original) scan numbers.
#' @param ts Vectors of MS1 retention times (for calculating area-under-a-peak).
#' @param df A data frame of MS1 and MS2 corresponding to \code{xs}.
#' @param gap_bf A preceding gap size.
#' @param gap_af A following gap size.
#' @param from The starting point in mass binning.
#' @param step_tr The bin size of MS1 tracing
#' @param tol_ms1 The tolerance in MS1 errors \code{ppm_ms1 * 1E-6}.
#' @param yco The cut-off in y values.
#' @param ytot_co A more permissive cut-off in peak area for finding all peaks.
#' @param y_perc The cut-off in intensity values in relative to the base peak.
#' @param min_y A more strict cut-off in peak area for finding tentatively the
#'   best peak.
#' @param min_n1 The first cut-off of a minimum number of points across an MS1
#'   peak.
#' @param min_n2 The second cut-off of a minimum number of points across an MS1
#'   peak.
#' @param min_n3 The third cut-off of a minimum number of points across an MS1
#'   peak.
#' @param sum_y Logical; sum Y values at near the same X values or not. Mostly
#'   FALSE with Thermo's data.
#' @param path_ms1 The path for saving the output of MS1 data.
#' @param out_name A file name of output.
#' @param out_cols The output column keys.
#' @inheritParams matchMS
htraceXY <- function (xs, ys, ss, ts, df, gap_bf = 256L, gap_af = 256L, 
                      out_name = NULL, n_dia_scans = 6L, from = 200L, 
                      step_tr = 6E-6, tol_ms1 = 1E-5, y_perc = .01, yco = 100, 
                      min_y = 2E6, ytot_co = 2E5, sum_y = FALSE, 
                      min_n1 = 10L, min_n2 = 20L, min_n3 = 15L, 
                      path_ms1 = NULL, 
                      out_cols = 
                        c("ms_level", "ms1_moverz", "ms1_int", "orig_scan", 
                          "apex_scan_num", "apex_xs", "apex_ys", "apex_ps", 
                          "apex_ts"))
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
      # orig_scan = null, # don't: need `character` not `list`
      apex_scan_num = null, 
      apex_xs = null, 
      apex_ys = null, 
      apex_ps = null, 
      apex_ts = null)

    return(df[, out_cols])
  }
  
  ### (1) Make MS1 X and Y matrices across adjacent scans
  ixs <- lapply(xs, index_mz, from = from, d = step_tr)
  
  # 1. remove duplicated ixs and collapse the Y values under the same ixs
  #    Y values are only for de-isotoping, not for precursor intensities
  if (sum_y) {
    for (i in seq_along(xs)) {
      ix <- ixs[[i]]
      x  <- xs[[i]]
      y  <- ys[[i]]
      ps <- .Internal(which(duplicated(ix)))
      
      if (l <- length(ps)) {
        for (j in 1:l) {
          okp <- .Internal(which(ix == ix[ps[[j]]]))
          y[okp[[1]]] <- sum(y[okp], na.rm = TRUE)
        }
        
        ixs[[i]] <- ix[-ps]
        xs[[i]]  <- x[-ps]
        ys[[i]]  <- y[-ps]
      }
    }
  }
  else {
    for (i in seq_along(xs)) {
      ix <- ixs[[i]]
      x  <- xs[[i]]
      y  <- ys[[i]]
      ps <- .Internal(which(duplicated(ix)))
      
      if (length(ps)) {
        ixs[[i]] <- ix[-ps]
        xs[[i]]  <- x[-ps]
        ys[[i]]  <- y[-ps]
      }
    }
  }
  # rm(list = c("x", "y", "ix", "ps"))
  
  ## maps ixs vectors to unv (presence or absence)
  unv  <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv  <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups  <- lapply(ixs, function (x) unv %fin% x)
  
  # note one-to-one correspondence between ixs and xs
  matx <- mapcoll_xyz(vals = xs, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)
  maty <- mapcoll_xyz(vals = ys, ups = ups, lenx = lenx, lenu = lenu, 
                      direct_out = TRUE)

  if (gap_bf) {
    if (gap_af) { # middle
      sta <- gap_bf + 1L
      end <- lenx - gap_af
    }
    else { # last
      sta <- gap_bf + 1L
      end <- lenx
    }
  }
  else { # first
    sta <- 1L
    end <- lenx - gap_af
  }
  
  # look both before and after scans
  ss <- as.integer(ss)
  
  args <- list(
    df = df, matx = matx, maty = maty, unv = unv, ss = ss, ts = ts, 
    row_sta = sta, row_end = end, from = from, step_tr = step_tr, 
    tol_ms1 = tol_ms1, min_y = min_y, yco = yco, ytot_co = ytot_co, 
    n_dia_scans = n_dia_scans, 
    min_n1 = min_n1, min_n2 = min_n2, min_n3 = min_n3)

  df <- do.call(updateMS1Int, args)

  attr(df, "row_sta") <- sta
  attr(df, "row_end") <- end
  attr(df, "from")    <- from
  attr(df, "step_tr") <- step_tr
  attr(df, "min_y")   <- min_y
  
  df <- df[, out_cols]
}


#' Updates MS1 intensity with apex values.
#'
#' Including forward and backward looking.
#'
#' @param df A data frame.
#' @param matx The matrix of m-over-z's Y: by masses; X: by LC scans.
#' @param maty The matrix of intensities. Y: by masses; X: by LC scans.
#' @param unv The universe of m-over-z bins.
#' @param ss A vector of scan numbers.
#' @param ts A vector of retention times.
#' @param row_sta The starting row of \code{matx}.
#' @param row_end The ending row of \code{matx}.
#' @param step_tr A step size for mass binning.
#' @param tol_ms1 The tolerance in MS1 errors \code{ppm_ms1 * 1E-6}.
#' @param from The starting point for mass binning.
#' @param min_y The cut-off in peak area.
#' @param yco An intensity (Y-value) cut-off.
#' @param ytot_co A more permissive cut-off in peak area for finding all peaks.
#' @param min_n1 The first cut-off of a minimum number of points across an MS1
#'   peak.
#' @param min_n2 The second cut-off of a minimum number of points across an MS1
#'   peak.
#' @param min_n3 The third cut-off of a minimum number of points across an MS1
#'   peak.
#' @inheritParams matchMS
#' @importFrom fastmatch %fin%
updateMS1Int <- function (df, matx, maty, unv, ss, ts, row_sta, row_end, 
                          from = 200L, step_tr = 6E-6, tol_ms1 = 1E-5, 
                          min_y = 2E6, yco = 100, ytot_co = 2E5, 
                          min_n1 = 10L, min_n2 = 20L, min_n3 = 15L, 
                          n_dia_scans = 6L)
  
{
  df <- init_lfq_cols(df)
  nr <- length(ss)
  nc <- length(unv)
  
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
  
  rt_apexs <- scan_apexs <- vector("list", nc)
  
  for (i in seq_along(ms1_stas)) {
    if (FALSE) {
      # look for penultimate scan number of psmQ.txt::pep_scan_num
      df$orig_scan[ms1_stas]
      i <- 246
      ms2sta <- ms2_stas[[i]]
      ms2end <- ms2_ends[[i]]
      ms2rng <- ms2sta:ms2end
      df2    <- df[ms2rng, ]
    }

    # i = 29; which(rownames(matx) == 78900) - gap -> i
    ms2sta <- ms2_stas[[i]]
    if (is.na(ms2sta)) { next }
    ms2end <- ms2_ends[[i]]
    ms2rng <- ms2sta:ms2end
    df2    <- df[ms2rng, ]
    rowi   <- gap + i    # row number in matx, maty
    sref   <- ss[[rowi]] # the MS1 scan number at the current row
    tref   <- ts[rowi]
    xs2    <- df2[["ms1_moverz"]]
    
    # by MS2 spectra
    for (j in 1:nrow(df2)) {
      # j <- 10L
      x1s  <- xs2[[j]] # MS1 masses associated with an MS2 scan
      nx   <- length(x1s) # nx > 1 at a chimeric spectrum
      if (!nx) { next }
      ix1s <- index_mz(x1s, from = from, d = step_tr)
      ks   <- lapply(ix1s, function (x) which(abs(x - unv) <= 1L))
      
      # by chimeric precursors
      for (m in 1:nx) {
        km   <- ks[[m]]
        nkm  <- length(km)
        if (!nkm) { next }
        xref <- x1s[[m]]
        
        if (nkm == 1L) {
          xvs  <- matx[, km, drop = TRUE]
          yvs  <- maty[, km, drop = TRUE]
          # later keep slightly greater than tol_ms1 but tangent entries...
          bads <- .Internal(which(abs(xvs / xref - 1.0) > tol_ms1))
          
          if (length(bads)) {
            xvs[bads] <- NA_real_
            yvs[bads] <- NA_real_
          }
        }
        else {
          mx   <- matx[, km, drop = FALSE]
          my   <- maty[, km, drop = FALSE]
          # later keep slightly greater than tol_ms1 but tangent entries...
          bads <- .Internal(which(abs(mx / xref - 1.0) > tol_ms1))
          
          if (length(bads)) {
            mx[bads] <- NA_real_ # can be all NAs after this
            my[bads] <- NA_real_
          }
          
          xvs <- rowMeans(mx, na.rm = TRUE)
          yvs <- rowSums(my, na.rm = TRUE) # all NA row: NA -> 0
          xvs[is.nan(xvs)] <- NA_real_
        }
        
        # yvs[yvs < yco] <- NA_real_
        # if (!length(oks <- .Internal(which(!is.na(yvs))))) { next }

        if (FALSE) {
          data.frame(x = ts/60, y = yvs) |>
            ggplot2::ggplot() + 
            ggplot2::geom_segment(mapping = aes(x = x, y = y, xend = x, yend = 0), 
                                  color = "gray", linewidth = .1)
        }
        
        gates <- find_lc_gates(
          xs = xvs, ys = yvs, ts = ts, 
          yco = yco, ytot_co = ytot_co, n_dia_scans = n_dia_scans)

        if (is.null(gates)) {
          next
        }
        
        rows  <- gates[["apex"]] # can be 0 with bad or unknown apexes
        # if (all(rows == 0L)) { next }
        rngs  <- gates[["ranges"]]
        yints <- gates[["yints"]]
        fwhms <- gates[["fwhm"]]
        
        # xbars <- lapply(rngs, function (rng) mean(xvs[rng], na.rm = TRUE))
        # xbars <- .Internal(unlist(xbars, recursive = FALSE, use.names = FALSE))
        if (length(rows)) {
          ans <- find_best_apex(
            xvs = xvs, yvs = yvs, ss = ss, 
            yints = yints, aps = ss[rows], rts = ts[rows], rngs = rngs, 
            fwhms = fwhms, xref = xref, ret_ms1 = tref, # scan_ms1 = sref, 
            # ??? only keep those apexes with tracing error <= step_tr ???
            step_tr = step_tr, min_y = min_y, 
            min_n1 = min_n1, min_n2 = min_n2, min_n3 = min_n3)
          if (is.null(ans)) { next }

          ap1 <- ans[["ap"]]
          y1  <- ans[["y"]]
          aps <- ans[["aps"]]
          rts <- ans[["rts"]]
          xs  <- ans[["xs"]]
          ys  <- ans[["ys"]]
          # fwhms <- ans[["fwhms"]]
        }
        else {
          next
        }
        
        # scalar
        df2[["apex_scan_num"]][[j]][[m]] <- ap1
        df2[["ms1_int"]][[j]][[m]]       <- y1
        
        # vector
        df2[["apex_ps"]][[j]][[m]] <- aps
        df2[["apex_ts"]][[j]][[m]] <- rts
        
        # list
        df2[["apex_xs"]][[j]][[m]]   <- xs
        df2[["apex_ys"]][[j]][[m]]   <- ys
      }
    }
    
    df[["apex_scan_num"]][ms2rng] <- df2[["apex_scan_num"]]
    df[["ms1_int"]][ms2rng]       <- df2[["ms1_int"]]
    df[["apex_ps"]][ms2rng]       <- df2[["apex_ps"]]
    df[["apex_ts"]][ms2rng]       <- df2[["apex_ts"]]
    df[["apex_xs"]][ms2rng]       <- df2[["apex_xs"]]
    df[["apex_ys"]][ms2rng]       <- df2[["apex_ys"]]
  }
  
  # contains df[["ms1_int"]] not zero (from deisotoping) and df2[["apex_ps"]] 
  # are NULLs, clear up in psmC.txt

  df
}


#' Helper to initiate LFQ output columns
#' 
#' @param df A data frame.
init_lfq_cols <- function (df)
{
  nr <- nrow(df)
  
  df[["apex_ys"]] <- df[["apex_xs"]] <- df[["apex_ts"]] <- 
    df[["apex_ps"]] <- df[["apex_scan_num"]] <- vector("list", nr)
  
  # replicate the lengths by the number of chimeric precursors
  rows2 <- which(df$ms_level == 2L)
  lens2 <- lengths(df$ms1_int[rows2])
  
  df[["apex_scan_num"]][rows2] <- lapply(lens2, function (x) rep_len(0L, x))
  
  df[["apex_ps"]][rows2] <- df[["apex_ts"]][rows2] <- 
    df[["apex_xs"]][rows2] <- df[["apex_ys"]][rows2] <- 
    lapply(lens2, function (x) rep_len(list(NULL), x))
  
  df
}


#' Find tentatively the apex under a mass column
#' 
#' Also remove low-quality apexes.
#' 
#' @param xvs X values.
#' @param yvs Y values.
#' @param yints Peak areas under under each gate.
#' @param aps All apexes scan numbers under a mass column.
#' @param rts All apexes retention times under a mass column.
#' @param rngs All apexes scan ranges under a column.
#' @param fwhms All FWHMs under a column.
#' @param xref The experimental X value for tracing.
#' @param ret_ms1 The preceding MS1 retention time of an MS2 event at
#'   \code{xref}. Or may be the MS2 retention time of \code{xref}.
#' @param ss All of the scan numbers.
#' @param rngs The scan ranges of apexes under column k.
#' @param step_tr A step size.
#' @param min_y The cut-off of intensity values.
#' @param max_scan_delta The maximum difference in scan numbers between an apex
#'   and \code{scan_ms1}.
#' @param max_rt_delta The maximum difference in retention times between an apex
#'   and \code{ret_ms1}.
#' @importFrom fastmatch %fin%
find_best_apex <- function (xvs, yvs, ss, 
                            yints, aps, rts, rngs, fwhms, xref, # scan_ms1, 
                            ret_ms1, step_tr = 6E-6, min_y = 0, 
                            min_n1 = 10L, min_n2 = 20L, min_n3 = 15L, 
                            max_scan_delta = 3500, max_rt_delta = 180)
{
  lens  <- lengths(rngs)
  naps  <- length(aps)
  xbars <- lapply(rngs, function (rng) mean(xvs[rng], na.rm = TRUE))
  xbars <- .Internal(unlist(xbars, recursive = FALSE, use.names = FALSE))

  # (1) first-pass spike removals: in case that a major peak is just out of 
  # the mass tolerance and the spike is within the bound
  if (noks1 <- length(oks1 <- .Internal(which(lens > min_n1)))) {
    if (!noks1) { return(NULL) }
    
    if (noks1 < naps) {
      aps   <- aps[oks1]
      rts   <- rts[oks1]
      rngs  <- rngs[oks1]
      xbars <- xbars[oks1]
      yints <- yints[oks1]
      fwhms <- fwhms[oks1]
      lens  <- lens[oks1]
      naps  <- noks1
    }
  }
  
  # (2) all apexes are < MS1 mass tolerance, but subset further those < step_tr?
  if (FALSE && length(ok_xs  <- .Internal(which(abs(xvs / xref - 1) <= step_tr)))) {
    if (noks2 <- length(ok_aps <- .Internal(which(aps %fin% ss[ok_xs])))) {
      # FALSE: allow large mass error which might be due to space charge effects
      if (FALSE) {
        if (!noks2) { return(NULL) }
        
        if (noks2 < naps) {
          
        }
      }

      if (noks2 && noks2 < naps) {
        aps  <- aps[ok_aps]
        rts  <- rts[ok_aps]
        rngs  <- rngs[ok_aps]
        xbars <- xbars[ok_aps]
        yints <- yints[ok_aps]
        fwhms <- fwhms[ok_aps]
        lens <- lens[ok_aps]
        naps <- noks2
      }
    }
  }
  
  # (3) remove one-hit-wonders and spikes
  if (noks3 <- length(oks3 <- .Internal(which(lens > min_n2)))) {
    if (noks3 && noks3 < naps) {
      aps   <- aps[oks3]
      rts   <- rts[oks3]
      rngs  <- rngs[oks3]
      xbars <- xbars[oks3]
      yints <- yints[oks3]
      fwhms <- fwhms[oks3]
      lens  <- lens[oks3]
      naps  <- noks3
    }
  }
  else if (noks4 <- length(oks4 <- .Internal(which(lens > min_n3)))) {
    if (FALSE) {
      if (!noks4) { return(NULL) }
      
      if (noks4 < naps) {
        aps   <- aps[oks4]
        rts   <- rts[oks4]
        rngs  <- rngs[oks4]
        xbars <- xbars[oks4]
        yints <- yints[oks4]
        fwhms <- fwhms[oks4]
        lens  <- lens[oks4]
        naps  <- noks4
      }
    }
    
    if (noks4 && noks4 < naps) {
      aps   <- aps[oks4]
      rts   <- rts[oks4]
      rngs  <- rngs[oks4]
      xbars <- xbars[oks4]
      yints <- yints[oks4]
      fwhms <- fwhms[oks4]
      lens  <- lens[oks4]
      naps  <- noks4
    }
  }
  
  if (FALSE) {
    if (length(lens) > 3L) { 
      return(NULL)
      
      oksn <- which_topx2(lens, 3L)
      aps  <- aps[oksn]
      rts  <- rts[oksn]
      rngs  <- rngs[oksn]
      xbars <- xbars[oksn]
      yints <- yints[oksn]
      fwhms <- fwhms[oks4]
      lens <- lens[oksn]
    }
  }
  
  # (4) subset apexes by distances between apex_scan and the triggering MS2 scan
  ds <- abs(rts - ret_ms1)
  p1 <- .Internal(which.min(ds))
  d1 <- ds[[p1]]
  if (d1 > max_rt_delta) { return(NULL) }
  px <- max(1L, p1 - 3L):min(naps, p1 + 3L)
  ds <- ds[px]
  
  oks_d  <- .Internal(which(ds <= max_rt_delta))
  noks_d <- length(oks_d)
  
  if (noks_d < naps) {
    px <- px[oks_d]
    ds <- ds[oks_d]
    
    aps   <- aps[px]
    rts   <- rts[px]
    rngs  <- rngs[px]
    xbars <- xbars[px]
    yints <- yints[px]
    fwhms <- fwhms[px]
    lens  <- lens[px]
    naps  <- noks_d
  }
  
  # at least one (the most centered) peak is within 3500 scans
  # arbitrary and may disable this; or a function of intensity
  if (FALSE) {
    if (d1 <= max_rt_delta) {
      oks_d  <- .Internal(which(ds <= max_rt_delta))
      noks_d <- length(oks_d)
      
      if (noks_d < naps) {
        px <- px[oks_d]
        ds <- ds[oks_d]
        
        aps   <- aps[px]
        rts   <- rts[px]
        rngs  <- rngs[px]
        xbars <- xbars[px]
        yints <- yints[px]
        fwhms <- fwhms[px]
        lens  <- lens[px]
        naps  <- noks_d
      }
    }
    else {
      noks_d <- length(px)
      
      if (noks_d < naps) {
        aps   <- aps[px]
        rts   <- rts[px]
        rngs  <- rngs[px]
        xbars <- xbars[px]
        yints <- yints[px]
        fwhms <- fwhms[px]
        lens  <- lens[px]
        naps  <- noks_d
      }
    }
  }
  
  # if (!noks_d) { return(NULL) } # should not occur

  # by intensity cut-off
  oks_y  <- .Internal(which(yints > min_y))
  noks_y <- length(oks_y)
  
  if (!noks_y) { return(NULL) }
  
  if (noks_y < naps) {
    aps   <- aps[oks_y]
    rts   <- rts[oks_y]
    rngs  <- rngs[oks_y]
    xbars <- xbars[oks_y]
    yints <- yints[oks_y]
    fwhms <- fwhms[oks_y]
    lens  <- lens[oks_y]
    naps  <- noks_y
    
    ds <- ds[oks_y]
  }
  
  if (naps == 1L) {
    x1   <- xbars[[1]]
    y1   <- yints[[1]]
    ap1  <- aps[[1]]
    rt1  <- rts[[1]]
    ssx  <- ss[rngs[[1]]]
    fw1  <- fwhms[[1]]
    len1 <- lens[[1]]
    return(list(x = x1, y = y1, ap = ap1, rt = rt1, n = len1, fwhm = fw1, 
                from = ssx[[1]], to = ssx[[len1]], rng = rngs[[1]], 
                xs = x1, ys = y1, aps = ap1, rts = rt1, ns = len1, 
                fwhms = fwhms, rngs = rngs))
  }
  
  # the closest may be a spike... 
  # look for more "regular" peak nearby...
  
  # ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, ds))
  # topa <- ds[[ord[[1]]]]
  
  topa <- .Internal(which.min(ds))
  xa   <- xbars[[topa]]
  ya   <- yints[[topa]]
  apa  <- aps[[topa]]
  da   <- ds[[topa]]
  rga  <- rngs[[topa]]
  ssa  <- ss[rga]
  rta  <- rts[[topa]]
  lena <- lens[[topa]]
  fwa  <- fwhms[[topa]]
  
  # oka  <- apa %in% ssa # must be
  
  topb <- .Internal(which.max(yints))
  xb   <- xbars[[topb]]
  yb   <- yints[[topb]]
  apb  <- aps[[topb]]
  db   <- ds[[topb]]
  rgb  <- rngs[[topb]]
  ssb  <- ss[rgb]
  rtb  <- rts[[topb]]
  lenb <- lens[[topb]]
  fwb  <- fwhms[[topb]]
  # yvs[rgb]
  # okb  <- apb %in% ssb # must be
  
  # for checking the continuum of values...
  if (FALSE) {
    ansa <- find_gate_edges(which(!is.na(yvs[rga])))
    ansb <- find_gate_edges(which(!is.na(yvs[rgb])))
  }

  if (topa == topb) {
    return(list(x = xa, y = ya, ap = apa, rt = rta, n = lena, fwhm = fwa, 
                # starting and ending MS1 scan numbers
                from = ssa[[1]], to = ssa[[lena]], rng = rga, 
                xs = xbars, ys = yints, aps = aps, rts = rts, ns = lens, 
                fwhms = fwhms, rngs = rngs))
  }
  
  # compare peak distances
  # arbitrary: 
  #  PSM1 closer to rta but less significant; 
  #  PSM2 closer to rtb but more significant; 
  #  compare local patterns
  if (abs(rta - rtb) > 30) {
    return(list(x = xa, y = ya, ap = apa, rt = rta, n = lena, fwhm = fwa, 
                from = ssa[[1]], to = ssa[[lena]], rng = rga, 
                xs = xbars, ys = yints, aps = aps, rts = rts, ns = lens, 
                fwhms = fwhms, rngs = rngs))
  }
  else {
    return(list(x = xb, y = yb, ap = apb, rt = rtb, n = lenb, fwhm = fwb, 
                from = ssb[[1]], to = ssb[[lenb]], rng = rgb, 
                xs = xbars, ys = yints, aps = aps, rts = rts, ns = lens, 
                fwhms = fwhms, rngs = rngs))
  }
}


#' Checks overlaps in start and end positions
#' 
#' @param sta_b Scalar; the start index of b.
#' @param end_b Scalar; the end index of b. Note that \code{end_b >= sta_b}.
#' @param stas_a Vector; the start indexes of a.
#' @param ends_a Vector; the end indexes of a.
#' @param margin The margin for considering being intersecting.
#' @examples
#' oks <- mzion:::are_ovlap_ranges(sta_b = 450, end_b = 774, stas_a = c(412, 775, 850), ends_a = c(455, 843, 900), margin = 5)
#' # stopifnot(identical(oks, c(TRUE, TRUE, FALSE)))
#' oks <- mzion:::are_ovlap_ranges(sta_b = 470, end_b = 774, stas_a = c(412, 775, 850), ends_a = c(455, 843, 900), margin = 5)
#' # stopifnot(identical(oks, c(FALSE, TRUE, FALSE)))
#' oks <- mzion:::are_ovlap_ranges(sta_b = 470, end_b = 600, stas_a = c(412, 775, 850), ends_a = c(455, 843, 900), margin = 5)
#' # stopifnot(identical(oks, c(FALSE, FALSE, FALSE)))
#' oks <- mzion:::are_ovlap_ranges(sta_b = 440, end_b = 600, stas_a = c(412, 775, 850), ends_a = c(455, 843, 900), margin = 5)
#' # stopifnot(identical(oks, c(TRUE, FALSE, FALSE)))
#' oks <- mzion:::are_ovlap_ranges(sta_b = 440, end_b = 850, stas_a = c(412, 775, 850), ends_a = c(455, 843, 900), margin = 5)
#' # stopifnot(identical(oks, c(TRUE, TRUE, TRUE)))
are_ovlap_ranges <- function (sta_b, end_b, stas_a, ends_a, margin = 0L)
{
  ok_lwrs <- sta_b < ends_a + margin & end_b > ends_a # b-start intersect with a-end
  ok_uprs <- sta_b < stas_a & end_b + margin > stas_a # b-end intersect with a-start
  
  ok_lwrs | ok_uprs
}


#' Merge adjacent gates between columns ka and kb
#' 
#' Allow overlaps.
#' 
#' @param xs_ka A vector of X values under column ka.
#' @param ys_ka A vector of Y values under column ka.
#' @param xs_kb A vector of X values under column kb.
#' @param ys_kb A vector of Y values under column kb.
#' @param rngs_ka A vector of scan range values under column ka.
#' @param rngs_kb A vector of scan range values under column kb.
#' @param rts_ka A vector of retention time values under column ka.
#' @param rts_kb A vector of retention time values under column kb.
#' @param scans_ka A vector of scan number values under column ka.
#' @param scans_kb A vector of scan number values under column kb.
#' @param margin The margin for considering being intersecting.
mergeAdjGates2 <- function (xs_ka, ys_ka, xs_kb, ys_kb, rngs_ka, rngs_kb, 
                            rts_ka, rts_kb, scans_ka, scans_kb)
{
  if (!(n_ranges <- length(scans_kb))) {
    return(list(xs_ka = xs_ka, ys_ka = ys_ka, xs_kb = xs_kb, ys_kb = ys_kb, 
                rngs_ka = rngs_ka, rngs_kb = rngs_kb, rts_ka = rts_ka, 
                rts_kb = rts_kb, scans_ka = scans_ka, scans_kb = scans_kb, 
                margin = 0L, merged = FALSE))
  }
  
  ends_ka <- 
    .Internal(unlist(lapply(rngs_ka, function (x) x[length(x)]), 
                     recursive = FALSE, use.names = FALSE))
  stas_ka <- 
    .Internal(unlist(lapply(rngs_ka, `[[`, 1L), 
                     recursive = FALSE, use.names = FALSE))
  stas_kb <- 
    .Internal(unlist(lapply(rngs_kb, `[[`, 1L), 
                     recursive = FALSE, use.names = FALSE))
  ends_kb <- 
    .Internal(unlist(lapply(rngs_kb, function (x) x[length(x)]), 
                     recursive = FALSE, use.names = FALSE))
  
  psa <- vector("list",    n_ranges)
  psb <- vector("integer", n_ranges)
  merged <- FALSE
  
  for (i in 1:n_ranges) {
    sta_bi <- stas_kb[[i]]
    end_bi <- ends_kb[[i]]
    ioks_a <- are_ovlap_ranges(sta_b = sta_bi, end_b = end_bi, stas_a = stas_ka, 
                               ends_a = ends_ka, margin = margin)
    ioks_a <- which(ioks_a)

    if (!(na <- length(ioks_a))) {
      next
    }
    
    stas_a <- stas_ka[ioks_a]
    ends_a <- ends_ka[ioks_a]
    sta_a1 <- stas_a[[1]]
    end_a1 <- ends_a[[1]]
    sta_an <- stas_a[[na]]
    end_an <- ends_a[[na]]
    
    rng_b <- rngs_kb[[i]]
    ib1 <- rng_b[[1]] 
    ibn <- rng_b[[length(rng_b)]]
    yb1 <- ys_kb[[ib1]]
    
    iok1 <- ioks_a[[1]]
    iokn <- ioks_a[[na]]
    ia1 <- rngs_ka[[iok1]][[1]]
    tmp <- rngs_ka[[iokn]]
    ian <- tmp[[length(tmp)]]
    
    # does this before updating, the first index of the OK gates
    rngs_a1 <- unlist(lapply(rngs_ka[ioks_a], `[[`, 1L))
    rngs_ka[[iokn]] <- rgx <- min(ia1, ib1):max(ian, ibn) 
    xsa <- xs_ka[rngs_a1]
    ysa <- ys_ka[rngs_a1]
    ys_ka[rgx] <- sum(yb1, ysa, na.rm = TRUE)
    xs_ka[rgx] <- mean(xsa, na.rm = TRUE) # sum(xsa) / na
    
    imax <- which.max(c(yb1, ysa))
    if (imax == 1L) {
      rts_ka[[iokn]] <- rts_kb[[i]]
      scans_ka[[iokn]] <- scans_kb[[i]]
    }
    else {
      i_best_a <- ioks_a[[imax - 1L]]
      rts_ka[[iokn]] <- rts_ka[[i_best_a]]
      scans_ka[[iokn]] <- scans_ka[[i_best_a]]
    }
    
    psa[[i]] <- if (na > 1L) iok1:(iokn - 1L) else integer()
    xs_kb[rng_b] <- NA_real_
    ys_kb[rng_b] <- NA_real_
    psb[[i]] <- i
  }
  
  if (any(psb)) {
    rts_kb   <- rts_kb[-psb]
    scans_kb <- scans_kb[-psb]
    rngs_kb  <- rngs_kb[-psb]
    merged <- TRUE
  }
  
  if (length(psa <- unlist(psa))) {
    rts_ka   <- rts_ka[-psa]
    scans_ka <- scans_ka[-psa]
    rngs_ka  <- rngs_ka[-psa]
    merged <- TRUE
  }
  
  list(xs_ka = xs_ka, ys_ka = ys_ka, xs_kb = xs_kb, ys_kb = ys_kb, 
       rngs_ka = rngs_ka, rngs_kb = rngs_kb, rts_ka = rts_ka, 
       rts_kb = rts_kb, scans_ka = scans_ka, scans_kb = scans_kb, 
       merged = merged)
}


#' Merge adjacent gates between columns ka and kb
#' 
#' Exact tangent gates for merging.
#' 
#' @param xs_ka A vector of X values under column ka.
#' @param ys_ka A vector of Y values under column ka.
#' @param xs_kb A vector of X values under column kb.
#' @param ys_kb A vector of Y values under column kb.
#' @param rngs_ka A vector of scan range values under column ka.
#' @param rngs_kb A vector of scan range values under column kb.
#' @param rts_ka A vector of retention time values under column ka.
#' @param rts_kb A vector of retention time values under column kb.
#' @param scans_ka A vector of scan number values under column ka.
#' @param scans_kb A vector of scan number values under column kb.
mergeAdjGates <- function (xs_ka, ys_ka, xs_kb, ys_kb, rngs_ka, rngs_kb, 
                           rts_ka, rts_kb, scans_ka, scans_kb)
{
  if (!(n_ranges <- length(scans_kb))) {
    return(list(xs_ka = xs_ka, ys_ka = ys_ka, xs_kb = xs_kb, ys_kb = ys_kb, 
                rngs_ka = rngs_ka, rngs_kb = rngs_kb, rts_ka = rts_ka, 
                rts_kb = rts_kb, scans_ka = scans_ka, scans_kb = scans_kb, 
                merged = FALSE))
  }
  
  ends_ka <- 
    .Internal(unlist(lapply(rngs_ka, function (x) x[length(x)]), 
                     recursive = FALSE, use.names = FALSE))
  stas_ka <- 
    .Internal(unlist(lapply(rngs_ka, `[[`, 1L), 
                     recursive = FALSE, use.names = FALSE))
  stas_kb <- 
    .Internal(unlist(lapply(rngs_kb, `[[`, 1L), 
                     recursive = FALSE, use.names = FALSE))
  ends_kb <- 
    .Internal(unlist(lapply(rngs_kb, function (x) x[length(x)]), 
                     recursive = FALSE, use.names = FALSE))
  
  mt_stas <- match(stas_kb, ends_ka + 1L)
  mt_ends <- match(ends_kb, stas_ka - 1L)
  
  psa <- psb <- vector("integer", n_ranges)
  merged <- FALSE
  
  for (r in 1:n_ranges) {
    mt_stai <- mt_stas[[r]]
    mt_endi <- mt_ends[[r]]
    ok_sta  <- !is.na(mt_stai)
    ok_end  <- !is.na(mt_endi)
    
    if (ok_sta && ok_end) {
      rng_b <- rngs_kb[[r]]
      ra_bf <- rngs_ka[[mt_stai]]
      ra_af <- rngs_ka[[mt_endi]]
      
      ib1 <- rng_b[[1]]
      ia1 <- ra_bf[[1]]
      ian <- ra_af[[length(ra_af)]]
      yb1 <- ys_kb[[ib1]]
      ya1 <- ys_ka[[ia1]]
      yan <- ys_ka[[ian]]
      
      imax <- which.max(c(ya1, yb1, yan))
      
      rngs_ka[[mt_endi]] <- rgx <- ia1:ian
      ys_ka[rgx] <- yb1 + ya1 + yan
      xs_ka[rgx] <- (xs_ka[[ia1]] + xs_ka[[ian]]) / 2
      
      if (imax == 2L) {
        rts_ka[[mt_endi]] <- rts_kb[[r]]
        scans_ka[[mt_endi]] <- scans_kb[[r]]
      }
      else if (imax == 1L) {
        rts_ka[[mt_endi]] <- rts_ka[[mt_stai]]
        scans_ka[[mt_endi]] <- scans_ka[[mt_stai]]
      }
      
      psa[[r]] <- mt_stai
    }
    else if (ok_sta) {
      rng_b <- rngs_kb[[r]]
      ra_bf <- rngs_ka[[mt_stai]]
      
      ib1 <- rng_b[[1]]
      ibn <- rng_b[[length(rng_b)]]
      
      if (mt_stai < length(rngs_ka)) {
        ra_af <- rngs_ka[[mt_stai + 1L]]
        ok_af <- ibn < ra_af[[1]]
      }
      else {
        ok_af <- TRUE
      }
      
      if (!ok_af) {
        next
      }

      ia1 <- ra_bf[[1]]
      yb1 <- ys_kb[[ib1]]
      ya1 <- ys_ka[[ia1]]
      
      imax <- which.max(c(ya1, yb1))
      
      rngs_ka[[mt_stai]] <- rgx <- ia1:ibn
      ys_ka[rgx] <- yb1 + ya1
      xs_ka[rgx] <- xs_ka[[ia1]]
      
      if (imax == 2L) {
        rts_ka[[mt_stai]] <- rts_kb[[r]]
        scans_ka[[mt_stai]] <- scans_kb[[r]]
      }
    }
    else if (ok_end) {
      rng_b <- rngs_kb[[r]]
      ra_af <- rngs_ka[[mt_endi]]
      
      ib1 <- rng_b[[1]]
      ibn <- rng_b[[length(rng_b)]]
      
      if (mt_endi > 1L) {
        ra_bf <- rngs_ka[[mt_endi - 1L]]
        ok_bf <- ib1 > ra_bf[[length(ra_bf)]]
      }
      else {
        ok_bf <- TRUE
      }
      
      if (!ok_bf) {
        next
      }

      ian <- ra_af[[length(ra_af)]]
      yb1 <- ys_kb[[ib1]]
      yan <- ys_ka[[ian]]
      
      imax <- which.max(c(yb1, yan))
      
      rngs_ka[[mt_endi]] <- rgx <- ib1:ian
      ys_ka[rgx] <- yb1 + yan
      xs_ka[rgx] <- xs_ka[[ian]]
      
      if (imax == 1L) {
        rts_ka[[mt_endi]] <- rts_kb[[r]]
        scans_ka[[mt_endi]] <- scans_kb[[r]]
      }
    }
    else {
      next
    }
    
    xs_kb[rng_b] <- NA_real_
    ys_kb[rng_b] <- NA_real_
    psb[[r]] <- r
  }
  
  if (any(psb)) {
    rts_kb   <- rts_kb[-psb]
    scans_kb <- scans_kb[-psb]
    rngs_kb  <- rngs_kb[-psb]
    merged <- TRUE
  }
  
  if (any(psa)) {
    rts_ka   <- rts_ka[-psa]
    scans_ka <- scans_ka[-psa]
    rngs_ka  <- rngs_ka[-psa]
    merged <- TRUE
  }
  
  list(xs_ka = xs_ka, ys_ka = ys_ka, xs_kb = xs_kb, ys_kb = ys_kb, 
       rngs_ka = rngs_ka, rngs_kb = rngs_kb, rts_ka = rts_ka, 
       rts_kb = rts_kb, scans_ka = scans_ka, scans_kb = scans_kb, 
       merged = merged)
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


