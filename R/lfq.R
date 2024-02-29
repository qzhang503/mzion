#' Subsets full MS by the universe of monoisotopic moverzs.
#' 
#' Some non-monoisotopic moverzs can be removed.
#' 
#' @param xs Vectors of moverzs values.
#' @param ys Vectors of intensity values.
#' @param unv The universe.
#' @param from The starting mozerz for calculating bin indexes.
#' @param step The size of bins.
#' @importFrom fastmatch %fin%
subMSfull <- function (xs, ys, unv, from = 200L, step = 1E-5)
{
  len <- length(xs)
  yout <- xout <- vector("list", len)
  
  unv <- unlist(unv, recursive = FALSE, use.names = FALSE)
  unv <- index_mz(unv, from, step)
  unv <- sort(unique(unv))
  
  # for each xi, keeps entries found in unv.
  for (i in 1:len) {
    xi <- xs[[i]]
    li <- length(xi)
    
    if (!li)
      next
    
    ix <- as.integer(ceiling(log(xi/from)/log(1+step)))
    ps0 <- fastmatch::fmatch(ix, unv)
    ps1 <- fastmatch::fmatch(ix + 1L, unv)
    ps2 <- fastmatch::fmatch(ix - 1L, unv)
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
  
  list(x = xout, y = yout)
}


#' Splits \code{df} for parallel tracing.
#' 
#' @param df A data frame of both MS1 and MS2.
#' @param from The starting point for mass binning.
#' @param step The step size for mass binning.
#' @param n_chunks The number of chunks.
#' @param gap The size of gapped entries.
#' @inheritParams matchMS
prep_traceXY <- function (df, from = 200L, step = 1e-5, n_chunks = 4L, 
                          gap = 128L, n_dia_scans = 4L)
{
  cols1 <- c("ms1_mass", "ms1_moverz", "ms1_int", "ms1_charge", 
             "msx_moverzs", "msx_ints", "msx_charges", "orig_scan")
  rows1 <- df[["ms_level"]] == 1L
  df1 <- df[rows1, cols1]
  rm(list = "cols1")
  
  ## Remove non-essential MS2 xyz values
  # try subset by the sub_unv of +/2 mins to reduce 13C peak interference...
  ans_bins <- subMSfull(df1$msx_moverzs, df1$msx_ints, df1$ms1_moverz, 
                        from = from, step = step)
  df1$msx_moverzs <- ans_bins$x
  df1$msx_ints <- ans_bins$y
  rm(list = "ans_bins")
  gc()
  
  # if (n_chunks > 1L) {}

  df1s <- chunksplit(df1, n_chunks, type = "row")
  rm(list = "df1")
  end1s <- cumsum(lapply(df1s, nrow))
  sta1s <- c(1L, end1s[1:(n_chunks-1L)] + 1L)

  ## Adds gaps
  gaps <- lapply(df1s, function (x) floor(min(gap, nrow(x)/2L)))
  df1s_bf <- df1s_af <- vector("list", n_chunks)
  
  for (i in 2:n_chunks) {
    df1s_bf[[i]] <- head(df1s[[i]], gaps[[i]])
  }
  for (i in 1:(n_chunks - 1L)) {
    df1s_af[[i]] <- tail(df1s[[i]], gaps[[i]])
  }
  for (i in 1:n_chunks) {
    df1s[[i]] <- dplyr::bind_rows(df1s_bf[[i]], df1s[[i]], df1s_af[[i]])
  }
  rm(list = c("df1s_bf", "df1s_af"))
  
  ## 
  if (n_chunks > 2L) {
    types <- c("first", rep("middle", n_chunks - 2L), "last")
  } else {
    types <- c("first", "last")
  }
  
  ##  Splits df, sta1s, end1s, n_chunks
  pos_levs <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")
  
  ms1_stax <- ms2_stax <- ms2_endx <- dfs <- vector("list", n_chunks)
  
  for (i in 1:n_chunks) {
    stai <- sta1s[[i]]
    endi <- end1s[[i]]
    ms1_stax[[i]] <- ms1_stas[stai:endi]
    ms2_stax[[i]] <- ms2_stas[stai:endi]
    ms2_endx[[i]] <- ms2_ends[stai:endi]
  }
  rm(list = c("stai", "endi"))
  
  cols <- c("ms_level", "ms1_moverz", "ms1_int")
  for (i in 1:(n_chunks - 1L)) {
    rowx <- ms1_stax[[i]][[1]]:(ms1_stax[[i+1]][[1]] - 1L) # ms2_endx may be NA
    dfs[[i]] <- df[rowx, cols]
  }
  dfs[[n_chunks]] <- df[ms1_stax[[n_chunks]][[1]]:nrow(df), cols]
  rm(list = c("rowx", "cols", "ms1_stax", "ms2_stax", "ms2_endx"))

  list(dfs = dfs, df1s = df1s, gaps = gaps, types = types)
}


#' Helper of \link{traceXY}.
#' 
#' @param xs Vectors of MS1 moverzs.
#' @param ys Vectors of MS1 intensities.
#' @param df A data frame of both MS1 and MS2.
#' @param gap A gap size.
#' @param type The type of data subtype.
#' @param from The starting point for mass binning.
#' @param step The step size for mass binning.
#' @inheritParams matchMS
htraceXY <- function (xs, ys, df, gap = 128L, type = c("first", "middle", "last"), 
                      n_dia_scans = 4L, from = 200L, step = 1E5)
{
  mat <- traceXY(xs = xs, ys = ys, 
                 n_dia_scans = n_dia_scans, from = from, step = step, 
                 reord = FALSE, cleanup = FALSE, replace_ms1_by_apex = TRUE)
  matx <- mat[["x"]]
  maty <- mat[["y"]]
  rm(list = "mat")
  gc()

  if (type == "first") {
    stai <- 1L
    endi <- nrow(matx) - gap
    matx <- matx[stai:endi, ]
    maty <- maty[stai:endi, ]
  } else if (type == "last") {
    stai <- gap + 1L
    endi <- nrow(matx)
    matx <- matx[stai:endi, ]
    maty <- maty[stai:endi, ]
  } else {
    stai <- gap + 1L
    endi <- nrow(matx) - gap
    matx <- matx[stai:endi, ]
    maty <- maty[stai:endi, ]
  }
  
  updateMS1Int(df = df, matx = matx, maty = maty, from = from, step = step)
}


#' Helper of MS1 tracing.
#'
#' @param xs moverzs.
#' @param ys intensities.
#' @param step Step size.
#' @param from The starting point for mass binning.
#' @param step A step size for mass binning.
#' @param reord Logical; re-order data or not.
#' @param cleanup Logical; cleans up xs, ys and zs or not. Set the value to
#'   FALSE to maintain one-to-one correspondence between input (data frame) and
#'   the outputs. This will help, e.g., keep track of scan numbers in the input.
#' @param replace_ms1_by_apex Logical; if TRUE, fill all entries within a gate
#'   by its apex values.
#' @inheritParams matchMS
traceXY <- function (xs, ys, n_dia_scans = 4L, from = 115L, step = 1E-5, 
                     reord = TRUE, cleanup = TRUE, replace_ms1_by_apex = FALSE)
{
  lens <- lengths(xs)
  
  if (reord) {
    ords <- lapply(xs, order)
    
    for (i in seq_along(xs)) {
      if (lens[[i]]) {
        ordi <- ords[[i]]
        xs[[i]] <- xs[[i]][ordi]
        ys[[i]] <- ys[[i]][ordi]
      }
    }
    rm(list = "ords", "ordi")
  }
  
  ## collapses MS data by the indexes of mass bins; 
  # matrix outputs; rows: scans; columns: masses at "x", intensities at "y" ...
  
  # xs can be numeric(0)?
  # RAM expensive step...
  ans <- collapse_mms1ints(xs = xs, ys = ys, lwr = from, step = step, 
                           reord = FALSE, coll = FALSE, cleanup = TRUE)
  ansx <- ans[["x"]]
  ansy <- ans[["y"]]
  rm(list = c("ans"))
  gc()
  
  ## traces MS data matrices across LC scans; rows: scans; columns: masses
  nrc <- dim(ansy)
  nr <- nrc[[1]]
  nc <- nrc[[2]]
  rm(list = "nrc")
  
  xmat <- ymat <- matrix(rep_len(NA_real_, nc * nr), ncol = nc)
  ranges <- apexes <- ns <- vector("list", nc)
  
  if (replace_ms1_by_apex) {
    for (i in 1:nc) {
      gates <- find_lc_gates(ansy[, i], n_dia_scans = n_dia_scans)
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
      gates <- find_lc_gates(ansy[, i], n_dia_scans = n_dia_scans)
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
  pos_levs <- getMSrowIndexes(df$ms_level, pad_nas = TRUE)
  ms1_stas <- pos_levs$ms1_stas
  ms2_stas <- pos_levs$ms2_stas
  ms2_ends <- pos_levs$ms2_ends
  rm(list = "pos_levs")

  for (i in seq_along(ms2_stas)) {
    ms2sta <- ms2_stas[[i]]
    
    if (is.na(ms2sta))
      next
    
    xs <- matx[i, ]
    ys <- maty[i, ]
    oks <- .Internal(which(!is.na(xs)))
    xs <- xs[oks]
    ys <- ys[oks]
    ixs <- as.integer(ceiling(log(xs/from)/log(1+step)))
    
    ms2end <- ms2_ends[[i]]
    df2 <- df[ms2sta:ms2end, ]
    
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


