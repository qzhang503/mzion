#' De-isotopes precursor masses.
#'
#' @param moverzs Mass-to-charge ratios.
#' @param msxints MS1 or MS2 peak intensities.
#' @param center The mass center of an isolation window (only for MS1).
#' @param ppm Allowance in mass error when deisotoping.
#' @param offset_upr A cardinal number of upper mass off-sets.
#' @param offset_lwr A cardinal number of lower mass off-sets.
#' @param bound Not yet used. Logical; if TRUE, removes precursors outside of
#'   the boundary of isolation window.
#' @param ms_lev MS level.
#' @param maxn_feats The maximum number of MS features.
#' @param max_charge The maximum charge state.
#' @param step Step size for mass binning.
#' @param order_mz Logical; if TRUE, orders peaks from low to high m-over-z's.
#' @param backward_mass_co A mass cut-off to initiate backward looking of
#'   an isotope envelop.
#' @examples
#' \donttest{
#' library(mzion)
#' moverzs <- c(881 + 1:10*.1, 882.0674, 882.0981, 882.4034, 882.60, 882.7372)
#' msxints <- c(1000 * 1:10, 1652869, 882.0981, 2043015, 2314111, 4314111)
#'
#' # out <- mzion:::deisotope(moverzs, msxints, ppm = 10L, ms_lev = 1L, 
#' #                          maxn_feats = 5L, max_charge = 4L, offset_upr = 8L,
#' #                          offset_lwr = 8L, order_mz = TRUE, bound = FALSE)
#' }
deisotope <- function (moverzs, msxints, center = 650.0, ppm = 6L, ms_lev = 1L, 
                       maxn_feats = 300L, max_charge = 4L, 
                       offset_upr = 30L, offset_lwr = 30L, order_mz = TRUE, 
                       step = ppm/1e6, backward_mass_co = 800/ms_lev, 
                       bound = FALSE)
{
  ###
  # if to apply intensity cut-offs, should note the difference intensity 
  # scales between Thermo's and Bruker's data
  ###
  
  # null_out <- list(masses = NA_real_, charges = NA_integer_, intensities = NA_real_)
  null_out <- list(masses = NULL, charges = NULL, intensities = NULL)
  
  if (!(len <- length(moverzs)))
    return(null_out)
  
  # if (len == 1L && is.na(moverzs))
  #   return(null_out)

  from <- moverzs[[1]] - .001
  ims  <- index_mz(moverzs, from, step)
  
  lenp <- min(len, maxn_feats)
  peaks <- intens <- rep(NA_real_, lenp)
  css <- rep(NA_integer_, lenp)
  p <- 1L
  
  while(len & p <= maxn_feats) {
    if (len == 1L) {
      intens[[p]] <- msxints
      peaks[[p]] <- moverzs
      css[[p]] <- 0L
      p <- p + 1L
      len <- len -1L
      next
    }
    
    imax <- if (p == 1L && ms_lev == 1L)
      .Internal(which.min(abs(moverzs - center)))
    else
      .Internal(which.max(msxints))

    mass <- moverzs[imax]
    mint <- msxints[imax]

    for (ch in max_charge:0) {
      if (ch == 0L) {
        peaks[[p]] <- mass
        intens[[p]] <- mint
        p <- p + 1L
        moverzs <- moverzs[-imax]
        msxints <- msxints[-imax]
        ims <- ims[-imax]
        len <- len - 1L
      }
      else {
        gap <- 1.003355/ch
        
        # forward looking up to 10 mass entries
        mx  <- mass + gap
        sta <- min(len, imax + 1L)
        end <- min(len, imax + 10L)
        oks <- abs((mx - moverzs[sta:end])/mx) * 1E6 <= ppm
        
        if (any(oks)) {
          # consider 13C off-sets at mass > backward_mass_co 
          if (ch * mass > backward_mass_co) {
            iths <- index_mz(mass + gap * (1L - ch):(1 + ch), from, step)
            if ((lwr <- imax - offset_lwr) < 1L) lwr <- 1L
            if ((upr <- imax + offset_upr) > len) upr <- len
            
            iexs <- ims[lwr:upr]
            oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
            hits <- .Internal(which(oks2)) + lwr - 1L
            
            # backward looking
            if ((idx <- .Internal(which(hits == imax))) > 1L) {
              for (i in (idx - 1L):1) {
                hi <- hits[[i]]
                
                if (msxints[[imax]]/msxints[[hi]] <= 3)
                  mass <- moverzs[[hi]]
                else
                  break
              }
              
              hits <- hits[i:length(hits)]
              lenh <- length(hits)
            }
            else {
              lenh <- length(hits)
            }
          }
          else {
            iths <- index_mz(mass + gap * 1:ch, from, step)
            iexs <- ims[sta:min(imax + offset_upr, len)]
            oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
            hits <- c(imax, .Internal(which(oks2)) + imax) 
            lenh <- length(hits)
          }
          
          intens[[p]] <- sum(msxints[hits])
          peaks[[p]] <- mass
          css[[p]] <- ch
          p <- p + 1L
          
          len <- len - lenh
          moverzs <- moverzs[-hits]
          msxints <- msxints[-hits]
          ims <- ims[-hits]
          
          break
        }
        else {
          # backward looking up to 10 mass entries
          mx  <- mass - gap
          end <- max(1L, imax - 1L)
          sta <- max(1L, imax - 10L)
          oks <- abs((mx - moverzs[sta:end])/mx) * 1E6 <= ppm
          
          if (any(oks)) {
            iths <- index_mz(mass + gap * -ch:0L, from, step)
            lwr <- max(imax - offset_lwr, 1L)
            upr <- imax
            # if ((lwr <- imax - offset_lwr) < 1L) lwr <- 1L
            # if ((upr <- imax + offset_upr) > len) upr <- len
            iexs <- ims[lwr:upr]
            oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
            hits <- .Internal(which(oks2)) + lwr - 1L
            lenh <- length(hits)
            
            intens[[p]] <- sum(msxints[hits])
            peaks[[p]] <- mass <- moverzs[hits[[1]]]
            css[[p]] <- ch
            p <- p + 1L
            
            len <- len - lenh
            moverzs <- moverzs[-hits]
            msxints <- msxints[-hits]
            ims <- ims[-hits]
            
            break
          }
        }
      }
    }
  }
  
  if (ms_lev == 1L)
    oks1 <- css > 1L & !is.na(css)
  else if (ms_lev == 2L)
    oks1 <- !is.na(peaks)
  else
    stop("Unhandled MS level = ", ms_lev)
  
  masses <- peaks[oks1]
  charges <- css[oks1]
  intensities <- intens[oks1]
  
  if (order_mz) {
    ord <- order(masses)
    masses <- masses[ord]
    charges <- charges[ord]
    intensities <- intensities[ord]
  }
  
  out <- list(masses = masses, charges = charges, intensities = intensities)
}


#' Detection of chromatographic peaks.
#' 
#' May be overkilling.
#' 
#' https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
#' https://stackoverflow.com/questions/58942623/isolating-peaks-in-from-time-series-data-in-r
#' 
#' @param y Values of y.
#' @param lag The lag length.
#' @param threshold The threshold.
#' @param influence The influence.
#' @param min_len The minimum length of y for peak-picking algorithm.
#' @param min_lag The minimum length of the lag.
#' @param max_lag The maximum length of the lag.
#' @param min_var The minimum variance.
#' @param static Is the time series static (or adaptive).
find_lcpeaks <- function (y, lag = 5L, min_lag = 3L, max_lag = 200L, 
                          min_len = 10L, min_var = 1E-5, static = TRUE,
                          threshold = 3, influence = 0) 
{
  if (min_len < min_lag)
    stop("`min_len` cannot be smaller than `min_lag`.")
  
  if (lag < min_lag)
    stop("`lag` cannot be smaller than `min_lag`.")
  
  if (lag > max_lag)
    stop("`lag` cannot be greater than `max_lag`.")
  
  if ((len <- length(y)) <= min_len) {
    # return(list(signals = NA_integer_, avgFilter = NA_real_, stdFilter = NA_real_))
    return(rep.int(1L, len))
  }
  
  if (lag >= (len2 <- len * .15)) 
    lag <- min(max_lag, ceiling(len2))
  else
    lag <- min(max_lag, lag)
  
  n <- lag - 1L
  n2 <- lag * n
  filteredY <- y
  
  # global (if static)
  # also need a default, static gl_va and gl_sd...
  gl_su <- sum(y, na.rm = TRUE)
  gl_me <- gl_su/len
  
  if ((gl_va <- sum(y^2, na.rm = TRUE)/len - gl_su^2/len/(len - 1L)) < 0) {
    gl_va <- var(y)
    gl_sd <- sqrt(gl_va)
  }
  else {
    gl_sd <- sqrt(gl_va) * .10
  }
  
  # initialization
  signals <- rep(0L, len)
  va_y <- su_y <- ss_y <- stdFilter <- avgFilter <- rep(NA_real_, len)
  
  i <- lag
  ys <- y[1:i]
  ss_y[i] <- ss_yi <- sum(ys^2, na.rm = TRUE)
  su_y[i] <- su_yi <- sum(ys, na.rm = TRUE)
  avgFilter[i] <- su_me <- su_yi/lag
  
  if ((va_yi <- ss_yi/n - (su_yi)^2/n2) < min_var) 
    va_yi <- gl_va
  
  if ((sd_yi <- sqrt(va_yi)) < gl_sd)
    sd_yi <- gl_sd
  
  va_y[i] <- va_yi
  stdFilter[i] <- sd_yi
  
  # rolling
  for (i in (lag + 1L):len) {
    p <- i - 1L
    l <- i - lag + 1L
    yi <- y[i]
    d <- yi - avgFilter[p]
    
    if (abs(d) > threshold * stdFilter[p]) {
      signals[i] <- if (d > 0) 1L else -1L
      filteredY[i] <- influence * yi + (1 - influence) * filteredY[p]
    }
    
    ss_y[i] <- ss_yi <- ss_y[p] + filteredY[i]^2 - filteredY[l]^2
    su_y[i] <- su_yi <- su_y[p] + filteredY[i] - filteredY[l]
    avgFilter[i] <- su_yi/lag
    
    if ((va_yi <- ss_yi/n - (su_yi)^2/n2) < 1E-5) 
      va_yi <- 0
    
    if ((sd_yi <- sqrt(va_yi)) < gl_sd)
      sd_yi <- gl_sd
    
    va_y[i] <- va_yi
    stdFilter[i] <- sd_yi
  }
  
  # list(signals = signals, avgFilter = avgFilter, stdFilter = stdFilter)
  signals
}


#' Chunks of precursor lists by continous retention-time bins.
#' 
#' Not used.
#' 
#' @param df A data frame of precursors.
sep_ms1rts <- function (df)
{
  dfs <- split(df, df$irt) 
  len <- length(dfs)
  frs <- as.integer(names(dfs))
  brs <- which(diff(frs) > 1L)
  
  ans <- vector("list", length(brs) + 1L)
  
  sta <- 0L
  
  for (i in seq_along(brs)) {
    end <- brs[i]
    ans[[i]] <- dfs[(sta + 1L):end]
    sta <- end
  }
  
  ans[[length(ans)]] <- dfs[(sta + 1L):len]
  
  ans
}


#' Three-frame binning of precursor retention times.
#' 
#' Not used.
#' 
#' @param df A data frame of precursors.
bin_ms1rts <- function (df)
{
  len <- length(df)
  
  if (len <= 3L) {
    return(dplyr::bind_rows(df))
  }
  
  ans <- vector("list", len_out <- len - 2L)
  
  for (i in 1:len_out) {
    ans[[i]] <- dplyr::bind_rows(df[[i]], df[[i+1L]], df[[i+2L]])
  }
  
  ans
}


#' Generate an isotope envelop
#' 
#' @param mass Precursor mass.
#' @param charge Charge state.
gen_isoenvlope <- function (mass = 560, charge = 2L)
{
  # averagine: C 4.9384; H 7.7583; N 1.3577; O 1.4773; S 0.0417
  # 4.9384*12 + 7.7583*1.007825 + 1.3577*14.003074 + 1.4773*15.99491462 + 0.0417*31.9720707
  
  p13c <- .0111   # .9889 + 0.0111
  p2h  <- .000156 # .999844 + .000156
  p18o <- .00187  # .99738 + .00187
  p15n <- .004    # .996 + .004
  p34s <- .043    # .952 + .0425 -> .957 + .0427

  size <- mass / 111
  nc <- round(size*5.1)
  nh <- round(size*7.7)
  no <- round(size*1.4773)
  nn <- round(size*1.3577)
  ns <- round(size*0.0417)
  
  nd <- mass - (nc * 12 + nh * 1.007825 + nn * 14.003074 + no * 15.99491462 + ns * 31.9720707)
  nh <- nh + round(nd/1.007825)
  
  # C19 H25 O5 N5 S3
  dc <- dbinom(0:5, nc, p13c)
  dh <- dbinom(0:5, nh, p2h)
  do <- dbinom(0:5, no, p18o)
  dn <- dbinom(0:5, nn, p15n)
  ds <- dbinom(0:5, ns, p34s)
  
  dc <- dc/dc[[1]]
  dh <- dh/dh[[1]]
  do <- do/do[[1]]
  dn <- dn/dn[[1]]
  ds <- ds/ds[[1]]
  
  ## +1
  p13c * nc # 0.2109; no collapsing at hi-res MS, good approximation
  p13c * nc + p2h * nh + p15n * nn # .235; no collapsing unless with low-res MS

  # no need to incur the computationally more involved version, 
  #   as MS isotope envelope from a single scan probably not have the accuracy
  xc <- nc * p13c * (1 - p13c)^(nc - 1L)/(1 - p13c)^nc
  xh <- nh * p2h  * (1 -  p2h)^(nh - 1L)/(1 - p2h)^nh
  xn <- nn * p15n * (1 - p15n)^(nn - 1L)/(1 - p15n)^nn
  xc + xh + xn # .237
  
  ## +2
  # Xcalibur 9.14E4/6.74E5 # .136
  
  if (ns) {
    # uses the bigger
    ds[[2]]
    dc[[3]]
  }
  else {
    dc[[3]]
  }

  ## +3
  # .111 * .889^18 * 19 # .253
  
  # + p34s * ns
  # p13c * (nc + 1L) + p2h * nh + p18o * no + p15n * nn + p34s * ns
  
  # +2 can be 2*13c, 13c + 2h, 
  
  
  # C19 H25 O5 N5 S3
  nc <- 19; nh <- 25; no <- 5; nn <- 5; ns <- 3
}


