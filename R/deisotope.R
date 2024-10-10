#' De-isotopes precursor masses.
#'
#' @param moverzs MS1 or MS2 Mass-to-charge ratios. The inputs are typically at
#'   weighted-mean statistics.
#' @param msxints MS1 or MS2 peak intensities. The inputs are typically at mean
#'   statistics.
#' @param n_ms1s The underlying counts that have contributed to \code{moverzs}
#'   and \code{maxints}.
#' @param center The mass center of an isolation window (only for MS1).
#' @param is_dda Logical; is DDA or not.
#' @param ppm Allowance in mass error when deisotoping.
#' @param offset_upr A cardinal number of upper mass off-sets in relative to a
#'   mass center when finding mono-isotopic species.
#' @param offset_lwr A cardinal number of lower mass off-sets in relative to a
#'   mass center when finding mono-isotopic species.
#' @param ms_lev MS level.
#' @param maxn_feats The maximum number of MS features.
#' @param max_charge The maximum charge state.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries when determining
#'   the charge state.
#' @param n_bwd Backward looking up to \code{n_bwd} mass entries when
#'   determining the charge state.
#' @param step Step size for mass binning.
#' @param backward_mass_co A mass cut-off to initiate backward looking of an
#'   isotope envelop.
#' @param fct_iso2 The multiplication factor for the second isotopic peak.
#' @param is_pasef Logical; is TIMS TOF data or not.
#' @param min_y Not yet used. The minimum intensity for consideration of
#'   deisotoping (guard against PASEF data). The setting should be MS platform
#'   dependent, e.g., a much smaller value with PASEF. Also note different
#'   settings between MS1 and MS2.
#' @param min_ratio A ratio threshold. Exit chimeric deisotoping if the Y_cr /
#'   Y_bf is smaller the threshold.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' moverzs <- c(881 + 1:10*.1, 882.0674, 882.0981, 882.4034, 882.60, 882.7372)
#' msxints <- c(1000 * 1:10, 1652869, 788368, 2043015, 2314111, 4314111)
#' n_ms1s <- rep_len(1L, length(msxints))
#'
#' if (FALSE) {
#'   out <- mzion:::find_ms1stat(moverzs, msxints, n_ms1s, ppm = 5L,
#'                               ms_lev = 1L, maxn_feats = 5L, max_charge = 4L,
#'                               offset_upr = 8L, offset_lwr = 8L)
#' }
#'
#' moverzs <- c(907.968018,908.136475,908.872711,909.073690,909.121,
#'              909.273670,909.280212,909.45,909.454981,909.459717,
#'              909.473072,909.675781,909.721191,909.78,909.79,
#'              909.791633,909.840088,909.969744,910.077278,910.091888,
#'              910.124076,910.130798,910.180115,910.218689,910.341797,
#'              910.456115,910.463650,910.474152,910.791669,911.126692,
#'              911.133219,911.471605)
#' msxints <- c(2971101,3013171,5774301,24508055,40267873,
#'              1E6,7340374,1E6,154459778,8312258,
#'              18375220,22909676,2863776,1.1E6,1.2E6,
#'              191655971,7273611,32089542,9125108,5008687,
#'              167021373,12912115,3242378,16278636,4731828,
#'              98601412,22828501,10539991,26255569,8499971,
#'              8305187,37899662)
#' n_ms1s <- rep_len(1L, length(msxints))
#'
#' if (FALSE) {
#'   out <- mzion:::find_ms1stat(moverzs, msxints, n_ms1s, center = 909.788696,
#'                       exclude_reporter_region = FALSE,
#'                       ppm = 5L, ms_lev = 1L, maxn_feats = 1L, max_charge = 4L,
#'                       step = 5/1e6)
#'   out <- mzion:::find_ms1stat(moverzs, msxints, n_ms1s, # center = 909.788696,
#'                       exclude_reporter_region = FALSE,
#'                       ppm = 5L, ms_lev = 1L, maxn_feats = 1L, max_charge = 4L,
#'                       step = 5/1e6)
#' }
#' }
find_ms1stat <- function (moverzs, msxints, n_ms1s = 1L, center = 0, 
                          exclude_reporter_region = FALSE, 
                          tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                          is_dda = TRUE, ppm = 8L, ms_lev = 1L, maxn_feats = 300L, 
                          max_charge = 4L, n_fwd = 20L, n_bwd = 20L, 
                          offset_upr = 30L, offset_lwr = 30L, step = ppm / 1e6, 
                          backward_mass_co = 800/ms_lev, grad_isotope = 1.6, 
                          fct_iso2 = 3.0, use_defpeaks = FALSE, 
                          is_pasef = FALSE, min_y = 10, min_ratio = .05)
{
  ###
  # if to apply intensity cut-offs, should note the difference intensity 
  # scales between Thermo's and Bruker's data
  ###
  
  null_out <- list(masses = NULL, charges = NULL, intensities = NULL)
  
  if (!(len_ms <- length(moverzs))) {
    return(null_out)
  }

  excl_rptrs <- exclude_reporter_region && ms_lev != 1L
  
  if (excl_rptrs) {
    if (FALSE) {
      tmts <- c(
        `126` = 126.127726, `127N` = 127.124761, `127C` = 127.131080,
        `128N` = 128.128115, `128C` = 128.134435, `129N` = 129.131470,
        `129C` = 129.137790, `130N` = 130.134825, `130C` = 130.141145,
        `131N` = 131.138180, `131C` = 131.144499, `132N` = 132.141535,
        `132C` = 132.147855, `133N` = 133.14489, `133C` = 133.15121,
        `134N` = 134.148245, `134C` = 134.155114, `135N` = 135.152149)
      
      tmt_lwr <- tmts - tmts * 1e-5
      tmt_upr <- tmts + tmts * 1e-5
      
      ok_rptrs <- mapply(function (x, y, m ) m >= x & m <= y, 
                         tmt_lwr, tmt_upr, MoreArgs = list(moverzs), 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
      ok_rptrs <- Reduce(`|`, oks_tmt)
    }
    
    ok_rptrs <- moverzs > tmt_reporter_lower & moverzs < tmt_reporter_upper
    no_rptrs <- .Internal(which(!ok_rptrs))
    ok_rptrs <- .Internal(which(ok_rptrs))
    
    rptr_moverzs <- moverzs[ok_rptrs]
    rptr_ints <- msxints[ok_rptrs]
    moverzs <- moverzs[no_rptrs]
    msxints <- msxints[no_rptrs]
    
    len_rptrs <- length(rptr_moverzs)
    len_ms <- len_ms - len_rptrs
    
    if (!len_ms) {
      return(null_out)
    }
  }
  
  grad_isotope2 <- fct_iso2 * grad_isotope
  from <- moverzs[[1]] - .001
  ims  <- index_mz(moverzs, from, step)
  
  tol <- ppm / 1E6 # the same as the default `step`; may remove arg `ppm`
  lenp <- min(len_ms, maxn_feats)
  peaks_fuz <- peaks <- intens_fuz <- intens <- rep_len(list(numeric(1)), lenp)
  css_fuz <- css <- rep_len(list(integer(1)), lenp)
  p_no_zero <- p_fuz <- p <- 1L
  ymean_fuz <- ymono_fuz <- mass_fuz <- NA_real_
  yref <- 0 # to keep track of Y values

  # Little gain with the addition guard of p_no_zero; 
  # At MS1 ch == 0 without p_no_zero, ch-0 can be picked up by chimeric searches
  while((len_ms > 0L) && (p_no_zero <= maxn_feats)) {
    if (len_ms == 1L) {
      intens[[p]] <- msxints[[1]]
      peaks[[p]] <- moverzs[[1]]
      css[[p]] <- NA_integer_
      break
    }
    
    if ((!is_pasef) && center > 0 && p == 1L && ms_lev == 1L) {
      imax <- .Internal(which.min(abs(moverzs - center)))
    }
    else {
      imax <- .Internal(which.max(msxints))
    }
    
    mass  <- moverzs[[imax]]
    mint  <- msxints[[imax]]
    n_ms1 <- n_ms1s[[imax]]
    
    # the first `mint` may be the one closest to the center with Thermo's and 
    # the first yref is 0
    if (yref > 0 && mint / yref < min_ratio) {
      break
    }
    
    ## find Z values
    ch <- find_charge_state(
      mass = mass, imax = imax, mint = mint, 
      moverzs = moverzs, msxints = msxints, n_ms1s = n_ms1s,
      max_charge = max_charge, n_fwd = n_fwd, n_bwd = n_bwd, 
      ms_lev = ms_lev, is_dda = is_dda, tol = tol)
    # if (mint < min_y) { ch <- NA_integer_ }

    if (is.na(ch)) {
      peaks[[p]] <- mass
      intens[[p]] <- mint
      len_ms <- len_ms - 1L
      moverzs <- moverzs[-imax] # copy-on-modify
      msxints <- msxints[-imax]
      ims <- ims[-imax]
      p <- p + 1L
      
      # Many MS2 charge states are undetermined (and move on)
      if (ms_lev == 2L) {
        p_no_zero <- p_no_zero + 1L
      }
      
      # yref <- mint # no update of reference intensity
      next
    }
    
    p_no_zero <- p_no_zero + 1L
    gap <- 1.003355 / ch
    sta1 <- min(len_ms, imax + 1L)
    end1 <- min(len_ms, imax + offset_upr)
    # yzero <- msxints[[imax]]
    ymono <- msxints[[imax]]
    ymean <- ymono/n_ms1

    ## find monoisotopic m/z
    if (ch * mass > backward_mass_co) {
      if ((lwr <- imax - offset_lwr) < 1L) lwr <- 1L
      if ((upr <- imax + offset_upr) > len_ms) upr <- len_ms
      
      iths <- index_mz(mass + gap * (-ch-1L):(1L+ch), from, step)
      iexs <- ims[lwr:upr]
      oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
      hits <- .Internal(which(oks2)) + lwr - 1L
      lenh <- length(hits)
      idx  <- .Internal(which(hits == imax))
      
      # only one preceding hit
      if (idx == 2L) {
        hi <- hits[[1]]
        
        if (abs((moverzs[[hi]] + gap)/mass - 1) <= tol) { # not satellite to `mass`
          yhi <- msxints[[hi]]
          ri  <- mint/yhi
          
          # the preceding intensity pass the 1st threshold
          if (ri <= grad_isotope) {
            mass  <- moverzs[[hi]]
            n_ms1 <- n_ms1s[hi]
            ymono <- yhi
            ymean <- ymono/n_ms1
          }
          else {
            # the preceding intensity fail on the 1st but pass the 2nd threshold
            if (ri <= grad_isotope2) {
              mass_fuz  <- moverzs[[hi]]
              n_ms1_fuz <- n_ms1s[hi]
              ymono_fuz <- yhi
              ymean_fuz <- ymono_fuz/n_ms1_fuz
            }
          }
        }
      }
      else if (idx > 2L) { # two or more preceding matches
        ix <- lenx <- idx - 1L
        hx <- NULL
        
        for (i in seq_len(lenx)) {
          mzsub <- moverzs[hits[1:ix]]
          ks <- .Internal(which(abs((mzsub + gap)/mass - 1) <= tol))
          
          # all satellite to `mass`
          if (!length(ks))
            break
          
          hsub <- hits[ks]
          isub <- .Internal(which.max(msxints[hsub]))
          hi <- hsub[[isub]]
          h <- ks[[isub]]
          gi <- grad_isotope * i
          yhi <- msxints[[hi]]
          ri <- mint / yhi
          
          if (ri <= gi) { # the preceding intensity pass the 1st threshold
            ymono <- yhi
            ymean <- ymono / n_ms1 # added on 2024-06-26
          }
          else {
            if (ri <= fct_iso2 * gi) { # the preceding pass the 2st threshold
              mass_fuz  <- moverzs[[hi]]
              ymono_fuz <- yhi
              n_ms1_fuz <- n_ms1s[[hi]]
              ymean_fuz <- ymono_fuz/n_ms1_fuz
              break
            }
            else {
              break
            }
          }
          
          mass <- moverzs[[hi]]
          hx <- c(hi, hx)
          
          if (h == 1L)
            break
          
          ix <- h - 1L
        }
        
        hits <- c(hx, hits[idx:lenh])
        lenh <- length(hits)
      }
    }
    else {
      iths <- index_mz(mass + gap * 1:ch, from, step)
      iexs <- ims[sta1:end1]
      oks1 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
      hits <- c(imax, .Internal(which(oks1)) + imax)
      lenh <- length(hits)
    }
    
    # 2% more with ymean than yzero;
    # refine later: e.g., ymean if the mass in the isolation window else yzero
    # or ymean more advantageous: e.g., at yzero with 13C1 > ymean, 
    #  but scrambling of m/z with 13C; ymean may average out the scrambling.
    # also possible both yzero (e.g., 13C1) and ymono both in the isoWindow 
    #  -> take the mean
    # or may be ymono the best...
    intens[[p]] <- ymean
    # intens[[p]] <- yzero
    peaks[[p]] <- mass
    css[[p]] <- ch
    
    if (!is.na(mass_fuz)) {
      intens_fuz[[p_fuz]] <- ymean_fuz
      peaks_fuz[[p_fuz]] <- mass_fuz
      css_fuz[[p_fuz]] <- ch
      p_fuz <- p_fuz + 1L
      mass_fuz <- NA_real_
    }
    
    len_ms <- len_ms - lenh
    moverzs <- moverzs[-hits] # copy-on-modify
    msxints <- msxints[-hits]
    ims <- ims[-hits]
    p <- p + 1L
    yref <- mint
  }
  
  css <- .Internal(unlist(css, recursive = FALSE, use.names = FALSE))
  css_fuz <- .Internal(unlist(css_fuz, recursive = FALSE, use.names = FALSE))
  intens <- .Internal(unlist(intens, recursive = FALSE, use.names = FALSE))
  intens_fuz <- .Internal(unlist(intens_fuz, recursive = FALSE, use.names = FALSE))
  peaks <- .Internal(unlist(peaks, recursive = FALSE, use.names = FALSE))
  peaks_fuz <- .Internal(unlist(peaks_fuz, recursive = FALSE, use.names = FALSE))

  # outputs
  if (ms_lev == 1L) {
    oks1 <- css >= 1L & !is.na(css)
  }
  else if (ms_lev == 2L) {
    oks1 <- !is.na(peaks)
  }
  else {
    stop("Unhandled MS level = ", ms_lev)
  }
  
  masses <- peaks[oks1]
  charges <- css[oks1]
  intensities <- intens[oks1]
  
  if (endf <- p_fuz - 1L) {
    rng_fuz <- 1:endf
    masses <- c(masses, peaks_fuz[rng_fuz])
    charges <- c(charges, css_fuz[rng_fuz])
    intensities <- c(intensities, intens_fuz[rng_fuz])
  }
  
  if (excl_rptrs && len_rptrs) {
    masses <- c(rptr_moverzs, masses)
    charges <- c(rep.int(NA_integer_, len_rptrs), charges)
    intensities <- c(rptr_ints, intensities)
  }
  
  # De-isotoping from high-to-low intensities
  # For MS2, the sequences need to be from low-to-high m-over-z values
  len_out <- length(masses)
  
  # MS2 need to be ordered from low to high m/z values
  if (len_out > 1L) {
    if (ms_lev == 1L) {
      ord <- .Internal(radixsort(na.last = TRUE, decreasing = TRUE, FALSE, TRUE, 
                                 intensities))
      masses <- masses[ord]
      charges <- charges[ord]
      intensities <- intensities[ord]
    }
    else {
      ord <- .Internal(radixsort(na.last = TRUE, decreasing = FALSE, FALSE, TRUE, masses))
      masses <- masses[ord]
      charges <- charges[ord]
      intensities <- intensities[ord]
    }
  }
  
  out <- list(masses = masses, charges = charges, intensities = intensities)
}


#' Finds charge state.
#'
#' @param mass The mass of a peak for questing its charge state.
#' @param imax The index of the peak, which is often the most intense.
#' @param mint Not yet used. The intensity of the peak at \code{imax}.
#' @param tol The tolerance in mass errors.
#' @inheritParams find_ms1stat
find_charge_state <- function (mass, imax, mint, moverzs, msxints, n_ms1s, 
                               max_charge = 4L, n_fwd = 15L, n_bwd = 20L, 
                               ms_lev = 1L, is_dda = TRUE, tol = 8E-6)
{
  # find results for all ch -> uses the best...
  # distinguish left, right or both left and right evidence
  # if left-only evidence -> monomass must < mass
  
  # stopifnot(max_charge >= 2L, tol >= 1E-6, lenm)
  lenm <- length(moverzs)
  gaps <- 1.003355 / 2:max_charge
  sta1 <- min(lenm, imax + 1L)
  sta2 <- max(1L, imax - n_bwd)
  end1 <- min(lenm, imax + n_fwd)
  end2 <- max(1L, imax - 1L)
  xsub1 <- moverzs[sta1:end1]
  xsub2 <- moverzs[sta2:end2]
  ysub1 <- msxints[sta1:end1]
  ysub2 <- msxints[sta2:end2]
  
  # if (ms_lev == 1L) {
  #   nsub1 <- n_ms1s[sta1:end1]
  #   nsub2 <- n_ms1s[sta2:end2]
  # }
  
  mxs1 <- mass + gaps
  ioks1 <- lapply(mxs1, function (x) .Internal(which(abs(xsub1/x - 1L) <= tol)))
  lens1 <- lengths(ioks1)
  chs1 <- .Internal(which(lens1 > 0L))
  l1 <- length(chs1)
  
  # l2 not use in the case
  if (l1 > 1L) {
    p1s <- lapply(ioks1[chs1], function (xs) {
      if (length(xs) == 1L)
        sta1 + xs - 1L
      else {
        ps <- sta1 + xs - 1L
        i <- .Internal(which.max(msxints[ps]))
        sta1 + xs[i] - 1L
      }
    })
    p1s <- .Internal(unlist(p1s, recursive = FALSE, use.names = FALSE))
    
    if (ms_lev == 1L) {
      ansduo <- check_chduo(duo = c(1L, 3L), chs1, p1s, msxints)
      chs1 <- ansduo$chs
      p1s <- ansduo$ps
      
      if (max_charge >= 6L) {
        ansduo <- check_chduo(duo = c(2L, 5L), chs1, p1s, msxints)
        chs1 <- ansduo$chs
        p1s <- ansduo$ps
      }
    }
    
    if (FALSE && ms_lev == 1L) {
      if (length(chs1) > 1L) {
        chx1 <- chs1 + 1L
        unqs <- !as.integer(chx1 * 2L) %in% chx1
        chs1 <- chs1[unqs]
        p1s <- p1s[unqs]
      }
      
      nsubs <- n_ms1s[p1s]
      
      if (max(nsubs) == min(nsubs)) {
        chs1 <- chs1[which.max(msxints[p1s])]
        return(chs1 + 1L)
      }
      
      chs1 <- chs1[which.max(n_ms1s[p1s])]
      return(chs1 + 1L)
    }
    else {
      chs1 <- chs1[which.max(msxints[p1s])]
      return(chs1 + 1L)
    }
  }
  
  mxs2 <- mass - gaps
  ioks2 <- lapply(mxs2, function (x) .Internal(which(abs(xsub2/x - 1L) <= tol)))
  lens2 <- lengths(ioks2)
  chs2 <- .Internal(which(lens2 > 0L))
  l2 <- length(chs2)
  
  if (l1 == 1L) {
    if (l2 == 0L) {
      return(chs1 + 1L)
    }
    
    if (l2 == 1L) {
      if (chs1 == chs2) {
        return(chs1 + 1L)
      }
      
      p1 <- sta1 + ioks1[[chs1]] - 1L
      p2 <- sta2 + ioks2[[chs2]] - 1L
      yp1 <- max(msxints[p1])
      yp2 <- max(msxints[p2])
      
      if (ms_lev == 1L && is_dda) {
        np1 <- max(n_ms1s[p1])
        np2 <- max(n_ms1s[p2])
        
        if (np1 == np2) {
          if (yp1 >= yp2) return(chs1 + 1L) else return(chs2 + 1L)
        }
        
        if (np1 > np2) return(chs1 + 1L) else return(chs2 + 1L)
      }
      else {
        if (yp1 >= yp2) return(chs1 + 1L) else return(chs2 + 1L)
      }
    }
    
    # l2 > 1L
    if (any(chs2 %in% chs1)) {
      return(chs1 + 1L)
    }
    
    p1 <- sta1 + ioks1[[chs1]] - 1L
    p2s <- lapply(ioks2[chs2], function (xs) {
      if (length(xs) == 1L)
        sta2 + xs - 1L
      else {
        ps <- sta2 + xs - 1L
        i <- which.max(msxints[ps])
        sta2 + xs[i] - 1L
      }
    })
    p2s <- .Internal(unlist(p2s, recursive = FALSE, use.names = FALSE))
    yp1 <- max(msxints[p1])
    yp2 <- max(msxints[p2s])
    
    if (FALSE & ms_lev == 1L) {
      np1 <- max(n_ms1s[p1])
      np2 <- max(ns2 <- n_ms1s[p2s])
      
      if (np1 == np2) {
        if (yp1 >= yp2) {
          return(chs1 + 1L)
        }
        else {
          chx2 <- chs2 + 1L
          unqs <- !as.integer(chx2 * 2L) %in% chx2
          chs2 <- chs2[unqs]
          p2s <- p2s[unqs]
          chs2 <- chs2[which.max(msxints[p2s])]
          return(chs2 + 1L)
        }
      }
      
      if (np1 > np2) {
        return(chs1 + 1L)
      }
      else {
        chs2 <- chs2[which.max(ns2)]
        return(chs2 + 1L)
      }
    }
    else {
      if (yp1 >= yp2) {
        return(chs1 + 1L)
      }
      else {
        if (ms_lev == 1L) {
          ansduo <- check_chduo(duo = c(1L, 3L), chs2, p2s, msxints)
          chs2 <- ansduo$chs
          p2s <- ansduo$ps
          
          if (max_charge >= 6L) {
            ansduo <- check_chduo(duo = c(2L, 5L), chs2, p2s, msxints)
            chs2 <- ansduo$chs
            p2s <- ansduo$ps
          }
        }
        
        if (FALSE && ms_lev == 1L) {
          if (length(chs2) > 1L) {
            chx2 <- chs2 + 1L
            unqs <- !as.integer(chx2 * 2L) %in% chx2
            chs2 <- chs2[unqs]
            p2s <- p2s[unqs]
          }
          
          nsubs <- n_ms2s[p2s]
          
          if (max(nsubs) == min(nsubs)) {
            chs2 <- chs2[which.max(msxints[p2s])]
            return(chs2 + 1L)
          }
          
          chs2 <- chs2[which.max(n_ms1s[p2s])]
          return(chs2 + 1L)
        }
        else {
          chs2 <- chs2[which.max(msxints[p2s])]
          return(chs2 + 1L)
        }
      }
    }
  }
  
  # z = 1
  # may only needed when min_change >= 2L 
  if (!(l1 || l2)) {
    gap <- 1.003355
    mx1 <- mass + gap
    mx2 <- mass - gap
    iok1 <- .Internal(which(abs(xsub1/mx1 - 1L) <= tol))
    iok2 <- .Internal(which(abs(xsub1/mx2 - 1L) <= tol))
    len1 <- length(iok1)
    len2 <- length(iok2)
    
    if (len1 || len2) return(1L) else return(NA_integer_)
  }
  
  if (l2 == 1L) { # && l1 == 0
    return(chs2 + 1L)
  }
  
  # l2 > 1 && l1 == 0
  p2s <- lapply(ioks2[chs2], function (xs) {
    if (length(xs) == 1L)
      sta2 + xs - 1L
    else {
      ps <- sta2 + xs - 1L
      i <- which.max(msxints[ps])
      sta2 + xs[i] - 1L
    }
  })
  p2s <- .Internal(unlist(p2s, recursive = FALSE, use.names = FALSE))
  
  if (ms_lev == 1L) {
    ansduo <- check_chduo(duo = c(1L, 3L), chs2, p2s, msxints)
    chs2 <- ansduo$chs
    p2s <- ansduo$ps
    
    if (max_charge >= 6L) {
      ansduo <- check_chduo(duo = c(2L, 5L), chs2, p2s, msxints)
      chs2 <- ansduo$chs
      p2s <- ansduo$ps
    }
  }
  
  if (FALSE && ms_lev == 1L) {
    chx2 <- chs2 + 1L
    unqs <- !as.integer(chx2 * 2L) %in% chx2
    chs2 <- chs2[unqs]
    p2s <- p2s[unqs]
    
    nsubs <- n_ms1s[p2s]
    
    if (max(nsubs) == min(nsubs)) {
      chs2 <- chs2[which.max(msxints[p2s])]
      return(chs2 + 1L)
    }
    
    chs2 <- chs2[which.max(n_ms1s[p2s])]
    return(chs2 + 1L)
  }
  else {
    chs2 <- chs2[which.max(msxints[p2s])]
    return(chs2 + 1L)
  }
}


#' Checks the duplicacy between charge states 2 and 4.
#'
#' @param duo The potentially duplicated charge states. The values are off by
#'   one less in that the input charge states are derived from searches
#'   beginning charge states 2 instead of 1.
#' @param chs The charge states for checking duplicacy.
#' @param ps The corresponding positions of the observed \code{chs} in a
#'   m-over-z sequence.
#' @param msxints A sequence of intensity values.
#' @param cutoff A discriminating cut-off threshold.
check_chduo <- function (duo = c(1L, 3L), chs, ps, msxints, cutoff = 1.5)
{
  if (length(chs) < 2L)
    return(list(chs = chs, ps = ps))
  
  idxs <- match(duo, chs)
  i1 <- idxs[[1]]
  i2 <- idxs[[2]]
  
  if (!(is.na(i1) || is.na(i2))) {
    p1 <- ps[[i1]]
    p2 <- ps[[i2]]
    
    if (msxints[p1]/msxints[p2] < cutoff) {
      chs <- chs[-i1]
      ps <- ps[-i1]
    }
    else {
      chs <- chs[-i2]
      ps <- ps[-i2]
    }
  }
  
  list(chs = chs, ps = ps)
}


#' Is logical one.
#' 
#' @param x A logical vector.
is_true <- function(x) `==`(x, 1L)


#' Searches for a possible doubled charge state.
#'
#' Not used.
#' 
#' @param ch The initial charge state.
#' @param p The position in \code{moverzs}.
#' @param sta The position of start.
#' @param mass The current m-over-z.
#' @param moverzs A vector of m-over-z values.
#' @param max_charge The maximum charge state for consideration.
#' @param ppm Mass error tolerance.
#' @param f A function of \code{+} or \code{-} for forward or backward
#'   searching.
find_dbl_z <- function(ch = 2L, p = 2L, sta, mass, moverzs, max_charge = 4L, 
                       ppm = 5L, f = `+`) 
{
  # stopifnot(p >= 1L) # no other peaks in between
  
  # no other peaks in between
  if (p <= 1L) 
    return(ch)

  ans <- ch

  while((ch <- ch * 2L) <= max_charge) {
    mx  <- f(mass, 1.003355/ch)
    sta <- f(sta, 1L)
    end <- f(sta, p - 1L)
    oks <- abs((mx - moverzs[sta:end])/mx) * 1E6 <= ppm
    ps  <- Position(is_true, oks)
    
    if (is.na(ps)) 
      break 
    else
      ans <- ch
  }
  
  ans
}


#' Detects chromatographic peaks.
#' 
#' May be over-killing.
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
  signals <- rep.int(0L, len)
  va_y <- su_y <- ss_y <- stdFilter <- avgFilter <- rep_len(NA_real_, len)
  
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


#' Chunks of precursor lists by continuous retention-time bins.
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
#' Not used.
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


