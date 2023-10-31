#' De-isotopes precursor masses.
#'
#' @param moverzs Mass-to-charge ratios.
#' @param msxints MS1 or MS2 peak intensities.
#' @param n_ms1s The underlying counts that have contributed to \code{moverzs}
#'   and \code{maxints}.
#' @param center The mass center of an isolation window (only for MS1).
#' @param ppm Allowance in mass error when deisotoping.
#' @param offset_upr A cardinal number of upper mass off-sets.
#' @param offset_lwr A cardinal number of lower mass off-sets.
#' @param bound Not yet used. Logical; if TRUE, removes precursors outside of
#'   the boundary of isolation window.
#' @param ms_lev MS level.
#' @param maxn_feats The maximum number of MS features.
#' @param max_charge The maximum charge state.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries. The default is
#'   20 for MS1 and 10 for MS2.
#' @param step Step size for mass binning.
#' @param order_mz Logical; if TRUE, orders peaks from low to high m-over-z's.
#' @param backward_mass_co A mass cut-off to initiate backward looking of an
#'   isotope envelop.
#' @param fct_iso2 The multiplication factor for the second isotopic peak.
#' @param fct_bg A ratio factor to discriminate background noises from isotopic
#'   features.
#' @param fct_bg2 A second ratio factor to discriminate background noises from
#'   isotopic features.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' library(mzion)
#' moverzs <- c(881 + 1:10*.1, 882.0674, 882.0981, 882.4034, 882.60, 882.7372)
#' msxints <- c(1000 * 1:10, 1652869, 788368, 2043015, 2314111, 4314111)
#' n_ms1s <- rep_len(1L, length(msxints))
#'
#' # out <- mzion:::find_ms1stat(moverzs, msxints, n_ms1s, ppm = 5L, ms_lev = 1L,
#' #                             maxn_feats = 5L, max_charge = 4L, offset_upr = 8L,
#' #                             offset_lwr = 8L, order_mz = TRUE, bound = FALSE)
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
#' # out <- find_ms1stat(moverzs, msxints, n_ms1s, center = 909.788696,
#' #                     exclude_reporter_region = FALSE,
#' #                     ppm = 5L, ms_lev = 1L, maxn_feats = 1L, max_charge = 4L,
#' #                     order_mz = TRUE, step = 5/1e6)
#' # out <- find_ms1stat(moverzs, msxints, n_ms1s, # center = 909.788696,
#' #                     exclude_reporter_region = FALSE,
#' #                     ppm = 5L, ms_lev = 1L, maxn_feats = 1L, max_charge = 4L,
#' #                     order_mz = TRUE, step = 5/1e6)
#' }
find_ms1stat <- function (moverzs, msxints, n_ms1s = 1L, center = 0, 
                          exclude_reporter_region = FALSE, 
                          tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                          ppm = 8L, ms_lev = 1L, maxn_feats = 300L, 
                          max_charge = 4L, n_fwd = 20L, offset_upr = 30L, 
                          offset_lwr = 30L, order_mz = TRUE, step = ppm/1e6, 
                          backward_mass_co = 800/ms_lev, grad_isotope = 2.5, 
                          fct_iso2 = 3.0, fct_bg = 15, fct_bg2 = 5, 
                          use_defpeaks = FALSE, bound = FALSE)
{
  ###
  # if to apply intensity cut-offs, should note the difference intensity 
  # scales between Thermo's and Bruker's data
  ###
  
  null_out <- list(masses = NULL, charges = NULL, intensities = NULL)
  
  if (!(len_ms <- length(moverzs)))
    return(null_out)
  
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
    rptr_moverzs <- moverzs[ok_rptrs]
    rptr_ints <- msxints[ok_rptrs]
    
    no_rptrs <- !ok_rptrs
    moverzs <- moverzs[no_rptrs]
    msxints <- msxints[no_rptrs]
    
    len_rptrs <- length(rptr_moverzs)
    len_ms <- len_ms - len_rptrs
    
    if (!len_ms)
      return(null_out)
  }
  
  from <- moverzs[[1]] - .001
  ims  <- index_mz(moverzs, from, step)
  
  lenp <- min(len_ms, maxn_feats)
  peaks_fuz <- intens_fuz <- peaks <- intens <- rep_len(NA_real_, lenp)
  css_fuz <- css <- rep.int(NA_integer_, lenp)
  p_fuz <- p <- 1L
  mass_fuz <- NA_real_
  
  while(len_ms & p <= maxn_feats) {
    if (len_ms == 1L) {
      intens[[p]] <- msxints
      peaks[[p]] <- moverzs
      css[[p]] <- 0L
      break
    }
    
    if (center > 0 && p == 1L && ms_lev == 1L) {
      imax <- .Internal(which.min(abs(moverzs - center)))
    } else {
      imax <- .Internal(which.max(msxints))
    }
    
    mass <- moverzs[imax]
    mint <- msxints[imax]

    ch <- find_charge_state(mass = mass, imax = imax, mint = mint, 
                            moverzs = moverzs, msxints = msxints, n_ms1s = n_ms1s, 
                            lenm = len_ms, max_charge = max_charge, n_fwd = n_fwd, 
                            ms_lev = ms_lev, ppm = ppm, fct_bg = fct_bg)

    if (ch == 0L) {
      if (FALSE) {
        # not ideal: e.g. at max_charge = 4, 
        #  a z = 5 peak will not be identified but assigned z = 2 here
        if (ms_lev == 1L && !use_defpeaks) {
          peaks[[p]] <- center
          css[[p]] <- 2L # arbitrary, better learn from MS2
        }
        else {
          peaks[[p]] <- mass
        }
      }

      peaks[[p]] <- mass
      intens[[p]] <- mint
      len_ms <- len_ms - 1L
      moverzs <- moverzs[-imax]
      msxints <- msxints[-imax]
      ims <- ims[-imax]
      p <- p + 1L
      next
    }
    
    gap <- 1.003355/ch
    err <- ppm/1e6
    sta1 <- min(len_ms, imax + 1L)
    end1 <- min(len_ms, imax + offset_upr)

    # find mono-isotopic m/z
    if (ch * mass > backward_mass_co) {
      iths <- index_mz(mass + gap * (-ch-1L):(1L+ch), from, step)
      if ((lwr <- imax - offset_lwr) < 1L) lwr <- 1L
      if ((upr <- imax + offset_upr) > len_ms) upr <- len_ms
      
      iexs <- ims[lwr:upr]
      oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
      hits <- .Internal(which(oks2)) + lwr - 1L
      lenh <- length(hits)
      idx <- .Internal(which(hits == imax))
      
      if (idx == 2L) {
        hi <- hits[[1]]
        
        if (abs((moverzs[[hi]] + gap)/mass - 1) <= err) {
          if ((ri <- mint/msxints[[hi]]) <= grad_isotope)
            mass <- moverzs[[hi]]
          else {
            if (ri <= fct_iso2 * grad_isotope) mass_fuz <- moverzs[[hi]]
          }
        }
      }
      else if (idx > 2L) {
        ix <- lenx <- idx - 1L
        hx <- NULL
        
        for (i in seq_len(lenx)) {
          mzsub <- moverzs[hits[1:ix]]
          ks <- .Internal(which(abs((mzsub + gap)/mass - 1) <= err))
          
          # all satellite to `mass`
          if (!length(ks)) {
            break
          }
          
          hsub <- hits[ks]
          isub <- which.max(msxints[hsub])
          hi <- hsub[isub]
          h <- ks[isub]
          gi <- grad_isotope * i
          ri <- mint/msxints[[hi]]
          
          if (ri > gi) {
            if (ri < fct_iso2 * gi) mass_fuz <- moverzs[[hi]]
            break
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
    
    intens[[p]] <- sum(msxints[hits])
    peaks[[p]] <- mass
    css[[p]] <- ch
    
    if (!is.na(mass_fuz)) {
      intens_fuz[[p_fuz]] <- intens[[p]]
      peaks_fuz[[p_fuz]] <- mass_fuz
      css_fuz[[p_fuz]] <- ch
      p_fuz <- p_fuz + 1L
      mass_fuz <- NA_real_
    }
    
    len_ms <- len_ms - lenh
    moverzs <- moverzs[-hits]
    msxints <- msxints[-hits]
    ims <- ims[-hits]
    p <- p + 1L
  }
  
  # outputs
  if (ms_lev == 1L)
    oks1 <- css > 1L & !is.na(css)
  else if (ms_lev == 2L)
    oks1 <- !is.na(peaks)
  else
    stop("Unhandled MS level = ", ms_lev)
  
  masses <- peaks[oks1]
  charges <- css[oks1]
  intensities <- intens[oks1]
  
  if (endf <- p_fuz - 1L) {
    masses <- c(masses, peaks_fuz[1:endf])
    charges <- c(charges, css_fuz[1:endf])
    intensities <- c(intensities, intens_fuz[1:endf])
  }
  
  if (excl_rptrs && len_rptrs) {
    masses <- c(rptr_moverzs, masses)
    charges <- c(rep.int(0L, len_rptrs), charges)
    intensities <- c(rptr_ints, intensities)
  }
  
  if (order_mz <- FALSE) {
    ord <- order(masses)
    masses <- masses[ord]
    charges <- charges[ord]
    intensities <- intensities[ord]
  }
  
  out <- list(masses = masses, charges = charges, intensities = intensities)
}


#' Finds charge state.
#' 
#' @param mass The current mass.
#' @param imax The index of the most intense peak.
#' @param mint The maximum intensity (corresponding to the peak at \code{imax}).
#' @param lenm The length of \code{moverzs}.
#' @inheritParams find_ms1stat 
find_charge_state <- function (mass, imax, mint, moverzs, msxints, n_ms1s, lenm, 
                               max_charge = 4L, n_fwd = 20L, 
                               ms_lev = 1L, ppm = 8, fct_bg = 10)
{
  # find results for all ch -> uses the best...
  # distinguish left, right or both left and right evidence
  # if left-only evidence -> monomass must < mass

  # stopifnot(max_charge >= 2L, ppm >= 1)
  
  gaps <- 1.003355/2:max_charge
  err <- ppm/1e6

  sta1 <- min(lenm, imax + 1L)
  sta2 <- max(1L, imax - n_fwd)
  end1 <- min(lenm, imax + n_fwd)
  end2 <- max(1L, imax - 1L)
  xsub1 <- moverzs[sta1:end1]
  xsub2 <- moverzs[sta2:end2]
  ysub1 <- msxints[sta1:end1]
  ysub2 <- msxints[sta2:end2]
  
  if (ms_lev == 1L) {
    nsub1 <- n_ms1s[sta1:end1]
    nsub2 <- n_ms1s[sta2:end2]
  }

  mxs1 <- mass + gaps
  mxs2 <- mass - gaps
  ioks1 <- lapply(mxs1, function (x) .Internal(which(abs(xsub1/x - 1L) <= err)))
  ioks2 <- lapply(mxs2, function (x) .Internal(which(abs(xsub2/x - 1L) <= err)))
  lens1 <- lengths(ioks1)
  lens2 <- lengths(ioks2)
  chs1 <- .Internal(which(lens1 > 0L))
  chs2 <- .Internal(which(lens2 > 0L))
  l1 <- length(chs1)
  l2 <- length(chs2)
  
  if (!(l1 || l2)) {
    gap <- 1.003355
    mx1 <- mass + gap
    mx2 <- mass - gap
    iok1 <- .Internal(which(abs(xsub1/mx1 - 1L) <= err))
    iok2 <- .Internal(which(abs(xsub1/mx2 - 1L) <= err))
    len1 <- length(iok1)
    len2 <- length(iok2)
    
    if (len1 || len2) return(1L) else return(0L)
  }
  
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
      
      if (ms_lev == 1L) {
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
  
  if (l1 > 1L) {
    p1s <- lapply(ioks1[chs1], function (xs) {
      if (length(xs) == 1L)
        sta1 + xs - 1L
      else {
        ps <- sta1 + xs - 1L
        i <- which.max(msxints[ps])
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


#' Checks the duplicacy detween charge states 2 and 4.
#'
#' @param duo The potentially duplicated charge states. The values are off by
#'   one less in that the input charge states are derived from searches
#'   beginning charge states 2 instead of 1.
#' @param chs The charge states for checking duplicacy.
#' @param ps The corresponding positions of the observed \codes{chs} in a
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


#' Finds charge state.
#' 
#' @inheritParams find_charge_state
find_charge_state_v1 <- function (mass, imax, mint, len_ms, moverzs, msxints, 
                                  n_ms1s, max_charge = 4L, n_fwd = 20L, 
                                  ms_lev = 1L, ppm = 8, fct_bg = 10)
{
  for (ch in max_charge:0) {
    if (ch == 0L) {
      break
    }
    
    gap <- 1.003355/ch
    mx1 <- mass + gap
    sta1 <- min(len_ms, imax + 1L)
    end1 <- min(len_ms, imax + n_fwd)
    xsub1 <- moverzs[sta1:end1]
    ioks1 <- .Internal(which(abs(1 - xsub1/mx1) * 1E6 <= ppm))
    
    mx2  <- mass - gap
    sta2 <- max(1L, imax - n_fwd)
    end2 <- max(1L, imax - 1L)
    xsub2 <- moverzs[sta2:end2]
    ioks2 <- .Internal(which(abs(1 - xsub2/mx2) * 1E6 <= ppm))
    
    if (!length(c(ioks2, ioks1))) {
      next
    }
    
    # noises
    isub1 <- sta1 + ioks1 - 1L
    isub2 <- sta2 + ioks2 - 1L
    ysub1 <- msxints[isub1]
    ysub2 <- msxints[isub2]
    
    if (all(c(mint/ysub1, mint/ysub2) > fct_bg)) {
      next
    }
    
    # one-hit wonders
    isub <- c(isub2, isub1)
    xsub1 <- moverzs[isub1]
    xsub2 <- moverzs[isub2]
    # xsub <- c(xsub2, xsub1)
    # ysub <- c(ysub2, ysub1)
    
    if (ms_lev == 1L && ch > 1L && all(n_ms1s[isub] == 1L)) {
      chx <- ch - 1L
      gapx <- 1.003355/chx
      mxx1 <- `+`(mass, gapx)
      mxx2 <- `-`(mass, gapx)
      
      ioksx1 <- .Internal(which(abs(1 - xsub1/mxx1) * 1E6 <= ppm))
      ioksx2 <- .Internal(which(abs(1 - xsub2/mxx2) * 1E6 <= ppm))
      isubx1 <- sta1 + ioksx1 - 1L
      isubx2 <- sta2 + ioksx2 - 1L
      isubx <- c(isubx2, isubx1)
      
      if (length(isubx) && any(n_ms1s[isubx] > 1L)) {
        ch <- chx
        gap <- gapx
        mx1 <- mxx1
        mx2 <- mxx2
        # ioks1 <- ioksx1
        # ioks2 <- ioksx2
        isub1 <- isubx1
        isub2 <- isubx2
        # isub <- isubx
        ysub1 <- msxints[isub1]
        ysub2 <- msxints[isub2]
        # ysub <- c(ysub2, ysub1)
        xsub1 <- moverzs[isub1]
        xsub2 <- moverzs[isub2]
        # xsub <- c(xsub2, xsub1)
      }
      # rm(list = c("mxx1", "mxx2", "oksx1", "oksx2", "ioksx1", "ioksx2"))
    }
    
    # check charge halving (for simplicity, no recursive halving)
    if (ch >= 4L && ch %% 2L == 0L) {
      ch1 <- check_half_cs(ch, mass, xsub1, ysub1, sta1, msxints, ppm, `+`)
      ch2 <- check_half_cs(ch, mass, xsub2, ysub2, sta2, msxints, ppm, `-`)
      chx <- if (ch1 <= ch2) ch1 else ch2
      if (chx != ch) ch <- chx
    }
    
    break
  }
  
  list(ch = ch, gap = gap, sta1 = sta1, end1 = end1)
}


#' De-isotopes precursor masses.
#' 
#' The right half first.
#'
#' @param moverzs Mass-to-charge ratios.
#' @param msxints MS1 or MS2 peak intensities.
#' @param n_ms1s The underlying counts that have contributed to \code{moverzs}
#'   and \code{maxints}.
#' @param center The mass center of an isolation window (only for MS1).
#' @param ppm Allowance in mass error when deisotoping.
#' @param offset_upr A cardinal number of upper mass off-sets.
#' @param offset_lwr A cardinal number of lower mass off-sets.
#' @param bound Not yet used. Logical; if TRUE, removes precursors outside of
#'   the boundary of isolation window.
#' @param ms_lev MS level.
#' @param maxn_feats The maximum number of MS features.
#' @param max_charge The maximum charge state.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries. The default is
#'   20 for MS1 and 10 for MS2.
#' @param step Step size for mass binning.
#' @param order_mz Logical; if TRUE, orders peaks from low to high m-over-z's.
#' @param backward_mass_co A mass cut-off to initiate backward looking of an
#'   isotope envelop.
#' @param fct_iso2 The multiplication factor for the second isotopic peak.
#' @param fct_bg A ratio factor to discriminate background noises from isotopic
#'   features.
#' @param fct_bg2 A second ratio factor to discriminate background noises from
#'   isotopic features.
#' @inheritParams matchMS
find_ms1stat_rf <- function (moverzs, msxints, n_ms1s = 1L, center = 0, 
                          exclude_reporter_region = FALSE, 
                          tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                          ppm = 5L, ms_lev = 1L, maxn_feats = 300L, 
                          max_charge = 4L, n_fwd = 20L, offset_upr = 30L, 
                          offset_lwr = 30L, order_mz = TRUE, step = ppm/1e6, 
                          backward_mass_co = 800/ms_lev, grad_isotope = 2.5, 
                          fct_iso2 = 3.0, fct_bg = 15, fct_bg2 = 5, bound = FALSE)
{
  ###
  # if to apply intensity cut-offs, should note the difference intensity 
  # scales between Thermo's and Bruker's data
  ###
  
  null_out <- list(masses = NULL, charges = NULL, intensities = NULL)
  
  if (!(len_ms <- length(moverzs)))
    return(null_out)
  
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
    rptr_moverzs <- moverzs[ok_rptrs]
    rptr_ints <- msxints[ok_rptrs]
    
    no_rptrs <- !ok_rptrs
    moverzs <- moverzs[no_rptrs]
    msxints <- msxints[no_rptrs]
    
    len_rptrs <- length(rptr_moverzs)
    len_ms <- len_ms - len_rptrs
    
    if (!len_ms)
      return(null_out)
  }
  
  from <- moverzs[[1]] - .001
  ims  <- index_mz(moverzs, from, step)
  
  lenp <- min(len_ms, maxn_feats)
  peaks_fuz <- intens_fuz <- peaks <- intens <- rep_len(NA_real_, lenp)
  css_fuz <- css <- rep.int(NA_integer_, lenp)
  p_fuz <- p <- 1L
  mass_fuz <- NA_real_

  while(len_ms & p <= maxn_feats) {
    if (len_ms == 1L) {
      intens[[p]] <- msxints
      peaks[[p]] <- moverzs
      css[[p]] <- 0L
      break
    }
    
    if (center > 0 && p == 1L && ms_lev == 1L) {
      imax <- .Internal(which.min(abs(moverzs - center)))
    } else {
      imax <- .Internal(which.max(msxints))
    }

    mass <- moverzs[imax]
    mint <- msxints[imax]
    
    for (ch in max_charge:0) {
      if (ch == 0L) {
        break
      }

      gap <- 1.003355/ch
      mx <- mass + gap
      sta <- min(len_ms, imax + 1L)
      end <- min(len_ms, imax + n_fwd)
      xsub <- moverzs[sta:end]
      oks <- abs((mx - xsub)/mx) * 1E6 <= ppm
      ioks <- .Internal(which(oks))
      
      if (!length(ioks)) {
        next
      }

      # noises
      rgsub <- sta + ioks - 1L
      ysub <- msxints[rgsub]
      rysub <- mint/ysub
      
      if (all(rysub > fct_bg)) {
        next
      }
      
      # one-hit wonders
      if (ms_lev == 1L && ch > 1L) {
        if (all(n_ms1s[rgsub] == 1L)) {
          ansch <- check_cs_next(ch, mass, xsub, ysub, sta, msxints, n_ms1s, 
                                 ppm, f = `+`)
          ch <- ansch$ch
          gap <- ansch$gap
          ysub <- ansch$ysub
        }
      }

      # check charge halving (for simplicity, no recursive halving)
      if (ch >= 4L && ch %% 2L == 0L) {
        ch <- check_half_cs(ch, mass, xsub, ysub, sta, msxints, ppm, `+`)
        gap <- 1.003355/ch
      }
      
      # find mono-isotopic m/z
      if (ch * mass > backward_mass_co) {
        iths <- index_mz(mass + gap * (1L - ch):(1L + ch), from, step)
        if ((lwr <- imax - offset_lwr) < 1L) lwr <- 1L
        if ((upr <- imax + offset_upr) > len_ms) upr <- len_ms
        
        iexs <- ims[lwr:upr]
        oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
        hits <- .Internal(which(oks2)) + lwr - 1L
        lenh <- length(hits)
        idx <- .Internal(which(hits == imax))
        
        if (idx == 2L) {
          hi <- hits[[1]]
          
          if (abs((moverzs[[hi]] + gap)/mass - 1) * 1E6 <= ppm) {
            if ((ri <- mint/msxints[[hi]]) <= grad_isotope)
              mass <- moverzs[[hi]]
            else {
              if (ri <= fct_iso2 * grad_isotope) mass_fuz <- moverzs[[hi]]
            }
          }
        }
        else if (idx > 2L) {
          ix <- lenx <- idx - 1L
          hx <- NULL
          
          for (i in seq_len(lenx)) {
            mzsub <- moverzs[hits[1:ix]]
            ks <- .Internal(which(abs((mzsub + gap)/mass -1) * 1E6 <= ppm))
            
            # all satellite to `mass`
            if (!length(ks))
              break
            
            hsub <- hits[ks]
            isub <- which.max(msxints[hsub])
            hi <- hsub[isub]
            h <- ks[isub]
            gi <- grad_isotope * i
            ri <- mint/msxints[[hi]]
            
            if (ri > gi) {
              if (ri < fct_iso2 * gi) mass_fuz <- moverzs[[hi]]
              break
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
        iexs <- ims[sta:min(imax + offset_upr, len_ms)]
        oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
        hits <- c(imax, .Internal(which(oks2)) + imax) 
        lenh <- length(hits)
      }
      
      intens[[p]] <- sum(msxints[hits])
      peaks[[p]] <- mass
      css[[p]] <- ch
      
      if (!is.na(mass_fuz)) {
        intens_fuz[[p_fuz]] <- intens[[p]]
        peaks_fuz[[p_fuz]] <- mass_fuz
        css_fuz[[p_fuz]] <- ch
        p_fuz <- p_fuz + 1L
        mass_fuz <- NA_real_
      }

      len_ms <- len_ms - lenh
      moverzs <- moverzs[-hits]
      msxints <- msxints[-hits]
      ims <- ims[-hits]

      break
    }
    
    if (FALSE) {
      # second chance for ch == 1L
      if (ch == 1L) {
        gap <- 1.003355/2
        mx  <- mass - gap
        end <- max(1L, imax - 1L)
        sta <- max(1L, imax - n_fwd)
        xsub <- moverzs[sta:end]
        ioks <- .Internal(which(abs(1 - xsub/mx) * 1E6 <= ppm))
        
        if (length(ioks)) {
          ysub <- msxints[sta + ioks - 1L]
          
          if (any((mint/ysub) < fct_bg)) {
            # -(2:1) since mass itself already used
            iths <- index_mz(mass + gap * -(2:1), from, step)
            lwr <- max(imax - offset_lwr, 1L)
            upr <- imax
            iexs <- ims[lwr:upr]
            oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
            hits <- .Internal(which(oks2)) + lwr - 1L
            lenh <- length(hits)
            
            intens[[p]] <- sum(msxints[hits])
            peaks[[p]] <- mass <- moverzs[hits[[1]]]
            css[[p]] <- ch <- 2L
            
            len_ms <- len_ms - lenh
            moverzs <- moverzs[-hits]
            msxints <- msxints[-hits]
            ims <- ims[-hits]
          }
        }
      }
    }
    
    # backward search
    if (ch > 0L)
      p <- p + 1L
    
    if (ch == 0L) {
      for (ch in max_charge:0) {
        if (ch == 0L) {
          break
        }

        gap <- 1.003355/ch
        mx  <- mass - gap
        end <- max(1L, imax - 1L)
        sta <- max(1L, imax - n_fwd)
        xsub <- moverzs[sta:end]
        oks <- abs(1 - xsub/mx) * 1E6 <= ppm
        ioks <- .Internal(which(oks))
        
        if (!length(ioks))
          next
        
        ysub <- msxints[sta + ioks - 1L]
        rysub <- mint/ysub
        
        if (all(rysub > fct_bg))
          next

        if (ch >= 4L && ch %% 2L == 0L) {
          gap2 <- gap * 2L
          mx2 <- mass - gap2
          oks2 <- abs((mx2 - xsub)/mx2) * 1E6 <= ppm
          ioks2 <- .Internal(which(oks2))
          
          if (length(ioks2) && (max(msxints[sta + ioks2 - 1L])/max(ysub) > 1.5)) {
            ch <- ch/2L
            gap <- gap2
          }
        }
        
        iths <- index_mz(mass + gap * -ch:0L, from, step)
        lwr <- max(imax - offset_lwr, 1L)
        upr <- imax
        iexs <- ims[lwr:upr]
        oks2 <- iexs %fin% iths | (iexs - 1L) %fin% iths | (iexs + 1L) %fin% iths
        hits <- .Internal(which(oks2)) + lwr - 1L
        lenh <- length(hits)
        
        intens[[p]] <- sum(msxints[hits])
        peaks[[p]] <- mass <- moverzs[hits[[1]]]
        css[[p]] <- ch

        len_ms <- len_ms - lenh
        moverzs <- moverzs[-hits]
        msxints <- msxints[-hits]
        ims <- ims[-hits]
        # p <- p + 1L
        
        break
      }
      
      if (ch > 0)
        p <- p + 1L
    }
    
    if (ch == 0L) {
      peaks[[p]] <- mass
      intens[[p]] <- mint
      moverzs <- moverzs[-imax]
      msxints <- msxints[-imax]
      ims <- ims[-imax]
      len_ms <- len_ms - 1L
      p <- p + 1L
    }
  }
  
  # outputs
  if (ms_lev == 1L)
    oks1 <- css > 1L & !is.na(css)
  else if (ms_lev == 2L)
    oks1 <- !is.na(peaks)
  else
    stop("Unhandled MS level = ", ms_lev)
  
  masses <- peaks[oks1]
  charges <- css[oks1]
  intensities <- intens[oks1]
  
  if (end2 <- p_fuz - 1L) {
    masses <- c(masses, peaks_fuz[1:end2])
    charges <- c(charges, css_fuz[1:end2])
    intensities <- c(intensities, intens_fuz[1:end2])
  }
  
  if (excl_rptrs && len_rptrs) {
    masses <- c(rptr_moverzs, masses)
    charges <- c(rep.int(0L, len_rptrs), charges)
    intensities <- c(rptr_ints, intensities)
  }
  
  if (order_mz <- FALSE) {
    ord <- order(masses)
    masses <- masses[ord]
    charges <- charges[ord]
    intensities <- intensities[ord]
  }
  
  out <- list(masses = masses, charges = charges, intensities = intensities)
}


#' Checks the goodness of a next charge state
#' 
#' Next charge state: ch - 1
#' 
#' @param ch The current charge state.
#' @param mass The current mass.
#' @param xsub A subset of moverzs.
#' @param ysub A subset of msxints.
#' @param msxints MS2 intensities.
#' @param sta The start position.
#' @param f A function of plus or minus.
#' @inheritParams find_ms1stat
check_cs_next <- function (ch, mass, xsub, ysub, sta, msxints, n_ms1s, 
                           ppm = 8, f = `+`)
{
  gap <- 1.003355/ch
  chx <- ch - 1L
  gapx <- 1.003355/chx
  mxx <- f(mass, gapx)
  oksx <- abs(1 - xsub/mxx) * 1E6 <= ppm
  ioksx <- .Internal(which(oksx))
  
  if (!length(ioksx))
    return(list(ch = ch, gap = gap, ysub = ysub))
  
  rgsubx <- sta + ioksx - 1L
  
  if (!any(n_ms1s[rgsubx] > 1L))
    return(list(ch = ch, gap = gap, ysub = ysub))
  
  list(ch = chx, gap = gapx, ysub = msxints[rgsubx])
}


#' Check the goodness of a half charge state.
#' 
#' @inheritParams check_cs_next
check_half_cs <- function (ch, mass, xsub, ysub, sta, msxints, ppm = 8, f = `+`)
{
  gap <- 1.003355/ch
  gapx <- gap * 2L
  mxx <- f(mass, gapx)
  oksx2 <- abs(1 - xsub/mxx) * 1E6 <= ppm
  ioksx2 <- .Internal(which(oksx2))
  
  if (length(ioksx2) && (max(msxints[sta + ioksx2 - 1L])/max(ysub) > 1.5)) {
    ch <- ch/2L
  }
  
  ch
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


#' Finds the charge state of a mass.
#' 
#' Not used.
#' 
#' @inheritParams find_dbl_z
#' @examples
#' \donttest{
#' library(mzion)
#' moverzs <- c(881 + 1:10*.1, 882.0674, 882.0981, 882.4034, 882.60, 882.7372)
#' mzion:::find_charage_state(1, moverzs[[1]], moverzs)
#' mzion:::find_charage_state(5, moverzs[[5]], moverzs)
#' mzion:::find_charage_state(7, moverzs[[7]], moverzs[-10])
#' }
find_charage_state <- function (imax, mass, moverzs, max_charge = 4L, ppm = 5L) 
{
  len <- length(moverzs)
  ans <- 0L
  
  # p <- function(f, u) function(x) f(x, u)
  # g <- function(x, u) Position(p(`>=`, u), x)
  # end <- g(moverzs, mx * (1 + ppm/1e6))
  # p <- function(u) function(x) `>=`(x, u)
  # g <- function(x, u) Position(p(u), x)
  
  sta <- min(len, `+`(imax, 1L))

  for (ch in 1:max_charge) {
    gap <- 1.003355/ch
    mx  <- `+`(mass, gap)
    rs  <- .Internal(which(moverzs[sta:len] >= `+`(mx, mx * ppm/1e6)))
    end <- if (length(rs)) rs[[1]] + sta - 1L else sta
    oks <- abs((mx - moverzs[sta:end])/mx) * 1E6 <= ppm
    ps1 <- Position(is_true, oks)
    
    if (!is.na(ps1)) {
      if (ch > 1L) {
        ans <- find_dbl_z(ch, ps1, sta, mass, moverzs, max_charge, ppm, `+`)
        break
      }
      else {
        ans <- ch
      }
    }
  }

  if (ans > 0L)
    return(ans)
  
  sta <- max(1L, `-`(imax, 1L))
  
  for (ch in 1:max_charge) {
    gap <- 1.003355/ch
    mx  <- `-`(mass, gap)
    rs  <- .Internal(which(moverzs[1:sta] >= `-`(mx, mx * ppm/1e6)))
    end <- if (length(rs)) rs[[1]] else sta
    oks <- abs((mx - moverzs[sta:end])/mx) * 1E6 <= ppm
    ps1 <- Position(is_true, oks)
    
    if (!is.na(ps1)) {
      if (ch > 1L) {
        ans <- find_dbl_z(ch, ps1, sta, mass, moverzs, max_charge, ppm, `-`)
        break
      }
      else {
        ans <- ch
      }
    }
  }
  
  ans
}


#' Detects chromatographic peaks.
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










