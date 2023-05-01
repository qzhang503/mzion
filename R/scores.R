#' Adds secondary ions of b0, y0 etc.
#' 
#' @param ms2s A vector of theoretical MS2 m-over-z values.
#' @inheritParams matchMS
add_seions <- function (ms2s, type_ms2ions = "by", digits = 4L) 
{
  len <- length(ms2s)
  
  if (type_ms2ions == "by") {
    proton <- 1.00727647
    h2o <- 18.010565
    nh3 <- 17.026549
    
    bs <- ms2s[1:(len/2)]
    ys <- ms2s[(len/2+1):len]
    
    b2s <- (bs + proton)/2
    bstars <- bs - nh3
    bstar2s <- (bstars + proton)/2
    b0s <- bs - h2o
    b02s <- (b0s + proton)/2
    
    y2s <- (ys + proton)/2
    ystars <- ys - nh3
    ystar2s <- (ystars + proton)/2
    y0s <- ys - h2o
    y02s <- (y0s + proton)/2
    
    round(c(b2s, bstars, bstar2s, b0s, b02s, y2s, ystars, ystar2s, y0s, y02s), 
          digits = digits)
  } 
  else if (type_ms2ions == "ax") {
    proton <- 1.00727647
    h2o <- 18.010565
    nh3 <- 17.026549
    
    as <- ms2s[1:(len/2)]
    xs <- ms2s[(len/2+1):len]
    
    a2s <- (as + proton)/2
    astars <- as - nh3
    astar2s <- (astars + proton)/2
    a0s <- as - h2o
    a02s <- (a0s + proton)/2
    
    x2s <- (xs + proton)/2
    
    round(c(a2s, astars, astar2s, a0s, a02s, x2s), digits = digits)
  } 
  else if (type_ms2ions == "cz") {
    proton <- 1.00727647
    
    cs <- ms2s[1:(len/2)]
    zs <- ms2s[(len/2+1):len]
    
    c2s <- (cs + proton)/2
    z2s <- (zs + proton)/2
    
    round(c(c2s, z2s), digits = digits)
  }
}


#' Matches two lists without making a data frame
#' 
#' Not currently used. 
#' 
#' @param a The left vector.
#' @param b The right vector.
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' a <- c(3, 4, 1, 2, 5)
#' b <- 2
#' 
#' mzion:::list_leftmatch(a, b)
#' }
list_leftmatch <- function (a, b) 
{
  ord <- order(a, decreasing = TRUE)
  a <- a[ord]
  
  oks <- a %in% b
  
  b2 <- rep(NA_real_, length(a))
  b2[oks] <- b
  
  b2
}


#' Helper for score calculations
#'
#' By the positions of variable modifications.
#' 
#' @section Model: 
#' N - the total number of features (white and black balls) \cr
#' k - the number of sampled features \cr
#' m - the numbers of theoretical features (white balls) \cr
#' n - the number of noise (black balls) \cr
#' 
#' * Subtracts \code{m} and the counts of secondary b0, y0 matches etc. from noise
#' ((N < m) -> (n < 0L); OK if n < 0L) \cr
#' * \code{ith} and \code{ith2} in ascending order, \code{iex} and \code{iex2} in adaptive order \cr
#' * One-to-one correspondence: \code{ith} <-> \code{iex}; \code{ith2} <-> \code{iex2} \cr
#' 
#' @section Check matches:
#' \code{abs(df2$theo[df2$ith] - expt_moverzs[df2$iex]) <= ppm_ms2} \cr
#' \code{identical(df2$int[df2$ith], expt_ints[df2$iex])}
#' 
#' @section y:
#' y$expt - experimental m/z's \cr
#' y$int  - experimental intensities (primary) \cr
#' y$int2 - experimental intensiteis (secondary) \cr
#' y$theo - matched theoretical m/z's \cr
#' y$idx  - values: the indexes in theoretical sequence (df$ith)
#'          positions: the position in expt_moverzs
#' 
#' @param df Two lists of \code{theo} and matched \code{expt} m-over-z.
#' @param nms The names (character strings indicating the names and position of
#'   variable modifications).
#' @param burn_ins The range of burn-ins where inputs will be excluded from
#'   probablity assessments.
#' @inheritParams calc_probi
#' @import dplyr
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' ##
#' pep <- "YGPQYGHPPPPPPPPDYGPHADSPVLMVYGLDQSK"
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' theos <- c(393.2335,450.2550,547.3078,675.3663,838.4297,
#'            895.4511,1032.5100,1129.5628,1226.6156,1323.6683,
#'            1420.7211,1517.7739,1614.8266,1711.8794,1808.9322,
#'            1923.9591,2087.0224,2144.0439,2241.0967,2378.1556,
#'            2449.1927,2564.2196,2651.2517,2748.3044,2847.3728,
#'            2960.4569,3091.4974,3190.5658,3353.6291,3410.6506,
#'            3523.7347,3638.7616,3766.8202,3853.8522,4211.1101,
#'            376.2757,463.3078,591.3663,706.3933,819.4773,
#'            876.4988,1039.5621,1138.6306,1269.6710,1382.7551,
#'            1481.8235,1578.8763,1665.9083,1780.9353,1851.9724,
#'            1989.0313,2086.0840,2143.1055,2306.1688,2421.1958,
#'            2518.2485,2615.3013,2712.3541,2809.4068,2906.4596,
#'            3003.5124,3100.5651,3197.6179,3334.6768,3391.6983,
#'            3554.7616,3682.8202,3779.8729,3836.8944,3999.9577)
#' 
#' expts <- c(NA,NA,NA,675.36646,838.42981,895.45056,1032.51025,
#'            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
#'            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,376.27603,463.30823,
#'            591.36676,706.39380,819.47552,876.50018,1039.56287,
#'            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
#'            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
#' 
#' names(expts) <- names(theos) <- c(nms, rev(nms))
#' 
#' ith <- c(4,5,6,7,36,37,38,39,40,41,42)
#' iex <- c(42,57,72,93,20,26,36,44,56,65,96)
#' ith <- as.integer(ith)
#' iex <- as.integer(iex)
#' m <- length(ith)
#' 
#' df <- list(theo = theos, expt = expts, ith = ith, iex = iex, m = m)
#' 
#' expt_moverzs <- c(126.05530,126.12798,127.12505,127.13139,128.12843,
#'                   128.13475,129.13177,129.13806,130.13516,130.14140,
#'                   131.13841,136.07588,143.11812,195.11308,221.09239,
#'                   230.17036,235.14416,244.09303,271.12891,376.27603,
#'                   377.28003,384.21320,400.54648,441.89435,442.22861,
#'                   463.30823,464.31158,502.76166,503.26221,516.75928,
#'                   517.26062,531.28296,565.28522,565.78705,574.33997,
#'                   591.36676,592.36957,613.81238,614.31451,648.34070,
#'                   669.80933,675.36646,696.86902,706.39380,707.39600,
#'                   713.32672,713.82715,714.33093,761.85156,762.33771,
#'                   762.85907,810.37610,810.88019,811.38434,811.88660,
#'                   819.47552,838.42981,846.49664,846.99695,855.74426,
#'                   859.91162,867.92706,868.43024,868.93225,876.50018,
#'                   877.50354,878.75262,884.42273,884.75665,885.09052,
#'                   885.42700,895.45056,896.45435,908.43860,908.94092,
#'                   909.44440,916.95612,949.79688,950.13049,950.46545,
#'                   950.79791,962.48212,964.98084,965.48273,987.49255,
#'                   987.82599,988.16028,988.50989,1004.51562,1005.52985,
#'                   1013.56659,1031.50378,1032.51025,1033.51355,1034.51636,
#'                   1039.56287,1040.56641,1062.48096,1063.48511,1161.55078)
#' 
#' expt_ints <- c(16921,29468,36904,28121,37829,23537,39307,36194,25192,33532,
#'                26551,91477,15182,24720,17471,50430,14282,14084,39681,99581,
#'                13900,52774,16289,17160,14127,115919,22728,29658,16358,77851,
#'                28386,18132,49422,32030,39155,78946,18520,25765,13728,19688,
#'                13347,12440,14183,95484,26643,20027,14759,13589,14318,18829,
#'                12227,19829,16826,31208,22442,18285,14859,31434,13299,16824,
#'                12843,21389,21024,12491,69633,25811,14946,28967,55172,29773,
#'                17955,40642,19888,35676,35222,14974,12466,46216,86923,48637,
#'                30824,12961,26805,31218,36352,49433,40495,31233,40643,18265,
#'                12316,25125,202241,90877,20903,40353,15008,31908,22554,13634)
#' 
#' mzion:::calc_probi_byvmods(df, nms = "0000000", expt_moverzs, expt_ints, N = 404)
#' 
#' ## 
#' pep <- "LFEEDEREK"
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' theos <- c(343.2543,490.3227,619.3653,748.4079,863.4348,992.4774,
#'            1148.5785,1277.6211,1634.8790,376.2757,505.3183,661.4194,
#'            790.4620,905.4890,1034.5316,1163.5742,1310.6426,1423.7266)
#' 
#' expts <- c(343.25455,490.32318,619.36487,NA,863.43542,NA,NA,
#'            NA,NA,376.27606,505.31924,NA,NA,NA,NA,NA,NA,NA)
#' 
#' names(expts) <- names(theos) <- c(nms, rev(nms))
#' 
#' ith <- c(1,2,3,5,10,11)
#' iex <- c(39,57,73,92,41,59)
#' ith <- as.integer(ith)
#' iex <- as.integer(iex)
#' m <- length(ith)
#' 
#' df <- list(theo = theos, expt = expts, ith = ith, iex = iex, m = m)
#' 
#' expt_moverzs <- c(115.08694,120.08116,126.12807,127.12512,127.13135,
#'                   128.12845,128.13469,129.13177,129.13800,130.13510,
#'                   130.14134,131.13840,132.14159,136.07581,158.09254,
#'                   173.15001,175.11917,175.15672,176.15985,186.15297,
#'                   188.15988,229.16661,230.17041,231.17406,248.18094,
#'                   249.18439,251.11345,269.12421,273.21259,287.19257,
#'                   291.17169,301.20755,301.24423,305.16913,315.25937,
#'                   322.70038,331.21384,331.71515,343.25455,361.71152,
#'                   376.27606,377.24023,377.27985,386.72955,387.22852,
#'                   395.23312,395.73502,396.23633,406.26596,407.27014,
#'                   433.20410,453.24899,470.30963,470.81094,472.27225,
#'                   487.30685,490.32318,499.81451,505.31924,517.77081,
#'                   519.27625,519.77734,551.32843,551.63458,551.82617,
#'                   551.86658,551.96796,552.32220,552.36621,591.37250,
#'                   592.34991,609.33093,619.36487,633.33594,634.33850,
#'                   635.34485,636.34497,644.39349,648.31799,655.82574,
#'                   706.39349,740.87518,741.37933,747.88177,748.39014,
#'                   748.88788,749.38904,749.89026,750.39099,791.45404,
#'                   819.47766,863.43542,864.43774,943.52075,944.51770,
#'                   945.51129,946.51874,947.52203,948.52100,1342.65393)
#' 
#' expt_ints <- c(8508,28501,194673,190890,160012,223742,147623,260388,198943,164580,
#'                149024,165036,7306,21035,8098,11560,35686,8953,15544,10771,
#'                17617,7660,227328,16612,88320,10075,8968,9123,20253,7351,
#'                49432,11202,8013,65251,21790,11174,48669,10941,40270,24887,
#'                81536,9488,9709,28670,8442,7151,185877,54537,9754,7176,
#'                10297,12091,42045,9021,14507,11520,27411,7991,45770,7927,
#'                18996,9520,111456,46109,49064,59457,41824,36185,9619,9065,
#'                8708,21105,26995,13041,8674,11169,14089,7697,9047,10710,
#'                12763,7125,7321,7483,31612,18050,26134,31474,11331,8789,
#'                8672,65913,25918,17819,34367,15767,15221,14225,7419,7284)
#' 
#' mzion:::calc_probi_byvmods(df, nms = "0000000", expt_moverzs, expt_ints, N = 434)
#' }
calc_probi_byvmods <- function (df, nms, expt_moverzs, expt_ints, 
                                N, type_ms2ions = "by", topn_ms2ions = 100L, 
                                ppm_ms2 = 20L, soft_secions = FALSE, 
                                burn_ins = c(1:2), min_ms2mass = 115L, 
                                d2 = 1E-5, index_mgf_ms2 = FALSE, 
                                tally_ms2ints = TRUE, digits = 4L) 
{
  df_theo <- df$theo
  m <- length(df_theo)

  ## df2
  tt2  <- add_seions(df_theo, type_ms2ions = type_ms2ions, digits = digits)
  df2  <- match_ex2th2(expt_moverzs, tt2, min_ms2mass, d2, index_mgf_ms2)
  ith2 <- df2[["ith"]]
  iex2 <- df2[["iex"]]

  ## 1. int2 (secondary intensities)
  len <- length(df2[["expt"]])
  df2[["int"]] <- rep(NA_real_, len)
  df2[["int"]][ith2] <- expt_ints[iex2] # works if iex2 contains NA
  
  facs <- rep(seq_len(len/m), each = m)
  int2 <- .Internal(split(df2[["int"]], as.factor(facs)))
  int2 <- Reduce(`%+%`, int2)
  
  # df2[["int"]]
  # NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
  # NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
  # NA     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA  48669 185877  12091   7927     NA  10710     NA
  # NA     NA   7697     NA     NA     NA     NA     NA     NA     NA     NA  11174   8442     NA     NA     NA     NA     NA
  # NA  11520     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA  28670     NA     NA     NA     NA     NA
  # 
  # int2
  #  0  11520   7697      0      0      0      0      0      0      0      0  59843 222989  12091   7927      0  10710      0


  ## 2. y 
  ith <- df[["ith"]]
  iex <- df[["iex"]]
  df[["int"]] <- rep(NA_integer_, m)
  df[["int"]][ith] <- expt_ints[iex]

  nudbl <- rep(NA_real_, topn_ms2ions)
  nuint <- rep(NA_integer_, topn_ms2ions)
  y <- list(expt = expt_moverzs, int = expt_ints, theo = nudbl, idx = nuint, int2 = nuint)
  y[["theo"]][iex] <- df_theo[ith]
  y[["idx"]][iex] <- ith
  
  # if (tally_ms2ints) {
    ## 3. join `int2` to `y`
    y_idx  <- y[["idx"]]
    ok_iex <- .Internal(which(!is.na(y_idx)))
    y_ith  <- y_idx[ok_iex]
    y[["int2"]][ok_iex] <- int2[y_ith]
    
    ## 4. collapses `int2` to `int`
    y[["int"]] <- y[["int"]] %+% y[["int2"]]
    y[["idx"]] <- y[["int2"]] <- NULL
  # }

  ###
  if (soft_secions) {
    ok_int2 <- .Internal(which(int2 > 0L & is.na(df[["int"]])))
    
    if (length(ok_int2)) {
      ok_iex2 <- iex2[match(ok_int2, ith2 %% m)]
      ok_iex2 <- ifelse(is.na(ok_iex2), m, ok_iex2)
      y[["int"]][ok_iex2] <- int2[ok_int2]
      # need to adjust the theoreticals, e.g., (df_theo[ok_int2] + 1.0078)/2; 
      # for simplicity use expt_moverzs, as long as they are not NA
      y[["theo"]][ok_iex2] <- expt_moverzs[ok_iex2]
    }
  }
  ###

  ## 5. arrange by "-int"
  ord_int <- order(y[["int"]], decreasing = TRUE, method = "radix", na.last = TRUE)
  y_theo  <- y[["theo"]][ord_int]
  maxi    <- .Internal(which(!is.na(y_theo)))
  maxi    <- maxi[length(maxi)]
  y_theo  <- y_theo[1:maxi]

  ## 6. mutate(k = row_number(), x = k - cumsum(is.na(theo)))
  k <- 1:maxi
  x <- k - cumsum(is.na(y_theo))
  
  ## 7. filter(!is.na(theo))
  # note: x <= k <= x + n
  ok_y <- !is.na(y_theo)
  k    <- k[ok_y]
  x    <- x[ok_y]

  ## 8. Probability
  # (to have sufficient counts of noise)
  # (also guaranteed n > 0L)
  # n <- max(n, topn_ms2ions + k[length(k)])
  
  # excludes unstable burn-in scores
  # burn_ins is different to min_n_ms2
  #   min_n_ms2 - guard against low-quality specta with low number of matches
  #   burn_ins  - guard against unstable scores
  x_ <- x[-burn_ins]
  k_ <- k[-burn_ins]
  
  if (length(x_)) {
    prs <- stats::dhyper(x = x_, m = m, n = N, k = k_)
    pr  <- min(prs, na.rm = TRUE)
  }
  else {
    pr <- .5
  }

  ## outputs
  list(pep_ivmod = nms, 
       pep_prob  = pr, 
       pri_matches = list(df), 
       sec_matches = list(df2))
}


#' Helper for score calculations
#'
#' By peptides.
#' 
#' @param nms The names (of peptides).
#' @inheritParams calc_probi
#' @import dplyr
#' @importFrom purrr map2
#' @importFrom tibble tibble
calc_probi_bypep <- function (mts, nms, expt_moverzs, expt_ints, 
                              N, type_ms2ions = "by", topn_ms2ions = 100L, 
                              ppm_ms2 = 20L, soft_secions = FALSE, 
                              min_ms2mass = 115L, d2 = 1E-5, 
                              index_mgf_ms2 = FALSE, tally_ms2ints = TRUE, 
                              digits = 4L) 
{
  ## for different positions: $TNLAMMR$`0000500`, $TNLAMMR$`0000050`
  #    the same `pep_seq`, `pep_calc_mr` for different mod positions
  #    different `pep_ivmod`, `pep_prob`, `pri_matches`, `sec_matches`
  
  # NAMES are hexcodes: 0000000
  res <- mapply(calc_probi_byvmods, 
                mts, names(mts), 
                MoreArgs = list(
                  expt_moverzs = expt_moverzs, 
                  expt_ints = expt_ints, 
                  N = N, 
                  type_ms2ions = type_ms2ions, 
                  topn_ms2ions = topn_ms2ions, 
                  ppm_ms2 = ppm_ms2, 
                  soft_secions = soft_secions, 
                  burn_ins = c(1:2),
                  min_ms2mass = min_ms2mass, 
                  d2 = d2, 
                  index_mgf_ms2 = index_mgf_ms2, 
                  tally_ms2ints = tally_ms2ints, 
                  digits = digits
                ), 
                SIMPLIFY = FALSE,
                USE.NAMES = TRUE)
  
  theo_ms1 <- attr(mts, "theo_ms1")
  
  len <- length(res)
  out <- vector("list", len)
  
  for (i in 1:len) {
    res_i <- res[[i]]
    
    out[[i]] <- list(
      pep_seq = nms,
      theo_ms1 = theo_ms1, 
      pep_ivmod = res_i$pep_ivmod, 
      pep_prob = res_i$pep_prob, 
      pri_matches = res_i["pri_matches"], 
      sec_matches = res_i["sec_matches"]
    )
  }
  
  out
}


#' Helper for score calculations
#'
#' @param mts Nested data frame of \code{theo} and matched \code{expt} m-over-z.
#' @param expt_moverzs Nested list of match and unmatched experimental m-over-z.
#' @param expt_ints Nested list of match and unmatched experimental intensity.
#' @param N Numeric; the number of MS2 features in an MGF query.
#' @param d2 Bin width in ppm divided by 1E6.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
#' @import dplyr
#' @importFrom purrr map
calc_probi <- function (mts, expt_moverzs, expt_ints, 
                        N, type_ms2ions = "by", topn_ms2ions = 100L, 
                        ppm_ms2 = 20L, soft_secions = FALSE, 
                        min_ms2mass = 115L, d2 = 1E-5, index_mgf_ms2 = FALSE, 
                        tally_ms2ints = TRUE, digits = 4L) 
{
  out <- mapply(
    calc_probi_bypep, 
    mts, names(mts), 
    MoreArgs = list(
      expt_moverzs = expt_moverzs, 
      expt_ints = expt_ints, 
      N = N, 
      type_ms2ions = type_ms2ions, 
      topn_ms2ions = topn_ms2ions, 
      soft_secions = soft_secions, 
      ppm_ms2 = ppm_ms2, 
      min_ms2mass = min_ms2mass, 
      d2 = d2, 
      index_mgf_ms2 = index_mgf_ms2, 
      tally_ms2ints = tally_ms2ints, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE)
  
  out <- .Internal(unlist(out, recursive = FALSE, use.names = FALSE))
}


#' Calculates peptide scores by single MGF entries (scan numbers).
#' 
#' Each entry corresponds to a row in \code{ion_matches.rds}.
#' 
#' @param entry A row of data.
#' @param d2 Bin width in ppm divided by 1E6.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
#' @import purrr
scalc_pepprobs <- function (entry, topn_ms2ions = 100L, type_ms2ions = "by", 
                            ppm_ms2 = 20L, soft_secions = FALSE, 
                            min_ms2mass = 115L, d2 = 1E-5, index_mgf_ms2 = FALSE, 
                            tally_ms2ints = TRUE, digits = 4L) 
{
  # only one experimental set of values and thus `[[1]]`
  expt_moverzs <- entry[["pep_ms2_moverzs"]][[1]]
  expt_ints <- entry[["pep_ms2_ints"]][[1]]

  ## matches between theoreticals and experimentals
  
  # [[1]] --- `entry$matches` (always at level-one and can be unlisted)
  # [[1]]$AMMASIGR --- `entry$matches[[1]]` (1:i peptides)
  # [[1]]$AMMASIGR$`00500000` --- `entry$matches[[1]][[1]]` (1:j positions)
  # A tibble: 18 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  301.   301.
  # 
  # [[1]]$TNLAMMR
  # [[1]]$TNLAMMR$`0000500`
  # A tibble: 16 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  331.   331.
  # 
  # [[1]]$TNLAMMR$`0000050`
  # A tibble: 16 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  331.   331.
  
  # (flattens by one level as is a list-column)
  mts <- entry[["matches"]][[1]]
  topn_ms2ions <- min(topn_ms2ions, entry$pep_n_ms2[[1]])
  N <- min(topn_ms2ions * 5L, 500L)
  
  out <- calc_probi(mts = mts, 
                    expt_moverzs = expt_moverzs, 
                    expt_ints = expt_ints, 
                    N = N, 
                    type_ms2ions = type_ms2ions, 
                    topn_ms2ions = topn_ms2ions, 
                    ppm_ms2 = ppm_ms2, 
                    soft_secions = soft_secions, 
                    min_ms2mass = min_ms2mass, 
                    d2 = d2, 
                    index_mgf_ms2 = index_mgf_ms2, 
                    tally_ms2ints = tally_ms2ints, 
                    digits = digits)

  uniq_id <- .Internal(unlist(entry$uniq_id, recursive = FALSE, use.names = FALSE))
  
  out <- lapply(out, function (x) {
    x[["uniq_id"]] <- uniq_id
    x
  })

  invisible(out)
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @param df Resulted data from ion matches.
#' @param d2 Bin width in ppm divided by 1E6.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
calc_pepprobs_i <- function (df, topn_ms2ions = 100L, type_ms2ions = "by", 
                             ppm_ms2 = 20L, soft_secions = FALSE, 
                             out_path = "~/mzion/outs", 
                             min_ms2mass = 115L, d2 = 1E-5, 
                             index_mgf_ms2 = FALSE, tally_ms2ints = TRUE, 
                             digits = 4L) 
{
  n_rows <- nrow(df)
  
  if (n_rows) {
    df <- split.data.frame(df, seq_len(n_rows)) 
    
    df <- lapply(df, scalc_pepprobs, 
                 topn_ms2ions = topn_ms2ions, 
                 type_ms2ions = type_ms2ions, 
                 ppm_ms2 = ppm_ms2, 
                 soft_secions = soft_secions, 
                 min_ms2mass = min_ms2mass, 
                 d2 = d2, 
                 index_mgf_ms2 = index_mgf_ms2, 
                 tally_ms2ints = tally_ms2ints, 
                 digits = digits)
    
    df <- .Internal(unlist(df, recursive = FALSE, use.names = FALSE))
    df <- dplyr::bind_rows(df)
  } 
  else {
    df <- data.frame(
      pep_seq = as.character(), 
      pep_ivmod = as.character(), 
      pep_prob = as.numeric(), 
      pri_matches = list(), 
      sec_matches = list(), 
      pep_scan_num = as.integer())
  }
  
  invisible(df)
}


#' Calculates the scores of peptides.
#' 
#' @param tally_ms2ints Logical; tally MS2 intensities or not.
#' @inheritParams matchMS
#' @import parallel
calc_pepscores <- function (topn_ms2ions = 100L, type_ms2ions = "by", 
                            target_fdr = 0.01, 
                            
                            # to be deleted later
                            fdr_type = "protein", 
                            
                            min_len = 7L, max_len = 40L, ppm_ms2 = 20L, 
                            soft_secions = FALSE, 
                            out_path = "~/mzion/outs", 
                            min_ms2mass = 115L, index_mgf_ms2 = FALSE, 
                            tally_ms2ints = TRUE, 
                            
                            mgf_path, maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                            maxn_vmods_sitescombi_per_pep = 64L, minn_ms2 = 6L, 
                            ppm_ms1 = 20L, quant = "none", 
                            ppm_reporters = 10L, fasta, acc_type, acc_pattern, 
                            fixedmods, varmods, 
                            enzyme = "trypsin_p", maxn_fasta_seqs = 200000L, 
                            maxn_vmods_setscombi = 64L, 
                            add_ms2theos = FALSE, add_ms2theos2 = FALSE, 
                            add_ms2moverzs = FALSE, add_ms2ints = FALSE,
                            digits = 4L) 
{
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    }, 
    add = TRUE
  )
  
  ## Check priors
  pat_i  <- "^ion_matches_"
  list_i <- find_targets(out_path, pattern = pat_i)$files
  len_i  <- length(list_i)
  
  if (!len_i) 
    stop("No target results with pattern '", pat_i, "'.")
  
  ## Check cached current
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  fml_nms <- names(formals(fun))
  
  args_except <- c("fdr_type")
  
  fml_incl <- if (length(args_except))
    fml_nms[!fml_nms %in% args_except]
  else
    fml_nms

  message("[x] For reprocessing (with new score function) ", 
          "delete cached 'pepscores_[...]' and 'calc_pepscores.rda'.\n")
  
  cache_pars <- find_callarg_vals(time = NULL, 
                                  path = file.path(out_path, "Calls"), 
                                  fun = paste0(fun, ".rda"), 
                                  args = fml_incl) 
  
  cache_pars <- cache_pars[sort(names(cache_pars))]
  call_pars  <- mget(fml_incl, envir = fun_env, inherits = FALSE) 
  call_pars  <- call_pars[sort(names(call_pars))]
  
  if (identical(cache_pars, call_pars)) {
    ok_scores <- find_targets(out_path, pattern = "^pepscores_")
    listsc <- ok_scores$files
    lensc  <- length(listsc)
    
    if (lensc == len_i) {
      message("Found cached 'pepscores_[...]'.")
      .savecall <- FALSE
      return(NULL)
    } 
    else {
      message("Recalculating peptide scores (not all 'pepscores_' found).")
      rm(list = c("listsc"))
    }
  } 
  else {
    message("Calculating peptide scores.")
    dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
    rm(list = c("cache_pars", "call_pars"))
  }
  
  d2 <- calc_threeframe_ppm(ppm_ms2) * 1E-6

  for (fi in list_i)
    calcpepsc(file = fi, 
              topn_ms2ions = topn_ms2ions, 
              type_ms2ions = type_ms2ions, 
              ppm_ms2 = ppm_ms2, 
              soft_secions = soft_secions, 
              out_path = out_path, 
              min_ms2mass = min_ms2mass, 
              d2 = d2, 
              index_mgf_ms2 = index_mgf_ms2, 
              tally_ms2ints = tally_ms2ints, 
              add_ms2theos = add_ms2theos, 
              add_ms2theos2 = add_ms2theos2, 
              add_ms2moverzs = add_ms2moverzs, 
              add_ms2ints = add_ms2ints,
              # slower with 48 cores
              n_cores = detect_cores(16L), 
              digits = digits)
    
  .savecall <- TRUE
  
  invisible(NULL)
}


#' Find the index or name of decoy results.
#' 
#' @param pattern The pattern of files.
#' @inheritParams matchMS
find_decoy <- function (out_path, pattern = "^ion_matches_") 
{
  pat <- paste0(pattern, "rev_[0-9]+\\.rds$")
  list_d <- list.files(path = file.path(out_path, "temp"), pattern = pat)

  if (length(list_d) > 1L) {
    warning("More than one decoy results found; ", 
            "'", list_d[1], "' will be used.", call. = FALSE)
    list_d <- list_d[1]
  }
  
  nms_d <- gsub(paste0(pattern, "(rev_[0-9]+)\\.rds$"), "\\1", list_d)
  
  list(idxes = nms_d, files = list_d)
}


#' Find the index or name of target results.
#' 
#' @param out_path An output path.
#' @param pattern The pattern of files.
find_targets <- function (out_path, pattern = "^ion_matches_") 
{
  pat <- paste0(pattern, "[0-9]+\\.rds$")
  list_t <- list.files(path = file.path(out_path, "temp"), pattern = pat)
  
  nms_t <- gsub(paste0(pattern, "(\\d+)\\.rds$"), "\\1", list_t)
  ord <- order(as.integer(nms_t))
  nms_t <- nms_t[ord]
  list_t <- list_t[ord]
  
  list(idxes = nms_t, files = list_t)
}


#' Helper of \link{calc_pepscores}.
#' 
#' @param file A file name of \code{ion_matches_}.
#' @param d2 Bin width in ppm divided by 1E6.
#' @param n_cores The number of CPU cores.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
calcpepsc <- function (file, topn_ms2ions = 100L, type_ms2ions = "by", 
                       ppm_ms2 = 20L, soft_secions = FALSE, out_path = NULL, 
                       min_ms2mass = 115L, d2 = 1E-5, index_mgf_ms2 = FALSE,
                       tally_ms2ints = TRUE, 
                       add_ms2theos = FALSE, add_ms2theos2 = FALSE, 
                       add_ms2moverzs = FALSE, add_ms2ints = FALSE, 
                       n_cores = 16L, digits = 4L) 
{
  message("\tModule: ", file)
  
  idx <- gsub("^ion_matches_(.*)\\.rds$", "\\1", file)
  file_lt <- file.path(out_path, "temp", paste0("list_table_", idx, ".rds"))
  file_sc <- file.path(out_path, "temp", paste0("pepscores_",  idx, ".rds"))

  cols_a  <- c("pep_scan_num", "raw_file")
  cols_b  <- c("pep_ms2_moverzs", "pep_ms2_ints", "pri_matches", "sec_matches")
  cols_lt <- c(cols_a, cols_b)
  
  cols_sc <- c("pep_seq", "pep_n_ms2", "pep_scan_title", "pep_exp_mz", "pep_exp_mr", 
               "pep_tot_int", "pep_exp_z", "pep_ret_range", "pep_scan_num", "raw_file", 
               "pep_mod_group", "pep_frame", "pep_fmod", "pep_vmod", "pep_isdecoy", 
               "pep_calc_mr", "pep_ivmod", "pep_prob", "pep_len", 
               "pep_ms2_moverzs", "pep_ms2_ints", 
               "pep_ms2_theos", "pep_ms2_theos2", 
               "pep_ms2_exptints", "pep_ms2_exptints2", 
               "pep_n_matches", "pep_n_matches2", "pep_ms2_deltas", 
               "pep_ms2_ideltas", "pep_ms2_deltas2", "pep_ms2_ideltas2", 
               "pep_ms2_deltas_mean", "pep_ms2_deltas_sd", 
               
               # for localization scores
               "pep_ms2_ideltas.")
  
  df <- qs::qread(file.path(out_path, "temp", file))
  n_rows <- nrow(df)
  
  df <- dplyr::rename(df, 
                      pep_ret_range = ret_time, 
                      pep_scan_title = scan_title,
                      pep_exp_mz = ms1_moverz, 
                      pep_n_ms2 = ms2_n, 
                      pep_exp_mr = ms1_mass, 
                      pep_tot_int = ms1_int, 
                      pep_scan_num = scan_num, 
                      pep_exp_z = ms1_charge, 
                      pep_ms2_moverzs = ms2_moverz, 
                      pep_ms2_ints = ms2_int, 
                      pep_frame = frame)

  if (!n_rows) {
    dfa <- data.frame(matrix(ncol = length(cols_lt), nrow = 0L))
    colnames(dfa) <- cols_lt
    qs::qsave(dfa, file_lt, preset = "fast")
    
    dfb <- data.frame(matrix(ncol = length(cols_sc), nrow = 0L))
    colnames(dfb) <- cols_sc
    qs::qsave(dfb, file_sc, preset = "fast")
    
    return (dfb)
  }

  tempdir <- create_dir(file.path(out_path, "sc_temp"))
  
  df[["uniq_id"]] <- paste(df[["pep_scan_num"]], df[["raw_file"]], sep = "@")
  esscols <- c("pep_ms2_moverzs", "pep_ms2_ints", "matches", "pep_n_ms2", "uniq_id")
  path_df2 <- file.path(tempdir, "df2_sc_temp.rda")
  df2 <- df[, -which(names(df) %in% esscols), drop = FALSE]
  qs::qsave(df2, path_df2, preset = "fast")
  df <- df[, esscols, drop = FALSE]
  rm(list = "df2")
  gc()

  # otherwise, chunksplit return NULL
  #   -> res[[i]] <- NULL 
  #   -> length(res) shortened by 1
  
  if (n_rows <= 5000L) {
    probs <- calc_pepprobs_i(
      df,
      topn_ms2ions = topn_ms2ions, 
      type_ms2ions = type_ms2ions, 
      ppm_ms2 = ppm_ms2,
      soft_secions = soft_secions, 
      out_path = out_path, 
      min_ms2mass = min_ms2mass, 
      d2 = d2, 
      index_mgf_ms2 = index_mgf_ms2, 
      tally_ms2ints = tally_ms2ints, 
      digits = digits)
  }
  else {
    path_df <- file.path(tempdir, "df_sc_temp.rda")
    max_rows <- 100000L

    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterExport(cl, list("calc_pepprobs_i", "scalc_pepprobs", 
                                     "calc_probi", "calc_probi_bypep", 
                                     "calc_probi_byvmods", "add_seions", 
                                     "find_ppm_outer_bycombi", "match_ex2th2", 
                                     "add_primatches"), 
                            envir = environment(mzion::matchMS))

    if (n_rows > max_rows) {
      dfs <- suppressWarnings(chunksplit(df, ceiling(n_rows/max_rows), "row"))
      len <- length(dfs)
      nms <- paste0("sc", 1:len, ".rds")
      mapply(qs::qsave, dfs, file.path(tempdir, nms), MoreArgs = list(preset = "fast"))
      rm(list = "dfs")
      gc()
      
      df <- df[, -which(names(df) == "matches"), drop = FALSE]
      qs::qsave(df, path_df, preset = "fast")
      rm(list = "df", envir = environment())
      gc()
      
      probs <- vector("list", len)
      
      for (i in seq_len(len)) {
        dfi <- suppressWarnings(
          # * 4L of smaller hashes for some slow phospho scoring
          chunksplit(qs::qread(file.path(tempdir, nms[[i]])), n_cores, "row"))

        probs[[i]] <- parallel::clusterApply(cl, dfi, 
                                             calc_pepprobs_i, 
                                             topn_ms2ions = topn_ms2ions, 
                                             type_ms2ions = type_ms2ions, 
                                             ppm_ms2 = ppm_ms2,
                                             soft_secions = soft_secions, 
                                             out_path = out_path, 
                                             min_ms2mass = min_ms2mass, 
                                             d2 = d2, 
                                             index_mgf_ms2 = index_mgf_ms2, 
                                             tally_ms2ints = tally_ms2ints, 
                                             digits = digits)
        probs[[i]] <- dplyr::bind_rows(probs[[i]])
      }
      
      parallel::stopCluster(cl)
      rm(list = c("dfi"))
      gc()
      
      probs <- dplyr::bind_rows(probs)
      
      df <- qs::qread(path_df)
    }
    else {
      dfs <- suppressWarnings(chunksplit(df, n_cores, "row"))
      
      # a case that `chunksplit` did not successfully split
      if (is.data.frame(dfs)) {
        probs <- calc_pepprobs_i(
          dfs,
          topn_ms2ions = topn_ms2ions, 
          type_ms2ions = type_ms2ions, 
          ppm_ms2 = ppm_ms2,
          soft_secions = soft_secions, 
          out_path = out_path, 
          min_ms2mass = min_ms2mass, 
          d2 = d2, 
          index_mgf_ms2 = index_mgf_ms2, 
          tally_ms2ints = tally_ms2ints, 
          digits = digits)
      }
      else {
        probs <- parallel::clusterApply(cl, dfs, 
                                        calc_pepprobs_i, 
                                        topn_ms2ions = topn_ms2ions, 
                                        type_ms2ions = type_ms2ions, 
                                        ppm_ms2 = ppm_ms2,
                                        soft_secions = soft_secions, 
                                        out_path = out_path, 
                                        min_ms2mass = min_ms2mass, 
                                        d2 = d2, 
                                        index_mgf_ms2 = index_mgf_ms2, 
                                        tally_ms2ints = tally_ms2ints, 
                                        digits = digits)
        
        parallel::stopCluster(cl)
        rm(list = c("dfs"))
        gc()
        
        probs <- dplyr::bind_rows(probs)
      }
    }
  }
  
  ## Reassemble `df`
  if ("matches" %in% names(df)) 
    df <- df[, -which(names(df) == "matches"), drop = FALSE]
  
  df2 <- qs::qread(path_df2)
  df <- dplyr::bind_cols(df, df2)
  rm(list = c("df2", "path_df2"))
  gc()
  
  df <- quick_rightjoin(df, probs, "uniq_id")
  rm(list = c("probs"))
  gc()
  
  df <- df[, -which(names(df) == "uniq_id"), drop = FALSE]
  df <- post_pepscores(df)
  qs::qsave(df[, cols_lt, drop = FALSE], file_lt, preset = "fast")
  
  ## scores
  n_rows <- nrow(df)
  max_rows <- 100000L
  n_chunks <- n_rows %/% max_rows + 1L

  if (n_chunks > 1L)
    mapply(function (x, i) qs::qsave(x, 
      file.path(tempdir, paste0("tempscores_", idx, "_", i, ".rds")), 
      preset = "fast"), chunksplit(df, n_chunks, "row"), 1:n_chunks)
  else
    qs::qsave(df, file.path(tempdir, paste0("tempscores_", idx, "_1", ".rds")), 
              preset = "fast")

  invisible(NULL)
}


#' Helper of \link{add_primatches}
#' 
#' @inheritParams matchMS
hadd_primatches <- function (out_path = NULL, 
                             add_ms2theos = FALSE, add_ms2theos2 = FALSE, 
                             add_ms2moverzs = FALSE, add_ms2ints = FALSE, 
                             index_mgf_ms2 = FALSE) 
{
  # the same as those in calcpepsc
  cols_sc <- c("pep_seq", "pep_n_ms2", "pep_scan_title", "pep_exp_mz", "pep_exp_mr", 
               "pep_tot_int", "pep_exp_z", "pep_ret_range", "pep_scan_num", "raw_file", 
               "pep_mod_group", "pep_frame", "pep_fmod", "pep_vmod", "pep_isdecoy", 
               "pep_calc_mr", "pep_ivmod", "pep_prob", "pep_len", 
               "pep_ms2_moverzs", "pep_ms2_ints", 
               "pep_ms2_theos", "pep_ms2_theos2", 
               "pep_ms2_exptints", "pep_ms2_exptints2", 
               "pep_n_matches", "pep_n_matches2", "pep_ms2_deltas", 
               "pep_ms2_ideltas", "pep_ms2_deltas2", "pep_ms2_ideltas2", 
               "pep_ms2_deltas_mean", "pep_ms2_deltas_sd", 
               
               # for localization scores
               "pep_ms2_ideltas.")
  
  tempdir <- file.path(out_path, "sc_temp")
  files <- list.files(path = tempdir, pattern = "^tempscores_\\d+_\\d+\\.rds$")
  n_cores <- min(detect_cores(48L), length(files))
  
  message("Adding theoretical MS2 m/z and intensity values: ", Sys.time())
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("add_primatches"), 
                          envir = environment(mzion::matchMS))
  parallel::clusterApplyLB(cl, files, 
                           add_primatches, 
                           tempdir = tempdir, 
                           add_ms2theos = add_ms2theos, 
                           add_ms2theos2 = add_ms2theos2, 
                           add_ms2moverzs = add_ms2moverzs, 
                           add_ms2ints = add_ms2ints, 
                           index_mgf_ms2 = index_mgf_ms2)
  parallel::stopCluster(cl)
  
  ms_files <- gsub("^tempscores", "tempms2info", files)
  # ms_files <- sort(ms_files)
  idxes <- as.integer(gsub("^tempms2info_(\\d+).*", "\\1", ms_files))
  fracs <- as.integer(gsub("^tempms2info_\\d+_(\\d+).*", "\\1", ms_files))
  fracs <- split(fracs, idxes)
  ms_files <- split(ms_files, idxes)
  idxes <- split(idxes, idxes)
  
  ords <- lapply(fracs, order)
  ms_files <- mapply(function (x, y) x[y], ms_files, ords)
  idxes <- names(ms_files)
  
  mapply(function (fis, idx) {
    df <- lapply(fis, function (x) qs::qread(file.path(tempdir, x)))
    df <- dplyr::bind_rows(df)
    df[["pep_isdecoy"]] <- ifelse(is.na(df[["pep_ivmod"]]), TRUE, FALSE)
    df <- dplyr::rename(df, pep_calc_mr = theo_ms1)
    
    if (!all(cols_sc %in% names(df)))
      stop("Developer needs to update the columns of peptide scores.")
    
    df <- df[, cols_sc, drop = FALSE]
    
    qs::qsave(df, file.path(out_path, "temp", paste0("pepscores_", idx, ".rds")), 
              preset = "fast")
  }, ms_files, idxes)
  
  message("Completed theoretical MS2 m/z and intensity values: ", Sys.time())
  unlink(tempdir, recursive = TRUE)
  
  invisible(NULL)
}


#' Adds sequences of primary and secondary matches.
#'
#' Applied to both targets and decoys as feature "pep_ms2_deltas_mean" may be
#' used in SVM-Percolator.
#'
#' @param file The file name of a temporary score file.
#' @param tempdir The temporary directory for score results.
#' @inheritParams matchMS
add_primatches <- function (file = NULL, tempdir = NULL, add_ms2theos = FALSE, 
                            add_ms2theos2 = FALSE, add_ms2moverzs = FALSE, 
                            add_ms2ints = FALSE, index_mgf_ms2 = FALSE) 
{
  df <- qs::qread(file.path(tempdir, file))
  
  df <- dplyr::mutate(df, 
                      pep_ms2_moverzs = NA_character_, 
                      pep_ms2_ints = NA_character_, 
                      pep_ms2_theos = NA_character_, 
                      pep_ms2_theos2 = NA_character_, 
                      
                      pep_ms2_exptints = NA_character_, 
                      pep_ms2_exptints2 = NA_character_, 
                      
                      pep_n_matches = NA_integer_, 
                      pep_n_matches2 = NA_integer_, 
                      
                      pep_ms2_deltas = NA_character_, 
                      pep_ms2_ideltas = NA_character_,
                      pep_ms2_deltas2 = NA_character_, 
                      pep_ms2_ideltas2 = NA_character_, 
                      
                      pep_ms2_deltas_mean = NA_real_, 
                      pep_ms2_deltas_sd = NA_real_, 
                      
                      pep_ms2_ideltas. = NA_integer_)

  # unlist from list table
  pris <- lapply(df$pri_matches, `[[`, 1)
  secs <- lapply(df$sec_matches, `[[`, 1)

  len <- length(pris)
  p1s. <- m2s <- m1s <- iys2 <- iys1 <- sd1s <- me1s <- p2s <- d2s <- p1s <- d1s <- 
    vector("list", len)

  for (i in 1:len) {
    mt1 <- pris[[i]]
    th1 <- mt1[["theo"]]
    ex1 <- mt1[["expt"]]
    iy1 <- mt1[["int"]]
    mt2 <- secs[[i]]
    th2 <- mt2[["theo"]]
    ex2 <- mt2[["expt"]]
    iy2 <- mt2[["int"]]

    ps1 <- mt1[["ith"]]
    ps2 <- mt2[["ith"]]
    
    ds1 <- (ex1[ps1] - th1[ps1]) * 1E3
    ds2 <- (ex2[ps2] - th2[ps2]) * 1E3
    me1 <- mean(ds1)
    sd1 <- stats::sd(ds1)
    
    # delayed rounding
    ds1 <- round(ds1, digits = 2L)
    ds2 <- round(ds2, digits = 2L)
    me1 <- round(me1, digits = 2L)
    sd1 <- round(sd1, digits = 2L)
    
    d1s[[i]] <- .Internal(paste0(list(ds1), collapse = ";", recycle0 = FALSE))
    d2s[[i]] <- .Internal(paste0(list(ds2), collapse = ";", recycle0 = FALSE))
    
    p1s[[i]]  <- .Internal(paste0(list(ps1), collapse = ";", recycle0 = FALSE))
    p2s[[i]]  <- .Internal(paste0(list(ps2), collapse = ";", recycle0 = FALSE))
    iys1[[i]] <- .Internal(paste0(list(iy1[ps1]), collapse = ";", recycle0 = FALSE))
    iys2[[i]] <- .Internal(paste0(list(iy2[ps2]), collapse = ";", recycle0 = FALSE))
    
    me1s[[i]] <- me1
    sd1s[[i]] <- sd1
    
    m1s[[i]]  <- mt1$m
    m2s[[i]]  <- mt2$m
    p1s.[[i]] <- ps1
  }
  
  if (index_mgf_ms2) {
    # need to convert theoretical m/z to integers; 
    # even so the resolution is limited by bin with (e.g. 10 ppm)
    pep_ms2_deltas <- NA_character_
    pep_ms2_deltas2 <- NA_character_
    pep_ms2_deltas_mean <- NA_real_
    pep_ms2_deltas_sd <- NA_real_
  }
  else {
    df[["pep_ms2_deltas"]] <- do.call(rbind, d1s)
    df[["pep_ms2_deltas2"]] <- do.call(rbind, d2s)
    df[["pep_ms2_deltas_mean"]] <- do.call(rbind, me1s)
    df[["pep_ms2_deltas_sd"]] <- do.call(rbind, sd1s)
  }

  df[["pep_ms2_ideltas"]] <- do.call(rbind, p1s)
  df[["pep_ms2_ideltas2"]] <- do.call(rbind, p2s)
  df[["pep_n_matches"]] <- do.call(rbind, m1s)
  df[["pep_n_matches2"]] <- do.call(rbind, m2s)
  df[["pep_ms2_exptints"]] <- do.call(rbind, iys1)
  df[["pep_ms2_exptints2"]] <- do.call(rbind, iys2)
  
  df[["pep_ms2_ideltas."]] <- p1s.
  
  if (add_ms2theos) df$pep_ms2_theos <- collapse_vecs(lapply(pris, `[[`, "theo"))
  if (add_ms2theos2) df$pep_ms2_theos2 <- collapse_vecs(lapply(secs, `[[`, "theo"))
  if (add_ms2moverzs) df$pep_ms2_moverzs <- collapse_vecs(df$ms2_moverz)
  if (add_ms2ints) df$pep_ms2_ints <- collapse_vecs(df$ms2_int)
  
  qs::qsave(df, file.path(tempdir, gsub("^tempscores", "tempms2info", file)), 
            preset = "fast")

  invisible(NULL)
}


#' Pastes vectors to character strings.
#'
#' @param vecs A list of vectors.
#' @param nm The name of sub list in \code{vecs}.
#' @param sep A separator.
collapse_vecs <- function (vecs, nm = "theo", sep = ";") 
{
  ans <- lapply(vecs, function (v) 
    .Internal(paste0(list(v), collapse = sep, recycle0 = FALSE))
  )
  
  do.call(rbind, ans)
}


#' Cleanups post \link{calc_pepscores}.
#' 
#' @param df A results after pep_scores.
post_pepscores <- function (df) 
{
  df[["pep_scan_num"]] <- as.character(df[["pep_scan_num"]])
  df[["pep_len"]] <- stringr::str_length(df[["pep_seq"]])
  df[["pep_scan_title"]] <- as.character(df[["pep_scan_title"]])
  df[["pep_exp_mz"]] <- as.numeric(df[["pep_exp_mz"]])
  df[["pep_exp_mz"]] <- as.numeric(df[["pep_exp_mz"]])
  df[["pep_tot_int"]] <- as.numeric(df[["pep_tot_int"]])
  df[["pep_exp_z"]] <- as.character(df[["pep_exp_z"]])
  df[["pep_ret_range"]] <- as.integer(df[["pep_ret_range"]])
  df[["pep_n_ms2"]] <- as.integer(df[["pep_n_ms2"]])
  df[["pep_fmod"]] <- as.character(df[["pep_fmod"]])
  df[["pep_vmod"]] <- as.character(df[["pep_vmod"]])
  
  invisible(df)
}


#' Finds the cut-off in peptide scores for a given \code{target_fdr}.
#' 
#' Assume normal distribution for log(decoy_score).
#' 
#' @param td A data frame of target-decoy results at a given peptide length.
#' @inheritParams matchMS
find_pepscore_co1 <- function (td, target_fdr = 0.01) 
{
  target <- dplyr::filter(td, !pep_isdecoy)
  decoy  <- dplyr::filter(td, pep_isdecoy)
  
  nt <- nrow(target)
  nd <- nrow(decoy)
  
  if (nd <= 5L) 
    return(NA_real_)
  
  n <- nt + nd
  lambt <- nt / n
  lambd <- 1 - lambt
  
  vecd <- log2(decoy$pep_score)
  sigmad <- sd(vecd, na.rm = TRUE)
  mud <- mean(vecd, na.rm = TRUE)
  
  if (is.na(sigmad)) 
    return(NA_real_)
  
  xs <- seq(mud + sigmad, mud + 3*sigmad, 0.014355293)
  
  for (i in seq_along(xs)) {
    y <- (1 - pnorm(xs[i], mud, sigmad)) * nd / n
    
    if (y <= target_fdr) 
      break
  }
  
  2^(xs[i])
}


#' Finds the cut-off in peptide scores for a given \code{target_fdr}.
#' 
#' Assume log-normal for decoy scores.
#' 
#' @param td A data frame of target-decoy results at a given peptide length.
#' @inheritParams matchMS
find_pepscore_co2 <- function (td, target_fdr = 0.01) 
{
  target <- dplyr::filter(td, !pep_isdecoy)
  decoy  <- dplyr::filter(td, pep_isdecoy)
  
  nt <- nrow(target)
  nd <- nrow(decoy)
  
  if (nd <= 5L) 
    return(NA_real_)
  
  n <- nt + nd
  lambt <- nt / n
  lambd <- 1 - lambt
  
  vecd <- decoy$pep_score
  sigmad <- sd(vecd, na.rm = TRUE)
  mud <- mean(vecd, na.rm = TRUE)
  
  if (is.na(sigmad)) 
    return(NA_real_)
  
  xs <- seq( mud + 4*sigmad, mud + sigmad, -.1)
  
  for (i in seq_along(xs)) {
    y <- (1 - plnorm(xs[i], mud, sigmad, lower.tail = FALSE)) * nd / n
    
    if (y >= target_fdr) 
      break
  }
  
  xs[i]
}


#' Helper of \link{calc_pepfdr}.
#'
#' Calculates the probability cut-off for target-decoy pairs at a given peptide
#' length.
#' 
#' @param td A target-decoy pair.
#' @param len Numeric; the length of peptides.
#' @inheritParams matchMS
probco_bypeplen <- function (len, td, fdr_type = "protein", target_fdr = 0.01, 
                             min_pepscores_co = 0, out_path) 
{
  td <- dplyr::filter(td, pep_len == len)
  td <- sub_td_byfdrtype(td, fdr_type)
  
  count <- nrow(td)
  
  if (count < (1 / target_fdr)) {
    if (count <= 10L)
      return(NA_real_)
    
    best_co <- tryCatch(
      (find_pepscore_co1(td, target_fdr) + find_pepscore_co2(td, target_fdr))/2,
      error = function(e) NA
    )
    
    prob_co <- 10^(-best_co/10)
    names(prob_co) <- count
    
    return(prob_co)
  }

  # ---
  rows <- which(td$fdr <= target_fdr)
  
  if (length(rows)) {
    max_pr <- max(td[["fdr"]], na.rm = TRUE)
    row <- max(rows, na.rm = TRUE)
    score_co <- td[["pep_score"]][row]

    ### guard against very low score_co
    if (FALSE) {
      lwr <- if (len <= 10L) {
        31 - len * target_fdr * 200
      }
      else if (len <= 15L) {
        29 - 11:15 * target_fdr * 160
      }
      else {
        2
      }

      lwr <- ifelse(lwr <= 5, 5, lwr)
      # if (len <= 15L) score_co <- max(score_co, lwr)
    }
    ####
    
    # fittings (the data range may affect the fitting)
    df <- data.frame(x = td[["pep_score"]], y = td[["fdr"]])
    
    if (max_pr <= target_fdr) {
      fit <- suppressWarnings(
        tryCatch(
          stats::nls(y ~ SSasymp(x, Asym, R0, lrc), data = df, 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA, 
          warning = function (w) NA))

      if (all(is.na(fit))) {
        fit <- suppressWarnings(
          tryCatch(
            stats::nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
                control = list(tol = 1e-03, warnOnly = TRUE), 
                algorithm = "port"), 
            error = function (e) NA, 
            warning = function (w) NA))
      }
    }
    else {
      fit <- suppressWarnings(
        tryCatch(
          stats::nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA, 
          warning = function (w) NA))

      if (all(is.na(fit))) {
        fit <- suppressWarnings(
          tryCatch(
            stats::nls(y ~ SSasymp(x, Asym, R0, lrc), data = df, 
                control = list(tol = 1e-03, warnOnly = TRUE), 
                algorithm = "port"), 
            error = function (e) NA, 
            warning = function (w) NA))
      }
    }

    if (!all(is.na(fit))) {
      min_x <- min(target_fdr, min(df$x, na.rm = TRUE))
      max_x <- max(target_fdr, max(df$x, na.rm = TRUE))
      del_x <- max_x - min_x
      
      step <- if (del_x > 10) 
        .25
      else if (del_x > 1)
        .02
      else
        .01

      newx <- seq(min_x, max_x, by = step)
      newy <- predict(fit, data.frame(x = newx)) %>% `names<-`(newx)
      
      # NA if not existed
      score_co2 <- as.numeric(names(which(newy <= target_fdr)[1]))
      score_co2 <- max(score_co2, min_pepscores_co)

      ### guard against vert low score_co2
      # if (len <= 15L) score_co2 <- max(score_co2, lwr)
      ###
      
      best_co <- min(score_co, score_co2, na.rm = TRUE)

      if (FALSE) {
        try(
          local({
            pdf(file.path(out_path, "temp", paste0("probco_peplen@", len, ".pdf"))) 
            plot(y ~ x, df, xlab = "pep_score", col = "blue", ylab = "FDR", pch = 19)
            title(main = paste0("Peptide length ", len))
            lines(newx, newy, col = "red", type = "b")
            abline(h = target_fdr, col = "green", lwd = 3, lty = 2)
            legend("topright", legend = c("Raw", "Smoothed"), 
                   col = c("blue", "red"), pch = 19, bty = "n")
            legend("center", pch = 19, 
                   legend = c(paste0("Best: ", round(best_co, 2)), 
                              paste0("Smoothed: ", round(score_co2, 2))))
            dev.off()
          })
        )
      }

      rm(list = c("newx", "newy", "score_co", "score_co2"))
    } 
    else {
      best_co <- score_co
    }
    
    rm(list = c("df", "fit"))
    
    prob_co <- 10^(-best_co/10)
  }
  else {
    best_co <- tryCatch(
      (find_pepscore_co1(td, target_fdr) + find_pepscore_co2(td, target_fdr))/2,
      error = function(e) NA)

    prob_co <- 10^(-best_co/10)
  }
  
  names(prob_co) <- count

  invisible(prob_co)
}


#' Subsets targets and decoys by fdr_type.
#' 
#' @param td A data frame of targets and decoys.
#' @inheritParams matchMS
sub_td_byfdrtype <- function (td, fdr_type)
{
  if (fdr_type %in% c("peptide", "protein")) {
    if (fdr_type == "protein") {
      td <- dplyr::arrange(td, pep_seq, pep_prob)
      td <- dplyr::group_by(td, pep_seq)
      td <- dplyr::filter(td, row_number() == 1L)
      td <- dplyr::ungroup(td)
    }
    else {
      # td[["pep_ivmod2"]] <- gsub(" [\\(\\[]\\d+[\\)\\[]$", "", td[["pep_ivmod"]])
      td[["pep_ivmod2"]] <- gsub(" .*", "", td[["pep_ivmod"]])
      td[["pep_seq_mod"]] <- ifelse(is.na(td[["pep_ivmod2"]]), td[["pep_seq"]], 
                                    paste0(td[["pep_seq"]], ".", td[["pep_ivmod2"]]))
      td <- dplyr::arrange(td, pep_seq_mod, pep_prob)
      td <- dplyr::group_by(td, pep_seq_mod)
      td <- dplyr::filter(td, row_number() == 1L)
      td <- dplyr::ungroup(td)
      td[["pep_ivmod2"]] <- NULL
      td[["pep_seq_mod"]] <- NULL
    }
  }
  
  td <- dplyr::select(td, pep_prob, pep_isdecoy)
  td <- dplyr::arrange(td, pep_prob)
  td <- dplyr::mutate(td, total = row_number())
  td <- dplyr::mutate(td, decoy = cumsum(pep_isdecoy))
  td <- dplyr::mutate(td, fdr = decoy/total)
  td <- dplyr::mutate(td, pep_score = -log10(pep_prob) * 10)
}


#' Find the suitable pep_len values for the fitting of probability cut-offs.
#'
#' Recursively decrease the value of \code{min_count} by half until some indexes
#' are found.
#'
#' @param all_lens A vector of all possible pep_len values found from data.
#' @param counts A vector of the counts of sequences at given pep_len values.
#' @param min_count The minimum counts for consideration in fitting.
find_optlens <- function (all_lens, counts, min_count = 128L) 
{
  idxes <- which(counts >= min_count)
  
  if (length(idxes)) 
    return(all_lens[idxes])
  else 
    find_optlens(all_lens, counts, min_count/2L)
}


#' Find the pep_len that yields the lowest probability.
#' 
#' Closest to 13L; favors smaller "left" if with a "right" tie.
#' 
#' @param prob_cos A vector of probability cut-offs.
#' @param guess An integer of guessed valley.
find_probco_valley <- function (prob_cos, guess = 12L) 
{
  len <- length(prob_cos)
  
  if (len == 1L) return(as.integer(names(prob_cos)[1]))

  min_len <- as.integer(names(prob_cos)[1])
  if (min_len > guess) return(min_len)
  
  deltas <- prob_cos[2:len] - prob_cos[1:(len-1)]
  idxes <- which(deltas > 0) 
  ups <- as.numeric(names(idxes))
  
  ups_left <- ups[ups <= guess & ups - guess >= -3L]
  lens_left <- length(ups_left)
  if (lens_left) return(ups_left[lens_left])
  
  ups_right <- ups[ups > guess & ups - guess <= 3L]
  lens_right <- length(ups_right)
  if (lens_right) return(ups_right[1])

  invisible(guess)
}


#' Prepares target-decoy data.
#' 
#' @param td A data frame of targets and decoys (for Percolator).
#' @inheritParams matchMS
prep_pepfdr_td <- function (td = NULL, out_path, enzyme = "trypsin_p", 
                            nes_fdr_group = "base", fdr_group = "base")
{
  files <- list.files(path = file.path(out_path, "temp"), 
                      pattern = "^pepscores_", full.names = TRUE)
  
  if (!length(files)) 
    stop("Score results not found.", call. = FALSE)
  
  top3s <- gsub(paste0("^.*pepscores_", "(\\d+)\\.rds$"), "\\1", 
                files[which_topx2(file.size(files), 3L)[1:3]])
  top3s <- top3s[!is.na(top3s)]
  
  max_i <- gsub(paste0("^.*pepscores_", "(\\d+)\\.rds$"), "\\1", 
                files[which.max(file.size(files))[[1]]])
  
  if (is.null(td)) {
    td <- lapply(files, qs::qread)
    td <- td[lapply(td, nrow) > 0L] # otherwise, error with bind_rows
    td <- dplyr::bind_rows(td)
  }

  enzyme <- tolower(enzyme)
  is_nes <- enzyme == "noenzyme" || grepl("^semi", enzyme)
  
  if (!is_nes)
    if (!nes_fdr_group %in% c("all", "base")) 
      nes_fdr_group <- fdr_group
  
  if (is_nes) {
    td  <- if (nes_fdr_group == "all")
      td
    else if (nes_fdr_group == "all_cterm_tryptic")
      td[grepl("[KR]$", td[["pep_seq"]]), ]
    else if (nes_fdr_group == "all_cterm_nontryptic")
      td[!grepl("[KR]$", td[["pep_seq"]]), ]
    else if (nes_fdr_group == "base")
      td[td[["pep_mod_group"]] == max_i, ]
    else if (nes_fdr_group == "base_cterm_tryptic")
      td[td[["pep_mod_group"]] == max_i & grepl("[KR]$", td[["pep_seq"]]), ]
    else if (nes_fdr_group == "base_cterm_nontryptic")
      td[td[["pep_mod_group"]] == max_i & !grepl("[KR]$", td[["pep_seq"]]), ]
    else if (nes_fdr_group == "top3")
      td[td[["pep_mod_group"]] %in% top3s, ]
    else if (nes_fdr_group == "top3_cterm_tryptic")
      td[td[["pep_mod_group"]] %in% top3s & grepl("[KR]$", td[["pep_seq"]]), ]
    else if (nes_fdr_group == "top3_cterm_nontryptic")
      td[td[["pep_mod_group"]] %in% top3s & !grepl("[KR]$", td[["pep_seq"]]), ]
    else
      stop("Invalid argument for \"nes_fdr_group\".")
  }
  else {
    td <- if (fdr_group == "all")
      td
    else if (fdr_group == "base")
      td[td[["pep_mod_group"]] == max_i, ]
    else if (fdr_group == "top3")
      td[td[["pep_mod_group"]] %in% top3s, ]
    else 
      stop("Invalid argument for \"fdr_group\".")
  }
  
  if (!nrow(td))
    stop("Found nothing: empty targets and decoys.")
  
  td
}


#' Keeps the best entries of targets and decoys.
#' 
#' @param td A data frame of targets and decoys.
#' @param cols Columns for grouping.
keep_pepfdr_best <- function (td, cols = c("pep_scan_num", "raw_file"))
{
  # keeps separated best hits for targets and decoys
  # two lists of "TRUE" and "FALSE"
  td <- split(td, td[["pep_isdecoy"]])
  
  td <- lapply(td, function (x) {
    x <- dplyr::group_by_at(x, cols)
    x <- dplyr::arrange(x, pep_prob)
    x <- dplyr::filter(x, row_number() == 1L)
    x <- dplyr::ungroup(x)
  })
  
  # decoy sequences may be present in targets
  ok <- is.na(fastmatch::fmatch(reverse_seqs(td[["FALSE"]][["pep_seq"]]), 
                                td[["TRUE"]][["pep_seq"]]))
  td[["FALSE"]] <- td[["FALSE"]][ok, ]
  
  #  keeps the best hit for each `scan_num`
  td <- if (nrow(td[["FALSE"]]))
    dplyr::bind_rows(td[c("TRUE", "FALSE")])
  else
    td[["TRUE"]]
  
  td <- dplyr::group_by_at(td, cols)
  td <- dplyr::arrange(td, pep_prob)
  td <- dplyr::filter(td, row_number() == 1L)
  td <- dplyr::ungroup(td)
}


#' Calculates the cut-off score at a peptide FDR.
#'
#' Needs \code{min_len} and \code{max_len} since the target-decoy pair may not
#' cover all \code{pep_len} values.
#' 
#' @param target_fdr Numeric; the levels of false-discovery rate (FDR).
#' @param fdr_type Character string; the type of FDR for controlling.
#' @param fct_score A trivial factor converting p-values to scores.
#' @inheritParams matchMS
#' @examples 
#' \donttest{
#' library(mzion)
#' 
#' if (FALSE) {
#'   prob_cos <- calc_pepfdr(target_fdr = .01, 
#'                           fdr_type = "protein", 
#'                           min_len = 7L, 
#'                           max_len = 50L, 
#'                           out_path = "~/mzion/bi_1")
#' }
#' 
#' }
calc_pepfdr <- function (target_fdr = .01, fdr_type = "protein", 
                         min_len = 7L, max_len = 40L, 
                         max_pepscores_co = 50, min_pepscores_co = 0, 
                         enzyme = "trypsin_p", 
                         fdr_group = "base", 
                         nes_fdr_group = "base", 
                         fct_score = 10, 
                         out_path) 
{
  message("Calculating peptide FDR.")
  
  td <- prep_pepfdr_td(out_path = out_path, 
                       enzyme = enzyme, 
                       nes_fdr_group = nes_fdr_group, 
                       fdr_group = fdr_group)
  
  # back-compatibility to new column keys (e.g. scan_num -> pep_scan_num)
  if (!"pep_scan_num" %in% names(td)) 
    stop("Seems like reprocessing of results from an earlier version.\n", 
         "  Please delete old \"temp\\pep_score[...].rds\".")
  
  if (!sum(td[["pep_isdecoy"]])) {
    warning("No decoys found.")
    seqs <- min_len:max(td[["pep_len"]], na.rm = TRUE)
    prob_cos <- rep(.5, length(seqs))
    
    return(data.frame(pep_len = seqs, pep_prob_co = prob_cos))
  }
  
  td <- keep_pepfdr_best(td)
  qs::qsave(td, file.path(out_path, "temp", "td_pepfdr.rds"), preset = "fast")

  # --- 
  all_lens <- sort(unique(td$pep_len))
  
  prob_cos <- lapply(all_lens, probco_bypeplen, 
                     td = td, 
                     fdr_type = fdr_type, 
                     target_fdr = target_fdr, 
                     min_pepscores_co = min_pepscores_co, 
                     out_path = out_path)
  prob_cos <- unlist(prob_cos)
  
  if (length(prob_cos) == 1L && !is.na(prob_cos))
    return(data.frame(pep_len = all_lens, pep_prob_co = prob_cos))
  else if (length(prob_cos) == 1L && is.na(prob_cos))
    return(data.frame(pep_len = all_lens, pep_prob_co = target_fdr))
  else if (all(is.na(prob_cos))) {
    seqs <- min_len:max(td$pep_len, na.rm = TRUE)
    prob_cos <- rep(target_fdr, length(seqs))
    
    return(data.frame(pep_len = seqs, pep_prob_co = prob_cos))
  } 
  
  counts <- as.integer(names(prob_cos))
  names(counts) <- all_lens
  names(prob_cos) <- all_lens
  
  prob_nas <- prob_cos[is.na(prob_cos)]
  prob_cos <- prob_cos[!is.na(prob_cos)]
  
  if (length(prob_cos) == 1L) {
    seqs <- min_len:max_len
    prob_cos <- rep(prob_cos, length(seqs))
    
    return(data.frame(pep_len = seqs, pep_prob_co = prob_cos))
  }
  
  #########################################
  # (At least two non-trivial prob_cos)
  #########################################
  
  lens <- find_optlens(all_lens, counts, 50L) # changed from 128L
  
  # no fittings
  if (length(lens) <= 3L) {
    prob_cos <- prob_cos[!is.na(names(prob_cos))]
    prob_no_uses <- NULL
    
    return(fill_probco_nas(prob_nas = prob_nas, prob_no_uses = prob_no_uses, 
                           prob_cos = prob_cos, target_fdr = target_fdr))
  }
  
  # quality ones for fittings
  ok_lens <- names(prob_cos) %in% lens
  prob_no_uses <- prob_cos[!ok_lens]
  count_no_uses <- counts[!ok_lens]
  prob_cos <- prob_cos[ok_lens]
  counts <- counts[ok_lens]
  
  valley <- find_probco_valley(prob_cos)
  best_co <- -fct_score * log10(prob_cos[as.character(valley)])
  
  prob_cos <- local({
    start <- which(names(prob_cos) == valley)
    end <- length(prob_cos)
    x <- -log10(prob_cos[start:end])
    
    # first pass
    fit_lm <- stats::lm(x ~ as.integer(names(x)))
    res_lm <- stats::residuals(fit_lm)
    bad_lm <- x[abs(res_lm) > 5]
    
    if (length(bad_lm)) {
      rm(list = c("fit_lm", "res_lm"))
      
      x_ok <- x[! names(x) %in% names(bad_lm)]
      fit_lm <- stats::lm(x_ok ~ as.integer(names(x_ok)))
      coefs <- stats::coef(fit_lm)
      slp <- coefs[[2]]
      
      if (is.na(slp)) slp <- 0
      
      x[names(bad_lm)] <- slp * as.integer(names(bad_lm)) + coefs[[1]]
      rm(list = c("x_ok", "fit_lm", "coefs", "bad_lm"))
    }
    
    # second pass
    ans <- tsoutliers(x)
    x[ans$index] <- ans$replacements
    prob_cos[start:end] <- 1/10^x
    
    prob_cos
  })
  
  ## fittings
  df <- data.frame(x = as.numeric(names(prob_cos)), 
                   y = -fct_score * log10(prob_cos))
  
  # valley left
  df_left <- df[df$x <= valley, ]
  nrow_left <- nrow(df_left)
  rank_left <- 4L # often results in: originals == fitted <-> no fitting
  # rank_left <- min(nrow_left, 4L)
  
  if (nrow_left == 0L) {
    newx_left <- NULL
    newy_left <- NULL
  }
  else if (nrow_left == 1L) {
    newx_left <- df_left$x[1]
    newy_left <- df_left$y[1]
  }
  else {
    fit_left <- if (nrow_left <= rank_left) {
      lm(y ~ x, df_left)
    } 
    else {
      local({
        fit_ns <- tryCatch(
          stats::lm(y ~ splines::ns(x, rank_left), df_left),
          error = function(e) NA
        )
        
        fit_bs <- tryCatch(
          stats::lm(y ~ splines::bs(x, rank_left), df_left),
          error = function(e) NA
        )
        
        res1 <- if (class(fit_ns) == "lm") sum(resid(fit_ns)^2) else Inf
        res2 <- if (class(fit_bs) == "lm") sum(resid(fit_bs)^2) else Inf
        
        if (res1 <= res2) fit_ns else fit_bs
      })
    }
    
    newx_left <- min(df_left$x, na.rm = TRUE):max(df_left$x, na.rm = TRUE)
    newy_left <- predict(fit_left, data.frame(x = newx_left))
    names(newy_left) <- newx_left
  }
  
  # valley right (small rank to down-weight wiggly high `pep_len` points)
  df_right <- df[df$x > valley, ]
  nrow_right <- nrow(df_right)
  rank_right <- 4L # changed from 3 to 4
  
  if (nrow_right == 0L) {
    newx_right <- NULL
    newy_right <- NULL
  } 
  else if (nrow_right == 1L) {
    newx_right <- df_right$x
    newy_right <- df_right$y
  } 
  else {
    if (nrow_right >= rank_right) {
      fit_right <- local({
        fit_ns <- tryCatch(
          stats::lm(y ~ splines::ns(x, rank_right), df_right),
          error = function(e) NA
        )
        
        fit_bs <- tryCatch(
          stats::lm(y ~ splines::bs(x, rank_right), df_right),
          error = function(e) NA
        )
        
        res1 <- if (class(fit_ns) == "lm") sum(resid(fit_ns)^2) else Inf
        res2 <- if (class(fit_bs) == "lm") sum(resid(fit_bs)^2) else Inf
        
        if (res1 <= res2) fit_ns else fit_bs
      })
      
      newx_right <- min(df_right$x, na.rm = TRUE):max(df_right$x, na.rm = TRUE)
      newy_right <- predict(fit_right, data.frame(x = newx_right))
      names(newy_right) <- newx_right
    } 
    else {
      fit_right <- tryCatch(
        stats::lm(y ~ x, df_right),
        error = function(e) NA
      )
      
      slope <- if (class(fit_right) == "lm") 
        unname(stats::coef(fit_right)[2])
      else 
        .1
      
      newx_right <- valley:max(df_right$x, na.rm = TRUE)
      newy_right <- best_co + slope * (newx_right - valley)
      names(newy_right) <- newx_right
      
      # excludes the `valley` itself (already in left fitting)
      newx_right <- newx_right[-1]
      newy_right <- newy_right[-1]
    }
  }
  
  # left + right
  newx <- c(newx_left, newx_right)
  newy <- c(newy_left, newy_right)
  newy[newy < min_pepscores_co] <- min_pepscores_co
  newy[newy > max_pepscores_co] <- max_pepscores_co
  
  local({
    n_row <- nrow(df)
    
    df_new <- rbind2(
      cbind(df, type = rep("Original", n_row)), 
      data.frame(x = newx, y = newy, type = "Fitted")
    )
    
    n_row2 <- nrow(df_new) - n_row
    
    try(
      local({
        pdf(file.path(out_path, "temp", "pepscore_len.pdf")) 
        plot(y ~ x, df_new, col = c(rep("blue", n_row), rep("red", n_row2)), 
             xlab = "pep_len", ylab = "score_co", pch = 19)
        legend("topright", legend = c("Raw", "Smoothed"), 
               col = c("blue", "red"), pch = 19, bty = "n")
        dev.off()
      })
    )
  })
  
  ## Outputs
  prob_cos <- local({
    nms <- names(prob_cos)
    prob_cos <- 10^(-newy/fct_score)
    names(prob_cos) <- nms
    
    # with `find_optlens`, some length in newy may not be in prob_cos
    prob_cos[!is.na(names(prob_cos))]
  })
  
  prob_cos <- fill_probco_nas(prob_nas, prob_no_uses, prob_cos, target_fdr) %T>% 
    qs::qsave(file.path(out_path, "temp", "pep_probco.rds"), preset = "fast")
}


#' Fills NA in probability cutoffs
#' 
#' @param prob_nas Named vector; peptide length in names and NA in values.
#' @param prob_no_uses Named vector; peptide length in names and NA in values.
#' @param prob_cos Named vector; non-NA probability cutoffs.
#' @inheritParams calc_pepfdr
fill_probco_nas <- function (prob_nas = NULL, prob_no_uses = NULL, prob_cos, 
                             target_fdr = .01)
{
  prob_nas <- fill_probs(prob_nas, prob_cos, target_fdr)
  prob_no_uses <- fill_probs(prob_no_uses, prob_cos, target_fdr)
  prob_cos <- c(prob_cos, prob_nas, prob_no_uses)
  
  prob_cos <- data.frame(pep_len = as.numeric(names(prob_cos)), 
                         pep_prob_co = prob_cos)
  
  prob_cos <- dplyr::arrange(prob_cos, pep_len)
}


#' Helper of \link{fill_probco_nas}
#'
#' @param nas A vector of NA or inferior probability cut-offs to be replaced.
#' @inheritParams fill_probco_nas
#' @return The replacements. Do not return the concatenated prob_cos and nas, to
#'   prevent the estimated values from affecting the consequent replacement of
#'   prob_no_uses.
fill_probs <- function (nas, prob_cos, target_fdr = .01) 
{
  len <- length(nas)
  
  if (!len)
    return(NULL)
  
  peplen_prs <- as.integer(names(prob_cos))
  peplen_nas <- as.integer(names(nas))
  nas <- nas[!peplen_nas %in% peplen_prs]
  peplen_nas <- as.integer(names(nas))

  co_nas <- if (all(peplen_nas > max(peplen_prs)))
    prob_cos[length(prob_cos)]
  else if (all(peplen_nas < min(peplen_prs)))
    prob_cos[1]
  else
    min(median(prob_cos), mean(prob_cos), target_fdr)
  
  val <- unname(co_nas)
  vals <-rep(val, len)
  names(vals) <- peplen_nas

  vals
}


#' Post processing after peptide FDR.
#' 
#' @param prob_cos Probability cut-offs (in data frame).
#' @param out_path An output path.
post_pepfdr <- function (prob_cos = NULL, out_path = NULL) 
{
  if (is.null(prob_cos)) {
    file_prob <- file.path(out_path, "temp", "pep_probco.rds")
    
    if (file.exists(file_prob))
      prob_cos <- qs::qread(file_prob)
    else 
      stop("File not found: ", file_prob)
  }
  
  fct_score <- 10
  
  ok_targets <- find_targets(out_path, "^pepscores_")
  files <- ok_targets$files
  
  if (!length(files)) 
    stop("Results of peptide scores not found.")
  
  td <- lapply(files, function (x) qs::qread(file.path(out_path, "temp", x)))
  names(td) <- ok_targets$idxes
  ok <- lapply(td, nrow) > 0L
  td <- dplyr::bind_rows(td[ok])
  rm(list = c("ok_targets", "files", "ok"))
  
  if (!nrow(td)) 
    stop("No PSM matches for scoring. Consider different search parameters.")
  
  # Adjusted p-values (just to moderate pep_score)
  td <- td %>% 
    dplyr::left_join(prob_cos, by = "pep_len") %>% 
    dplyr::mutate(pep_issig = ifelse(pep_prob <= pep_prob_co, TRUE, FALSE), 
                  pep_adjp = p.adjust(pep_prob, "BH"))
  
  adjp_cos <- purrr::map_dbl(prob_cos$pep_prob_co, function (x) {
    row <- which.min(abs(log10(td$pep_prob/x)))
    td[row, ]$pep_adjp
  }) 
  
  prob_cos <- dplyr::bind_cols(prob_cos, pep_adjp_co = adjp_cos)
  
  td <- td %>% 
    dplyr::left_join(prob_cos[, c("pep_len", "pep_adjp_co")], by = "pep_len") %>% 
    dplyr::mutate(pep_score = -log10(pep_adjp) * fct_score, 
                  pep_score = ifelse(pep_score > 250, 250, pep_score), 
                  pep_score_co = -log10(pep_adjp_co) * fct_score) %>% 
    dplyr::select(-c("pep_prob", "pep_adjp", "pep_prob_co", "pep_adjp_co"))
  
  qs::qsave(td, file.path(out_path, "temp", "pepfdr.rds"), preset = "fast")
  
  invisible(td)
}


#' Calculates the cut-offs of protein scores.
#'
#' @param df An output from upstream steps.
#' @param out_path An output path.
#' @inheritParams calc_pepfdr
#' @inheritParams matchMS
calc_protfdr <- function (df = NULL, target_fdr = .01, max_protscores_co = Inf, 
                          max_protnpep_co = 10L, method_prot_es_co = "median", 
                          out_path = NULL) 
{
  message("Calculating peptide-protein FDR.")
  
  # score cut-offs as a function of prot_n_pep
  max_n_pep <- max(df$prot_n_pep, na.rm = TRUE)
  all_n_peps <- unique(df$prot_n_pep)
  
  # protein enrichment score cut-offs at each `prot_n_pep`
  td <- df[with(df, pep_issig), ]
  
  score_co <- split(td, td[["prot_n_pep"]])
  score_co <- lapply(score_co, calc_protfdr_i, 
                     target_fdr = target_fdr, 
                     max_protnpep_co = max_protnpep_co, 
                     method_prot_es_co = method_prot_es_co, 
                     out_path = out_path)
  score_co <- unlist(score_co, recursive = FALSE, use.names = TRUE)
  score_co[score_co > max_protscores_co] <- max_protscores_co
  
  # fitted score cut-offs
  score_co_bf <- data.frame(prot_n_pep = as.integer(names(score_co)), 
                            prot_score_co_bf = score_co)
  
  score_co <- score_co %>% 
    fit_protfdr(max_n_pep, out_path) %>% 
    dplyr::filter(prot_n_pep %in% all_n_peps) %>% 
    dplyr::left_join(score_co_bf, by = "prot_n_pep") %>% 
    dplyr::mutate(prot_es_co = ifelse(prot_score_co <= prot_score_co_bf, 
                                      prot_score_co, prot_score_co_bf)) %>% 
    dplyr::select(-c("prot_score_co", "prot_score_co_bf"))
  
  # add protein enrichment score
  prot_es <- df %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(pep_issig) %>% 
    dplyr::mutate(pep_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_es = max(pep_es, na.rm = TRUE))
  
  # puts together
  df <- df %>% 
    dplyr::left_join(prot_es, by = "prot_acc")
  
  df <- df %>% 
    dplyr::left_join(score_co, by = "prot_n_pep") %>% 
    dplyr::mutate(prot_issig = ifelse(prot_es >= prot_es_co, TRUE, FALSE)) %>% 
    dplyr::mutate(pep_score = round(pep_score, digits = 1L), 
                  pep_score_co = round(pep_score_co, digits = 1L), 
                  prot_es = round(prot_es, digits = 1L), 
                  prot_es_co = round(prot_es_co, digits = 1L))
}


#' Function factories for descriptive statistics.
#'
#' @param f One of \code{max, mean, median, min}. Data frame \code{df} is a PSM
#'   table.
aggr_prot_es <- function(f) 
  function (df, ...) dplyr::summarise(df, prot_es = f(prot_es, ...))


#' Helper of \link{calc_protfdr}.
#' 
#' For prot_n_pep at value i.
#' 
#' @param td A data frame with paired target-decoys at prot_n_pep = i.
#' @param n_burnin The minimum number of burn-ins in protein enrichment scores. 
#' @param out_path An output path.
#' @inheritParams calc_protfdr
#' @return A score cut-off at a given prot_n_pep.
calc_protfdr_i <- function (td, target_fdr = .01, max_protnpep_co = 10L, 
                            method_prot_es_co = "median", n_burnin = 3L, out_path) 
{
  options(digits = 9L)
  
  if (td[["prot_n_pep"]][[1]] > max_protnpep_co)
    return(0L)

  td <- td %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup()

  ## no decoys
  if (sum(td$pep_isdecoy) == 0L) 
    return(0L)
  
  ## all decoys
  if (sum(!td$pep_isdecoy) == 0L) {
    if (nrow(td) <= 5L) return(0L) else return(20L)
  }

  ## both targets and decoys
  if (nrow(td) <= 20L) 
    return(1L)
  
  # (NOT to use `local` for hard `return (0L)`)
  prot_scores <- td %>% 
    dplyr::mutate(prot_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) 
  
  prot_scores <- switch(method_prot_es_co, 
                        max = aggr_prot_es(max)(prot_scores, na.rm = TRUE), 
                        mean = aggr_prot_es(mean)(prot_scores, na.rm = TRUE), 
                        median = aggr_prot_es(median)(prot_scores, na.rm = TRUE),
                        min = aggr_prot_es(min)(prot_scores, na.rm = TRUE), 
                        aggr_prot_es(max)(prot_scores, na.rm = TRUE))

  # no decoys
  if (!any(grepl("^-", prot_scores$prot_acc))) 
    return (0L)
  
  # Multiple dipping of -NP_003310, -NP_597676 with the same set of 
  #  identifying pep_seqs and identical `prot_es`
  # 
  # prot_acc     prot_es
  # <chr>          <dbl>
  # 1 -NP_003310      43.1
  # 2 -NP_597676      43.1
  # 3 -NP_597681      43.1
  # 4 NP_001004067    96.3
  #
  # Not yet to parse out same-set, sub-set proteins. For simplicity and performance, 
  # keep one unique `prot_es` (by roughly assuming `prot_es` is an indicator for 
  # the redundancy). 

  prot_scores <- prot_scores %>% 
    dplyr::mutate(prot_es = round(prot_es, digits = 5L)) %>% 
    dplyr::filter(!duplicated(prot_es))
  
  # Removes the last three in case only one or two decoy spikes near the end 
  # the list 
  # 
  # prot_acc      prot_es
  # <chr>           <dbl>
  # NP_001028448     53.4
  # -NP_001157012    54.4
  # NP_077772        62.7

  prot_scores <- prot_scores %>% 
    dplyr::arrange(prot_es) %>% 
    dplyr::filter(row_number() > n_burnin)

  # no decoys
  if (!any(grepl("^-", prot_scores$prot_acc))) 
    return (0L)

  # both targets and decoys
  td <- td %>% 
    dplyr::right_join(prot_scores, by = "prot_acc") %>% 
    dplyr::arrange(-prot_es) %>% 
    dplyr::mutate(total = row_number()) %>% 
    dplyr::mutate(decoy = cumsum(pep_isdecoy)) %>% 
    dplyr::mutate(fdr = decoy/total)
  
  rm(list = "prot_scores")
  
  rows <- which(td$fdr <= target_fdr)
  row <- if (length(rows)) max(rows, na.rm = TRUE) else -Inf

  if (row == -Inf) {
    score_co <- 0L
    score_co2 <- score_co
  } 
  else {
    score_co <- td$prot_es[row]

    score_co2 <- local({
      data <- data.frame(x = td[["prot_es"]], y = td[["fdr"]])
      
      fit <- suppressWarnings(
        tryCatch(
          stats::nls(y ~ SSlogis(x, Asym, xmid, scal), data = data, 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      )
      
      if (all(is.na(fit))) {
        score_co2 <- score_co
      } 
      else {
        min_score <- min(data$x, na.rm = TRUE)
        max_score <- max(data$x, na.rm = TRUE)
        newx <- min_score:max_score
        newy <- predict(fit, data.frame(x = newx)) %>% `names<-`(newx)

        # NA if not existed
        score_co2 <- which(newy <= target_fdr)[1] %>% 
          names() %>% 
          as.numeric()
        
        local({
          title <- paste0("prot_n_pep@", td$prot_n_pep[1])
          
          if (FALSE) {
            try(
              local({
                pdf(file.path(out_path, "temp", paste0(title, ".pdf"))) 
                plot(y ~ x, data, xlab = "Enrichment score", col = "blue", 
                     ylab = "FDR", pch = 19)
                title(main = title)
                lines(newx, newy, col = "red", type = "b")
                abline(h = target_fdr, col = "green", lwd = 3, lty = 2)
                legend("topright", legend = c("Raw", "Smoothed"), 
                       col = c("blue", "red"), pch = 19, bty = "n")
                legend("center", pch = 19, 
                       legend = c(paste0("Raw: ", round(score_co, 2)), 
                                  paste0("Smoothed: ", round(score_co2, 2))))
                dev.off()
              })
            )
          }
        })
      }

      score_co2
    })
  }
  
  invisible(min(score_co, score_co2, na.rm = TRUE))
}


#' Fits the raw cut-offs of protein scores.
#'
#' Assumed a sigmoidal function.
#'
#' @param vec Named numeric vector. The values are the score cut-offs as a
#'   function of \code{prot_n_pep}. The names correspond to the number of
#'   peptides being identified.
#' @param max_n_pep Integer; the maximum value of \code{prot_n_pep} for
#'   prediction.
#' @param out_path An output path.
fit_protfdr <- function (vec, max_n_pep = 1000L, out_path) 
{
  if (length(vec) <= 10L) 
    return(data.frame(prot_n_pep = as.numeric(names(vec)), prot_score_co = vec))

  rv <- rev(vec)
  df <- data.frame(x = as.numeric(names(rv)), y = rv)
  elbow <- min(df[which(df$y == min(df$y, na.rm = TRUE)), "x"], na.rm = TRUE)
  amp <- max(df$y, na.rm = TRUE) * .8
  sca <- 0.5
  
  ## nls
  f <- function (x, m = 0, s = 1, a = 1) { a - a / (1 + exp(-(x-m)/s)) }
  
  fit <- suppressWarnings(
    tryCatch(
      stats::nls(y ~ f(x, m, s, a), data = df, 
          start = list(a = amp, m = elbow, s = sca), 
          control = list(tol = 1e-03, warnOnly = TRUE), 
          algorithm = "port"), 
      error = function (e) NA)
  )
  
  # should not occur
  if (all(is.na(fit))) {
    fits <- suppressWarnings(
      purrr::map(seq_len(elbow-1), ~ {
        tryCatch(
          stats::nls(y ~ f(x, m, s, a), data = df, 
              start = list(a = amp, m = .x, s = sca), 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      })
    )
    
    fit <- if (all(is.na(fits))) {
      NA
    } 
    else {
      fits %>% 
        .[!is.na(.)] %>% 
        .[[length(.)]]
    }
  }

  ## lm
  df_left <- df[df$x <= elbow, ]
  df_right <- df[df$x > elbow, ]
  
  fit_left <- tryCatch(
    lm(y ~ x, data = df_left), 
    error = function (e) NA)
  
  fit_right <- tryCatch(
    lm(y ~ x, data = df_right), 
    error = function (e) NA)
  
  res_left <- if (class(fit_left) == "lm") sum(resid(fit_left)^2) else Inf
  res_right <- if (class(fit_right) == "lm") sum(resid(fit_right)^2) else Inf
  res_lm <- sum(res_left + res_right)
  
  res_nls <- if (class(fit) == "nls") sum(resid(fit)^2) else Inf

  ## nls vs. lm
  if (res_nls <= res_lm) {
    newx <- 1:max_n_pep
    newy <- predict(fit, data.frame(x = newx))
    
    out <- data.frame(
      prot_n_pep = newx, 
      prot_score_co = newy) %>% 
      dplyr::mutate(prot_score_co = ifelse(prot_n_pep >= elbow, 0, prot_score_co))
  } 
  else {
    newx_left <- 1:elbow
    newy_left <- predict(fit_left, data.frame(x = newx_left))
    newx_right <- (elbow + 1L):max_n_pep
    newy_right <- predict(fit_right, data.frame(x = newx_right))
    
    newx <- c(newx_left, newx_right)
    newy <- c(newy_left, newy_right)

    out <- data.frame(
      prot_n_pep = newx, 
      prot_score_co = newy) %>% 
      dplyr::mutate(prot_score_co = ifelse(prot_n_pep >= elbow, 0, prot_score_co))
    
  }
  
  try(
    local({
      pdf(file.path(out_path, "temp", "protein_score_co.pdf")) 
      plot(y ~ x, df, xlab = "prot_n_pep", col = "blue", 
           ylab = "Protein score", pch = 19)
      title(main = "Protein score CO")
      # lines(newx, newy, col = "red", type = "b")
      # lines(newx, newy, col = "red")
      lines(out$prot_n_pep , out$prot_score_co, col = "red", type = "b")
      lines(out$prot_n_pep , out$prot_score_co, col = "red")
      legend("topright", legend = c("Raw", "Smoothed"), 
             col = c("blue", "red"), pch = 19, bty = "n")
      dev.off()
      
      qs::qsave(df, file.path(out_path, "temp", "protein_score_co.rds"))
    })
  )
  
  invisible(out)
}


#' Finds the the outer products for a vector of MS2 ions at a given ion series.
#'
#' \emph{Experiment} values go first and theoretical seconded.
#'
#' @param X Numeric vector; one series experimental MS2s.
#' @param Y Numeric vector; one series of theoretical MS2s.
#' @inheritParams matchMS
#' 
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' expts <- c(101.0714, 102.0554, 110.0717, 115.0505, 126.1279, 127.0504,
#'            127.1249, 127.1312, 128.1283, 128.1346, 129.0660, 129.1316,
#'            129.1379, 130.1350, 130.1412, 131.1383, 133.0431, 133.0609,
#'            134.0448, 136.0757, 145.0608, 158.0924, 173.1496, 175.1190,
#'            176.1594, 191.0663, 193.0971, 198.0873, 201.0869, 201.1233,
#'            202.0821, 230.1702, 248.0876, 248.1806, 257.1527, 298.6711,
#'            312.6687, 335.1190, 361.1768, 361.6685, 367.2292, 369.6903,
#'            370.1821, 370.6830, 376.2756, 377.2790, 384.6952, 412.7264,
#'            420.2143, 423.2597, 475.2341, 475.3442, 476.3469, 484.2382,
#'            484.7389, 487.7817, 488.2842, 496.2717, 523.3182, 537.3159,
#'            537.8178, 572.8343, 573.3373, 623.3586, 624.3298, 636.8528,
#'            637.3558, 637.8575, 638.4089, 645.3453, 645.8500, 646.3501,
#'            687.8795, 688.3799, 692.8696, 695.4283, 701.3759, 701.8768,
#'            702.3782, 702.8797, 737.3956, 737.8962, 738.3671, 739.3572,
#'            740.3575, 745.9085, 746.4121, 752.4502, 801.4272, 809.4722,
#'            810.4739, 852.4418, 967.4690, 968.4714, 969.5046, 1084.5293,
#'            1132.5676, 1133.5686, 1134.5735, 1135.5786)
#'
#' theos <- c(116.0342, 203.0662, 290.0983, 418.1569, 475.1783, 572.2311,
#'            732.2617, 861.3043, 958.3571, 1071.4412, 1168.4939, 1225.5154,
#'            1322.5681, 1435.6522, 1536.6999, 1664.7585, 1761.8112, 1917.9123,
#'            175.1190, 272.1717, 400.2303, 501.2780, 614.3620, 711.4148,
#'            768.4363, 865.4890, 978.5731, 1075.6259, 1204.6684, 1364.6991,
#'            1461.7519, 1518.7733, 1646.8319, 1733.8639, 1820.8960, 1935.9229)
#'
#' names(theos) <- c("D", "S", "S", "Q", "G", "P", "C", "E", "P",
#'                   "L", "P", "G", "P", "L", "T", "Q", "P", "R",
#'                   "R", "P", "Q", "T", "L", "P", "G", "P", "L",
#'                   "P", "E", "C", "P", "G", "Q", "S", "S", "D")
#'
#' mzion:::find_ppm_outer_bycombi(theos, expts, ppm_ms2 = 25L)
#' 
#' # No secondary matches
#' theos <- c(56.5233, 235.1522, 284.6864, 326.2050, 361.7235, 
#'            418.2656, 474.8076, 653.4366, 95.0128, 452.2707, 
#'            551.3391, 634.3762, 705.4133, 818.4974, 931.5815, 
#'            1288.8394, 48.0100, 226.6390, 276.1732, 317.6917, 
#'            353.2103, 409.7523, 466.2944, 644.9233, 94.0287, 
#'            451.2866, 550.3550, 633.3921, 704.4292, 817.5133, 
#'            930.5974, 1287.8553, 47.5180, 226.1470, 275.6812, 
#'            317.1997, 352.7183, 409.2603, 465.8024, 644.4313,
#'            188.6415, 245.1835, 301.7256, 337.2441, 378.7627, 
#'            428.2969, 606.9258, 670.9551, 359.2492, 472.3333, 
#'            585.4174, 656.4545, 739.4916, 838.5600, 1195.8179, 
#'            1323.8765, 180.1282, 236.6703, 293.2123, 328.7309, 
#'            370.2494, 419.7836, 598.4126, 662.4419, 358.2651, 
#'            471.3492, 584.4333, 655.4704, 738.5075, 837.5759, 
#'            1194.8338, 1322.8924, 179.6362, 236.1783, 292.7203, 
#'            328.2389, 369.7574, 419.2916, 597.9206, 661.9499)
#' 
#' a <- c("Q", "K", "V", "M", "A", "I", "I", "K")
#' b <- rev(a)
#' names(theos) <- c(rep(a, 5), rep(b, 5))
#' 
#' expts <- c(101.01134,110.07162,112.03947,116.01675,118.83968,
#'            120.08091,126.12779,127.12479,127.13109,128.12816,
#'            128.13440,129.13150,129.13777,130.13483,130.14053,
#'            131.13824,150.26431,155.08147,156.07658,159.07658,
#'            173.14934,175.05345,175.15630,176.13777,176.15942,
#'            176.72833,187.07147,188.15984,203.04829,203.35188,
#'            204.05164,219.14928,227.06596,229.16634,230.16991,
#'            231.17354,232.95274,248.17976,268.16522,315.25906,
#'            318.07504,319.07828,333.21524,361.11679,361.20996,
#'            372.28027,376.27521,377.27872,389.11212,390.11490,
#'            404.25201,432.24750,489.26880,489.35962,504.09521,
#'            510.33450,515.28448,527.35297,533.29462,542.31720,
#'            560.34607,577.31445,588.33716,594.34283,599.39276,
#'            602.44226,604.36749,614.35211,632.36328,659.46362,
#'            694.44006,694.93970,695.39063,707.20782,734.39813,
#'            752.41211,753.41321,756.51660,757.52002,830.91998,
#'            851.48267,852.48303,855.58301,899.52069,956.62952,
#'            1027.66833,1230.73962,1232.75513,1233.75024)
#' 
#' mzion:::find_ppm_outer_bycombi(theos, expts, ppm_ms2 = 25)
#' }
find_ppm_outer_bycombi <- function (X, Y, ppm_ms2 = 20L) 
{
  d <- outer(X, Y, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)
  ix <- row_cols[, 1]
  iy <- row_cols[, 2]

  es <- rep(NA_real_, length(Y))
  names(es) <- names(Y)
  es[iy] <- X[ix]
  
  list(theo = Y, expt = es, ith = iy, iex = ix, m = length(ix))
}


#' Matches between secondary experimentals and theoreticals.
#' 
#' At \code{index_mgf_ms2 = FALSE}.
#' 
#' @param expt A vector of experimental m-over-z values.
#' @param theo A vector of theoretical m-over-z values.
#' @param d Bin size, e.g., \eqn{20 ppm / 2 * 1E-6}.
#' @inheritParams matchMS
match_ex2th2 <- function (expt, theo, min_ms2mass = 115L, d = 1E-5, 
                          index_mgf_ms2 = FALSE) 
{
  th <- index_mz(theo, from = min_ms2mass, d = d)
  ex <- if (index_mgf_ms2) expt else index_mz(expt, from = min_ms2mass, d = d)
  ith <- .Internal(which(th %fin% ex | (th - 1L) %fin% ex | (th + 1L) %fin% ex))

  # if: e.g. th[ith+1] = th[ith] + 1 -> can have NA in iex:
  #   th[ith+1] not in ex but th[ith+1] - 1
  # OK to keep the NA:
  #   es initiated as all NA, OK assign NA <- NA during es[ith] <- expt[iex]
  # in intensity tally of experimental intensity: 
  #   `%+%` default with na.rm = TRUE

  thi <- th[ith]
  iex <- fastmatch::fmatch(thi, ex)
  
  # indexes before and after
  nas <- .Internal(which(is.na(iex)))
  
  if (length(nas)) {
    bf <- fastmatch::fmatch(thi - 1L, ex)
    
    if (all(is.na(bf))) {
      af <- fastmatch::fmatch(thi + 1L, ex)
      iex[nas] <- af[nas]
    }
    else {
      iex[nas] <- bf[nas]
      nas <- .Internal(which(is.na(iex)))

      if (length(nas)) {
        af <- fastmatch::fmatch(thi + 1L, ex)
        iex[nas] <- af[nas]
      }
    }
  }
  
  es <- rep(NA_real_, length(th))
  names(es) <- names(th)
  es[ith] <- expt[iex]
  
  list(theo = theo, expt = es, ith = ith, iex = iex, m = length(iex))
}


#' Calculates the delta scores of \code{pep_seq}.
#'
#' A score delta between the best and the second best.
#'
#' There should not be any duplicated rows at the combination of
#' c("pep_isdecoy", "scan_num", "raw_file", "pep_seq", "pep_ivmod")
#'
#' @param x The results from \link{calc_pepscores}.
#' @param out_path An output path.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @inheritParams matchMS
#' @rawNamespace import(data.table, except = c(last, first, between, transpose,
#'   melt, dcast))
calc_peploc <- function (x = NULL, out_path = NULL, mod_indexes = NULL, 
                         locmods = c("Phospho (S)", "Phospho (T)", "Phospho (Y)"), 
                         topn_mods_per_seq = 3L, topn_seqs_per_query = 3L) 
{
  message("Calculating peptide localization scores and deltas.")

  # some shallow copy warnings from data.table
  old_opts <- options()
  options(warn = -1L)
  on.exit(options(old_opts), add = TRUE)
  
  if (is.null(x)) {
    file <- file.path(out_path, "temp", "pepfdr.rds")
    
    if (file.exists(file)) 
      x <- qs::qread(file) 
    else 
      stop("File not found: ", file)
    
    rm(list = "file")
  }
  
  x <- data.table::data.table(x)
  gc()
  
  n_cores <- detect_cores(16L)
  para <- nrow(x) > 10000L
  
  # For simplicity `pep_seq` uses interchangeably with `uniq_id` and 
  # `pep_seq_mod` with `uniq_id2` where everything is on top of the same 
  # pep_isdecoy, pep_scan_num, raw_file.
  # 
  # uniq_id --- can differentiate the same `pep_seq` at different mod locations & NL 
  # uniq_id2 --- can differentiate the same `pep_seq_mod` at different NLs
  # uniq_id3 --- can differentiate the same *query* at different `pep_seq`s
  # 
  # pep_rank --- `pep_seq_mod`s under the same `pep_seq`
  #   (note that each `pep_seq_mod` is represented by its best NL)
  # ppe_rank2 --- NLs under the same `pep_seq_mod`
  # pep_rank at output --- `pep_seq`s under the same `query`
  
  
  ## 1. compile `uniq_id`, `uniq_id2` and `pep_rank2`
  message("\tRank peptides by neutral losses.")
  
  x[, pep_isdecoy := as.integer(pep_isdecoy)]
  x[, uniq_id := paste(pep_isdecoy, pep_scan_num, raw_file, pep_seq, sep = ".")]
  x[, "pep_ivmod2" := gsub(" [\\(\\[]\\d+[\\)\\[]$", "", pep_ivmod)]
  x[, uniq_id2 := paste(uniq_id, pep_ivmod2, sep = ".")]
  
  if (para) {
    x <- x[order(x[["uniq_id2"]]), ] # for group split

    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    xs <- parallel::clusterApply(
      cl, split(x, find_chunkbreaks(x[["uniq_id2"]], n_cores)), calcpeprank_1)
    parallel::stopCluster(cl)
    
    x <- data.table::rbindlist(xs, use.names = FALSE)
    rm(list = c("xs"))
  }
  else {
    x[, pep_rank2 := data.table::frank(-pep_score, ties.method = "min"), 
      by = list(uniq_id2)]
  }
  

  ## 2 separate the best NL (x0) from the rest (y0, z0) at the same `pep_seq_mod`
  # 
  # the same pep_seq_mod only differ by NLs:
  # 
  # uniq_id2                                            pep_rank      pep_ivmod
  # 0.14332.1.ENGGTEDMFVMYLGNKDASK.00000007007000500000    1   00000007007000500000 (1)
  # 0.14332.1.ENGGTEDMFVMYLGNKDASK.00000007007000500000    1   00000007007000500000 (2)
  # 0.14332.1.ENGGTEDMFVMYLGNKDASK.00000007007000500000    1   00000007007000500000 (3)
  # 0.14332.1.ENGGTEDMFVMYLGNKDASK.00000007007000500000    1   00000007007000500000 (4)
  # one pep_seq_mod at different NLs can results in multiple rows in x0 
  # 
  # x0 --- best NL at a pep_seq_mod and only the first one if with ties in NLs
  # y0 --- not the best NL at a pep_seq_mod
  # z0 --- still the best NL at a pep_seq_mod but tied with the first one
  
  # need to order after ranking (better NLs first)
  x <- x[order(-pep_score), by = list(uniq_id2)] 
  x[, nl_id := seq_len(.N), by = uniq_id2]
  x0 <- x[x$pep_rank2 == 1L & x$nl_id == 1L, ]
  y0 <- x[x$pep_rank2 > 1L, ]
  z0 <- x[pep_rank2 == 1L & x$nl_id > 1L, ]
  x0[["pep_rank2"]] <- NULL
  y0[["pep_rank2"]] <- NULL
  z0[["pep_rank2"]] <- NULL
  x0[["nl_id"]] <- NULL
  y0[["nl_id"]] <- NULL
  z0[["nl_id"]] <- NULL
  rm(list = c("x"))
  gc()

  
  ## 3. keep the top-3 `pep_seq_mod`
  # the `pep_rank` here is after the "collapse" of NLs by only using the best;
  # net effect: the top-3 pep_seq_mod's, each represented by its best NL;
  # nevertheless, there can be ties in NL.
  message("\tSubset peptides by modifications: \"topn_mods_per_seq <= ", 
          topn_mods_per_seq, "\".")

  if (para) {
    x0 <- x0[order(x0[["uniq_id"]]), ] # for group split

    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    x0s <- parallel::clusterApply(
      cl, split(x0, find_chunkbreaks(x0[["uniq_id"]], n_cores)), calcpeprank_2)
    parallel::stopCluster(cl)
    
    x0 <- data.table::rbindlist(x0s, use.names = FALSE)
    rm(list = c("x0s"))
  }
  else {
    x0[, pep_rank := data.table::frank(-pep_score, ties.method = "min"), 
       by = list(uniq_id)]
  }
  
  x0 <- x0[pep_rank <= topn_mods_per_seq, ]
  x0[["pep_rank"]] <- NULL
  gc()

  
  ## 4 `pep_locprob` (the same `pep_seq`, different `pep_seq_mod`)
  # 4.1 separations into ambiguous x0 and non-ambiguous x1
  x0 <- x0[ , n_pep_seq_mod := .N, by = .(uniq_id)]
  x1 <- x0[n_pep_seq_mod == 1L, ] # single pep_seq_mod, no location ambiguity
  x0 <- x0[n_pep_seq_mod > 1L, ]
  x0[["n_pep_seq_mod"]] <- NULL
  x1[["n_pep_seq_mod"]] <- NULL
  
  # 4.2 probability and delta
  if (nrow(x0)) {
    locmod_indexes <- if (length(locmods))
      unname(mod_indexes[names(mod_indexes) %in% locmods])
    else 
      NULL
    
    if (length(locmod_indexes)) {
      message("\tCalculates peptide localization scores.")
      
      us <- split(x0[, c("uniq_id2", "pep_ivmod2", "pep_ms2_ideltas.")], 
                  x0[["uniq_id"]])
      gc()
      
      probs <- lapply(us, findLocFracsDF, locmod_indexes)
      
      if (FALSE) {
        deltas <- lapply(probs, function (x) {
          if (all(is.na(x)))
            NA_real_
          else {
            # may be more NA values, e.g. 3L, than 2L
            # topx <- x[which_topx2(x, 2L)]
            topx <- x[which_topx2(x, 2L, na.last = FALSE)]
            abs(topx[1] - topx[2])
          }
        })
      }

      deltas <- lapply(probs, function (x) {
        x <- x[!is.na(x)]
        len <- length(x)
        
        if (len == 1L)
          x
        else if (len > 1L) {
          topx <- x[which_topx2(x, 2L)]
          abs(topx[1] - topx[2])
        }
        else 
          NA_real_
      })
      
      us <- mapply(function (x, y, z) {
        x[["pep_locprob"]] <- y
        x[["pep_locdiff"]] <- z
        x
      }, us, probs, deltas, 
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
      
      us <- data.table::rbindlist(us)
      
      x0 <- quick_leftjoin(x0, us[, c("uniq_id2", "pep_locprob", "pep_locdiff")], 
                           by = "uniq_id2")
      
      rm(list = c("probs", "deltas", "us"))
      gc()
    }
    else {
      x0[["pep_locprob"]] <- NA_real_
      x0[["pep_locdiff"]] <- NA_real_
    }
  }
  else {
    col_nms <- c(names(x0), "pep_locprob", "pep_locdiff")
    x0 <- data.table::data.table(matrix(ncol = length(col_nms), nrow = 0L))
    colnames(x0) <- col_nms
    rm(list = "col_nms")
  }
  
  # 4.3 adds back x1
  if (nrow(x1)) {
    x1[["pep_locprob"]] <- 1.0
    x1[["pep_locdiff"]] <- 1.0
    x0 <- data.table::rbindlist(list(x0, x1), use.names = FALSE)
  }
  rm(list = c("x1"))
  
  # 4.4 adds back z0
  # ok to replace dplyr::left_join here:
  #   the same number of rows between the new and the old z0
  z0 <- quick_leftjoin(z0, x0[, c("uniq_id2", "pep_locprob", "pep_locdiff")], 
                       by = "uniq_id2")
  x0 <- data.table::rbindlist(list(x0, z0), use.names = FALSE)
  rm(list = "z0")
  
  # 4.5 adds back y0
  if (nrow(y0)) {
    y0[["pep_locprob"]] <- NA_real_
    y0[["pep_locdiff"]] <- NA_real_
    x0 <- data.table::rbindlist(list(x0, y0), use.names = FALSE)
  }
  rm(list = "y0")
  
  ## 5. clean-ups
  x0 <- x0[, -c("uniq_id", "uniq_id2")]
  x0[ , "pep_score" := round(pep_score, 2L)]
  x0[ , "pep_locprob" := round(pep_locprob, 2L)]
  x0[ , "pep_locdiff" := round(pep_locdiff, 2L)]
  gc()
  
  # 5.1 NEW `pep_rank`s across different `pep_seq`s under the same `query`
  # (this is different to the earlier `uniq_id` to differentiate LOCATIONS)
  # 
  # e.g., if MS evidence is equally feasible for : 
  #   pep_seq_1: EVEEDSEDEEMSEDE[E]D[D]S[SG]EEVVIPQKK
  #   pep_seq_2: EVEEDSEDEEMSEDE[D]D[S]S[GE]EEVVIPQKK
  # both will be kept (at the same rank)
  
  message("\tSubset peptide sequences by query: \"topn_seqs_per_query <= ", 
          topn_seqs_per_query, "\".")
  
  x0[, uniq_id3 := paste(pep_isdecoy, pep_scan_num, raw_file, sep = ".")]
  
  if (para) {
    x0 <- x0[order(x0[["uniq_id3"]]), ]

    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    x0s <- parallel::clusterApply(cl, 
      split(x0, find_chunkbreaks(x0[["uniq_id3"]], n_cores)), calcpeprank_3)
    parallel::stopCluster(cl)
    
    x0 <- data.table::rbindlist(x0s, use.names = FALSE)
    rm(list = c("x0s"))
  }
  else {
    x0[, pep_rank := data.table::frank(-pep_score, ties.method = "min"), 
       by = list(uniq_id3)]
  }
  
  x0 <- x0[pep_rank <= topn_seqs_per_query, ]
  data.table::setorder(x0, uniq_id3, -pep_score)
  x0$uniq_id3 <- NULL
  
  x0 <- x0[, pep_isdecoy := as.logical(pep_isdecoy)]
  x0[["pep_ms2_ideltas."]] <- NULL
  x0[["pep_ivmod2"]] <- NULL
  qs::qsave(x0, file.path(out_path, "temp", "peploc.rds"), preset = "fast")
  
  invisible(x0)
}


#' Helper for parallel frank
#' 
#' @param x A data.table object
calcpeprank_1 <- function (x)
{
  x[, pep_rank2 := data.table::frank(-pep_score, ties.method = "min"), 
    by = "uniq_id2"]
}


#' Helper for parallel frank
#' 
#' @param x0 A data.table object
calcpeprank_2 <- function (x0)
{
  x0[, pep_rank := data.table::frank(-pep_score, ties.method = "min"), 
     by = "uniq_id"]
}


#' Helper for parallel frank
#' 
#' @param x0 A data.table object
calcpeprank_3 <- function (x0) 
{
  x0[, pep_rank := data.table::frank(-pep_score, ties.method = "min"), 
     by = "uniq_id3"]
}


#' Finds the break points
#' 
#' @param vals An order vector.
#' @param n_chunks The number of chunks.
#' 
#' @examples 
#' \donttest{
#' library(mzion)
#' 
#' vals <- c(rep("B", 6), rep("A", 6), rep("E", 5), rep("D", 5), rep("C", 7))
#' brs <- mzion:::find_chunkbreaks(vals, 3L)
#' split(vals, brs)
#' }
find_chunkbreaks <- function (vals, n_chunks) 
{
  lenv <- length(vals)
  n_chunks <- min(length(unique(vals)), n_chunks)
  
  if (n_chunks == 1L)
    return(rep(1L, lenv))
  
  pos <- floor(lenv/n_chunks) * 1:n_chunks
  len <- length(pos)
  pos[len] <- lenv
  
  for (i in 1:(len - 1L)) {
    icr <- pos[i]
    inx <- icr + 1L
    vcr <- vals[icr]
    vnx <- vals[inx]
    
    while(vcr == vnx) {
      icr <- inx
      inx <- inx + 1L
      vcr <- vals[icr]
      vnx <- vals[inx]
    }
    
    pos[i] <- icr
  }
  
  ds <- c(pos[1], diff(pos))
  rep(seq_along(pos), ds)
}


#' Finds the localization fractions of STY.
#'
#' Counting statistics.
#'
#' @param df A data frame containing columns \code{pep_ivmod2} and
#'   \code{pep_ms2_ideltas.}.
#' @param locmod_indexes A vector to the modification indexes.
#'
#' @examples
#' library(mzion)
#' 
#' locmod_indexes <- c("8", "9", "a")
#' df <- data.frame(pep_ivmod2 = c("0080000", "0009000"), pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(1,2,3,4,5,6,8,9,10,11,12,13,14))
#' df$pep_ms2_ideltas.[2] <- list(c(1,2,3,4,5,6,8,9,10,12,13,14))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
#'
#' # Variable Acetyl (K) and fixed TMT6plex (K)
#' locmod_indexes <- "4"
#' df <- data.frame(pep_ivmod2 = c("04000000", "00000000"), pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(2,5,6,7,9,11,14,15,16))
#' df$pep_ms2_ideltas.[2] <- list(c(1,2,5,6,7,9,11,14))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
#'
#' df <- data.frame(pep_ivmod2 = c("0004000", "4000000"), pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(2,3,4,5,6,8,10,11,12,14))
#' df$pep_ms2_ideltas.[2] <- list(c(4,5,6,8,10,14))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
#'
#' df <- data.frame(pep_ivmod2 = c("00040402", "00040204", "00020202"), 
#'   pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(1,4,6,7,9,10,11,12,15,16))
#' df$pep_ms2_ideltas.[2] <- list(c(1,4,11,12,15,16))
#' df$pep_ms2_ideltas.[3] <- list(c(1,4,11,12,15,16))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
#'
#' locmod_indexes <- c("2", "4")
#' df <- data.frame(pep_ivmod2 = c("0004004000402", "0002004000402", "0004004000204"), 
#'   pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(1,3,4,5,6,7,8,10,14,15,16,17,18,19))
#' df$pep_ms2_ideltas.[2] <- list(c(4,5,6,7,8,10,14,15,16,17,18,19))
#' df$pep_ms2_ideltas.[3] <- list(c(1,3,4,5,6,7,8,10,14,16,17,18,19))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
#'
#' # Variable Digly (K), TMT6plex+Digly (K) fixed TMT6plex (K)
#' # 9 - TMT6plex+Digly (K)
#' # 6 - Digly (K)
#' # 2 - TMT6plex (K)
#' # different lengths but equal mass: 26 and 9
#' locmod_indexes <- c(2, 6, 9)
#' df <- data.frame(pep_ivmod2 = c("000000020000000602", "000000090000000002"), 
#'   pep_ms2_ideltas. = NA)
#' df$pep_ms2_ideltas.[1] <- list(c(1,2,3,4,5,6,19,20,21,22,23))
#' df$pep_ms2_ideltas.[2] <- list(c(1,2,3,4,5,6,19,20,      23))
#' ans <- mzion:::findLocFracsDF(df, locmod_indexes)
findLocFracsDF <- function (df, locmod_indexes = NULL) 
{
  ivms <- df[["pep_ivmod2"]]
  seqs <- df[["pep_ms2_ideltas."]]
  ivms <- .Internal(strsplit(ivms, "", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  naas <- length(ivms[[1]])
  
  # if (is.null(locmod_indexes))
  #   ps <- lapply(ivms, function (x) which(x != "0"))
  # else
  #   ps <- lapply(ivms, function (x) which(x %in% locmod_indexes))
  
  ps <- lapply(ivms, function (x) which(x %in% locmod_indexes))
  ns <- .Internal(unlist(lapply(ps, length), recursive = FALSE, use.names = FALSE))
  oks <- ns > 0L

  # (1) set aside entries without the target modifications, e.g. No STY sites
  ans <- numeric(length(ps))
  ans[!oks] <- NA_real_
  ps <- ps[oks]
  seqs <- seqs[oks]
  ns <- ns[oks]
  
  # (2) single non-trivial entry and thus unambiguous
  lenp <- length(ps)
  
  if (!lenp)
    return(ans)

  if (lenp == 1L) {
    ans[oks] <- 1.00
    return(ans)
  }

  # (3) multiple entries (lenp > 1L)
  
  len <- lenp - 1L
  ncr <- nnx <- vector("integer", len)
  
  for (i in 1:len) {
    inx <- i + 1L
    
    pcr <- ps[[i]] # all(pcr <= naas)
    pnx <- ps[[inx]]
    seqcr <- seqs[[i]]
    seqnx <- seqs[[inx]]
    lencr <- ns[[i]]
    lennx <- ns[[inx]]
    
    # ps different lengths: see the Digly (K) and TMT (K) example
    bmin <- min(pcr[[1]], pnx[[1]])
    bmax <- max(pcr[[lencr]], pnx[[lennx]]) # bmax <= naas
    
    # b-ions
    bcr <- seqcr <= naas
    bnx <- seqnx <= naas
    bseqcr <- seqcr[bcr]
    bseqnx <- seqnx[bnx]
    
    # y-ions
    yseqcr <- seqcr[!bcr] - naas
    yseqnx <- seqnx[!bnx] - naas
    ymax <- naas - bmin
    ymin <- naas - bmax
    
    # b & y
    ncr[[i]] <- sum(bseqcr >= bmin & bseqcr < bmax) + sum(yseqcr > ymin & yseqcr <= ymax)
    nnx[[i]] <- sum(bseqnx >= bmin & bseqnx < bmax) + sum(yseqnx > ymin & yseqnx <= ymax)
  }
  
  ans[oks] <- concatFracs(ncr, nnx)
  
  ans
}


#' Concatenates localization fractions.
#'
#' The probability of the second localization with \eqn{y} in \eqn{x = 0; y = 6}
#' is not any more probable than that in \eqn{x = 0; y = 1}. This is different
#' to the binomial model in A-score, which will lead to the finding that the
#' probability of localization 2 at \eqn{x = 0; y = 6} is much more probable.
#'
#' @param x A vector of counts of MS2 matches between two adjacent STY sites
#'   (for a preceding match: pep_ivmod).
#' @param y A vector of counts of MS2 matches between two adjacent STY sites
#'   (for a next match: pep_ivmod).
#' @param d A small positive number to handle value \eqn{0}.
#'
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' mzion:::concatFracs(c(3, 3, 5), c(2, 1, 3))
#' mzion:::concatFracs(c(0, 6), c(1, 3))
#' mzion:::concatFracs(c(0, 6), c(0, 3))
#' mzion:::concatFracs(c(0, 2), c(0, 2))
#' mzion:::concatFracs(c(1, 2), c(1, 2))
#' 
#' mzion:::concatFracs(c(1, 2), c(1, 4))
#' mzion:::concatFracs(c(0, 2), c(0, 4))
#' mzion:::concatFracs(c(1, 2), c(0, 4)) # near 1 0 0
#' 
#' mzion:::concatFracs(0, 2)
#' mzion:::concatFracs(2, 0)
#' mzion:::concatFracs(0, 0) # should not occur
#' }
concatFracs <- function (x, y, d = .001) 
{
  len <- length(x)
  
  if (identical(x, y))
    return(rep(1/(len + 1L), (len + 1L)))
  
  x <- x + d
  y <- y + d
  
  if (len <= 1L)
    return(c(x, y)/sum(x, y))
  
  for (i in 2:len) {
    fct <- y[i-1]/x[i]
    x[i] <- x[i] * fct
    y[i] <- y[i] * fct
  }
  
  ans <- c(x[1], y)
  
  ans/sum(ans)
}


#' Interpolates missing values in a time series.
#'
#' From \code{forecast}. 
#' 
#' @param linear Linear.
#' @inheritParams tsoutliers
na.interp <- function (x, lambda = NULL, 
                       linear = (frequency(x) <= 1 | 
                                   sum(!is.na(x)) <= 2 * frequency(x))) 
{
  missng <- is.na(x)
  
  if (sum(missng) == 0L) 
    return(x)

  origx <- x
  rangex <- range(x, na.rm = TRUE)
  drangex <- rangex[2L] - rangex[1L]
  
  if (is.null(tsp(x))) 
    x <- ts(x)

  if (length(dim(x)) > 1) {
    if (NCOL(x) == 1) 
      x <- x[, 1]
    else 
      stop("The time series is not univariate.")
  }
  
  if (!is.null(lambda)) {
    x <- BoxCox(x, lambda = lambda)
    lambda <- attr(x, "lambda")
  }
  
  freq <- frequency(x)
  tspx <- tsp(x)
  n <- length(x)
  tt <- 1:n
  idx <- tt[!missng]
  
  if (linear) 
    x <- ts(approx(idx, x[idx], tt, rule = 2)$y)
  else {
    if ("msts" %in% class(x)) 
      K <- pmin(trunc(attributes(x)$msts/2), 20L)
    else 
      K <- min(trunc(freq/2), 5)

    X <- cbind(fourier(x, K), 
               poly(tt, degree = pmin(pmax(trunc(n/10), 1), 6L)))

    fit <- lm(x ~ X, na.action = na.exclude)
    pred <- predict(fit, newdata = data.frame(X))
    x[missng] <- pred[missng]
    fit <- mstl(x, robust = TRUE)
    sa <- seasadj(fit)
    sa <- approx(idx, sa[idx], 1:n, rule = 2)$y
    seas <- seasonal(fit)
    
    if (NCOL(seas) > 1) 
      seas <- rowSums(seas)
    
    x[missng] <- sa[missng] + seas[missng]
  }
  
  if (!is.null(lambda)) 
    x <- InvBoxCox(x, lambda = lambda)
  
  tsp(x) <- tspx
  
  if (!linear & (max(x) > rangex[2L] + 0.5 * drangex | min(x) < 
                 rangex[1L] - 0.5 * drangex)) 
    return(na.interp(origx, lambda = lambda, linear = TRUE))
  else 
    return(x)
}


#' Checks if data are constant.
#' 
#' @param x A vector.
is.constant <- function (x) 
{
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  
  return(isTRUE(all.equal(x, y)))
}


#' Identifies and replace outliers in a time series.
#' 
#' From \code{forecast}.
#' 
#' @param x Time series
#' @param lambda Box-Cox transformation parameter. If lambda="auto", then a
#'   transformation is automatically selected using BoxCox.lambda. The
#'   transformation is ignored if NULL. Otherwise, data transformed before model
#'   is estimated.
#' @param iterate The number of iterations.
tsoutliers <- function (x, iterate = 2, lambda = NULL) 
{
  n <- length(x)
  freq <- frequency(x)
  missng <- is.na(x)
  nmiss <- sum(missng)
  
  xx <- if (nmiss > 0L) 
    na.interp(x, lambda = lambda)
  else 
    x
  
  if (is.constant(xx)) 
    return(list(index = integer(0), replacements = numeric(0)))
  
  tt <- 1:n
  mod <- supsmu(tt, xx)
  resid <- xx - mod$y
  
  if (nmiss) resid[missng] <- NA

  resid.q <- quantile(resid, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(resid.q)
  limits <- resid.q + 3 * iqr * c(-1, 1)
  
  outliers <- if ((limits[2] - limits[1]) > 1e-14) 
    which((resid < limits[1]) | (resid > limits[2]))
  else 
    numeric(0)
  
  x[outliers] <- NA
  x <- na.interp(x, lambda = lambda)
  
  if (iterate > 1) {
    tmp <- tsoutliers(x, iterate = 1, lambda = lambda)
    
    if (length(tmp$index)) {
      outliers <- sort(unique(c(outliers, tmp$index)))
      x[outliers] <- NA
      
      if (sum(!is.na(x)) == 1L) 
        x[is.na(x)] <- x[!is.na(x)]
      else 
        x <- na.interp(x, lambda = lambda)
    }
  }
  
  invisible(list(index = outliers, replacements = x[outliers]))
}


