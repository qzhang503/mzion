#' Adds secondary ions of b0, y0 etc.
#' 
#' @param ms2s A vector of theoretical MS2 m-over-z values.
#' @inheritParams matchMS
add_seions <- function (ms2s, type_ms2ions = "by", digits = 5L) 
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
  } else if (type_ms2ions == "ax") {
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
  } else if (type_ms2ions == "cz") {
    proton <- 1.00727647
    
    cs <- ms2s[1:(len/2)]
    zs <- ms2s[(len/2+1):len]
    
    c2s <- (cs + proton)/2
    z2s <- (zs + proton)/2
    
    round(c(c2s, z2s), digits = digits)
  }
}


#' Matches two lists.
#' 
#' Not currently used. Without making a data frame.
#' 
#' @param a The left vector.
#' @param b The right vector.
#' @examples
#' \donttest{
#' a <- c(3, 4, 1, 2, 5)
#' b <- 2
#' 
#' list_leftmatch(a, b)
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
#' @param df Two lists of \code{theo} and matched \code{expt} m-over-z.
#' @param nms The names (character strings indicating the names and position of
#'   variable modifications).
#' @inheritParams calc_probi
#' @import dplyr
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @examples
#' \donttest{
#' df <- list(theo = c(390.2009, 550.2315, 710.2622, 809.3306, 880.3677, 
#'                     995.3946, 1151.4957, 175.1190, 290.1459, 361.1830, 
#'                     460.2514, 620.2821, 780.3127, 940.3434), 
#'            expt = c(390.2008, 550.2323, 710.2624, 809.3301, 880.3662, 
#'                     995.3970, NA, 175.1191, 290.1458, 361.1832, 
#'                     460.2517, 620.2880, 780.3126, 940.3438))
#' 
#' expt_moverzs <- c(110.0717, 112.0509, 112.0873, 113.0713, 115.0869, 
#'                   116.0167, 116.0709, 126.1280, 127.1250, 127.1313, 
#'                   128.1284, 128.1346, 129.1317, 129.1380, 130.0978, 
#'                   130.1350, 130.1413, 131.1384, 133.0432, 136.0619, 
#'                   139.0504, 157.1083, 158.0563, 158.0925, 159.0763, 
#'                   173.1498, 175.1191, 176.1223, 176.1600, 178.0646, 
#'                   186.1531, 188.1597, 212.1031, 230.1144, 230.1702, 
#'                   232.1115, 248.1808, 255.1082, 261.6220, 264.6251, 
#'                   273.1193, 275.6198, 284.1326, 284.6348, 290.1458, 
#'                   310.1302, 321.0684, 329.1280, 344.1566, 359.6652, 
#'                   361.1832, 362.2062, 390.2008, 407.2268, 420.1369, 
#'                   459.2220, 460.2517, 481.0993, 491.1726, 522.2371, 
#'                   539.7526, 540.2550, 550.2323, 551.2338, 567.2598, 
#'                   576.7457, 584.7561, 585.2570, 585.7582, 586.2587, 
#'                   619.2524, 620.2880, 682.2661, 683.2708, 710.2624, 
#'                   711.2637, 718.3199, 780.3126, 781.3342, 782.3379, 
#'                   809.3301, 810.3351, 880.3662, 881.3693, 921.3688, 
#'                   922.3726, 923.3927, 924.3849, 940.3438, 941.3491, 
#'                   995.3970, 996.3967, 997.3690, 998.3657, 1011.3803, 
#'                   1012.3803, 1013.3842, 1014.3911, 1015.3893, 1016.3904)
#' 
#' expt_ints <- c(12810.80, 14142.40, 58754.70, 12451.00, 29055.70, 
#'                45291.00, 63865.00, 250674.00, 261949.00, 179089.00, 
#'                253049.00, 190448.00, 240766.00, 275813.00, 28354.30, 
#'                219360.00, 189991.00, 229268.00, 60450.10, 12415.10, 
#'                11351.30, 17766.50, 29119.50, 105925.00, 10832.10, 
#'                15792.60, 707208.00, 29073.60, 12632.30, 18499.40, 
#'                18826.60, 33715.50, 12418.70, 18046.80, 164112.00, 
#'                15920.80, 13090.70, 24475.10, 14995.90, 49102.20, 
#'                56960.70, 17143.40, 93462.60, 15536.50, 23416.20, 
#'                15584.90, 30465.80, 12715.50, 31551.40, 12031.20, 
#'                22367.40, 55729.80, 77277.60, 19092.70, 29280.50, 
#'                17574.20, 21419.20, 10681.80, 8850.03, 27071.20, 
#'                25317.90, 13518.70, 55593.70, 14856.20, 30853.00, 
#'                11551.30, 22966.50, 515863.00, 243235.00, 11709.10, 
#'                47070.90, 19336.10, 57334.40, 11747.80, 147943.00, 
#'                42007.30, 18287.80, 51689.70, 62069.30, 13403.00, 
#'                84499.10, 24180.00, 47260.80, 13985.00, 14132.90, 
#'                10097.10, 12578.40, 13326.10, 49003.80, 12951.10, 
#'                98873.10, 38704.30, 12971.00, 8924.51, 89986.00, 
#'                162963.00, 119860.00, 114223.00, 125885.00, 32305.60)
#' 
#' calc_probi_byvmods(df, nms = "0000000", expt_moverzs, expt_ints, N = 190)
#' 
#' # 
#' df2 <- df
#' df2$expt[8] <- NA
#' calc_probi_byvmods(df2, nms = "0000000", expt_moverzs, expt_ints, N = 190)
#' }
calc_probi_byvmods <- function (df, nms, expt_moverzs, expt_ints, 
                                N, type_ms2ions = "by", topn_ms2ions = 100L, 
                                penalize_sions = TRUE, ppm_ms2 = 25L, 
                                digits = 5L) 
{
  # N - the total number of features (white and black balls)
  # k - the number of sampled features
  # m - the numbers of theoretical features (white balls)
  # n - the number of noise (black balls)
  
  df_theo <- df$theo
  m <- length(df$theo)
  
  # OK: (N < m) -> (n < 0L)
  # if (N < m) N <- m
  
  ## matches additionally against secondary ions
  df2 <- add_seions(df_theo, type_ms2ions = type_ms2ions, digits = digits)
  df2 <- find_ppm_outer_bycombi(df2, expt_moverzs, ppm_ms2) # 132 us
  df2$theo <- round(df2$theo, digits = digits)
  
  # subtracts `m` and the counts of secondary b0, y0 matches etc. from noise
  # (OK if n < 0L)
  
  n <- N - m - sum(!is.na(df2$expt))
  
  ## step 0: the original tidyverse approach
  
  # m2 <- nrow(df2)
  # 
  # y2 <- df2 %>% 
  #   left_join(expts, by = "expt") %>% 
  #   `[[`("int") %>% 
  #   split(rep(seq_len(m2/m), each = m)) %>% 
  #   Reduce(`%+%`, .) %>% 
  #   data.frame(idx = seq_len(m), int2 = .)
  # 
  # y <- left_join(expts, df %>% mutate(idx = row_number()), by = "expt") %>% 
  #   dplyr::left_join(y2, by = "idx") %>% 
  #   mutate(int = ifelse(is.na(int2), int, int + int2)) %>% 
  #   select(-c("int2", "idx")) %>% 
  #   arrange(-int) %>% 
  #   mutate(k = row_number(), x = k - cumsum(is.na(theo))) %>% 
  #   filter(!is.na(theo))
  
  ## step 1: compiles secondary intensities
  
  # (1.1)
  len <- length(df2$expt)
  
  i_se <- match(df2$expt, expt_moverzs) # secondary ions (df2) in expts
  i_s <- which(!is.na(i_se)) # indexes in df2
  i_e <- i_se[i_s] # indexes of the matches in expts
  
  df2$int <- rep(NA, len) # matched intensities; .7 us
  df2$int[i_s] <- expt_ints[i_e] # the corresponding intensities from expts
  
  # (1.2) collapse b0, b*, b2 etc.
  f <- rep(seq_len(len/m), each = m)
  int2 <- .Internal(split(df2$int, as.factor(f)))
  int2 <- Reduce(`%+%`, int2)

  y2 <- list(idx = 1:m, int2 = int2)
  
  ## step 2 (join expts and "theo", "idx" from the primary df): 
  # (y <- left_join(expts, df %>% mutate(idx = row_number()), by = "expt"))
  # 
  #   expts: "expt" (expt_moverzs), "int" (expt_ints)
  #   df: "theo" (m/z), "expt" (m/z), "idx" (indexes)
  #   -> y: "theo", "expt", "int" "idx"
  
  # (2.1)
  df$idx <- 1:m
  
  i_ep <- match(expt_moverzs, df$expt) # expts in primary ions (df)
  i_e <- which(!is.na(i_ep)) # indexes in expts
  i_p <- i_ep[i_e] # indexes of matches in primary df
  
  # (2.2)
  nu <- rep(NA, topn_ms2ions)
  y <- list(expt = expt_moverzs, int = expt_ints)
  
  y$theo <- nu
  y$idx <- nu
  
  y$theo[i_e] <- df$theo[i_p]
  y$idx[i_e] <- df$idx[i_p]
  
  ## step 3: add `int2` from `y2`
  #  (dplyr::left_join(y, y2, by = "idx"))
  
  i_yy2 <- match(y$idx, y2$idx)
  i_y <- which(!is.na(i_yy2))
  i_y2 <- i_yy2[i_y]
  
  y$int2 <- nu
  y$int2[i_y] <- y2$int2[i_y2]
  
  ## step 4: collapses `int2` to `int`
  #  (mutate(y, int = ifelse(is.na(int2), int, int + int2)))

  y$int <- y$int %+% y$int2
  y$int2 <- NULL
  y$idx <- NULL
  
  ## step 5: arrange(-int)
  
  idx <- order(y$int, decreasing = TRUE, method = "radix", na.last = TRUE) # 16.4 us
  y$expt <- y$expt[idx]
  y$int <- y$int[idx]
  y$theo <- y$theo[idx]
  
  ## step 6: mutate(k = row_number(), x = k - cumsum(is.na(theo)))
  
  k <- 1:topn_ms2ions
  x <- k - cumsum(is.na(y$theo))
  
  y$k <- k
  y$x <- x
  
  ## step 7: filter(!is.na(theo))
  
  idx <- !is.na(y$theo)
  
  y$expt <- y$expt[idx]
  y$int <- y$int[idx]
  y$theo <- y$theo[idx]
  y$k <- y$k[idx]
  y$x <- y$x[idx]
  
  ## Probability
  # note: x <= k <= x + n
  
  x <- y$x
  k <- y$k
  
  # (to have sufficient counts of noise)
  # (also guaranteed n > 0L)
  
  n <- max(n, topn_ms2ions + k[length(k)])
  
  prs <- mapply(dhyper, x[-c(1:2)], m, n, k[-c(1:2)])
  pr <- min(prs, na.rm = TRUE)

  list(pep_ivmod = nms, 
       pep_prob = pr, 
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
                              N, type_ms2ions, topn_ms2ions, 
                              penalize_sions, ppm_ms2, digits) 
{
  ## for different positions: $TNLAMMR$`0000500`, $TNLAMMR$`0000050`
  #    the same `pep_seq`, `theo_ms1` for different mod positions
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
                  penalize_sions = penalize_sions, 
                  ppm_ms2 = ppm_ms2, 
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
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
#' @import dplyr
#' @importFrom purrr map
calc_probi <- function (mts, expt_moverzs, expt_ints, 
                        N, type_ms2ions = "by", topn_ms2ions = 100L, 
                        penalize_sions = TRUE, ppm_ms2 = 25L, digits = 5L) 
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
      penalize_sions = penalize_sions, 
      ppm_ms2 = ppm_ms2, 
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
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
#' @import purrr
scalc_pepprobs <- function (entry, topn_ms2ions = 100L, type_ms2ions = "by", 
                            penalize_sions = FALSE, ppm_ms2 = 25L, digits = 4L) 
{
  # only one experimental set of values and thus `[[1]]`
  expt_moverzs <- entry$ms2_moverz[[1]]
  expt_ints <- entry[["ms2_int"]][[1]]

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
  mts <- entry$matches[[1]]
  
  N <- entry$ms2_n[[1]]
  topn_ms2ions <- min(topn_ms2ions, N)
  
  out <- calc_probi(mts = mts, 
                    expt_moverzs = expt_moverzs, 
                    expt_ints = expt_ints, 
                    N = N, 
                    type_ms2ions = type_ms2ions, 
                    topn_ms2ions = topn_ms2ions, 
                    penalize_sions = penalize_sions, 
                    ppm_ms2 = ppm_ms2, 
                    digits = digits)

  uniq_id <- .Internal(unlist(entry$uniq_id, recursive = FALSE, use.names = FALSE))
  
  out <- lapply(out, function (x) {
    x$uniq_id <- uniq_id
    x
  })

  invisible(out)
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
calc_pepprobs_i <- function (res, topn_ms2ions = 100L, type_ms2ions = "by", 
                             penalize_sions = FALSE, ppm_ms2 = 25L, 
                             out_path = "~/proteoM/outs", digits = 5L) 
{
  n_rows <- nrow(res)
  
  if (n_rows) { # 3 ms
    probs <- split.data.frame(res, seq_len(n_rows)) 
    
    probs <- lapply(probs, scalc_pepprobs, 
                    topn_ms2ions = topn_ms2ions, 
                    type_ms2ions = type_ms2ions, 
                    penalize_sions = penalize_sions, 
                    ppm_ms2 = ppm_ms2, 
                    digits = digits)
    
    probs <- .Internal(unlist(probs, recursive = FALSE, use.names = FALSE))
    
    probs <- dplyr::bind_rows(probs) # 276 us
  } else {
    probs <- data.frame(
      pep_seq = as.character(), 
      pep_ivmod = as.character(), 
      pep_prob = as.numeric(), 
      pri_matches = list(), 
      sec_matches = list(), 
      scan_num = as.integer())
  }
  
  invisible(probs)
}


#' Calculates the scores of peptides.
#'
#' @param penalize_sions Logical; if TRUE, penalizes secondary ions of b0, y0
#'   etc. with lower weights in peptide scoring.
#' @inheritParams matchMS
#' @import parallel
calc_pepscores <- function (topn_ms2ions = 100L, type_ms2ions = "by", 
                            target_fdr = 0.01, fdr_type = "psm", 
                            min_len = 7L, max_len = 50L, 
                            penalize_sions = FALSE, ppm_ms2 = 25L, 
                            out_path = "~/proteoM/outs", digits = 5L) 
{
  ok_res <- file.path(out_path, "scores.rds")
  
  if (file.exists(ok_res)) {
    message("Found cached peptide scores.")
    return(readRDS(ok_res))
  } 
  
  
  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  message("Calculating peptide scores.")
  
  ## --- Target ---
  list_t <- list.files(path = file.path(out_path, "temp"), 
                       pattern = "^ion_matches_\\d+\\.rds$")
  
  if (!length(list_t)) 
    stop("Target matches not found.", call. = FALSE)
  
  nms_t <- gsub("^ion_matches_(\\d+)\\.rds$", "\\1", list_t)
  
  ord <- order(as.integer(nms_t))
  nms_t <- nms_t[ord]
  list_t <- list_t[ord]
  
  for (i in seq_along(list_t)) {
    calcpepsc(file = list_t[i], 
              topn_ms2ions = topn_ms2ions, 
              type_ms2ions = type_ms2ions, 
              penalize_sions = penalize_sions, 
              ppm_ms2 = ppm_ms2, 
              out_path = out_path, 
              digits = digits)
    
    gc()
  }
  
  ## --- Decoy ---
  list_d <- list.files(path = file.path(out_path, "temp"), 
                       pattern = "^ion_matches_rev_[0-9]+\\.rds$")
  
  if (!length(list_d)) 
    stop("Decoy matches not found.", call. = FALSE)
  
  if (length(list_d) > 1L) {
    warning("More than one decoy matches found: '", decoy, "'.\n", 
            "The first will be used.", call. = FALSE)
    list_d <- list_d[[1]]
  }
  
  nms_d <- gsub("^ion_matches_(rev_[0-9]+)\\.rds$", "\\1", list_d)
  
  rev <- calcpepsc(file = list_d, 
                   topn_ms2ions = topn_ms2ions, 
                   type_ms2ions = type_ms2ions, 
                   penalize_sions = penalize_sions, 
                   ppm_ms2 = ppm_ms2, 
                   out_path = out_path, 
                   digits = digits)
  
  if (is.null(rev)) 
    rev <- readRDS(file.path(out_path, "temp", paste0("pepscores_", nms_d, ".rds")))
  
  rev <- rev %>% 
    list() %>% 
    `names<-`(nms_d)

  ## --- FDR (target + decoy) --- 
  # Reload targets
  out <- local({
    files <- list.files(path = file.path(out_path, "temp"), 
                        pattern = "^pepscores_\\d+\\.rds$")
    
    if (!length(files)) 
      stop("Target scores not found.", call. = FALSE)
    
    nms <- gsub("^pepscores_(\\d+)\\.rds$", "\\1", files)
    
    ord <- order(as.integer(nms))
    nms <- nms[ord]
    files <- files[ord]
    
    out <- lapply(files, function (x) readRDS(file.path(out_path, "temp", x)))
    names(out) <- nms
    
    out <- c(out, rev)
  })
  
  rm(list = c("rev"))
  gc()
  
  ## Probability cut-offs
  prob_cos <- calc_pepfdr(out, 
                          nms = nms_d, 
                          target_fdr = target_fdr, 
                          fdr_type = fdr_type, 
                          min_len = min_len, 
                          max_len = max_len, 
                          out_path = out_path)
  
  prob_cos <- local({
    idxes <- prob_cos %>% .[. <= target_fdr]
    
    fct_homol <- if (length(idxes)) 
      target_fdr/max(idxes, na.rm = TRUE)
    else 
      1L

    prob_cos <- prob_cos * fct_homol
  })
  
  ## Outputs
  prob_cos <- prob_cos%>% 
    data.frame(pep_len = as.numeric(names(.)), pep_prob_co = .) 

  fct_score <- 10
  
  oks <- purrr::map_lgl(out, function (x) nrow(x) > 0L)
  out <- dplyr::bind_rows(out[oks])
  
  if (!nrow(out)) 
    stop("No PSM matches for scoring. Consider different search parameters.", 
         call. = FALSE)

  # Adjusted p-values
  out <- out %>% 
    dplyr::left_join(prob_cos, by = "pep_len") %>% 
    dplyr::mutate(pep_issig = ifelse(pep_prob <= pep_prob_co, TRUE, FALSE), 
                  pep_adjp = p.adjust(pep_prob, "BH"))
  
  prob_cos <- purrr::map_dbl(prob_cos$pep_prob_co, function (x) {
    row <- which.min(abs(log10(out$pep_prob/x)))
    out[row, ]$pep_adjp
  }) %>% 
    dplyr::bind_cols(prob_cos, pep_adjp_co = .) %T>% 
    saveRDS(file.path(out_path, "temp", "pep_prob_cos.rds"))
  
  out <- out %>% 
    dplyr::left_join(prob_cos[, c("pep_len", "pep_adjp_co")], by = "pep_len") %>% 
    dplyr::mutate(pep_score = -log10(pep_adjp) * fct_score, 
                  pep_score = ifelse(pep_score > 250, 250, pep_score), 
                  pep_score_co = -log10(pep_adjp_co) * fct_score) %>% 
    dplyr::select(-c("pep_prob", "pep_adjp", "pep_prob_co", "pep_adjp_co")) %T>%
    saveRDS(ok_res)
}


#' Helper of \link{calc_pepscores}.
#' 
#' @param file A file name of \code{ion_matches_}.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
calcpepsc <- function (file, topn_ms2ions = 100L, type_ms2ions = "by", 
                       penalize_sions = FALSE, ppm_ms2 = 25L, out_path = NULL, 
                       digits = 4L) 
{
  df <- readRDS(file.path(out_path, "temp", file))
  
  if (!nrow(df)) {
    idx <- gsub("^ion_matches_(.*)\\.rds$", "\\1", file)
    
    # list tables
    cols_lt <- c("raw_file", "pep_mod_group", "scan_num", 
                 "ms2_moverz", "ms2_int", "pri_matches", "sec_matches")
    dfa <- data.frame(matrix(ncol = length(cols_lt), nrow = 0L))
    colnames(dfa) <- cols_lt
    saveRDS(dfa, file.path(out_path, "temp", paste0("list_table_", idx, ".rds")))

    # scores
    cols_scores <- c("ms2_n", "scan_title", "ms1_moverz", "ms1_mass", 
                     "ms1_int", "ms1_charge", "ret_time", "scan_num", "raw_file", 
                     "pep_mod_group", "frame", "pep_fmod", "pep_vmod", "pep_isdecoy", 
                     "pep_seq", "theo_ms1", "pep_ivmod", "pep_prob", "pep_len")
    dfb <- data.frame(matrix(ncol = length(cols_scores), nrow = 0L))
    colnames(dfb) <- cols_scores
    saveRDS(dfb, file.path(out_path, "temp", paste0("pepscores_", idx, ".rds")))
    
    return (NULL)
  }
  
  df$uniq_id <- paste(df$scan_num, df$raw_file, sep = "@")
  
  esscols <- c("ms2_moverz", "ms2_int", "matches", "ms2_n", "uniq_id")
  df2 <- df[, -which(names(df) %in% esscols), drop = FALSE]
  df <- df[, esscols, drop = FALSE]
  
  # otherwise, chunksplit return NULL
  #   -> res[[i]] <- NULL 
  #   -> length(res) shortened by 1

  n_chunks <- detect_cores(16L)^2

  if (!is.null(df)) {
    df <- suppressWarnings(chunksplit(df, n_chunks, "row"))
    gc()
  }
  
  #  8L: 1.00 mins
  # 16L: 57.04 secs
  # 32L: 1.04 mins
  
  if (length(df) >= n_chunks) {
    
    n_cores <- detect_cores(16L)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(cl, list("%>%"), 
                            envir = environment(magrittr::`%>%`))
    
    parallel::clusterExport(cl, list("scalc_pepprobs"), 
                            envir = environment(proteoM:::scalc_pepprobs))

    # parallel::clusterExport(cl, 
    #                         list("calc_pepprobs_i", 
    #                              "scalc_pepprobs", 
    #                              "calc_probi", 
    #                              "calc_probi_bypep", 
    #                              "calc_probi_byvmods", 
    #                              "add_seions", 
    #                              "find_ppm_outer_bycombi"), 
    #                         envir = environment(proteoM:::scalc_pepprobs))
    
    probs <- parallel::clusterApplyLB(cl, df, 
                                      calc_pepprobs_i, 
                                      topn_ms2ions = topn_ms2ions, 
                                      type_ms2ions = type_ms2ions, 
                                      penalize_sions = penalize_sions, 
                                      ppm_ms2 = ppm_ms2,
                                      out_path = out_path, 
                                      digits = digits)
    
    parallel::stopCluster(cl)
    gc()
    
    probs <- dplyr::bind_rows(probs)
  } else {
    if (is.data.frame(df)) df <- list(df)
    
    probs <- lapply(df, calc_pepprobs_i, 
                    topn_ms2ions = topn_ms2ions, 
                    type_ms2ions = type_ms2ions, 
                    penalize_sions = penalize_sions, 
                    ppm_ms2 = ppm_ms2,
                    out_path = out_path, 
                    digits = digits) %>% 
      dplyr::bind_rows()
  }
  
  gc()
  
  ## Reassemble `df`
  df <- dplyr::bind_rows(df)
  df <- df[, -which(names(df) == "matches"), drop = FALSE]
  
  df <- dplyr::bind_cols(df, df2)
  df <- quick_rightjoin(df, probs, "uniq_id")
  df <- df[, -which(names(df) == "uniq_id"), drop = FALSE]
  df <- post_pepscores(df)
  
  rm(list = c("df2", "probs"))
  gc()
  
  ## Outputs 
  # (can be decoy => .*, not \\d+)
  idx <- gsub("^ion_matches_(.*)\\.rds$", "\\1", file)
  
  cols_a <- c("raw_file", "pep_mod_group", "scan_num")
  cols_b <- c("ms2_moverz", "ms2_int", "pri_matches", "sec_matches")
  
  saveRDS(df[, c(cols_a, cols_b), drop = FALSE], 
          file.path(out_path, "temp", paste0("list_table_", idx, ".rds")))
  
  df <- df[, -which(names(df) %in% cols_b), drop = FALSE]
  saveRDS(df, file.path(out_path, "temp", paste0("pepscores_", idx, ".rds")))

  invisible(df)
}


#' Cleanups post \link{calc_pepscores}.
#' 
#' @param df A results after pep_scores.
post_pepscores <- function (df) 
{
  df[["scan_num"]] <- as.numeric(df[["scan_num"]]) # 70 us
  df[["pep_len"]] <- stringr::str_length(df[["pep_seq"]]) # 6 us
  df[["scan_title"]] <- as.character(df[["scan_title"]]) # 69 us
  df[["ms1_moverz"]] <- as.numeric(df[["ms1_moverz"]]) # 4.5 us
  df[["ms1_mass"]] <- as.numeric(df[["ms1_mass"]]) # 5.4 us
  df[["ms1_int"]] <- as.numeric(df[["ms1_int"]]) # 5.5 us
  df[["ms1_charge"]] <- as.character(df[["ms1_charge"]]) # 5.9 us
  df[["ret_time"]] <- as.integer(df[["ret_time"]]) # 5.6 us
  df[["ms2_n"]] <- as.integer(df[["ms2_n"]]) # 5.3 us
  df[["pep_fmod"]] <- as.character(df[["pep_fmod"]]) # 4.6 us
  df[["pep_vmod"]] <- as.character(df[["pep_vmod"]]) # 4.8 us
  
  invisible(df)
}


#' Finds the cut-off in peptide scores for a given \code{target_fdr}.
#' 
#' Assume normal distribution for log(decoy_score).
#' 
#' @param td A data frame of target-decoy results at a given peptide length.
#' @inheritParams matchMS
find_pepscore_co1 <- function (td, target_fdr) 
{
  target <- dplyr::filter(td, !pep_isdecoy)
  decoy <- dplyr::filter(td, pep_isdecoy)
  
  nt <- nrow(target)
  nd <- nrow(decoy)
  
  if (nd <= 5L) return(NA)
  
  n <- nt + nd
  lambt <- nt/(n)
  lambd <- 1 - lambt
  
  vecd <- log2(decoy$pep_score)
  sigmad <- sd(vecd, na.rm = TRUE)
  mud <- mean(vecd, na.rm = TRUE)
  
  if (is.na(sigmad)) return(NA)
  
  xs <- seq(mud + sigmad, mud + 3*sigmad, 0.014355293)
  
  for (i in seq_along(xs)) {
    y <- (1 - pnorm(xs[i], mud, sigmad)) * nd / n
    
    if (y <= target_fdr) break
  }
  
  2^(xs[i])
}


#' Finds the cut-off in peptide scores for a given \code{target_fdr}.
#' 
#' Assume log-normal for decoy scores.
#' 
#' @param td A data frame of target-decoy results at a given peptide length.
#' @inheritParams matchMS
find_pepscore_co2 <- function (td, target_fdr) 
{
  target <- dplyr::filter(td, !pep_isdecoy)
  decoy <- dplyr::filter(td, pep_isdecoy)
  
  nt <- nrow(target)
  nd <- nrow(decoy)
  
  if (nd <= 5L) return(NA)
  
  n <- nt + nd
  lambt <- nt/(n)
  lambd <- 1 - lambt
  
  vecd <- decoy$pep_score
  sigmad <- sd(vecd, na.rm = TRUE)
  mud <- mean(vecd, na.rm = TRUE)
  
  if (is.na(sigmad)) return(NA)
  
  xs <- seq( mud + 4*sigmad, mud + sigmad, -.1)
  
  for (i in seq_along(xs)) {
    y <- (1 - plnorm(xs[i], mud, sigmad, lower.tail = FALSE)) * nd / n
    
    if (y >= target_fdr) break
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
probco_bypeplen <- function (len, td, fdr_type, target_fdr, out_path) 
{
  td <- dplyr::filter(td, pep_len == len)
  
  if (fdr_type %in% c("peptide", "protein")) {
    td <- dplyr::arrange(td, pep_seq, pep_prob)
    td <- dplyr::group_by(td, pep_seq)
    td <- dplyr::filter(td, row_number() == 1L)
    td <- dplyr::ungroup(td)
  }

  td <- td %>% 
    dplyr::select(pep_prob, pep_isdecoy) %>% 
    dplyr::arrange(pep_prob) %>% 
    dplyr::mutate(total = row_number()) %>% 
    dplyr::mutate(decoy = cumsum(pep_isdecoy)) %>% 
    dplyr::mutate(fdr = decoy/total) %>% 
    dplyr::mutate(pep_score = -log10(pep_prob) * 10)
  
  # ---
  # len 7:19
  # ans_param <- c(46.0, 42.8, 21.5, 20.3, 16.0, 19.25, 20.5, 22.1, 24.7, 
  #                26.4, 29.8, 33.3, 38.5)
  # ans_nonparam <- c(45.6, 37.3, 26.2, 21.4, 16.6, 19.5, 21.6, 22.4, 24.8, 
  #                   27.0, 31.1, 33.4, 37.0)

  # ---
  count <- nrow(td)
  
  if (count < (1 / target_fdr)) {
    # return(NA)
    if (count <= 20L) 
      return(NA)
    
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
    row <- max(rows, na.rm = TRUE)
    score_co <- td$pep_score[row]
    
    # fittings (the data range may affect the fitting)
    df <- data.frame(x = td[["pep_score"]], y = td[["fdr"]])

    fit <- suppressWarnings(
      tryCatch(
        nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
            control = list(tol = 1e-03, warnOnly = TRUE), 
            algorithm = "port"), 
        error = function (e) NA)
    )
    
    # p <- suppressWarnings(
    #   ggplot(df, aes(x = x, y = y)) + 
    #     geom_point() + 
    #     stat_smooth(method = "nls", formula = y ~ SSlogis(x, Asym, xmid, scal), 
    #                 se = FALSE) + 
    #     geom_hline(yintercept = target_fdr, size = .5) + 
    #     scale_x_continuous(minor_breaks = seq(1, 100, 10)) + 
    #     labs(title = paste("pep_len = ", len), x = "Score", y = "Prob")
    # )
    # 
    # try(suppressMessages(
    #   ggsave(filename = file.path(out_path, "temp", paste0("peplen_sco_", len, ".png")))
    # ))
    
    if (!all(is.na(fit))) {
      newx <- min(df$x, na.rm = TRUE):max(df$x, na.rm = TRUE)
      newy <- predict(fit, data.frame(x = newx)) %>% `names<-`(newx)
      
      # NA if not existed
      score_co2 <- which(newy <= target_fdr)[1] %>% names() %>% as.numeric()
      best_co <- min(score_co, score_co2, na.rm = TRUE)
      
      rm(list = c("newx", "newy", "score_co", "score_co2"))
    } else {
      best_co <- score_co
    }
    
    rm(list = c("df", "fit"))
    
    prob_co <- 10^(-best_co/10)
    
  } else {
    best_co <- tryCatch(
      (find_pepscore_co1(td, target_fdr) + find_pepscore_co2(td, target_fdr))/2,
      error = function(e) NA
    )
    
    prob_co <- 10^(-best_co/10)
    # prob_co <- NA
  }
  
  names(prob_co) <- count

  invisible(prob_co)
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
    find_optlens(all_lens, min_count/2L)
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


#' Calculates the cut-off score at a peptide FDR.
#'
#' Needs \code{min_len} and \code{max_len} since the target-decoy pair may not
#' cover all \code{pep_len} values.
#'
#' @param out The output from \link{calc_pepscores}.
#' @param nms The name(s) of \code{out} that correspond(s) to decoy results.
#' @param target_fdr Numeric; the levels of false-discovery rate (FDR).
#' @param fdr_type Character string; the type of FDR for controlling.
#' @inheritParams matchMS
#' @examples 
#' \donttest{
#' out_path <- "~/proteoM/bi_1"
#' 
#' # Reload decoy
#' rev <- local({
#'   files <- list.files(path = file.path(out_path, "temp"), 
#'                       pattern = "^pepscores_rev_\\d+\\.rds$")
#'   
#'   if (!length(files)) stop("Decoy scores not found.", call. = FALSE)
#'   
#'   nms <- gsub("^pepscores_rev_(\\d+)\\.rds$", "\\1", files)
#'   
#'   ord <- order(as.integer(nms))
#'   nms <- nms[ord]
#'   files <- files[ord]
#'   
#'   out <- lapply(files, function (x) readRDS(file.path(out_path, "temp", x)))
#'   names(out) <- paste0("rev_", nms)
#'   
#'   out
#' })
#' 
#' 
#' # Reload targets and concatenate decoy
#' out <- local({
#'   files <- list.files(path = file.path(out_path, "temp"), 
#'                       pattern = "^pepscores_\\d+\\.rds$")
#'   
#'   if (!length(files)) stop("Target scores not found.", call. = FALSE)
#'   
#'   nms <- gsub("^pepscores_(\\d+)\\.rds$", "\\1", files)
#'   
#'   ord <- order(as.integer(nms))
#'   nms <- nms[ord]
#'   files <- files[ord]
#'   
#'   out <- lapply(files, function (x) readRDS(file.path(out_path, "temp", x)))
#'   names(out) <- nms
#'   
#'   out <- c(out, rev)
#' })
#' 
#' prob_cos <- calc_pepfdr(out, 
#'                         nms = "rev_1", 
#'                         target_fdr = .01, 
#'                         fdr_type = "protein", 
#'                         min_len = 7L, 
#'                         max_len = 50L, 
#'                         out_path = "~/proteoM/bi_1")
#' 
#' }
calc_pepfdr <- function (out, nms, target_fdr = .01, fdr_type = "psm", 
                         min_len = 7L, max_len = 50L, out_path) 
{
  if (!is.null(nms)) {
    nms_t <- gsub("^rev_", "", nms)
    
    # keeps the best hit for each `scan_num`
    # (separated bests for targets and decoys)
    nrow_t <- nrow(out[[nms_t]])
    
    if (!nrow_t) {
      lens <- min_len:max_len
      prob_cos <- rep(.5, length(lens))
      names(prob_cos) <- lens
      
      return(prob_cos)
    }
    
    out[c(nms_t, nms)] <- lapply(out[c(nms_t, nms)], function (x) {
      x %>% 
        dplyr::group_by(scan_num, raw_file) %>% 
        dplyr::arrange(pep_prob) %>% 
        dplyr::filter(row_number() == 1L) %>% 
        dplyr::ungroup()
    })
    
    # decoy sequences may be present in targets
    out[[nms]] <- purge_decoys(target = out[[nms_t]], decoy = out[[nms]])
    
    #  keeps the best hit for each `scan_num`
    td <- dplyr::bind_rows(out[c(nms_t, nms)]) %>% 
      dplyr::group_by(scan_num, raw_file) %>% 
      dplyr::arrange(pep_prob) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
    
    rm(list = c("out", "nrow_t"))
    gc()

    # ---
    all_lens <- sort(unique(td$pep_len))
    
    prob_cos <- lapply(all_lens, probco_bypeplen, 
                       td, fdr_type, target_fdr, out_path)
    prob_cos <- unlist(prob_cos)

    if (length(prob_cos) == 1L && !is.na(prob_cos)) {
      names(prob_cos) <- all_lens
      return(prob_cos)
    }

    if (all(is.na(prob_cos))) {
      newx <- min_len:max_len
      prob_cos <- rep(target_fdr, length(newx))
      names(prob_cos) <- newx

      return(prob_cos)
    }
    
    counts <- as.numeric(names(prob_cos))
    names(counts) <- all_lens
    names(prob_cos) <- all_lens
    prob_cos <- prob_cos[!is.na(prob_cos)]
    
    # pdf(file.path(out_path, "pepprob_len.pdf")) 
    # plot(prob_cos ~ names(prob_cos), xlab = "pep_len", ylab = "prob_co")
    # dev.off()
    
    lens <- find_optlens(all_lens, counts, 128L)
    
    # no fittings
    if (length(lens) <= 3L) return(prob_cos)
    
    # quality ones for fittings
    prob_cos <- prob_cos %>% .[names(.) %in% lens]
    counts <- counts %>% .[names(.) %in% lens]
    
    valley <- find_probco_valley(prob_cos)
    best_co <- -log10(prob_cos[as.character(valley)])
    
    # low cut-offs at high pep_len may be by chance
    prob_cos <- local({
      len <- length(prob_cos)
      prs <- prob_cos[valley:len]
      prs[prs > .02] <- NA
      prob_cos[valley:len] <- prs
      
      prob_cos[!is.na(prob_cos)]
    })
    
    ## fittings
    df <- data.frame(x = as.numeric(names(prob_cos)), y = -log10(prob_cos))

    # valley left
    df_left <- df[df$x <= valley, ]
    rank_left <- 4L
    
    fit_left <- if (nrow(df_left) <= rank_left) {
      lm(y ~ x, df_left)
    }
    else {
      local({
        fit_ns <- tryCatch(
          lm(y ~ splines::ns(x, rank_left), df_left),
          error = function(e) NA
        )
        
        fit_bs <- tryCatch(
          lm(y ~ splines::bs(x, rank_left), df_left),
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

    # valley right (small rank to down-weight wiggly high `pep_len` points)
    df_right <- df[df$x > valley, ]
    rank_right <- 2L
    
    if (nrow(df_right) > rank_right) {
      fit_right <- lm(y ~ splines::ns(x, rank_right), df_right)

      newx_right <- min(df_right$x, na.rm = TRUE):max(df_right$x, na.rm = TRUE)
      newy_right <- predict(fit_right, data.frame(x = newx_right))
      names(newy_right) <- newx_right
    } 
    else {
      slope <- .28
      newx_right <- valley:max(df_right$x, na.rm = TRUE)
      newy_right <- best_co + slope * (newx_right - valley)
      names(newy_right) <- newx_right

      # excludes the `valley` itself (already in left fitting)
      newx_right <- newx_right[-1]
      newy_right <- newy_right[-1]
    }
    
    # left + right
    newx <- c(newx_left, newx_right)
    newy <- c(newy_left, newy_right)
    
    local({
      n_row <- nrow(df)
      
      df_new <- rbind2(
        cbind(df, type = rep("Original", n_row)), 
        data.frame(x = newx, y = newy, type = "Fitted")
      )

      n_row2 <- nrow(df_new) - n_row
      
      # `*10` only for plots
      df_new$y <- df_new$y * 10 

      pdf(file.path(out_path, "pepscore_len.pdf")) 
      plot(y ~ x, df_new, col = c(rep("blue", n_row), rep("red", n_row2)), 
           xlab = "pep_len", ylab = "score_co", pch = 19)
      legend("bottomright", legend = c("Original", "Fitted"), 
             col = c("blue", "red"), pch = 19, bty = "n")
      dev.off()
      
      # ggplot2::ggplot(df, aes(x = x, y = y)) + geom_point() +
      #   stat_smooth(method = "lm", formula = y ~ splines::bs(x, rank_left), se = FALSE)
    })

    prob_cos <- 10^-newy
    
  } else {
    seqs <- min_len : max_len
    prob_cos <- rep(.05, length(seqs))
    names(prob_cos) <- seqs
  }
  
  invisible(prob_cos)
}


#' Calculates the cut-offs of protein scores.
#'
#' @param out An output from upstream steps.
#' @inheritParams calc_pepfdr
calc_protfdr <- function (out, target_fdr = .01) 
{
  message("Calculating peptide-protein FDR.")
  
  # target-decoy pair
  nms_d <- unique(out$pep_mod_group) %>% 
    .[grepl("^rev_\\d+", .)]
  
  nms_t <- gsub("^rev_", "", nms_d)
  
  out2 <- dplyr::filter(out, pep_mod_group %in% c(nms_t, nms_d), pep_issig)

  # score cut-offs as a function of prot_n_pep
  max_n_pep <- max(out$prot_n_pep, na.rm = TRUE)
  all_n_peps <- unique(out$prot_n_pep)

  # protein enrichment score cut-offs at each `prot_n_pep`
  score_co <- out2 %>% 
    split(.$prot_n_pep) %>% 
    purrr::map_dbl(calc_protfdr_i, target_fdr) 

  # fitted score cut-offs
  score_co <- score_co %>% 
    fit_protfdr(max_n_pep) %>% 
    dplyr::filter(prot_n_pep %in% all_n_peps) %>% 
    dplyr::rename(prot_es_co = prot_score_co)
  
  # add protein enrichment score
  prot_es <- out %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(pep_issig) %>% 
    dplyr::mutate(pep_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_es = max(pep_es, na.rm = TRUE))
  
  # puts together
  out <- out %>% 
    dplyr::left_join(prot_es, by = "prot_acc")

  out <- out %>% 
    dplyr::left_join(score_co, by = "prot_n_pep") %>% 
    dplyr::mutate(prot_issig = ifelse(prot_es >= prot_es_co, TRUE, FALSE)) %>% 
    dplyr::mutate(pep_score = round(pep_score, digits = 1), 
                  pep_score_co = round(pep_score_co, digits = 1), 
                  prot_es = round(prot_es, digits = 1), 
                  prot_es_co = round(prot_es_co, digits = 1))
  
  rm(list = c("out2", "prot_es"))
  gc()

  invisible(out)
}


#' Helper of \link{calc_protfdr}.
#' 
#' For prot_n_pep at value i.
#' 
#' @param td A data frame with paired target-decoys at prot_n_pep = i.
#' @inheritParams calc_pepfdr
#' @return A score cut-off at a given prot_n_pep.
calc_protfdr_i <- function (td, target_fdr = .01) 
{
  options(digits = 9L)
  
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
    if (nrow(td) <= 5L) 
      return(0L) 
    else 
      return(20L)
  }

  ## both targets and decoys
  if (nrow(td) <= 20L) 
    return(1L)
  
  # (NOT to use `local` for hard `return (0L)`)
  prot_scores <- td %>% 
    dplyr::mutate(prot_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_es = max(prot_es, na.rm = TRUE))
  
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
    dplyr::filter(row_number() > 3L)

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
  } else {
    score_co <- td$prot_es[row]

    score_co2 <- local({
      df <- data.frame(x = td[["prot_es"]], y = td[["fdr"]])
      
      fit <- suppressWarnings(
        tryCatch(
          nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      )
      
      # ggplot(df, aes(x = x, y = y)) + geom_point() + 
      #  stat_smooth(method = "nls", formula = y ~ SSlogis(x, Asym, xmid, scal), 
      #  se = FALSE)

      if (all(is.na(fit))) {
        score_co2 <- score_co
      } else {
        min_score <- min(df$x, na.rm = TRUE)
        max_score <- max(df$x, na.rm = TRUE)
        newx <- min_score:max_score
        newy <- predict(fit, data.frame(x = newx)) %>% 
          `names<-`(newx)
        
        # NA if not existed
        score_co2 <- which(newy <= target_fdr)[1] %>% 
          names() %>% 
          as.numeric()
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
fit_protfdr <- function (vec, max_n_pep = 1000L) 
{
  if (length(vec) <= 10L) 
    return(data.frame(prot_n_pep = as.numeric(names(vec)), prot_score_co = vec))

  rv <- rev(vec)
  df <- data.frame(x = as.numeric(names(rv)), y = rv)
  elbow <- min(df[which(df$y == min(df$y, na.rm = TRUE)), "x"], na.rm = TRUE)
  amp <- max(df$y, na.rm = TRUE) * .8
  sca <- 0.5
  
  f <- function (x, m = 0, s = 1, a = 1) { a - a / (1 + exp(-(x-m)/s)) }
  
  fit <- suppressWarnings(
    tryCatch(
      nls(y ~ f(x, m, s, a), data = df, 
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
          nls(y ~ f(x, m, s, a), data = df, 
              start = list(a = amp, m = .x, s = sca), 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      })
    )
    
    fit <- if (all(is.na(fits))) {
      NA
    } else {
      fits %>% 
        .[!is.na(.)] %>% 
        .[[length(.)]]
    }
  }
  
  # ggplot(df, aes(x = x, y = y)) + 
  #   geom_point() + 
  #   stat_smooth(method='nls',formula = y ~ f(x, m, s, a), se = FALSE, 
  #               method.args = list(algorithm = "port", 
  #                                  start = c(m = elbow, s = .5,a = amp))) + 
  #   labs(title = "Protein score cut-offs", x = "prot_n_pep", y = "Score")
  
  
  # --- Outputs
  newx <- seq(1, max_n_pep, by = 1)
  
  out <- data.frame(
    prot_n_pep = newx, 
    prot_score_co = predict(fit, data.frame(x = newx))) %>% 
    dplyr::mutate(prot_score_co = ifelse(prot_n_pep >= elbow, 0, prot_score_co))

  invisible(out)
}


#' Finds the the outer products for a vector of MS2 ions at a given ion series.
#'
#' A theoretical peptide at a given position permutation.
#'
#' @param expts Numeric vector; one series experimental MS2s.
#' @param theos Numeric vector; one series of theoretical MS2s.
#' @inheritParams matchMS
#' @examples
#' \donttest{
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
#' find_ppm_outer_bycombi(theos, expts)
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
#' find_ppm_outer_bycombi(theos, expts)
#' }
find_ppm_outer_bycombi <- function (theos, expts, ppm_ms2 = 25L) 
{
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)
  
  e1 <- expts[row_cols[, 2]]
  names(e1) <- row_cols[, 1]
  
  len <- length(theos)
  
  es <- rep(NA, len)
  names(es) <- seq_len(len)
  es[names(e1)] <- e1
  
  # the first half are b-ions and the second half are y-ions
  list(theo = theos, expt = es)
}


#' Helper of \link{calc_peploc}.
#'
#' Not yet currently used.
#'
#' @inheritParams psmC2Q
#' @importFrom magrittr %>% %T>%
try_calc_peploc <- function (out = NULL) 
{
  out <- tryCatch(
    calc_peploc(out),
    error = function(e) NA
  )
  
  if (is.na(out)) {
    message("Retry with a new R session: \n\n",
            "proteoM:::calc_peploc(\n",
            "  out = \"", file.path(out_path, "scores.rds"), "\" \n",
            ")")
    
    fileConn <- file(file.path("~/calc_peploc.R"))
    
    lines <- c(
      "library(proteoM)\n",
      "proteoM:::calc_peploc(",
      paste0("  out = readRDS(", "\"", file.path(out_path, "scores.rds"), "\"", ")"),
      ")\n",
      "unlink(\"~/calc_peploc.R\")"
    )
    
    writeLines(lines, fileConn)
    close(fileConn)
    
    rstudioapi::restartSession(command='source("~/calc_peploc.R")')
  }
  
  invisible(out)
}


#' Calculates the delta scores of `pep_seq`.
#'
#' A score delta between the best and the second best.
#'
#' There should not be any duplicated rows by the combination of
#' c("pep_isdecoy", "scan_num", "raw_file", "pep_seq", "pep_ivmod")
#'
#' @param x The result from \link{calc_pepscores}.
#' @rawNamespace import(data.table, except = c(last, first, between, transpose,
#'   melt, dcast))
calc_peploc <- function (x) 
{
  ## memory inefficient
  # x <- x %>%
  #   dplyr::mutate(pep_score = round(pep_score, 2)) %>%
  #   dplyr::filter(!duplicated(.[, c("pep_isdecoy", "scan_num", "raw_file",
  #                                   "pep_seq", "pep_score")])) %>%
  #   dplyr::group_by(pep_isdecoy, scan_num, raw_file) %>%
  #   # does addl' row ordering by pep_isdecoy, scan_num, raw_file
  #   dplyr::arrange(-pep_score, .by_group = TRUE) %>%
  #   dplyr::mutate(pep_rank = row_number()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::filter(pep_rank <= 3L)
  #
  # gc()
  
  ## memory more efficient
  x <- data.table::data.table(x)
  gc()
  
  ## keeps the best for ties in each pep_seq (0000500, 0005000) -> top.3
  ## (`pep_score` as a surrogate for tied `pep_ivmod`)
  
  # x[ , "pep_score" := round(pep_score, 2)]
  # x <- x[!duplicated(x[, c("pep_isdecoy", "scan_num", "raw_file", "pep_seq", "pep_score")]), ]
  # x[order(-pep_score), pep_rank := seq_len(.N), by = list(pep_isdecoy, scan_num, raw_file)]
  # x <- x[x[, pep_rank <= 3L], ]
  #
  # gc()
  
  # uniq_id --- across different mod locations under the same *pep_seq*
  # uniq_id2 --- the same *uniq_id (~ pep_seq_mod)* but different NLs
  # uniq_id3 --- the same *query* across pep_seqs
  # 
  # pep_rank2 --- within NLs by uniq_id + pep_ivmod2
  # pep_rank --- across uniq_ids
  # pep_rank at output --- the query across pep_seqs
  
  # separates the best score across different NLs under the same `uniq_id2`
  # (`x0` corresponds to the best NL under each `uniq_id2`)
  x[, uniq_id := paste(pep_isdecoy, scan_num, raw_file, pep_seq, sep = ".")]
  x[, "pep_ivmod2" := gsub(" [\\(\\[]\\d+[\\)\\[]$", "", pep_ivmod)]
  x[, uniq_id2 := paste(uniq_id, pep_ivmod2, sep = ".")]
  x[order(-pep_score), pep_rank2 := seq_len(.N), by = list(uniq_id2)]
  
  x0 <- x[x$pep_rank2 == 1L, ] 
  y0 <- x[x$pep_rank2 > 1L, ]
  
  x0 <- x0[, -c("pep_rank2")]
  y0 <- y0[, -c("pep_rank2")]
  
  rm(list = c("x"))
  gc()
  
  # keep the top-3 mod LOCATIONS under each `uniq_id`
  # (the pep_rank here if after the "collapse" of NLs;
  # net effect: best NLs -> top-3 mod)
  topn_mods <- 3L
  x0[order(-pep_score), pep_rank := seq_len(.N), by = list(uniq_id)]
  x0 <- x0[pep_rank <= topn_mods, ]
  
  # `pep_locprob`
  x0[, "sscore" := sum(pep_score, na.rm = TRUE), by = list(uniq_id)]
  x0[, "pep_locprob" := (pep_score/sscore)]
  x0 <- x0[, -c("sscore")]
  gc()
  
  # `pep_locdiff`
  x1 <- x0[pep_rank == 1L, ]
  x2 <- x0[pep_rank == 2L, ]
  gc()
  
  delta <- dplyr::left_join(x1[, c("uniq_id", "pep_locprob")],
                            x2[, c("uniq_id", "pep_locprob")],
                            by = "uniq_id")
  
  delta$pep_locprob.y <- ifelse(is.na(delta$pep_locprob.y), 0, delta$pep_locprob.y)
  
  delta[["pep_locdiff"]] <- delta[["pep_locprob.x"]] - delta[["pep_locprob.y"]]
  delta <- delta[, c("uniq_id", "pep_locdiff")]
  
  rm(list = c("x1", "x2"))
  gc()
  
  # joins `delta`
  x0 <- dplyr::left_join(x0, delta, by = "uniq_id")
  
  rm(list = c("delta"))
  gc()
  
  # adds back `y0`
  # 
  # 1. some NA pep_ranks in y0 (inferior ones at ranks > 3 and no matches in x0);
  # 2. NLs at <= pep_score under the same uniq_id2 (y0) share the same best pep_rank
  #   (share the best pep_score across different NLs for the same pep_seq_mod);
  # 3. the NL entries will be used for proteoQ quantitation at `exprs(pep_rank == 1)`, 
  #   and can be excluded by `exprs(pep_rank_nl == 1)`.

  y0 <- dplyr::left_join(y0, x0[, c("uniq_id2", "pep_rank", "pep_locprob", 
                                    "pep_locdiff")], by = "uniq_id2")
  y0 <- y0[!is.na(y0$pep_rank), ]
  
  x0 <- data.table::rbindlist(list(x0, y0), use.names = FALSE)
  x0[, pep_rank_nl := seq_len(.N), by = "uniq_id2"]

  rm(list = c("y0"))
  
  # assumes NLS under the same c("uniq_id2", "pep_score") contain 
  #   no additional information and excluded
  x0[ , "pep_score" := round(pep_score, 2)]
  x0 <- unique(x0, by = c("uniq_id2", "pep_score"))
  
  # cleanups
  x0 <- x0[, -c("uniq_id", "uniq_id2", "pep_ivmod2")]
  x0[ , "pep_locprob" := round(pep_locprob, 2)]
  x0[ , "pep_locdiff" := round(pep_locdiff, 2)]
  
  gc()
  
  # new `pep_rank`s across different `pep_seq`s
  # (this is different to the earlier uniq_id at different LOCATIONS)
  x0[, uniq_id3 := paste(pep_isdecoy, scan_num, raw_file, sep = ".")]
  x0[order(-pep_score), pep_rank := seq_len(.N), by = list(uniq_id3)]
  topn_seqs <- 3L
  x0 <- x0[pep_rank <= topn_seqs, ]
  x0 <- x0[, -c("uniq_id3")]

  invisible(x0)
}


