#' Reporter-ion quantitation.
#'
#' @param data An upstream result from \link{matchMS}.
#' @param quant A quantitation method. The default is "none". Additional choices
#'   include \code{tmt6} etc. For other multiplicities of \code{tmt}, use the
#'   compatible higher plexes, for example, \code{tmt16} for \code{tmt12} etc.
#'   and \code{tmt10} for \code{tmt8} etc.
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
calc_tmtint <- function (data = NULL,
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"),
                         ppm_reporters = 10) 
{
  if (quant == "none") {
    out <- data
  } 
  else {
    nms_tmt6 <- c("126", "127N", "128N", "129N", "130N", "131N")
    
    nms_tmt10 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                   "130N", "130C", "131N")
    
    nms_tmt11 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                   "130N", "130C", "131N", "131C")
    
    nms_tmtpro <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                    "130N", "130C", "131N", "131C", "132N", "132C",
                    "133N", "133C", "134N", "134C", "135N")
    
    tmts <- c(
      `126` = 126.127726, `127N` = 127.124761, `127C` = 127.131080,
      `128N` = 128.128115, `128C` = 128.134435, `129N` = 129.131470,
      `129C` = 129.137790, `130N` = 130.134825, `130C` = 130.141145,
      `131N` = 131.138180, `131C` = 131.144499, `132N` = 132.141535,
      `132C` = 132.147855, `133N` = 133.14489, `133C` = 133.15121,
      `134N` = 134.148245, `134C` = 134.154565, `135N` = 135.15160)
    
    theos <- switch(quant,
                    tmt6 = tmts %>% .[names(.) %in% nms_tmt6],
                    tmt10 = tmts %>% .[names(.) %in% nms_tmt10],
                    tmt11 = tmts %>% .[names(.) %in% nms_tmt11],
                    tmt16 = tmts %>% .[names(.) %in% nms_tmtpro],
                    tmt18 = tmts %>% .[names(.) %in% nms_tmtpro],
                    stop("Unknown TMt type.", call. = FALSE))
    
    ul <- switch(quant,
                 tmt6 = c(126.1, 131.2),
                 tmt10 = c(126.1, 131.2),
                 tmt11 = c(126.1, 131.2),
                 tmt16 = c(126.1, 134.2),
                 tmt18 = c(126.1, 135.2),
                 stop("Unknown TMt type.", call. = FALSE))
    
    # stopifnot(all(c("ms2_moverz", "ms2_int") %in% names(data)))
    
    out <- purrr::map2(data$ms2_moverz, data$ms2_int,
                       find_reporter_ints,
                       theos = theos,
                       ul = ul,
                       ppm_reporters = ppm_reporters,
                       len = length(theos),
                       nms = names(theos)) %>%
      dplyr::bind_rows() 
    
    if (!nrow(out)) {
      out <- data.frame(matrix(ncol = length(theos), nrow = 0L))
      colnames(out) <- theos
      
      for (i in seq_along(out)) 
        out[[i]] <- as.numeric(out[[i]])
    }

    out <- dplyr::bind_cols(data, out)
  }
  
  names(out)[grep("^([0-9]{3}[NC]{0,1})", names(out))] <-
    find_int_cols(length(theos))
  
  invisible(out)
}


#' Adds back reporter-ion intensities
#'
#' @param df Results from \link{calc_protfdr}.
#' @inheritParams matchMS
add_rptrs <- function (df = NULL, quant = "none", out_path = NULL) 
{
  if (grepl("^tmt[0-9]+$", quant)) {
    files <- list.files(path = file.path(out_path, "temp"),
                        pattern = "^reporters_\\d+\\.rds$")
    
    idxes <- gsub("^reporters_([0-9]+)\\.rds$", "\\1", files) %>%
      as.integer() %>%
      order()
    
    files <- files[idxes]
    
    reporters <- lapply(files, function (x) 
      readRDS(file.path(out_path, "temp", x))) %>%
      dplyr::bind_rows()

    rm(list = c("idxes", "files"))
    
    df <- df %>%
      tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".",
                   remove = FALSE) %>%
      dplyr::left_join(reporters, by = "uniq_id") %>%
      dplyr::select(-uniq_id)
  }
  
  invisible(df)
}


#' Finds the intensities of reporter-ions.
#'
#' @param ms2_moverzs Numeric vector; a series of experimental MS2 m-over-z's
#'   (in the region of reporter ions).
#' @param ms2_ints Numeric vector; a series of experimental MS2 intensities (in
#'   the region of reporter ions).
#' @param theos The theoretical m-over-z of reporter ions.
#' @param ul The upper and lower bound for reporter-ion m-over-z's.
#' @param len The length of reporter-ion plexes.
#' @param nms The names of reporter-ion channels.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' ms2_moverzs <- c(112.0873, 126.1280, 127.1251, 127.1313, 128.1250,
#'                  128.1284, 128.1347, 129.1317, 129.1380, 130.0654,
#'                  130.1351, 130.1413, 131.1384)
#'
#' ms2_ints <- c(5113.79, 135569.00, 120048.00, 122599.00, 3397.98,
#'               140551.00, 144712.00, 103166.00, 145452.00, 3851.82,
#'               148218.00, 135393.00, 131215.00)
#'
#' theos <- c(126.1277, 127.1248, 127.1311, 128.1281, 128.1344,
#'            129.1315, 129.1378, 130.1348, 130.1411, 131.1382)
#' names(theos) <- c("126", "127N", "127C", "128N", "128C",
#'                   "129N", "129C", "130N", "130C", "131N")
#'
#' ppm_reporters <- 10
#' ul <- c(126.1, 131.2)
#' len <- 10
#' nms <- names(theos)
#'
#' x <- find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 10,
#'                         len , nms)
#'
#' x <- find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 25,
#'                         len , nms)
#'
#' # Two `129C`, no `127N` etc.
#' ms2_moverzs <- c(105.1503, 107.0428, 111.7716, 120.0811, 126.1281, 127.1312,
#'                  128.1282, 128.1349, 129.1317, 129.1365, 129.1382, 230.1694,
#'                  233.4857, 233.4964, 337.3533, 352.1844, 376.2764, 463.3083,
#'                  525.2150, 562.3732, 569.3899, 591.2545, 596.0308, 632.3300,
#'                  636.3959, 703.3637, 789.0423, 816.4487, 817.4516, 839.9531,
#'                  864.3056, 914.7645, 921.5302, 1479.9816)
#'
#' ms2_ints <- c(1201.79, 1319.32, 1603.45, 1595.34, 2148.66, 1785.74, 1254.24,
#'               1986.43, 10127.40, 1522.60, 1562.71, 2926.01, 1590.48, 1692.17,
#'               1347.88, 1412.64, 3050.10, 3231.10, 1355.21, 2424.18, 1783.26,
#'               1365.32, 1727.12, 2661.72, 1660.05, 5525.95, 1399.96, 4654.03,
#'               1990.57, 1758.72, 1655.09, 1460.68, 1641.39, 1721.33)
#'
#' x <- find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 25,
#'                         len , nms)
#' }
find_reporter_ints <- function (ms2_moverzs, ms2_ints, theos, ul,
                                ppm_reporters = 10, len, nms) 
{
  range <- findInterval(ul, ms2_moverzs)
  
  ms <-ms2_moverzs[range[1]:range[2]]
  is <-ms2_ints[range[1]:range[2]]
  
  idxes <- find_reporters_ppm(theos, ms, ppm_reporters, len, nms)
  
  if (!length(idxes)) 
    return(rep(NA, len) %>% `names<-`(nms))
  
  # 126      127N      127C      128N      128N      128C
  # 135569.00 120048.00 122599.00   3397.98 140551.00 144712.00
  
  if (anyDuplicated(names(idxes))) {
    idxes <- idxes %>%
      split(., names(.)) %>%
      purrr::imap_int(~ {
        if (length(.x) > 1L) {
          p <- which.min(abs(ms[.x] - theos[.y]))
          .x <- .x[p]
        }
        
        .x
      }) %>%
      .[nms]
  }
  
  # missing channels:
  # 126 <NA> 127C 128N 128C 129N 129C <NA> <NA> <NA>
  #  2   NA    3    4    5    6    8   NA   NA   NA
  
  if (anyNA(names(idxes))) names(idxes) <- nms
  
  rptr_ints <- is[idxes] %>%
    `names<-`(names(idxes))
  
  if (length(rptr_ints) < len) {
    es <- rep(NA, len)
    names(es) <- nms
    es[names(rptr_ints)] <- rptr_ints
  } 
  else {
    es <- rptr_ints
  }
  
  es
}


#' Finds the indexes of reporter ions.
#'
#' @param expts Numeric vector; a series of experimental MS2s (in the region of
#'   reporter ions).
#' @inheritParams find_reporter_ints
#' @return A vector of indexes
find_reporters_ppm <- function (theos, expts, ppm_reporters = 10, len, nms) 
{
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_reporters, arr.ind = TRUE)
  
  row_cols[, 2]
}


#' Adds prot_acc to a peptide table
#'
#' @param out_path An output path.
#' @param df The results after scoring.
add_prot_acc <- function (df, out_path = "~/proteoM/outs") 
{
  message("Adding protein accessions.")
  
  uniq_peps <- unique(df$pep_seq)
  
  # Targets, theoretical
  .path_ms1masses <- get(".path_ms1masses", envir = .GlobalEnv, inherits = FALSE)
  
  fwd_prps <- readRDS(file.path(.path_ms1masses, .time_stamp, "prot_pep_annots.rds"))
  fwd_prps <- fwd_prps[fwd_prps$pep_seq %in% uniq_peps, ]
  
  # Decoys, theoretical
  rev_prps <- readRDS(file.path(.path_ms1masses, .time_stamp, "prot_pep_annots_rev.rds"))
  rev_prps <- rev_prps[rev_prps$pep_seq %in% uniq_peps, ]
  rev_prps <- purge_decoys(target = fwd_prps, decoy = rev_prps)
  
  # Adds `prot_acc` (with decoys being kept)
  rm(list = c("uniq_peps"))
  gc()
  
  hadd_prot_acc(df, fwd_prps, rev_prps)
}


#' Adds prot_acc to a peptide table
#'
#' Cached results are under sub dirs.
#' 
#' @param out_path An output path.
#' @param df The results after scoring.
add_prot_acc2 <- function (df, out_path) 
{
  message("Adding protein accessions.")
  
  sub_dirs <- dir(out_path, pattern = "^sub[0-9]+_[0-9]_[0-9]+$", full.names = TRUE)
  len_dirs <- length(sub_dirs)
  
  if (!len_dirs) {
    warnning("Cached results not found for protein annotation.")
    return (df)
  }

  uniq_peps <- unique(df$pep_seq)
  
  # Targets, theoretical
  fwd_prps <- lapply(sub_dirs, function (x) {
    file <- file.path(x, "prot_pep_annots.rds")
    fwd <- if (file.exists(file)) readRDS(file) else NULL
  }) %>% 
    dplyr::bind_rows()
  
  fwd_prps <- fwd_prps[fwd_prps$pep_seq %in% uniq_peps, ]
  
  # Decoys, theoretical
  rev_prps <- lapply(sub_dirs, function (x) {
    file <- file.path(x, "prot_pep_annots_rev.rds")
    fwd <- if (file.exists(file)) readRDS(file) else NULL
  }) %>% 
    dplyr::bind_rows()
  
  rev_prps <- rev_prps[rev_prps$pep_seq %in% uniq_peps, ]
  rev_prps <- purge_decoys(target = fwd_prps, decoy = rev_prps)
  
  rm(list = c("uniq_peps"))
  gc()
  
  hadd_prot_acc(df, fwd_prps, rev_prps)
}


#' Helper of \link{add_prot_acc}.
#' 
#' @param df A data frame.
#' @param fwd_prps The look-ups of forward protein and peptides.
#' @param rev_prps The look-ups of reversed protein and peptides.
hadd_prot_acc <- function (df, fwd_prps, rev_prps) 
{
  # Adds `prot_acc` (with decoys being kept)
  out <- dplyr::bind_rows(fwd_prps, rev_prps) %>%
    dplyr::right_join(df, by = "pep_seq")
  
  rm(list = c("fwd_prps", "rev_prps"))
  gc()
  
  # Adds prot_n_psm, prot_n_pep for protein FDR
  x <- out[out$pep_issig, ]
  
  prot_n_psm <- x %>%
    dplyr::select(prot_acc) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_psm = n())
  
  prot_n_pep <- x %>%
    dplyr::select(pep_seq, prot_acc) %>% 
    tidyr::unite(pep_prot, c("pep_seq", "prot_acc"), sep = ".", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_prot)) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_pep = n())
  
  # inconsistent Protein[NC]-term
  out <- local({
    pnt_nots <- grepl("Protein N-term", out$pep_vmod) & !out$is_pnt
    pct_nots <- grepl("Protein C-term", out$pep_vmod) & !out$is_pct
    
    if (sum(pnt_nots) > 0L) 
      out <- out[!pnt_nots, ]
    
    if (sum(pct_nots) > 0L) 
      out <- out[!pct_nots, ]
    
    # `is_pnt` and `is_pct` are not EXACT facts
    # but prefer terminal over interiror matches
    # so remove them to avoid misleading uses or interpretations
    out$is_pnt <- NULL
    out$is_pct <- NULL
    
    out
  })

  out <- list(out, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = "prot_acc") %>%
    dplyr::arrange(-prot_n_pep, -prot_n_psm)
  
  invisible(out)
}




#' Helper of \link{groupProts}.
#'
#' @param df A data frame contains proteins and peptides.
#' @param out_path The output path.
grp_prots <- function (df, out_path = NULL) 
{
  if (!nrow(df))
    stop("Zero row of data for protein groupings.")
  
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # Significant (df1) and trivial (df0) entries
  df <- df[with(df, order(pep_seq)), ]
  rows <- (df$pep_issig & (!df$pep_isdecoy) & (df$pep_rank <= 3L))
  
  df1 <- df[rows, ]
  df0 <- df[!rows, ]
  rm(list = c("df"))
  gc()
  
  # (passes the whole `df1` to avoid expensive left_join)
  if (nrow(df1) > 1L) 
    df1 <- groupProts(df1, out_path)
  else {
    df1 <- dplyr::mutate(df1, 
                         prot_isess = TRUE,
                         prot_hit_num = 1L,
                         prot_family_member = 1L, 
                         pep_literal_unique = TRUE, 
                         pep_razor_unique = TRUE)
  }

  # Trivial entries
  ess_prots <- df1 %>%
    dplyr::filter(!duplicated(prot_acc), prot_isess) %>%
    `[[`("prot_acc")

  hits_n_fams <- df1[, c("prot_acc", "prot_hit_num", "prot_family_member")] %>%
    dplyr::filter(!duplicated(prot_acc))
  
  df2 <- df0 %>%
    dplyr::select(prot_acc) %>%
    dplyr::filter(!duplicated(prot_acc)) %>%
    dplyr::mutate(prot_isess = ifelse(prot_acc %in% ess_prots, TRUE, FALSE)) %>%
    dplyr::left_join(hits_n_fams, by = "prot_acc")

  df2 <- df2 %>%
    dplyr::right_join(df0, by = "prot_acc") %>%
    dplyr::mutate(pep_literal_unique = NA, pep_razor_unique = NA) %>%
    dplyr::select(names(df1))
  
  if (nrow(df2))
    df1 <- dplyr::bind_rows(df1, df2)
  
  df1 %>% 
    dplyr::select(-which(names(.) %in% c("prot_n_psm", "prot_n_pep")))
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
#' @param out_name The output filename.
#' 
#' @examples
#' \donttest{
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ]
#'
#' out <- proteoM:::groupProts(df)
#'
#' # One peptide, multiple proteins
#' df <- data.frame(prot_acc = LETTERS[1:3], pep_seq = rep("X", 3))
#' out <- proteoM:::groupProts(df)
#' stopifnot(nrow(out) == 3L)
#'
#' # One peptide, one proteins
#' df <- data.frame(prot_acc = "A", pep_seq = "X")
#' out <- proteoM:::groupProts(df)
#' stopifnot(nrow(out) == 1L)
#'
#' # One proteins
#' df <- data.frame(prot_acc = rep("A", 3), pep_seq = LETTERS[24:26])
#' out <- proteoM:::groupProts(df)
#' stopifnot(nrow(out) == 3L)
#' }
groupProts <- function (df, out_path = NULL, out_name = "prot_pep_setcover.rds") 
{
  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  ## (1) builds protein ~ peptide map
  Mats <- map_pepprot(df[, c("prot_acc", "pep_seq")], out_path)
  
  Mat_upr_left <- Mats$upr_left
  Mat_lwr_left <- Mats$lwr_left
  Mat_lwr_right <- Mats$lwr_right
  
  peps_shared <- rownames(Mat_upr_left)
  prots_upr_right <- colnames(Mat_lwr_right)
  # peps_unique <- rownames(Mat_lwr_left)
  
  if (is.null(Mat_lwr_right))
    Mat_upr_right <- NULL
  else {
    # works for zero-column matrix
    Mat_upr_right <- Matrix::sparseMatrix(
      dims = c(nrow(Mat_upr_left), ncol(Mat_lwr_right)), 
      i={}, j={}
    )
    colnames(Mat_upr_right) <- prots_upr_right
    rownames(Mat_upr_right) <- peps_shared
  }
  
  # Empty `Mat_upr_left` is a zero-row data.frame, not NULL

  rm(list = c("Mats"))
  gc()
  
  # Mat_left (`Mat_upr_left` + `Mat_lwr_left`)
  # - `Mat_upr_left`
  #   * column names: proteins with shared peptides
  #   * row names: shared peptides 
  # - `Mat_lwr_left`
  #   * column names: the same as `Mat_upr_left`
  #   * row names: unique peptides
  # 
  # `Mat_right` (`Mat_upr_right` + `Mat_lwr_right`)
  # - `Mat_upr_right`
  #   * column names: protein with only unique peptides
  #   * row names: shared peptides (as in `Mat_upr_left`)
  #   * (the values are all zeros)
  # - `Mat_lwr_right`
  #   * column names: the same as `Mat_upr_right`
  #   * row names: unique peptides (as in `Mat_lwr_left`)

  ## (2) establishes protein groups
  # each row in cbind(Mat_lwr_left, Mat_lwr_right) has only one "1"
  #   -> Mat_lwr_left does not affect logical distance of 0/1
  
  if (nrow(Mat_upr_left)) {
    prot_grps <- local({
      grps_1 <- cut_proteinGroups(Mat_upr_left, out_path)
      gc()
      
      max <- max(grps_1$prot_hit_num, na.rm = TRUE)
      idxes <- seq_along(prots_upr_right) + max
      
      if (is.null(prots_upr_right) || !length(prots_upr_right))
        grps_2 <- NULL
      else 
        grps_2 <- data.frame(prot_acc = prots_upr_right, 
                             prot_hit_num = idxes, 
                             prot_family_member = 1L)

      rbind2(grps_1, grps_2)
    })
  }
  else {
    # no shared peptides
    prot_grps <- data.frame(prot_acc = prots_upr_right, 
                            prot_hit_num = seq_along(prots_upr_right), 
                            prot_family_member = 1L)
  }

  ## (3) finds essential protein entries
  ess_prots <- local({
    if (is.null(Mat_lwr_left))
      df_shared <- greedysetcover3(Mat_upr_left)
    else {
      rows_lwr_left_is_one <- Matrix::rowSums(Mat_lwr_left) > 0
      Mat_lwr_left_is_one <- Mat_lwr_left[rows_lwr_left_is_one, , drop = FALSE]
      df_shared <- greedysetcover3(rbind(Mat_upr_left, Mat_lwr_left_is_one))
      gc()
    }

    df_uniq <- df[, c("prot_acc", "pep_seq")]
    df_uniq <- df_uniq[df_uniq$prot_acc %in% prots_upr_right, , drop = FALSE]
    df_uniq <- df_uniq[!duplicated.data.frame(df_uniq), , drop = FALSE]
    df_uniq <- df_uniq[with(df_uniq, order(prot_acc, pep_seq)), , drop = FALSE]
    
    sets <- rbind2(df_shared, df_uniq)
    
    if (!is.null(out_path)) 
      saveRDS(sets, file.path(out_path, out_name))
    
    unique(sets$prot_acc)
  })

  ## Parses literal or razor uniqueness of peptides
  # sets aside df0
  df <- dplyr::mutate(df, prot_isess = prot_acc %in% ess_prots)
  df0 <- dplyr::filter(df, !prot_isess)
  df <- dplyr::filter(df, prot_isess)
  gc()
  
  # combines four quadrants
  M4 <- rbind(
    cbind(Mat_upr_left, Mat_upr_right), 
    cbind(Mat_lwr_left, Mat_lwr_right)
  )
  
  rm(list = c("Mat_upr_left", "Mat_upr_right", 
              "Mat_lwr_left", "Mat_lwr_right"))
  gc()
  
  M4_ess <- if (nrow(M4) == 1L) 
    M4
  else 
    M4[, colnames(M4) %in% ess_prots, drop = FALSE]
  
  # literal or razor
  peps_uniq <- local({
    rsums <- Matrix::rowSums(M4)
    rsums2 <- Matrix::rowSums(M4_ess)
    
    peps <- data.frame(pep_seq = rownames(M4)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })
  
  rm(list = c("M4", "M4_ess"))
  gc()
  
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)
  
  df <- df %>%
    dplyr::left_join(prot_grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")

  invisible(df)
}


#' Helper of \link{groupProts}.
#'
#' Builds the logical map between peptide (in rows) and proteins (in columns).
#'
#' The \code{lwr_left} and \code{lwr_right} can be NULL with early exit or
#' "empty" sparse matrix with end return. Maybe uniform later to "empty" matrix.
#'
#' @param df The data frame from upstream steps. It must contains the two
#'   columns of \code{prot_acc} and \code{pep_seq}.
#' @param out_path An output path.
#' @examples
#' \donttest{
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ]
#'
#' out <- proteoM:::map_pepprot(df)
#'
#' # One peptide, multiple proteins
#' df <- data.frame(prot_acc = LETTERS[1:3], pep_seq = rep("X", 3))
#' out <- proteoM:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == "X", colnames(out[[1]]) == LETTERS[1:3])
#'
#' # One peptide, one proteins
#' df <- data.frame(prot_acc = "A", pep_seq = "X")
#' out <- proteoM:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == "X", colnames(out[[1]]) == "A")
#'
#' # One proteins
#' df <- data.frame(prot_acc = rep("A", 3), pep_seq = LETTERS[24:26])
#' out <- proteoM:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == LETTERS[24:26], colnames(out[[1]]) == "A")
#' }
map_pepprot <- function (df, out_path = NULL) 
{
  # collapse rows of the same pep_seq
  #
  #      pep_seq prot_acc
  # 1       A        X
  # 2       A        Y
  # 3       B        X
  # 4       C        Y
  #
  #   X Y
  # 1 1 0
  # 2 0 1
  # 3 1 0
  # 4 0 1
  #
  # pep_seq   X     Y
  # 1 A       TRUE  TRUE
  # 2 B       TRUE  FALSE
  # 3 C       FALSE TRUE
  
  if (!identical(names(df), c("prot_acc", "pep_seq")))
    stop("The two columns of `df` need to be in the order of ", "
         \"prot_acc\" and \"pep_seq\".")
  
  # FIRST ordered by `pep_seq`, SECOND by `prot_acc` 
  # (for continuity of the same `pep_seq`)
  df <- df[!duplicated.data.frame(df), ]
  df <- df[with(df, order(pep_seq, prot_acc)), ]

  peps <- df$pep_seq

  # (one peptide, ONE protein)
  if (length(peps) == 1L) {
    out <- matrix(1)
    colnames(out) <- df$prot_acc
    rownames(out) <- peps
    Mat <- Matrix::Matrix(out, sparse = TRUE)
    
    return(list(upr_left = Mat, lwr_left = NULL, lwr_right = NULL))
  }

  ## Separates into Mat0 and Mat1
  uniq_prots <- unique(df$prot_acc)
  
  # (One protein, multiple peptides)
  if (length(uniq_prots) == 1L) {
    Mat <- Matrix::Matrix(matrix(rep(1L, length(peps))), sparse = TRUE)
    colnames(Mat) <- uniq_prots
    rownames(Mat) <- peps
    
    return(list(upr_left = Mat, lwr_left = NULL, lwr_right = NULL))
  }
  
  Mat <- Matrix::sparse.model.matrix(~ -1 + prot_acc, df)
  prots <- stringi::stri_replace_first_fixed(colnames(Mat), "prot_acc", "")
  colnames(Mat) <- prots
  
  Mat <- Mat == 1L
  rownames(Mat) <- peps
  gc()
  
  dpeps <- peps[duplicated.default(peps)]
  drows <- peps %in% dpeps
  Mat0 <- Mat[!drows, ]
  Mat1 <- Mat[drows, ]

  rm(list = c("dpeps", "drows", "Mat", "peps"))
  gc()
  
  ## Mat1: (a) pre sparse-matrix vector
  ncol <- as.numeric(ncol(Mat1))
  peps1 <- rownames(Mat1)
  vec <- pcollapse_sortpeps(Mat1, ncol, peps1)
  rm(list = c("Mat1"))
  gc()

  ## Mat1: (b) vector -> sparse matrix
  # (not to use function to avoid copying large vector)
  upeps <- unique(peps1)
  n_upeps <- as.numeric(length(upeps))
  llen <- n_upeps * ncol
  
  if (object.size(vec)/1024^3 > 5) {
    rows_per_chunk <- 10000
    size_chunk <- ncol*rows_per_chunk
    n_chunks <- ceiling(n_upeps/rows_per_chunk)

    out <- NULL
    
    for (i in 1:n_chunks) {
      start <- size_chunk*(i-1)+1
      end <- min(llen, size_chunk*i)
      
      vsub <- vec[start:end]
      msub <- Matrix::Matrix(vsub, ncol = ncol, byrow = TRUE, sparse = TRUE)
      out <- rbind2(out, msub)
      
      rm(list = c("vsub", "msub", "start", "end"))
      gc()
    }
    
    rm(list = c("rows_per_chunk", "size_chunk", "n_chunks"))
  } 
  else {
    out <- Matrix::Matrix(vec, ncol = ncol, byrow = TRUE, sparse = TRUE)
  }
  
  rm(list = c("vec"))
  gc()

  colnames(out) <- prots
  rownames(out) <- upeps
  gc()

  ## To logical sparse matrix
  out <- out == 1L
  gc()

  if (FALSE) {
    # grp_prots -> groupProts -> map_pepprot called multiple times
    # need additional "tier" information to prevent overwrites
    if (!is.null(out_path)) {
      Matrix::writeMM(out, file = file.path(out_path, "prot_pep_map.mtx"))
      saveRDS(colnames(out), file.path(out_path, "prot_pep_map_col.rds"))
      saveRDS(rownames(out), file.path(out_path, "prot_pep_map_row.rds"))
    }
  }

  ## Cleans up
  cols_1 <- Matrix::colSums(out) > 0
  upr_left <- out[, cols_1, drop = FALSE]
  lwr_left <- Mat0[, cols_1, drop = FALSE]
  lwr_right <- Mat0[, !cols_1, drop = FALSE]

  invisible(list(upr_left = upr_left, 
                 lwr_left = lwr_left, 
                 lwr_right = lwr_right))
}


#' Helper of \link{map_pepprot}.
#'
#' Collapses the counts the number of peptide sequences under proteins.
#'
#' The row names in the input matrix need to be \emph{sorted}. The output is a
#' vector and will be later wrapped into a sparse matrix.
#'
#' @param mat A dgCMatrix object. Column names are SORTED protein accessions.
#'   Rownames are SORTED peptide sequences.
#' @param ncol The number of columns in \code{mat}.
#' @param peps Peptide sequences as the row names of \code{mat}.
collapse_sortpeps <- function (mat, ncol = NULL, peps = NULL) 
{
  if (is.null(peps)) 
    peps <- rownames(mat)
  
  if (is.null(ncol))
    ncol <- as.numeric(ncol(mat))

  # !!! `peps` must be SORTED !!!

  cts <- cumsum(table(peps))
  llen <- as.numeric(length(cts)) * ncol

  out <- rep.int(0L, llen)
  
  start <- 1
  end <- ncol
  
  r1 <- 1
  r2 <- 0
  
  for (i in seq_along(cts)) {
    r1 <- r2 + 1
    r2 <- cts[i]
    
    out[start:end] <- Matrix::colSums(mat[r1:r2, ])
    
    start <- start + ncol
    end <- end + ncol
    
    # if (i %% 100 == 0) gc()
  }
  
  # rm(list = c("mat"))

  invisible(out)
}


#' Helper of \link{map_pepprot}.
#' 
#' Parallel version of \link{collapse_sortpeps}.
#' 
#' @param n_cores The number of CPU cores.
#' @param Mat A sparse matrix.
#' @inheritParams collapse_sortpeps
pcollapse_sortpeps <- function (Mat, ncol = NULL, peps = NULL, n_cores = NULL) 
{
  if (is.null(peps)) 
    peps <- rownames(Mat)
  
  if (is.null(ncol))
    ncol <- as.numeric(ncol(Mat))
  
  if (is.null(n_cores))
    n_cores <- detect_cores(16L)

  # !!! `peps` must be SORTED !!!
  
  size <- local({
    dim <- dim(Mat)
    as.numeric(dim[1]) * as.numeric(dim[2])
  })

  if (size <= 100000000 || n_cores <= 1L) {
    vec <- collapse_sortpeps(Mat, ncol, peps)
  }
  else {
    Mats <- chunksplit_spmat(Mat, peps, n_cores)
    rm(list = "Mat")
    gc()
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterEvalQ(cl, library(Matrix))
    vecs <- parallel::clusterApply(cl, Mats, collapse_sortpeps)
    parallel::stopCluster(cl)

    vec <- NULL
    
    for (i in seq_along(vecs)) {
      vec <- c(vec, vecs[[i]])
      vecs[i] <- list(NULL)
      gc()
    }
  }
  
  invisible(vec)
}


#' Splits sparse matrix by chunks.
#'
#' The same peptide sequence spans multiple consecutive rows will stay in the
#' same chunk.
#' 
#' @param Mat A sparse matrix.
#' @param peps The names of peptide sequences.
#' @param n_chunks The number of chunks.
chunksplit_spmat <- function (Mat, peps = NULL, n_chunks = 4L) 
{
  if (is.null(peps))
    peps <- rownames(Mat)

  breaks <- find_group_breaks(peps, n_chunks)
  breaks <- c(unname(breaks), length(peps))
  
  Mats <- vector("list", n_chunks)
  
  start <- 1L
  
  for (i in seq_along(Mats)) {
    end <- breaks[i]
    Mats[[i]] <- Mat[start:end, ]
    start <- end + 1
  }
  
  Mats
}


#' Chunksplits by groups.
#' 
#' @param vec A sorted vector.
#' @param n_chunks The number of chunks
#' @examples 
#' \donttest{
#' vec <- rep(LETTERS[1:5], 1:5)
#' vec <- sort(vec)
#' 
#' find_group_breaks(vec, 3)
#' }
find_group_breaks <- function (vec, n_chunks = 5L) 
{
  # !!! vec must be sorted !!!

  if (n_chunks <= 1L)
    return (vec)
  
  tv <- table(vec)
  
  if (n_chunks >= length(tv)) 
    return(split(vec, vec))
  
  len <- length(vec)
  clens <- cumsum(tv)
  labs <- levels(cut(1:len, n_chunks))
  
  x <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
             upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
  
  x1 <- x[, 1]
  
  inds <- lapply(x1[2:length(x1)], function (x) max(which(clens <= x)))
  inds <- unlist(inds)
  
  clens[inds] 
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides.
#'
#' @param M A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file path to outputs.
cut_proteinGroups <- function (M = NULL, out_path = NULL) 
{
  prots <- colnames(M)
  n_prots <- length(prots)
  
  if (ncol(M) == 1L) {
    D <- matrix(1.0)
    colnames(D) <- prots
    rownames(D) <- prots
  } else {
    D <- proxyC::simil(M, margin = 2) # dsTMatrix
  }
  
  rm(list = c("M"))
  gc()
  
  if (n_prots > 10000) {
    dm <- matrix(nrow = n_prots, ncol = n_prots)
    colnames(dm) <- prots
    rownames(dm) <- prots
    gc()
    
    # max_chunksize <- find_free_mem() * 3600000 / 1024
    max_chunksize <- 2500 * 40000
    nrows_per_chunk <- ceiling(max_chunksize/n_prots)
    n_chunks <- ceiling(n_prots/nrows_per_chunk)
    
    cols <- 1:n_prots
    
    for (i in 1:n_chunks) {
      start <- (i-1) * nrows_per_chunk + 1
      end <- min(i * nrows_per_chunk, n_prots)
      rows <- start:end
      
      X <- D[rows, cols] # dsTMatrix
      gc()
      X <- (X == 0) # lsyMatrix
      gc()
      
      x <- as.matrix(X)
      dm[rows, cols] <- x
      rm(list = c("X", "x"))
      gc()
    }
    
    rm(list = c("D", "start", "end", "rows", "max_chunksize", "nrows_per_chunk", 
                "n_chunks", "cols"))
    gc()
  } 
  else {
    D <- (D == 0) # lsyMatrix
    gc()
    
    dm <- as.matrix(D)
    rm(list = c("D"))
    gc()
  }
  
  # Diagonal values are `FALSE`
  # TRUE - orthogonal (without shared peptides)
  # FALSE - with shared peptides
  # 
  #             KKA1_ECOLX NP_000005 NP_000007
  # KKA1_ECOLX      FALSE      TRUE      TRUE
  # NP_000005        TRUE     FALSE      TRUE
  # NP_000007        TRUE      TRUE     FALSE
  
  # --- finds protein groups
  d <- as_lgldist(dm, diag = FALSE, upper = FALSE) # logical distance
  rm(list = "dm")
  gc()
  
  if (length(d)) {
    hc <- hclust(d, method = "single")
    gc()
    
    grps <- data.frame(prot_hit_num = cutree(hc, h = .9))
  } 
  else {
    hc <- NULL
    grps <- data.frame(prot_hit_num = 1L)
    rownames(grps) <- prots
  }
  
  grps <- grps %>% 
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::group_by(prot_hit_num) %>%
    dplyr::mutate(prot_family_member = dplyr::row_number()) %>%
    dplyr::ungroup()

  if (!is.null(out_path)) 
    saveRDS(grps, file.path(out_path, "prot_grps.rds"))
  
  # stopifnot(identical(grps$prot_acc, prots))
  
  invisible(grps)
}


#' Builds manually distance sparse matrix.
#'
#' Not yet used.
#' 
#' @param M_ul The upper-left matrix.
#' @param ncols_ur The number of columns for the matrix block on the upper
#'   right.
#' @examples 
#' m <- sparseD_fourquad(ul, 6)
sparseD_fourquad <- function (M_ul, ncols_ur = 0) 
{
  nrows_ul <- ncols_ul <- ncol(M_ul)
  nrows_l <- ncols_ur
  ncols <- ncols_ul + ncols_ur
  
  M_ur <- sparseMatrix(dims = c(nrows_ul, ncols_ur), i={}, j={})
  m_lwr <- sparseMatrix(i = 1:nrows_l, j = (ncols_ul+1):ncols, x = 1)
  
  list(upr = cbind2(M_ul, M_ur), lwr = m_lwr)
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Not yet used; assumed the input is already a symmetric matrix.
#' 
#' @inheritParams stats::as.dist
as_dist <- function (m, diag = FALSE, upper = FALSE) 
{
  p <- nrow(m)
  
  ans <- m[row(m) > col(m)]
  gc()
  attributes(ans) <- NULL
  
  if (!is.null(rownames(m))) 
    attr(ans, "Labels") <- rownames(m)
  else if (!is.null(colnames(m))) 
    attr(ans, "Labels") <- colnames(m)
  
  attr(ans, "Size") <- p
  attr(ans, "call") <- match.call()
  class(ans) <- "dist"
  
  if (is.null(attr(ans, "Diag")) || !missing(diag)) 
    attr(ans, "Diag") <- diag
    
  if (is.null(attr(ans, "Upper")) || !missing(upper)) 
    attr(ans, "Upper") <- upper

  ans
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Assumed the input is already a symmetric matrix.
#' 
#' @inheritParams stats::as.dist
as_lgldist <- function(m, diag = FALSE, upper = FALSE) 
{
  d = proteoCpp::to_lgldistC(m)
  
  if (!is.null(rownames(m))) 
    attr(d, "Labels") <- rownames(m)
  else if (!is.null(colnames(m))) 
    attr(d, "Labels") <- colnames(m)

  attr(d, "class") = "dist"
  attr(d, "Size") = nrow(m)
  attr(d, "call") = match.call()
  attr(d, "Diag") = diag
  attr(d, "Upper") = upper

  d
}


#' Greedy set cover.
#' 
#' A bool matrix input. Output both essential sets and elements.
#' 
#' @param mat A bool matrix of protein (cols)-peptide (rows) map. 
#' 
#' @return A two-column data frame of prot_acc and pep_seq. 
greedysetcover3 <- function (mat) 
{
  if (is.matrix(mat)) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    gc()
  }
  
  if (nrow(mat) == 1L || ncol(mat) == 1L) 
    return(data.frame(prot_acc = colnames(mat), pep_seq = rownames(mat)))

  prot_acc <- NULL
  pep_seq <- NULL
  
  while(nrow(mat)) {
    max <- which.max(Matrix::colSums(mat, na.rm = TRUE))
    
    if (max == 0L) break
    
    prot <- names(max)
    rows <- which(mat[, max])
    # peps <- names(rows) # name dropped if only one row
    peps <- rownames(mat)[rows]
    
    prot_acc <- c(prot_acc, rep(prot, length(peps)))
    pep_seq <- c(pep_seq, peps)
    
    mat <- mat[-rows, -max, drop = FALSE]
  }
  
  rm(list = c("mat"))
  gc()
  
  dplyr::bind_cols(prot_acc = prot_acc, pep_seq = pep_seq)
}


