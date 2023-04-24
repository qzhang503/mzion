#' Reporter-ion quantitation.
#' 
#' Not yet used: \code{`134C` = 134.154565}, \code{`135N` = 135.15160}
#'  
#' @param data An upstream result from \link{matchMS}.
#' @param quant A quantitation method. The default is "none". Additional choices
#'   include \code{tmt6} etc. For other multiplicities of \code{tmt}, use the
#'   compatible higher plexes, for example, \code{tmt16} for \code{tmt12} etc.
#'   and \code{tmt10} for \code{tmt8} etc.
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
#' @inheritParams matchMS
calc_tmtint <- function (data = NULL,
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"),
                         ppm_reporters = 10L, index_mgf_ms2 = FALSE) 
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
                    "133N", "133C", "134N")
    
    tmts <- c(
      `126` = 126.127726, `127N` = 127.124761, `127C` = 127.131080,
      `128N` = 128.128115, `128C` = 128.134435, `129N` = 129.131470,
      `129C` = 129.137790, `130N` = 130.134825, `130C` = 130.141145,
      `131N` = 131.138180, `131C` = 131.144499, `132N` = 132.141535,
      `132C` = 132.147855, `133N` = 133.14489, `133C` = 133.15121,
      `134N` = 134.148245)

    theos <- switch(quant,
                    tmt6  = tmts[names(tmts) %in% nms_tmt6],
                    tmt10 = tmts[names(tmts) %in% nms_tmt10],
                    tmt11 = tmts[names(tmts) %in% nms_tmt11],
                    tmt16 = tmts[names(tmts) %in% nms_tmtpro],
                    stop("Unknown TMt type.", call. = FALSE))
    
    ul <- switch(quant,
                 tmt6 = c(126.1, 131.2),
                 tmt10 = c(126.1, 131.2),
                 tmt11 = c(126.1, 131.2),
                 tmt16 = c(126.1, 134.2),
                 stop("Unknown TMt type.", call. = FALSE))
    
    # stopifnot(all(c("rptr_moverz", "rptr_int") %in% names(data)))
    
    col_rptr_mzs <- "rptr_moverz"
    col_rptr_int <- "rptr_int"

    out <- mapply(find_reporter_ints, data[["rptr_moverz"]], data[["rptr_int"]], 
                  MoreArgs = list(
                    theos = theos,
                    ul = ul,
                    ppm_reporters = ppm_reporters,
                    len = length(theos),
                    nms = names(theos)
                  ), USE.NAMES = FALSE, SIMPLIFY = FALSE)
      
    out <- dplyr::bind_rows(out)

    if (!nrow(out)) {
      out <- data.frame(matrix(ncol = length(theos), nrow = 0L))
      colnames(out) <- theos
      
      for (i in seq_along(out)) 
        out[[i]] <- as.numeric(out[[i]])
    }

    data[["rptr_moverz"]] <- data[["rptr_int"]] <- NULL

    out <- dplyr::bind_cols(data, out)
  }
  
  cols <- grep("^([0-9]{3}[NC]{0,1})", names(out))
  names(out)[cols] <- find_int_cols(length(theos))
  
  invisible(out)
}


#' Adds back reporter-ion intensities
#'
#' @param df Results from \link{calc_protfdr}.
#' @inheritParams matchMS
add_rptrs <- function (df = NULL, quant = "none", out_path = NULL) 
{
  if (!grepl("^tmt[0-9]+$", quant))
    return(df)
  
  pattern <- "^reporters_\\d+\\.rds$"
  files <- list.files(path = file.path(out_path, "temp"), pattern = pattern)
  idxes <- order(as.integer(gsub("^reporters_([0-9]+)\\.rds$", "\\1", files)))
  files <- files[idxes]
  
  reporters <- lapply(files, function (x) qs::qread(file.path(out_path, "temp", x)))
  reporters <- dplyr::bind_rows(reporters)
  
  df <- df %>%
    tidyr::unite(uniq_id, raw_file, pep_mod_group, pep_scan_num, sep = ".",
                 remove = FALSE) %>%
    dplyr::left_join(reporters, by = "uniq_id") %>%
    dplyr::select(-uniq_id)
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
#' library(mzion)
#' 
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
#' x <- mzion:::find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 10,
#'                         len , nms)
#'
#' x <- mzion:::find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 25,
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
#' x <- mzion:::find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 25,
#'                         len , nms)
#' }
find_reporter_ints <- function (ms2_moverzs, ms2_ints, theos, ul,
                                ppm_reporters = 10L, len, nms) 
{
  range <- findInterval(ul, ms2_moverzs)
  
  ms <- ms2_moverzs[range[1]:range[2]]
  is <- ms2_ints[range[1]:range[2]]
  
  idxes <- find_reporters_ppm(theos, ms, ppm_reporters, len, nms)
  
  if (!length(idxes)) 
    return(rep(NA, len) %>% `names<-`(nms))
  
  # 126      127N      127C      128N      128N      128C
  # 135569.00 120048.00 122599.00   3397.98 140551.00 144712.00
  
  if (anyDuplicated(names(idxes))) {
    idxes <- idxes %>%
      split(names(.)) %>%
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
  
  if (anyNA(names(idxes))) 
    names(idxes) <- nms
  
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
find_reporters_ppm <- function (theos, expts, ppm_reporters = 10L, len, nms) 
{
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_reporters, arr.ind = TRUE)
  
  row_cols[, 2]
}


#' Helper of \link{sub_protpep}
#'
#' Subsets the protein-peptide lookup (prps) against the uniq_peps list from
#' search results.
#'
#' @param prps A protein-peptide lookup.
#' @param uniq_peps A non-redundant list of peptides from search results.
msub_protpep <- function (prps, uniq_peps)
{
  mts <- lapply(prps, fastmatch::fmatch, uniq_peps)
  
  ans <- mapply(sub_protpep, prps, mts, names(mts),
                SIMPLIFY = FALSE, USE.NAMES = FALSE)

  dplyr::bind_rows(ans)
}


#' Subsets the protein-peptide lookups.
#' 
#' Duplicated peptide sequences with a protein also removed.
#'
#' @param prp A vector of peptides under a protein.
#' @param mt A vector of matches. The values are NA if no matches. The matched
#'   integers are the indexes in the list of uniq_peps for protein annotations.
#' @param nm The names of \code{prp}.
sub_protpep <- function (prp, mt, nm)
{
  nas <- is.na(mt)
  
  if (all(nas))
    return(NULL)
  
  oks <- !nas
  ps  <- which(oks)
  
  len <- length(prp)
  vc <- vn <- vector("logical", len)

  if (length(ps)) {
    pnt <- ps[ps %fin% attr(prp, "pnt_idxes", exact = TRUE)]
    pct <- ps[ps %fin% attr(prp, "pct_idxes", exact = TRUE)]
  }
  else {
    pct <- pnt <- numeric()
  }
  
  vn[pnt] <- TRUE
  vc[pct] <- TRUE

  u <- !duplicated.default(mt[ps])
  # pep_id = mt[ps][u]
  
  list(prot_acc = rep(nm, sum(u)), 
       pep_seq = prp[ps][u], 
       is_pnt = vn[ps][u], 
       is_pct = vc[ps][u])
}


#' Adds prot_acc to a peptide table
#'
#' Cached results are under sub dirs.
#' 
#' @param out_path An output path.
#' @param df The results after scoring.
#' @inheritParams matchMS
add_protacc2 <- function (df = NULL, out_path = NULL, .path_cache = NULL, 
                          .path_fasta = NULL) 
{
  message("Adding protein accessions (no enzyme specificity).")
  
  if (is.null(df)) {
    file <- file.path(out_path, "temp", "peploc.rds")
    
    if (file.exists(file))
      df <- qs::qread(file)
    else
      stop("File not found: ", file)
    
    rm(list = c("file"))
  }
  
  sub_dirs <- dir(out_path, pattern = "^sub[0-9]+_[0-9]+_[0-9]+$", 
                  full.names = TRUE)
  len_dirs <- length(sub_dirs)
  
  if (!len_dirs) {
    warning("Cached results not found for protein annotation.")
    return (df)
  }
  
  prps <- vector("list", len_dirs)
  
  for (i in seq_along(sub_dirs)) {
    cache_file <- file.path(sub_dirs[[i]], "Calls/.cache_info.rds")
    
    if (!file.exists(cache_file))
      stop("Cached file not found: ", cache_file)
    
    cache_info <- load_cache_info(cache_file)
    .time_stamp <- cache_info[[".time_stamp"]]
    .path_ms1masses <- cache_info[[".path_ms1masses"]]
    
    file <- file.path(.path_ms1masses, .time_stamp, "simple_prot_pep.rds")
    prps[[i]] <- if (file.exists(file)) qs::qread(file) else NULL
  }
  
  upeps <- unique(df$pep_seq)
  prps <- lapply(prps, msub_protpep, upeps)
  prps <- dplyr::bind_rows(prps)
  
  df <- hannot_decoys(df, prps)
}


#' Adds prot_acc to a peptide table.
#'
#' Decoys being kept.
#'
#' @param out_path An output path.
#' @param df The results after scoring.
#' @inheritParams matchMS
#' @importFrom fastmatch fmatch %fin% 
add_protacc <- function (df = NULL, out_path = NULL, .path_cache = NULL, 
                         .path_fasta = NULL) 
{
  message("Adding protein accessions.")
  
  if (is.null(df)) {
    file <- file.path(out_path, "temp", "peploc.rds")
    
    if (file.exists(file))
      df <- qs::qread(file)
    else
      stop("File not found: ", file)
    
    rm(list = c("file"))
  }
  
  # need to inverse the sequence at NA pep_ivmod: 
  # decoy MS2 are appended to targets and share the same pep_seq name as targets 
  #   where the pep_ivmod values are NA with the decoys
  
  .path_ms1masses <- create_dir(file.path(.path_fasta, "ms1masses"))
  .time_stamp <- find_ms1_times(out_path)
  
  if (length(.time_stamp) > 1L)
    stop("Multiple matches in time stamps in annotating protein acccessions.")
  
  prps <- qs::qread(file.path(.path_ms1masses, .time_stamp, "simple_prot_pep.rds"))
  prps <- msub_protpep(prps, unique(df[["pep_seq"]]))
  
  df <- hannot_decoys(df, prps)
}


#' Helper of annotating decoy peptides.
#' 
#' @param df A data frame.
#' @param prps The look-ups of protein and peptides.
hannot_decoys <- function (df, prps)
{
  # keep prot_acc be the first column
  df <- dplyr::right_join(prps, df, by = "pep_seq")
  rows <- is.na(df[["pep_ivmod"]])
  tars <- df[!rows, ]
  
  decs <- df[rows, ]
  decs[["pep_seq"]] <- reverse_seqs(decs[["pep_seq"]])
  oks <- fastmatch::fmatch(decs[["pep_seq"]], tars[["pep_seq"]])
  decs <- decs[is.na(oks), ]
  decs[["prot_acc"]] <- paste0("-", decs[["prot_acc"]])
  df <- dplyr::bind_rows(tars, decs)
  
  # Adds prot_n_psm, prot_n_pep for protein FDR estimates
  x <- df[df[["pep_issig"]], ]
  
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
  pnt_nots <- grepl("Protein N-term", df[["pep_vmod"]], fixed = TRUE) & !df[["is_pnt"]]
  pct_nots <- grepl("Protein C-term", df[["pep_vmod"]], fixed = TRUE) & !df[["is_pct"]]
  df <- df[(!pnt_nots) & (!pct_nots), ]
  
  # peptide `is_pnt` and `is_pct` are not EXACT facts
  # but preference of terminal over interior matches
  # so remove them to avoid misleading uses or interpretations
  df[["is_pnt"]] <- NULL
  df[["is_pct"]] <- NULL
  
  list(df, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = "prot_acc") %>%
    dplyr::arrange(-prot_n_pep, -prot_n_psm)
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} etc. to
#' \code{psmQ.txt}.
#'
#' In general, non-significant and decoy peptides should have been removed from
#' the input \code{df}, as well as decoy proteins.
#'
#' In addition, there is no duplicated entries.
#'
#' @param df A two-column data frame contains \code{prot_acc} and
#'   \code{pep_seq}.
#' @param out_path The output path.
#' @param out_name The output filename.
#' @param fct A factor for data splitting into chunks.
#' 
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ]
#'
#' out <- mzion:::groupProts(df, "~")
#'
#' # One peptide, multiple proteins
#' df <- data.frame(prot_acc = LETTERS[1:3], pep_seq = rep("X", 3))
#' out <- mzion:::groupProts(df, "~")
#' stopifnot(nrow(out) == 3L)
#'
#' # One peptide, one proteins
#' df <- data.frame(prot_acc = "A", pep_seq = "X")
#' out <- mzion:::groupProts(df, "~")
#' stopifnot(nrow(out) == 1L)
#'
#' # One proteins
#' df <- data.frame(prot_acc = rep("A", 3), pep_seq = LETTERS[24:26])
#' out <- mzion:::groupProts(df, "~")
#' stopifnot(nrow(out) == 3L)
#' }
groupProts <- function (df, out_path = NULL, fct = 4L, 
                        out_name = "prot_pep_setcover.rds") 
{
  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  if (!identical(names(df), c("prot_acc", "pep_seq")))
    stop("The two columns of `df` need to be in the order of ", "
         \"prot_acc\" and \"pep_seq\".")
  
  if (!nrow(df))
    stop("Zero row of data for protein groupings.")
  
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  
  df <- df[with(df, order(pep_seq)), ]
  
  if (nrow(df) <= 1L) {
    df <- dplyr::mutate(df, 
                        prot_isess = TRUE,
                        prot_hit_num = 1L,
                        prot_family_member = 1L, 
                        pep_literal_unique = TRUE, 
                        pep_razor_unique = TRUE)
    
    return(df)
  }

  ## (1) builds protein ~ peptide map
  Mats <- map_pepprot(df, out_path = out_path, fct = fct)
  message("Completed protein-peptide maps.")
  
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
  message("Completed protein grouping.")

  ## (3) finds essential protein entries
  ess_prots <- local({
    if (is.null(Mat_lwr_left))
      df_shared <- greedysetcover3(Mat_upr_left)
    else {
      # unique peptides of proteins with shared peptides
      rows_lwr_left_is_one <- Matrix::rowSums(Mat_lwr_left) > 0
      Mat_lwr_left_is_one <- Mat_lwr_left[rows_lwr_left_is_one, , drop = FALSE]
      
      # set covers of shared + unique peptides under shared proteins
      df_shared <- greedysetcover3(rbind2(Mat_upr_left, Mat_lwr_left_is_one))
      gc()
    }

    # proteins with exclusive unique peptides
    df_uniq <- df[, c("prot_acc", "pep_seq")]
    df_uniq <- df_uniq[df_uniq$prot_acc %in% prots_upr_right, , drop = FALSE]
    df_uniq <- df_uniq[!duplicated.data.frame(df_uniq), , drop = FALSE]
    df_uniq <- df_uniq[with(df_uniq, order(prot_acc, pep_seq)), , drop = FALSE]
    
    sets <- rbind2(df_shared, df_uniq)
    
    if (!is.null(out_path)) {
      qs::qsave(sets, file.path(out_path, out_name), preset = "fast")
    }

    unique(sets$prot_acc)
  })
  message("Established essential proteins.")

  ## Parses literal or razor uniqueness of peptides
  # sets aside df0
  df <- dplyr::mutate(df, prot_isess = prot_acc %in% ess_prots)
  df0 <- dplyr::filter(df, !prot_isess)
  df <- dplyr::filter(df, prot_isess)
  gc()
  
  # combines four quadrants
  M4 <- rbind2(
    cbind2(Mat_upr_left, Mat_upr_right), 
    cbind2(Mat_lwr_left, Mat_lwr_right)
  )
  
  rm(list = c("Mat_upr_left", "Mat_upr_right", 
              "Mat_lwr_left", "Mat_lwr_right"))
  gc()
  
  M4_ess <- if (nrow(M4) == 1L) 
    M4
  else 
    M4[, colnames(M4) %in% ess_prots, drop = FALSE]
  
  # literal: unique in M4
  # razor: unique in M4_ess
  peps_uniq <- local({
    rsums <- Matrix::rowSums(M4)
    rsums2 <- Matrix::rowSums(M4_ess)
    
    peps <- data.frame(pep_seq = rownames(M4)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })
  
  rm(list = c("M4", "M4_ess"))
  gc()
  message("Parsed unique versus shared peptides.")
  
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
#'   columns of \code{prot_acc} and \code{pep_seq}. It should be TRUE that
#'   \code{df} is identical to \code{unique(df)}.
#' @param out_path An output path.
#' @param fct A factor for data splitting into chunks.
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ]
#'
#' out <- mzion:::map_pepprot(df)
#'
#' # One peptide, multiple proteins
#' df <- data.frame(prot_acc = LETTERS[1:3], pep_seq = rep("X", 3))
#' out <- mzion:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == "X", colnames(out[[1]]) == LETTERS[1:3])
#'
#' # One peptide, one proteins
#' df <- data.frame(prot_acc = "A", pep_seq = "X")
#' out <- mzion:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == "X", colnames(out[[1]]) == "A")
#'
#' # One proteins
#' df <- data.frame(prot_acc = rep("A", 3), pep_seq = LETTERS[24:26])
#' out <- mzion:::map_pepprot(df)
#' stopifnot(rownames(out[[1]]) == LETTERS[24:26], colnames(out[[1]]) == "A")
#' }
map_pepprot <- function (df, out_path = NULL, fct = 4L) 
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
  
  # df <- df[!duplicated.data.frame(df), ] # should not contain duplicated entries
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
  Mat0  <- Mat[!drows, ]
  Mat1  <- Mat[drows, ]

  rm(list = c("dpeps", "drows", "Mat", "peps"))
  gc()
  
  ## Mat1: (a) pre sparse-matrix vector
  ncol  <- as.numeric(ncol(Mat1))
  peps1 <- rownames(Mat1)
  vec   <- pcollapse_sortpeps(Mat = Mat1, ncol = ncol, peps = peps1, fct = fct)
  
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
  
  rm(list = c("peps"))
  gc()

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
    
    if (i %% 100 == 0) gc()
  }
  
  rm(list = c("mat"))
  gc()

  invisible(out)
}


#' Helper of \link{map_pepprot}.
#' 
#' Parallel version of \link{collapse_sortpeps}.
#' 
#' @param fct A factor for data splitting into chunks.
#' @param Mat A sparse matrix.
#' @inheritParams collapse_sortpeps
pcollapse_sortpeps <- function (Mat, ncol = NULL, peps = NULL, fct = 4L) 
{
  if (is.null(peps)) 
    peps <- rownames(Mat)
  
  if (is.null(ncol))
    ncol <- as.numeric(ncol(Mat))
  
  size <- local({
    dim <- dim(Mat)
    as.numeric(dim[1]) * as.numeric(dim[2])
  })
  
  n_cores <- detect_cores(16L)
  n_cores <- min(n_cores, floor(n_cores * 7E9 /size))

  # !!! `peps` must be SORTED !!!

  if (size <= 1E8 || n_cores <= 1L) {
    vec <- collapse_sortpeps(Mat, ncol, peps)
  }
  else {
    Mats <- chunksplit_spmat(Mat, peps, n_cores * fct)
    rm(list = c("Mat", "peps"))
    gc()
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    # don't delete, otherwise return 0L
    parallel::clusterEvalQ(cl, library(Matrix))
    vecs <- parallel::clusterApplyLB(cl, Mats, collapse_sortpeps)
    parallel::stopCluster(cl)
    gc()
    
    vec <- integer(sum(unlist(lapply(vecs, length))))
    len <- length(vecs)
    sta <- 0L
    
    while(len) {
      vec_1 <- vecs[[1]]
      len_1 <- length(vec_1)
      end <- sta + len_1
      
      if (is.na(end)) {
        message("Handling integer overflow.")
        
        len_1 <- as.numeric(len_1)
        vec[(sta + 1):(sta + len_1)] <- vec_1
        sta <- sta + len_1
      }
      else {
        vec[(sta + 1L):end] <- vec_1
        ok_sta <- sta + len_1
        sta <- if (is.na(sta)) sta + as.numeric(len_1) else ok_sta
      }

      vecs[1] <- NULL
      len <- len - 1L
      rm(list = c("vec_1", "len_1"))
      gc()
    }
  }
  
  message("\tCompleted matrix-to-vector conversion.")
  
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
#' library(mzion)
#' 
#' vec <- rep(LETTERS[1:5], 1:5)
#' vec <- sort(vec)
#' 
#' mzion:::find_group_breaks(vec, 3)
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
  } 
  else {
    D <- proxyC::simil(M, margin = 2) # dsTMatrix
  }
  
  rm(list = c("M"))
  gc()
  
  # n_prots > 10000L
  if (FALSE) {
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
  # d <- as_lgldist(dm, diag = FALSE, upper = FALSE) # logical distance
  d <- as.dist(dm)
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

  if (!is.null(out_path)) {
    qs::qsave(grps, file.path(out_path, "prot_grps.rds"), preset = "fast")
  }

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
#' \donttest{
#' library(mzion)
#' # m <- mzion:::sparseD_fourquad(ul, 6)
#' }
sparseD_fourquad <- function (M_ul, ncols_ur = 0L) 
{
  nrows_ul <- ncols_ul <- ncol(M_ul)
  nrows_l <- ncols_ur
  ncols <- ncols_ul + ncols_ur
  
  M_ur <- Matrix::sparseMatrix(dims = c(nrows_ul, ncols_ur), i={}, j={})
  m_lwr <- Matrix::sparseMatrix(i = 1:nrows_l, j = (ncols_ul+1):ncols, x = 1)
  
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
  d <- proteoCpp::to_lgldistC(m)
  
  if (!is.null(rownames(m))) 
    attr(d, "Labels") <- rownames(m)
  else if (!is.null(colnames(m))) 
    attr(d, "Labels") <- colnames(m)

  attr(d, "class") <- "dist"
  attr(d, "Size") <- nrow(m)
  attr(d, "call") <- match.call()
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper

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
  # dense matrix to sparse matrix
  if (is.matrix(mat)) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
    gc()
  }
  
  if (nrow(mat) == 1L || ncol(mat) == 1L) 
    return(data.frame(prot_acc = colnames(mat), pep_seq = rownames(mat)))

  prot_acc <- NULL
  pep_seq <- NULL
  
  while(nrow(mat)) {
    max <- which.max(Matrix::colSums(mat, na.rm = TRUE))
    
    if (max == 0L) 
      break
    
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


