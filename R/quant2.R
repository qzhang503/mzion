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
                         ppm_reporters = 10) {
  
  if (quant == "none") {
    out <- data
  } else {
    # message("Calculating reporter-ion intensities.")
    
    nms_tmt6 <- c("126", "127N", "128N", "129N", "130N", "131N")
    
    nms_tmt10 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                   "130N", "130C", "131N")
    
    nms_tmt11 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                   "130N", "130C", "131N", "131C")
    
    nms_tmt16 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
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
                    tmt6 = tmts %>% .[names(.) %in% nms_tmt6],
                    tmt10 = tmts %>% .[names(.) %in% nms_tmt10],
                    tmt11 = tmts %>% .[names(.) %in% nms_tmt11],
                    tmt16 = tmts %>% .[names(.) %in% nms_tmt16],
                    stop("Unknown TMt type.", call. = FALSE))
    
    ul <- switch(quant,
                 tmt6 = c(126.1, 131.2),
                 tmt10 = c(126.1, 131.2),
                 tmt11 = c(126.1, 131.2),
                 tmt16 = c(126.1, 134.2),
                 stop("Unknown TMt type.", call. = FALSE))
    
    # stopifnot(all(c("ms2_moverz", "ms2_int") %in% names(data)))
    
    out <- purrr::map2(data$ms2_moverz, data$ms2_int,
                       find_reporter_ints,
                       theos = theos,
                       ul = ul,
                       ppm_reporters = ppm_reporters,
                       len = length(theos),
                       nms = names(theos)) %>%
      dplyr::bind_rows() %>%
      dplyr::bind_cols(data, .)
  }
  
  names(out)[grep("^([0-9]{3}[NC]{0,1})", names(out))] <-
    find_int_cols(length(theos))
  
  invisible(out)
}


#' Adds back reporter-ion intensities
#'
#' @param df Results from \link{calc_protfdr}.
#' @inheritParams matchMS
add_rptrs <- function (df = NULL, quant = "none", out_path = NULL) {
  
  if (grepl("^tmt[0-9]+$", quant)) {
    files <- list.files(path = file.path(out_path, "temp"),
                        pattern = "^reporters_\\d+\\.rds$")
    
    idxes <- gsub("^reporters_([0-9]+)\\.rds$", "\\1", files) %>%
      as.integer() %>%
      order()
    
    files <- files[idxes]
    
    reporters <- lapply(files, function (x) {
      readRDS(file.path(out_path, "temp", x))
    }) %>%
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
                                ppm_reporters = 10, len, nms) {
  
  range <- findInterval(ul, ms2_moverzs)
  
  ms <-ms2_moverzs[range[1]:range[2]]
  is <-ms2_ints[range[1]:range[2]]
  
  idxes <- find_reporters_ppm(theos, ms, ppm_reporters, len, nms)
  
  if (!length(idxes)) {
    return(rep(NA, len) %>% `names<-`(nms))
  }
  
  # 126      127N      127C      128N      128N      128C
  # 135569.00 120048.00 122599.00   3397.98 140551.00 144712.00
  
  if (anyDuplicated(names(idxes))) {
    idxes <- idxes %>%
      split(., names(.)) %>%
      imap_int(~ {
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
  
  if (anyNA(names(idxes))) {
    names(idxes) <- nms
  }
  
  rptr_ints <- is[idxes] %>%
    `names<-`(names(idxes))
  
  if (length(rptr_ints) < len) {
    es <- rep(NA, len) %>%
      `names<-`(nms)
    
    es[names(rptr_ints)] <- rptr_ints
  } else {
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
find_reporters_ppm <- function (theos, expts, ppm_reporters = 10, len, nms) {
  
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_reporters, arr.ind = TRUE)
  
  row_cols[, 2]
}


#' Adds prot_acc to a peptide table
#'
#' @param out_path An output path.
#' @param df The results after scoring.
add_prot_acc <- function (df, out_path = "~/proteoM/outs") {
  
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
  out <- dplyr::bind_rows(fwd_prps, rev_prps) %>%
    dplyr::right_join(df, by = "pep_seq")
  
  rm(list = c("fwd_prps", "rev_prps", "uniq_peps"))
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
  
  out <- list(out, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = "prot_acc") %>%
    dplyr::arrange(-prot_n_pep, -prot_n_psm)
  
  rm(list = c("x", "prot_n_psm", "prot_n_pep"))
  gc()
  
  invisible(out)
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides. parDist <- cut_protgrps
#'
#' @param mat A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file pth to outputs.
cut_protgrps <- function (mat, out_path = NULL) {
  
  cns <- colnames(mat)
  
  mat <- as.list(mat)
  len <- length(mat)
  
  if (len <= 200L) {
    out <- vector("list", len)
    
    for (i in seq_len(len)) {
      out[[i]] <- map_dbl(mat[i:len], ~ sum(.x & mat[[i]]))
      out[[i]] <- c(out[seq_len(i-1)] %>% map_dbl(`[[`, i), out[[i]])
    }
  } else {
    out <- parDist(mat)
  }
  
  out <- do.call(rbind, out)
  rownames(out) <- colnames(out)
  
  # stopifnot(identical(out, t(out)))
  
  # --- finds protein groups
  out[out == 0L] <- 1000000
  out[out < 1000000] <- 0
  out <- out %>% as.dist(diag = TRUE, upper = TRUE)
  
  hc <- hclust(out, method = "single")
  
  grps <- data.frame(prot_hit_num = cutree(hc, h = 1)) %>%
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::group_by(prot_hit_num) %>%
    dplyr::mutate(prot_family_member = row_number()) %>%
    dplyr::ungroup() 
  
  # stopifnot(identical(grps$prot_acc, cns))
  
  invisible(grps)
}


#' Parallel distance calculations.
#'
#' @param mat A bool matrix.
#' @import parallel
parDist <- function (mat) {
  
  message("Calculating distance matrix.")
  
  gc()
  
  size <- object.size(mat)/1024^3
  mem <- find_free_mem() *.45/1024
  n_cores <- floor(min(mem/size, detect_cores(16L)))
  
  if (n_cores <= 1L) {
    stop("Not enough memory for parallel distance calculation.", 
         call. = FALSE)
  }
  
  idxes <- chunksplit(seq_along(mat), 2 * n_cores, "list")
  len <- length(mat)
  nms <- names(mat)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  
  out <- parallel::clusterApplyLB(cl, idxes, proteoCpp::par_distC, mat) %>%
    purrr::flatten()
  
  parallel::stopCluster(cl)
  rm(list = c("mat"))
  gc()
  
  if (len > 1L) {
    for (i in 2:len) {
      out[[i]] <- c(out[1:(i-1)] %>% purrr::map_dbl(`[[`, i), out[[i]])
    }
  }
  
  lapply(out, function (x) {
    names(x) <- nms
    x
  })
}


#' Helper of \link{groupProts}.
#'
#' @param out The data frame from upstream steps.
#' @param out_path The output path.
grp_prots <- function (out, out_path = NULL) {
  
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  
  out <- dplyr::arrange(out, pep_seq)
  
  # essential entries
  rows <- (out$pep_issig & (!out$pep_isdecoy) & (out$pep_rank <= 3L))
  
  df <- out[rows, ]
  
  if (nrow(df) > 1L) {
    df <- groupProts2(df, out_path)
  } else {
    df <- df %>%
      dplyr::mutate(prot_isess = TRUE,
                    prot_hit_num = 1L,
                    prot_family_member = 1L)
  }
  
  # non-essential entries
  prot_accs <- df %>%
    dplyr::filter(!duplicated(prot_acc), prot_isess) %>%
    `[[`("prot_acc")
  
  df2 <- out[!rows, ] %>%
    dplyr::filter(!duplicated(prot_acc)) %>%
    dplyr::select(prot_acc) %>%
    dplyr::mutate(prot_isess = ifelse(prot_acc %in% prot_accs, TRUE, FALSE)) %>%
    dplyr::left_join(df[, c("prot_acc", "prot_hit_num", "prot_family_member")] %>%
                       dplyr::filter(!duplicated(prot_acc)),
                     by = "prot_acc") %>%
    dplyr::right_join(out[!rows, ], by = "prot_acc") %>%
    dplyr::mutate(pep_literal_unique = NA, pep_razor_unique = NA) %>%
    dplyr::select(names(df))
  
  rm(list = c("prot_accs"))
  
  out <- dplyr::bind_rows(df, df2) %>%
    dplyr::select(-which(names(.) %in% c("prot_n_psm", "prot_n_pep")))
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts2 <- function (df, out_path = NULL) {
  
  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  # --- (1) protein ~ peptide map ---
  mat <- map_pepprot2(df[, c("prot_acc", "pep_seq")], out_path)
  gc()
  
  # --- (2) protein ~ protein groups by distance map ---
  grps <- cut_protgrps2(mat, out_path)
  gc()
  
  # --- (3) set covers by groups ---
  sets <- greedysetcover3(mat)
  gc()
  
  if (!is.null(out_path)) {
    saveRDS(sets, file.path(out_path, "prot_pep_setcover.rds")) 
  }
  
  sets <- sets %>% 
    `[[`("prot_acc") %>%
    unique()
  gc()
  
  # --- set aside df0 ---
  df <- dplyr::mutate(df, prot_isess = prot_acc %in% sets)
  df0 <- dplyr::filter(df, !prot_isess)
  df <- dplyr::filter(df, prot_isess)

  mat_ess <- mat[, colnames(mat) %in% unique(df$prot_acc)]
  
  peps_uniq <- local({
    rsums <- Matrix::rowSums(mat)
    rsums2 <- Matrix::rowSums(mat_ess)
    
    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })
  
  # --- put together ---
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)
  
  df <- df %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")
  
  gc()
  
  invisible(df)
}


#' Helper of \link{groupProts2}.
#'
#' Builds the logical map between peptide (in rows) and proteins (in columns).
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
#' out <- proteoM:::map_pepprot2(df)
#' }
map_pepprot2 <- function (df, out_path = NULL) {

  df <- df[, c("prot_acc", "pep_seq")] 
  df <- unique(df)
  
  gc()
  
  peps <- df$pep_seq
  
  mat <- Matrix::sparse.model.matrix(~ -1 + prot_acc, df)
  colnames(mat) <- gsub("prot_acc", "", colnames(mat))
  mat <- mat == 1L
  rownames(mat) <- peps
  gc()
  
  # ---
  dpeps <- peps[duplicated(peps)]
  drows <- (rownames(mat) %in% dpeps)
  mat0 <- mat[!drows, ]
  mat <- mat[drows, ]
  
  # ---
  mpeps <- unique(rownames(mat))
  
  len <- as.numeric(length(mpeps))
  ncol <- as.numeric(ncol(mat))
  len2 <- len * ncol
  out <- rep(0L, len2)
  
  start <- 1
  end <- ncol

  for (i in seq_len(len)) {
    pep <- rownames(mat)[[1]]
    rows <- rownames(mat) == pep
    
    mati <- mat[rows, ]
    out[start:end] <- Matrix::colSums(mati)

    mat <- mat[!rows, ]
    
    start <- start + ncol
    end <- end + ncol
    
    if (i %% 100 == 0) gc()
  }
  
  rm(list = c("mat", "mati"))
  gc()
  
  # ---
  if (object.size(out)/1024^3 > 5) {
    size <- 10000
    n_chunks <- ceiling(len/size)
    
    x0 <- NULL
    
    for (i in 1:n_chunks) {
      x <- out[(ncol*(size*(i-1))+1):min(len2, (ncol*(size*i)))]
      x <- Matrix::Matrix(x, ncol = ncol, byrow = TRUE, sparse = TRUE)
      gc()
      x0 <- rbind2(x0, x)
    }
    
    rownames(x0) <- mpeps
    out <- x0
    
    rm(list = c("x", "x0"))
    gc()
  } else {
    out <- Matrix::Matrix(out, ncol = ncol, byrow = TRUE, sparse = TRUE)
    rownames(out) <- mpeps
    gc()
  }

  out <- out == 1L
  gc()
  
  out <- rbind2(out, mat0)
  gc()
  
  # collapse rows of the same pep_seq; may use `sum`
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
  # A tibble: 3 x 3
  # pep_seq X     Y
  # <chr>   <lgl> <lgl>
  # 1 A       TRUE  TRUE
  # 2 B       TRUE  FALSE
  # 3 C       FALSE TRUE
  
  if (!is.null(out_path)) {
    Matrix::writeMM(out, file = file.path(out_path, "prot_pep_map.mtx"))
    saveRDS(colnames(out), file.path(out_path, "prot_pep_map_col.rds"))
    saveRDS(rownames(out), file.path(out_path, "prot_pep_map_row.rds"))
  }
  
  invisible(out)
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides.
#'
#' @param mat A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file pth to outputs.
cut_protgrps2 <- function (mat = NULL, out_path = NULL) {

  dista = proxyC::simil(mat, margin = 2) # sparse distance matrix
  cns <- colnames(mat)
  rm(list = c("mat"))
  gc()
  
  # ---
  ncol <- ncol(dista)
  cols <- 1:ncol

  if (ncol > 10000) {
    mat2 <- matrix(nrow = ncol, ncol = ncol)
    colnames(mat2) <- cns
    rownames(mat2) <- cns
    
    max_rc <- 2500 * 40000 # 5000 out of memory
    gc()
    
    size <- ceiling(max_rc/ncol)
    n_chunks <- ceiling(ncol/size)
    
    for (i in 1:n_chunks) {
      rows <- (1+(i-1)*size):min(i*size, ncol)

      x <- dista[rows, cols] 
      gc()
      
      x <- (x == 0) # sparse logical matrix
      gc()
      
      mat2[rows, cols] <- as.matrix(x) # regular matrix
      gc()
    }
    
    rm(list = c("x"))
    gc()
  } else {
    dista <- (dista == 0) # sparse logical matrix
    gc()
    
    dista <- as.matrix(dista) # regular matrix
    gc()
    
    mat2 <- dista
    
    rm(list = c("dista"))
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
  mat2 <- as_lgldist(mat2, diag = FALSE, upper = FALSE)
  gc()
  
  hc <- hclust(mat2, method = "single")
  gc()
  
  grps <- data.frame(prot_hit_num = cutree(hc, h = .9)) %>% 
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::group_by(prot_hit_num) %>%
    dplyr::mutate(prot_family_member = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  if (!is.null(out_path)) {
    saveRDS(grps, file.path(out_path, "prot_grps.rds"))
  }
  
  stopifnot(identical(grps$prot_acc, cns))
  
  invisible(grps)
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Assumed the input is already a symmetric matrix (not yet being used).
#' 
#' @inheritParams stats::as.dist
as_dist <- function (m, diag = FALSE, upper = FALSE) {
  
  p <- nrow(m)
  
  ans <- m[row(m) > col(m)]
  gc()
  attributes(ans) <- NULL
  
  if (!is.null(rownames(m))) {
    attr(ans, "Labels") <- rownames(m)
  } else if (!is.null(colnames(m))) {
    attr(ans, "Labels") <- colnames(m)
  }
  
  attr(ans, "Size") <- p
  attr(ans, "call") <- match.call()
  class(ans) <- "dist"
  
  if (is.null(attr(ans, "Diag")) || !missing(diag)) {
    attr(ans, "Diag") <- diag
  }
    
  if (is.null(attr(ans, "Upper")) || !missing(upper)) {
    attr(ans, "Upper") <- upper
  }

  ans
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Assumed the input is already a symmetric matrix.
#' 
#' @inheritParams stats::as.dist
as_lgldist <- function(m, diag = FALSE, upper = FALSE) {
  
  d = proteoCpp::to_lgldistC(m)
  
  if (!is.null(rownames(m))) {
    attr(d, "Labels") <- rownames(m)
  } else if (!is.null(colnames(m))) {
    attr(d, "Labels") <- colnames(m)
  }

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
greedysetcover3 <- function (mat) {
  
  if (is.matrix(mat)) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    gc()
  }
  
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


