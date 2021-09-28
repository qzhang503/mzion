#' Helper in binning precursor masses.
#' 
#' For target peptides.
#' 
#' @param res The results from \link{calc_pepmasses2}.
#' @param min_mass A minimum mass of precursors.
#' @param max_mass A maximum mass of precursors.
#' @inheritParams matchMS
bin_ms1masses <- function (res = NULL, min_mass = 500L,
                           max_mass = 10000L, ppm_ms1 = 20L) {
  
  .path_fasta <- get(".path_fasta", envir = .GlobalEnv)
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv)

  out_path <- file.path(.path_fasta, "pepmasses", .time_stamp)
  bins <- list.files(path = out_path, pattern = "binned_theopeps_\\d+\\.rds$")

  bin_ms1masses_td(bins = bins,
                   type = "target",
                   res = res,
                   min_mass = min_mass,
                   max_mass = max_mass,
                   ppm_ms1 = ppm_ms1,
                   out_path = out_path)

  invisible(NULL)
}


#' Helper of \link{bin_ms1masses} by the separate sets of \code{target}s and
#' \code{decoy}s.
#' 
#' For either target or decoy peptides.
#'
#' @param bins The names of bins in the format of
#'   "binned_theopeps_(rev_)*\\d+\\.rds$".
#' @param type The type of data: "target" or "decoy".
#' @inheritParams binTheoSeqs
bin_ms1masses_td <- function (bins = NULL, type = c("target", "decoy"),
                              res = NULL, min_mass = 500L, max_mass = 10000L,
                              ppm_ms1 = 20L, out_path = NULL) {

  len_b <- length(bins)

  type2 <- switch(type,
                  target = "",
                  "decoy" = "rev_",
                  stop("Unknown `type = ", type, "`."))

  masses <- list.files(
    path = out_path,
    pattern = paste0("^pepmasses_", type2, "\\d+\\.rds$"))

  len_m <- length(masses)

  if (!len_m) {
    stop("File not found: ",
         file.path(out_path, paste0("pepmasses_", type2, "[...].rds")))
  }

  ## already binned
  if (len_b == len_m) {
    message("Loading bins of MS1 masses from cache.")

    return(NULL)
  }

  ## to be binned
  message("Binning MS1 masses (theoretical ", type, ").")

  # (a) process directly
  if (!is.null(res)) {
    binTheoSeqs(idxes = NULL,
                res = res,
                min_mass = min_mass,
                max_mass = max_mass,
                ppm_ms1 = ppm_ms1,
                out_path = file.path(out_path, "binned_theopeps.rds"))

    return(NULL)
  }

  # (b) reload
  idxes <- local({
    idxes <- gsub("^(pepmasses_|pepmasses_rev_)(\\d+)\\.rds$", "\\2", masses)
    idxes <- as.integer(idxes)
    idxes <- idxes[order(idxes)]
    idxes <- paste0(type2, idxes)
  })

  n_cores <- local({
    max_n_cores <- 8L
    
    fct <- 20 
    free_mem <- find_free_mem()
    max_sz <- max(file.size(file.path(out_path, masses)))/1024^2

    n_cores <- min(max_n_cores, 
                   floor(free_mem/max_sz/fct), 
                   detect_cores())

    if (n_cores < 1L) {
      warning("May be out of memory with large peptide tables.")
      n_cores <- 1L
    }

    n_cores
  })

  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(
      cl,
      c("binTheoSeqs_i", 
        "binTheoSeqs2", 
        "bin_theoseqs", 
        "find_ms1_cutpoints"), 
      envir = environment(proteoM:::binTheoSeqs_i))
    
    # No need of purrr::flatten() as saveRDS by INDIVIDUAL idx (and return NULL)
    parallel::clusterApplyLB(
      cl = cl, 
      x = chunksplit(idxes, n_cores, "list"), 
      fun = lapply, 
      FUN = "binTheoSeqs_i", 
      min_mass = min_mass, 
      max_mass = max_mass, 
      ppm_ms1 = ppm_ms1, 
      out_path = out_path
    )

    parallel::stopCluster(cl)
  } else {
    lapply(idxes, binTheoSeqs_i, min_mass, max_mass, ppm_ms1, out_path)
  }

  invisible(NULL)
}


#' Helper of \link{binTheoSeqs2}.
#'
#' @inheritParams binTheoSeqs2
binTheoSeqs_i <- function (idx = 1L, min_mass = 500L,max_mass = 10000L,
                           ppm_ms1 = 20L, out_path = NULL) {
  
  message("\tSet: ", idx)
  
  in_nm <- paste0("pepmasses_", idx, ".rds")
  res <- s_readRDS(in_nm, out_path)
  
  binTheoSeqs2(idx = idx,
               res = res,
               min_mass = min_mass,
               max_mass = max_mass,
               ppm_ms1 = ppm_ms1,
               out_path = file.path(out_path, "binned_theopeps.rds"))
  
}


#' Separates theoretical peptides into mass groups.
#'
#' @param idx An index, e.g. "1" for \code{pepmasses_1.rds} and "rev_1" for
#'   \code{pepmasses_rev_1.rds}.
#' @inheritParams binTheoSeqs
binTheoSeqs2 <- function (idx = 1L, res = NULL, min_mass = 500L,
                          max_mass = 10000L, ppm_ms1 = 20L, out_path = NULL) {
  
  if (is.null(res)) {
    stop("`res` cannot be NULL.")
  }
  
  if (is.null(out_path)) {
    stop("`out_path` cannot be NULL.")
  }
  
  out_dir <- create_dir(gsub("(^.*/).*$", "\\1", out_path))
  
  out_nm <- gsub("^.*/(.*)\\.[^\\.].*$", "\\1", out_path) %>%
    paste(idx, sep = "_") %>%
    paste0(".rds")
  
  res <- attr(res, "data")
  gc()
  
  bin_theoseqs(res, out_nm, min_mass, max_mass, ppm_ms1)
  
  rm(list = c("res"))
  gc()
  
  invisible(NULL)
}


#' Separates theoretical peptides into mass groups.
#'
#' @param peps A list of theoretical peptides with masses.
#' @param out_nm A output name.
#' @inheritParams binTheoSeqs
bin_theoseqs <- function (peps = NULL, out_nm = NULL, min_mass = 500L, 
                          max_mass = 10000L, ppm_ms1 = 20L) {

  ps <- find_ms1_cutpoints(min_mass, max_mass, ppm_ms1)
  frames <- findInterval(peps, ps)
  
  out <- data.frame(pep_seq = names(peps),
                    mass = peps,
                    frame = frames, 
                    row.names = NULL)
  
  out <- dplyr::arrange(out, frame, pep_seq)
  out <- split(out, out$frame, drop = FALSE)
  
  saveRDS(out, out_nm)
  
  invisible(NULL)
}


#' Helper of \link{bin_theoseqs}.
#' 
#' @param idxes A set of indexes, e.g. "1" for \code{pepmasses_1.rds} and
#'   "rev_1" for \code{pepmasses_rev_1.rds}.
#' @param res Lists of data containing theoretical peptides and masses from
#'   \link{readRDS}.
#' @param min_mass Numeric; the minimum MS1 mass.
#' @param max_mass Numeric; the maximum MS1 mass.
#' @param ppm_ms1 Numeric; the error tolerance of MS1 mass in ppm.
#' @param out_path The output path.
#' @examples
#' \donttest{
#' res <- readRDS("~/proteoM/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
#' theopeps <- proteoM:::binTheoSeqs(res)
#' }
#' @return Lists of theoretical peptides binned by MS1 masses. The lists
#'   correspond to the lists of \code{res}.
#' @import parallel
binTheoSeqs <- function (idxes = NULL, res = NULL, min_mass = 500L,
                         max_mass = 10000L, ppm_ms1 = 20L, out_path = NULL) {
  
  if (is.null(res)) {
    stop("`res` cannot be NULL.")
  }
  
  if (is.null(out_path)) {
    stop("`out_path` cannot be NULL.")
  }
  
  if (is.null(idxes)) {
    idxes <- seq_along(res)
  }
  
  out_dir <- create_dir(gsub("(^.*/).*$", "\\1", out_path))
  
  out_nms <- gsub("^.*/(.*)\\.[^\\.].*$", "\\1", out_path) %>%
    paste(idxes, sep = "_") %>%
    paste0(".rds")
  
  res <- res %>%
    lapply(attributes) %>%
    lapply(`[[`, "data")
  
  gc()
  
  n_cores <- local({
    len <- length(res)
    n_cores <- detect_cores()
    
    if (len > n_cores) {
      n_cores <- min(floor(n_cores/2L), len)
    } else {
      n_cores <- min(n_cores, len)
    }
    
    n_cores <- max(1L, n_cores)
  })

  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(
    cl,
    c("bin_theoseqs", 
      "find_ms1_cutpoints"), 
    envir = environment(proteoM:::bin_theoseqs))
  
  out <- parallel::clusterMap(cl, bin_theoseqs, 
                              res, file.path(out_dir, out_nms), 
                              MoreArgs = list(min_mass = min_mass, 
                                              max_mass = max_mass, 
                                              ppm_ms1 = ppm_ms1), 
                              SIMPLIFY = FALSE, USE.NAMES = FALSE, 
                              .scheduling = "dynamic")
  
  parallel::stopCluster(cl)
  rm(list = c("res"))
  gc()
  
  invisible(NULL)
}


#' Finds the cut-points of MS1 masses for binning.
#'
#' The cut-points will be used as the boundary in data binning. Note that the
#' upper bound is open.
#'
#' @param from Numeric; the starting MS1 mass.
#' @param to Numeric; the ending MS1 mass.
#' @param ppm Numeric; the ppm for data binning.
#' @return Cut points.
#' @seealso find_ms1_interval
find_ms1_cutpoints <- function (from = 500L, to = 10000L, ppm = 20L) {

  d <- ppm/1e6
  n <- ceiling(log(to/from)/log(1+d))

  x <- vector("numeric", n)
  x[1] <- from

  for (i in seq_len(n-1)) {
    x[i+1] <- x[i] * (1 + d)
  }

  x
}


#' Helper.
#'
#' Reads single rds.
#'
#' @param file A file name.
#' @param out_path An output path.
s_readRDS <- function (file, out_path) {
  
  ans <- tryCatch(
    readRDS(file = file.path(out_path, file)),
    error = function (e) NULL
  )
  
  if (is.null(ans)) {
    stop("Files found: ", file.path(out_path, file), call. = FALSE)
  }
  
  ans
}

