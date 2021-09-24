#' Helper in binning precursor masses
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
                   res = res$fwd,
                   min_mass = min_mass,
                   max_mass = max_mass,
                   ppm_ms1 = ppm_ms1,
                   out_path = out_path)

  invisible(NULL)
}


#' Helper of \link{bin_ms1masses} by the separate sets of \code{target}s and
#' \code{decoy}s.
#'
#' @param bins The names of bins in the format of
#'   "binned_theopeps_(rev_)*\\d+\\.rds$".
#' @param type The type of data: "target" or "decoy".
#' @inheritParams binTheoPeps
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

  # (a) small data and process directly
  if (!is.null(res)) {
    binTheoPeps(idxes = NULL,
                res = res,
                min_mass = min_mass,
                max_mass = max_mass,
                ppm_ms1 = ppm_ms1,
                out_path = file.path(out_path, "binned_theopeps.rds"))

    return(NULL)
  }

  # (b) large data set; `res` loaded from disk
  idxes <- local({
    idxes <- gsub("^(pepmasses_|pepmasses_rev_)(\\d+)\\.rds$", "\\2", masses)
    idxes <- as.integer(idxes)
    idxes <- idxes[order(idxes)]
    idxes <- paste0(type2, idxes)
  })

  n_cores <- local({
    max_n_cores <- 8L
    
    fct <- 26 
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
      c("binTheoPeps_i", 
        "binTheoPeps2", 
        "bin_theopeps", 
        "find_ms1_cutpoints", 
        "cbind_theopepes"), 
      envir = environment(proteoM:::binTheoPeps_i))
    
    # No need of purrr::flatten() as saveRDS by INDIVIDUAL idx (and return NULL)
    parallel::clusterApplyLB(
      cl = cl, 
      x = chunksplit(idxes, n_cores, "list"), 
      fun = lapply, 
      FUN = "binTheoPeps_i", 
      min_mass = min_mass, 
      max_mass = max_mass, 
      ppm_ms1 = ppm_ms1, 
      out_path = out_path
    )

    parallel::stopCluster(cl)
  } else {
    lapply(idxes, binTheoPeps_i, min_mass, max_mass, ppm_ms1, out_path)
  }

  invisible(NULL)
}


#' Helper of \link{binTheoPeps2}.
#'
#' @inheritParams binTheoPeps2
binTheoPeps_i <- function (idx = 1L, min_mass = 500L,max_mass = 10000L,
                           ppm_ms1 = 20L, out_path = NULL) {

  message("\tSet: ", idx)

  in_nm <- paste0("pepmasses_", idx, ".rds")
  res <- s_readRDS(in_nm, out_path)

  binTheoPeps2(idx = idx,
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
#' @inheritParams binTheoPeps
binTheoPeps2 <- function (idx = 1L, res = NULL, min_mass = 500L,
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

  out <- bin_theopeps(res, min_mass, max_mass, ppm_ms1)
  rm(list = c("res"))
  gc()

  out <- cbind_theopepes(out, file.path(out_dir, out_nm))

  invisible(NULL)
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


#' Helper of \link{bin_ms1masses_td}.
#'
#' Not used. Bin peptide masses by sets of data.
#'
#' @inheritParams binTheoPeps
binTheoPeps_bysets <- function (idxes = NULL, min_mass = 500L, max_mass = 10000L,
                                ppm_ms1 = 20L, out_path = NULL) {

  message("\tSets: ", paste(idxes, collapse = ", "), "...")

  in_nms <- paste0("pepmasses_", idxes, ".rds")

  n_cores <- min(detect_cores(), length(idxes))
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("s_readRDS"),
                          envir = environment(proteoM:::s_readRDS))
  
  res <- parallel::clusterApply(cl, in_nms, s_readRDS, out_path)
  
  parallel::stopCluster(cl)
  gc()

  binTheoPeps(idxes = idxes,
              res = res,
              min_mass = min_mass,
              max_mass = max_mass,
              ppm_ms1 = ppm_ms1,
              out_path = file.path(out_path, "binned_theopeps.rds"))
}


#' Separates theoretical peptides into mass groups.
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
#' theopeps <- proteoM:::binTheoPeps(res)
#' }
#' @return Lists of theoretical peptides binned by MS1 masses. The lists
#'   correspond to the lists of \code{res}.
#' @import parallel
binTheoPeps <- function (idxes = NULL, res = NULL, min_mass = 500L,
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

  n_cores <- detect_cores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(
    cl,
    c("bin_theopeps", 
      "find_ms1_cutpoints", 
      "cbind_theopepes"), 
    envir = environment(proteoM:::bin_theopeps))
  
  out <- lapply(res, function (x) {
    x <- chunksplit(x, n_cores)
    
    parallel::clusterApply(cl, x, bin_theopeps, min_mass, max_mass, ppm_ms1) %>%
      purrr::flatten()
  })
  
  parallel::stopCluster(cl)
  rm(list = c("res"))
  gc()

  ## Split by frames (involves memory-demanding `split`)
  # (and note that `res` has to be by paralleled across lists)
  n_cores <- min(detect_cores(), length(out), floor(memory.limit()/memory.size()))

  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out <- parallel::clusterMap(cl, cbind_theopepes, out, file.path(out_dir, out_nms))
    parallel::stopCluster(cl)
    gc()
  } else if (n_cores == 1L) {
    for (i in seq_along(out)) {
      out[[i]] <- cbind_theopepes(out[[i]], file.path(out_dir, out_nms))
      gc()
    }
  } else {
    stop("Not enough memory.")
  }

  invisible(NULL)
}


#' Helper: separates theoretical peptides into mass groups.
#'
#' @param peps A list of theoretical peptides with masses.
#' @inheritParams binTheoPeps
bin_theopeps <- function (peps = NULL, min_mass = 500L, max_mass = 10000L,
                          ppm_ms1 = 20L) {

  ps <- find_ms1_cutpoints(min_mass, max_mass, ppm_ms1)
  
  mapply(function (prot_peps, prot_accs, ps) {
    frames <- findInterval(prot_peps, ps)
    
    list(pep_seq = names(prot_peps),
         mass = prot_peps,
         frames = frames,
         prot_acc = prot_accs)
  }, peps, names(peps), MoreArgs = list(ps = ps), SIMPLIFY = FALSE)
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


#' Combines theoretical peptides after binning.
#'
#' @param out A list of binned theoretical peptides.
#' @param out_nm The output file path and name.
#' @importFrom dplyr arrange
cbind_theopepes <- function (out, out_nm) {
  
  prot_acc <- mapply(function (x, y) rep(y, length(x$pep_seq)), out, names(out), 
                     SIMPLIFY = FALSE) %>%
    do.call(`c`, .) %>%
    unname()
  
  pep_seq <- lapply(out, `[[`, "pep_seq") %>%
    do.call(`c`, .) %>%
    unname()

  mass <- lapply(out, `[[`, "mass") %>%
    do.call(`c`, .) %>%
    unname()

  frame <- lapply(out, `[[`, "frames") %>%
    do.call(`c`, .) %>%
    unname()

  out <- data.frame(pep_seq = pep_seq,
                    mass = mass,
                    frame = frame,
                    prot_acc = prot_acc) %>%
    dplyr::arrange(frame, pep_seq, prot_acc) %>%
    split(., .$frame, drop = FALSE)

  saveRDS(out, out_nm)

  rm(list = c("out", "prot_acc", "pep_seq", "mass", "frame"))
  gc()

  invisible(NULL)
}

