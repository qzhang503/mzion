#' Helper in binning precursor masses.
#'
#' For target peptides.
#'
#' @param res The results from \link{calc_pepmasses2}.
#' @param min_mass A minimum mass of precursors.
#' @param max_mass A maximum mass of precursors.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @inheritParams calc_pepmasses2
bin_ms1masses <- function (res = NULL, min_mass = 500L, max_mass = 6000L, 
                           ppm_ms1 = 20L, use_ms1_cache = TRUE, 
                           .path_cache = NULL, .path_ms1masses = NULL, 
                           is_ms1_three_frame = TRUE) 
{
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        res <- NULL
        save_call2(path = file.path(.path_cache, "calc_pepmasses2", .time_stamp), 
                   fun = fun, time = .time_bin)
      }
    },
    add = TRUE
  )
  
  ## Initial setups
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  
  ppm_ms1_new <- if (is_ms1_three_frame) 
    as.integer(ceiling(ppm_ms1 * .5))
  else 
    ppm_ms1

  # checks pre-existed precursor masses
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv, inherits = FALSE)
  .path_mass <- file.path(.path_ms1masses, .time_stamp)
  
  masses <- list.files(path = .path_mass, 
                       pattern = paste0("^pepmasses_", "\\d+\\.rds$"))
  
  len_m <- length(masses)
  
  if (!len_m) 
    stop("File not found: ", 
         file.path(.path_mass, paste0("pepmasses_", "[...].rds")))

  # checks pre-existed, binned precursor masses
  .time_bin <- match_calltime(path = file.path(.path_cache, 
                                               "calc_pepmasses2", 
                                               .time_stamp), 
                              fun = fun,
                              nms = c("min_mass", "max_mass", "ppm_ms1")) 
  
  # already binned
  len_bts <- length(.time_bin)
  
  if (len_bts > 1L) 
    stop("More than one cached results found: \n\n", 
         paste(file.path(.path_ms1masses, .time_stamp, fun), collapse = "\n"), 
         "\n\nDelete the caches and start over.", 
         call. = FALSE)
  
  if (len_bts && use_ms1_cache) {
    .path_bin <- file.path(.path_ms1masses, .time_stamp, fun, .time_bin)
    
    bins <- list.files(path = .path_bin, pattern = "binned_theopeps_\\d+\\.rds$")
    len_b <- length(bins)

    if (len_b == len_m) {
      message("Loading bins of MS1 masses from cache.")
      
      .savecall <- FALSE
      
      # no need of global `.time_bin`
      assign(".path_bin", .path_bin, envir = .GlobalEnv)

      return(NULL)
    }
  }
  
  
  # to be binned
  message("Binning MS1 masses...")
  
  .time_bin <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")
  .path_bin <- create_dir(file.path(.path_ms1masses, .time_stamp, fun, .time_bin))

  if (!is.null(res)) {
    # (a) process directly
    binTheoSeqs(idxes = NULL,
                res = res,
                min_mass = min_mass,
                max_mass = max_mass,
                ppm_ms1 = ppm_ms1_new,
                out_path = file.path(.path_bin, "binned_theopeps.rds"))
  } else {
    # (b) reload
    idxes <- local({
      idxes <- gsub("^pepmasses_(\\d+)\\.rds$", "\\1", masses)
      idxes <- as.integer(idxes)
      idxes <- idxes[order(idxes)]
    })
    
    n_cores <- local({
      fct <- 20 
      free_mem <- find_free_mem()
      max_sz <- max(file.size(file.path(.path_mass, masses)))/1024^2
      
      n_cores <- min(floor(free_mem/max_sz/fct), detect_cores(8L))

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
      
      # No need of flatten() as saveRDS by INDIVIDUAL idx (and return NULL)
      parallel::clusterApplyLB(
        cl = cl, 
        x = chunksplit(idxes, n_cores, "list"), 
        fun = lapply, 
        FUN = "binTheoSeqs_i", 
        min_mass = min_mass, 
        max_mass = max_mass, 
        ppm_ms1 = ppm_ms1_new, 
        in_path = .path_mass,
        out_path = .path_bin
      )
      
      parallel::stopCluster(cl)
    } else {
      lapply(idxes, binTheoSeqs_i, min_mass, max_mass, ppm_ms1_new, 
             .path_mass, .path_bin)
    }
  }
  
  .savecall <- TRUE

  assign(".time_bin", .time_bin, envir = .GlobalEnv)
  assign(".path_bin", .path_bin, envir = .GlobalEnv)

  invisible(NULL)
}


#' Helper of \link{binTheoSeqs2}.
#' 
#' @param in_path An input path of \code{pepmasses_}.
#' @inheritParams binTheoSeqs2
binTheoSeqs_i <- function (idx = 1L, min_mass = 500L,max_mass = 6000L,
                           ppm_ms1 = 20L, in_path = NULL, out_path = NULL) 
{
  if (is.null(in_path)) 
    stop("`in_path` cannot be NULL.")
  
  message("\tSet: ", idx)
  
  in_nm <- paste0("pepmasses_", idx, ".rds")
  res <- s_readRDS(in_nm, in_path)
  
  binTheoSeqs2(idx = idx,
               res = res,
               min_mass = min_mass,
               max_mass = max_mass,
               ppm_ms1 = ppm_ms1,
               out_path = out_path)
}


#' Separates theoretical peptides into mass groups.
#'
#' @param idx An index, e.g. "1" for \code{pepmasses_1.rds} and "rev_1" for
#'   \code{pepmasses_rev_1.rds}.
#' @inheritParams binTheoSeqs
binTheoSeqs2 <- function (idx = 1L, res = NULL, min_mass = 500L,
                          max_mass = 6000L, ppm_ms1 = 20L, 
                          out_path = NULL) 
{
  if (is.null(res)) 
    stop("`res` cannot be NULL.")
  
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.")
  
  out_dir <- create_dir(gsub("(^.*/).*$", "\\1", out_path))
  out_nm <- paste0(paste0("binned_theopeps_", idx), ".rds")
  
  res <- attr(res, "data")
  gc()
  
  bin_theoseqs(res, file.path(out_path, out_nm), min_mass, max_mass, ppm_ms1)
  
  rm(list = c("res"))
  gc()
  
  invisible(NULL)
}


#' Separates theoretical peptides into mass groups.
#'
#' @param peps A list of theoretical peptides with masses.
#' @param out_nm A output name with prepending file path.
#' @inheritParams binTheoSeqs
bin_theoseqs <- function (peps = NULL, out_nm = NULL, min_mass = 500L, 
                          max_mass = 6000L, ppm_ms1 = 20L) 
{
  if (!length(peps)) {
    # out <- data.frame(pep_seq = character(), mass = numeric(), frame = integer(), row.names = NULL)
    out <- NULL                
    saveRDS(out, out_nm)
    
    return(NULL)
  }

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
                         max_mass = 6000L, ppm_ms1 = 20L, out_path = NULL) 
{
  if (is.null(res)) 
    stop("`res` cannot be NULL.")
  
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.")
  
  if (is.null(idxes)) 
    idxes <- seq_along(res)
  
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
    n_cores <- detect_cores(16L)
    
    if (len > n_cores) 
      n_cores <- min(floor(n_cores/2L), len)
    else 
      n_cores <- min(n_cores, len)

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
find_ms1_cutpoints <- function (from = 500L, to = 10000L, ppm = 20L) 
{
  d <- ppm/1e6
  n <- ceiling(log(to/from)/log(1+d))

  x <- vector("numeric", n)
  x[1] <- from

  for (i in seq_len(n-1)) x[i+1] <- x[i] * (1 + d)

  x
}


#' Helper.
#'
#' Reads single rds.
#'
#' @param file A file name.
#' @param out_path An output path.
s_readRDS <- function (file, out_path) 
{
  ans <- tryCatch(
    readRDS(file = file.path(out_path, file)),
    error = function (e) NULL
  )
  
  if (is.null(ans)) 
    stop("Files found: ", file.path(out_path, file), call. = FALSE)
  
  ans
}


