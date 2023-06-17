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
bin_ms1masses <- function (res = NULL, min_mass = 200L, max_mass = 4500L, 
                           min_len = 7L, max_len = 40L, ppm_ms1 = 20L, 
                           use_ms1_cache = TRUE, .path_cache = NULL, 
                           .path_ms1masses = NULL, is_ms1_three_frame = TRUE, 
                           out_path = NULL, enzyme = "trypsin_p") 
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
  
  ppm_ms1_bin <- calc_threeframe_ppm(ppm_ms1)

  # checks pre-existed precursor masses
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv, inherits = FALSE)
  .path_mass  <- file.path(.path_ms1masses, .time_stamp)
  masses <- list.files(path = .path_mass, pattern = paste0("^pepmasses_", "\\d+\\.rds$"))
  
  if (!(len_m  <- length(masses))) 
    stop("File not found: ", 
         file.path(.path_mass, paste0("pepmasses_", "[...].rds")))

  # checks pre-existed, binned precursor masses
  .time_bin <- match_calltime(path = file.path(.path_cache, 
                                               "calc_pepmasses2", 
                                               .time_stamp), 
                              fun = fun,
                              nms = c("min_mass", "max_mass", "min_len", 
                                      "max_len", "ppm_ms1")) 
  
  # already binned
  if ((len_bts <- length(.time_bin)) > 1L) 
    stop("More than one cached results found: \n\n", 
         paste(file.path(.path_ms1masses, .time_stamp, fun), collapse = "\n"), 
         "\n\nDelete the caches and start over.")
  
  if (len_bts && use_ms1_cache) {
    .path_bin <- file.path(.path_ms1masses, .time_stamp, fun, .time_bin)
    bins <- list.files(path = .path_bin, pattern = "binned_theopeps_\\d+\\.rds$")

    if (length(bins) == len_m) {
      message("Loading bins of MS1 masses from cache.")
      .savecall <- FALSE
      return(.path_bin)
    }
  }
  
  message("Binning MS1 masses...")
  
  .time_bin <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")
  .path_bin <- create_dir(file.path(.path_ms1masses, .time_stamp, fun, .time_bin))
  
  if (is.null(res)) {
    # (b) reload
    idxes <- local({
      idxes <- gsub("^pepmasses_(\\d+)\\.rds$", "\\1", masses)
      idxes <- as.integer(idxes)
      idxes <- idxes[order(idxes)]
    })
    
    n_cores <- set_bin_ncores(len_m, enzyme)

    if (n_cores > 1L) {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(cl, list("qread", "qsave"), envir = environment(qs::qsave))
      
      parallel::clusterExport(
        cl,
        c("binTheoSeqs_i", 
          "binTheoSeqs2", 
          "bin_theoseqs", 
          "s_readRDS", 
          "find_ms1_cutpoints"), 
        envir = environment(mzion::matchMS))
      
      # No need of flatten() as saveRDS by INDIVIDUAL idx (and return NULL)
      parallel::clusterApplyLB(
        cl = cl, 
        x = chunksplit(idxes, n_cores, "list"), 
        fun = lapply, 
        FUN = "binTheoSeqs_i", 
        min_mass = min_mass, 
        max_mass = max_mass, 
        ppm_ms1 = ppm_ms1_bin, 
        in_path = .path_mass,
        out_path = .path_bin
      )
      
      parallel::stopCluster(cl)
    } 
    else
      lapply(idxes, binTheoSeqs_i, min_mass, max_mass, ppm_ms1_bin, 
             .path_mass, .path_bin)
  }
  else {
    # (a) process directly
    binTheoSeqs(idxes = NULL,
                res = res,
                min_mass = min_mass,
                max_mass = max_mass,
                ppm_ms1 = ppm_ms1_bin,
                enzyme = enzyme, 
                out_path = file.path(.path_bin, "binned_theopeps.rds"))
  }
  
  pat_b <- "^binned_theopeps_[0-9]+\\.rds"
  len_b <- length(list.files(.path_bin, pattern = pat_b))
  
  if (len_b != len_m)
    stop("May need more RAM: expect ", len_m, " \"", pat_b, "\" files, ", 
         "but found ", len_b, " files under \n\"", .path_bin, "\"\n")
  
  .savecall <- TRUE

  local({
    file <- file.path(out_path, "Calls", ".cache_info.rds")
    
    if (file.exists(file))
      .cache_info <- qs::qread(file)

    .cache_info$.time_bin <- .time_bin
    .cache_info$.path_bin <- .path_bin
    
    qs::qsave(.cache_info, file, preset = "fast")
  })
  
  message("Completed precusor bins at: ", Sys.time())
  
  invisible(.path_bin)
}


#' Helper of \link{binTheoSeqs2}.
#' 
#' @param in_path An input path of \code{pepmasses_}.
#' @inheritParams binTheoSeqs2
binTheoSeqs_i <- function (idx = 1L, min_mass = 200L, max_mass = 4500L,
                           ppm_ms1 = 10L, in_path = NULL, out_path = NULL) 
{
  if (is.null(in_path)) 
    stop("`in_path` cannot be NULL.")
  
  message("\tSet: ", idx)
  
  res <- s_readRDS(paste0("pepmasses_", idx, ".rds"), in_path)
  
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
binTheoSeqs2 <- function (idx = 1L, res = NULL, min_mass = 200L,
                          max_mass = 4500L, ppm_ms1 = 10L, out_path = NULL) 
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
bin_theoseqs <- function (peps = NULL, out_nm = NULL, min_mass = 200L, 
                          max_mass = 4500L, ppm_ms1 = 10L) 
{
  if (!length(peps)) {
    out <- NULL   
    qs::qsave(out, out_nm, preset = "fast")
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
  out <- lapply(out, function (x) { x[["frame"]] <- NULL; x })
  qs::qsave(out, out_nm, preset = "fast")

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
#' @param ppm_ms1 Numeric; (half of) the error tolerance of MS1 mass in ppm.
#' @param out_path The output path.
#' @param enzyme The assume enzyme activity.
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' # res <- readRDS("~/mzion/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
#' # theopeps <- mzion:::binTheoSeqs(res)
#' }
#' @return Lists of theoretical peptides binned by MS1 masses. The lists
#'   correspond to the lists of \code{res}.
binTheoSeqs <- function (idxes = NULL, res = NULL, min_mass = 200L,
                         max_mass = 4500L, ppm_ms1 = 10L, enzyme = "trypsin_p", 
                         out_path = NULL) 
{
  if (is.null(res)) 
    stop("`res` cannot be NULL.")
  
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.")
  
  if (is.null(idxes)) 
    idxes <- seq_along(res)
  
  out_dir <- create_dir(gsub("(^.*/).*$", "\\1", out_path))
  
  out_nms <- gsub("^.*/(.*)\\.[^\\.].*$", "\\1", out_path)
  out_nms <- paste0(out_nms, "_", idxes, ".rds")
  
  res <- lapply(res, attributes)
  res <- lapply(res, `[[`, "data")
  gc()
  
  n_cores <- set_bin_ncores(length(res), enzyme)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("qread", "qsave"), 
                          envir = environment(qs::qsave))
  parallel::clusterExport(cl, c("bin_theoseqs", "find_ms1_cutpoints"), 
                          envir = environment(mzion::matchMS))
  
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
find_ms1_cutpoints <- function (from = 200L, to = 4500L, ppm = 10L) 
{
  d <- ppm/1e6
  n <- ceiling(log(to/from)/log(1+d))

  x <- vector("numeric", n)
  x[1] <- from

  for (i in seq_len(n-1)) 
    x[i+1] <- x[i] * (1 + d)

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
    qs::qread(file = file.path(out_path, file)),
    error = function (e) NULL
  )
  
  if (is.null(ans)) 
    stop("Files not found: ", file.path(out_path, file))
  
  ans
}


#' Sets the number of CPU cores for precursor mass binning.
#' 
#' @param len_m The number of \code{aa_masses} modules.
#' @param enzyme The assume enzymatic activity.
set_bin_ncores <- function (len_m, enzyme)
{
  n_cores <- detect_cores(15L)
  
  n_cores <- if (len_m > n_cores) 
    min(floor(n_cores/2L), len_m)
  else 
    min(n_cores, len_m)
  
  if (enzyme == "noenzyme")
    n_cores <- floor(n_cores/2L)
  
  max(1L, n_cores)
}


