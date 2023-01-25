#' Finds the indexes of top-n entries without re-ordering.
#'
#' At length(x) >= n, the length of output may be shorter than n with ties.
#'
#' @param x A numeric vector.
#' @param n The number of top entries to keep.
#' @param ... Additional arguments for base function sort.
#' @return The indexes of the top-n entries.
#' @examples
#' \donttest{
#' which_topx(c(1:5), 50)
#'
#' length(which_topx(sample(5000, 500), 100))
#'
#' length(which_topx(sample(100, 100, replace = TRUE), 100))
#' }
which_topx <- function(x, n = 50L, ...) 
{
  len <- length(x)
  p <- len - n
  
  if (p  <= 0L) 
    return(seq_along(x))
  
  xp <- sort(x, partial = p, ...)[p]
  
  which(x > xp)
}


#' Finds the indexes of top-n entries without re-ordering.
#' 
#' Handles ties.
#' 
#' @inheritParams which_topx
#' @return The indexes of the top-n entries.
#' @examples 
#' p <- 100
#' set.seed(1)
#' x <- sample(1:150, replace = T)
#' 
#' # 103, not 100
#' xp <- sort(x, partial = p)[p]
#' 
#' # multiple ties
#' x <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
#' p <- 2L
#' 
#' x <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5)
#' set.seed <- (3)
#' x <- sample(x)
#' p <- 4L
#' 
#' ## more NA values (3) than n (2)
#' x <- c(0.11, 0.11, NA, rep(0.11, 7), NA, NA)
#' ans <- which_topx2(x, 2L, na.last = FALSE)
#' 
#' # OK
#' ans <- which_topx2(x, 8L, na.last = FALSE)
#' 
#' # Bad
#' ans <- which_topx2(x, 9L, na.last = FALSE)
#' # OK
#' ans <- which_topx2(x, 9L, na.last = TRUE)
#'  
#' \donttest{
#' which_topx2(5000, NA_integer_)
#' }
which_topx2 <- function(x, n = 50L, ...) 
{
  if (is.na(n)) 
    return(NULL)
  
  len <- length(x)
  p <- len - n

  if (p  <= 0L)
    return(seq_along(x))

  # not yet handle xp is NA
  # or:  xp <- sort(x, partial = p, na.last = TRUE)[p]
  xp <- sort(x, partial = p, ...)[p]
  
  inds <- which(x > xp)
  
  # in case of ties -> length(inds) < n
  # detrimental e.g. ms2_n = 500 and n = 100
  #   -> expect 100 `ms2_moverzs` guaranteed but may be only 99
  #
  # MGF `ms2_moverzs` is increasing
  # `inds2` goes first to ensure non-decreasing index for `ms2_moverzs`
  
  d <- n - length(inds)
  
  if (d) {
    # must exist and length(ties) >= length(d)
    ties <- which(x == xp)
    for (i in seq_len(d)) inds <- insVal(ties[i], inds)
  }
  
  invisible(inds)
}


#' Gets top-n values
#' 
#' @param vals a vector of values
#' @param n The top-n values to be retained
get_topn_vals <- function (vals, n) vals[which_topx2(vals, n)]


#' Inserts a value.
#' 
#' @param x A value to be inserted.
#' @param sv A sorted vector of indexes(thus without ties).
#' 
#' @examples
#' sv <- c(4, 10)
#' 
#' sv <- insVal(2, sv)
#' sv <- insVal(6, sv)
#' sv <- insVal(20, sv)
insVal <- function (x, sv) 
{
  grs <- sv > x
  c(sv[!grs], x, sv[grs])
}


#' Finds the top-n entries without re-ordering.
#'
#' @inheritParams which_topx
#' @return The top-n entries.
topx <- function(x, n = 50L, ...) 
{
  len <- length(x)
  p <- len - n
  
  if (p  <= 0L) 
    return(x)
  
  xp <- sort(x, partial = p, ...)[p]
  
  x[x > xp]
}


#' Finds the numeric difference in ppm.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The difference between \eqn{x} and \eqn{y} in ppm.
find_ppm_error <- function (x = 1000, y = 1000.01) (y - x)/y * 1E6


#' Finds the error range of a number.
#'
#' Assumes \eqn{x} is positive without checking.
#'
#' @param x A numeric value.
#' @param ppm Numeric; the ppm allowed from \code{x}.
#' @return The lower and the upper bound to \eqn{x} by \eqn{ppm}.
find_mass_error_range <- function (x = 500L, ppm = 20L) 
{
  d <- x * ppm/1E6
  c(x - d, x + d)
}


#' Sums elements across lists.
#'
#' Each list has the same length. NA values are removed.
#'
#' @param x A numeric value.
#' @param y A numeric value.
`%+%` <- function(x, y) mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))


#' Post processing after ms2match.
#'
#' @param df An output from various ms2match(es).
#' @inheritParams ms2match_base
post_ms2match <- function (df, i, aa_masses, out_path) 
{
  create_dir(file.path(out_path, "temp"))
  
  if (is.null(df)) {
    qs::qsave(df, file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")), 
              preset = "fast")
    return(NULL)
  }

  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
  is_decoy <- if (grepl("^rev", i)) TRUE else FALSE

  df <- df %>%
    dplyr::mutate(pep_fmod = nm_fmods,
                  pep_vmod = nm_vmods,
                  pep_mod_group = as.character(i), 
                  scan_num = as.character(scan_num)) %>%
    { if (is_decoy) dplyr::mutate(., pep_isdecoy = TRUE) else
      dplyr::mutate(., pep_isdecoy = FALSE) }

  df %>%
    reloc_col_after("raw_file", "scan_num") %>%
    reloc_col_after("pep_mod_group", "raw_file") %T>%
    qs::qsave(file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")), 
              preset = "fast")
}


#' Post frame advancing.
#'
#' Flattens mgfs within each frame (the number of entries equals to the number
#' of mgfs).
#' 
#' @param res Results from frame-advanced searches.
#' @param mgf_frames Data of MGF frames.
post_frame_adv <- function (res, mgf_frames) 
{
  res <- unlist(res, recursive = FALSE)

  lens <- lapply(res, length)
  lens <- unlist(lens, recursive = FALSE, use.names = FALSE)
  empties <- !lens

  out <- do.call(rbind, mgf_frames)
  out <- dplyr::mutate(out, matches = res)
  
  out[!empties, ]
}


#' Subsets the search space.
#' 
#' @param n_cores The number of CPU cores.
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' @inheritParams post_ms2match
purge_search_space <- function (i, aa_masses, mgf_path, n_cores, ppm_ms1 = 10L,
                                fmods_nl = NULL) 
{
  # loads freshly mgfs (as will be modified)
  mgf_frames <- 
    qs::qread(file.path(mgf_path, "mgf_queries.rds")) %>% 
    dplyr::group_by(frame) %>%
    dplyr::group_split() %>%
    setNames(purrr::map_dbl(., function (x) x$frame[1]))

  mgf_frames <- local({
    ranges <- seq_along(mgf_frames)
    labs   <- levels(cut(ranges, n_cores^2))
    lower  <- floor(as.numeric( sub("\\((.+),.*", "\\1", labs)))
    grps   <- findInterval(ranges, lower)

    split(mgf_frames, grps)
  })

  # parses aa_masses
  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
  msg_end  <- if (grepl("^rev_", i)) " (decoy)." else "."

  message("Matching against: ",
          paste0(nm_fmods,
                 nm_vmods %>% { if (nchar(.) > 0L) paste0(" | ", .) else . },
                 msg_end))

  # reads theoretical peptide data
  .path_bin <- get(".path_bin", envir = .GlobalEnv, inherits = FALSE)
  theopeps  <- qs::qread(file.path(.path_bin, paste0("binned_theopeps_", i, ".rds")))
  
  if (is.null(theopeps)) {
    return(list(mgf_frames = mgf_frames, 
                theopeps = NULL, 
                theopeps2 = NULL))
  }
  
  theopeps <- lapply(theopeps, function (x) x[, c("pep_seq", "mass")])

  # (1) for a given aa_masses_all[[i]], some mgf_frames[[i]]
  #     may not be found in theopeps[[i]]
  frames_theo <- names(theopeps)
  
  mgf_frames <- lapply(mgf_frames, function (x) {
    oks <- names(x) %in% frames_theo
    x <- x[oks]
    
    empties <- purrr::map_lgl(x, purrr::is_empty)
    x[!empties]
  })
  
  rm(list = "frames_theo")
  
  # (2) splits `theopeps` in accordance to `mgf_frames` with
  #     preceding and following frames: (o)|range of mgf_frames[[1]]|(o)
  frames_mgf <- lapply(mgf_frames, function (x) as.integer(names(x)))
  
  mins <- purrr::map_int(frames_mgf, function (x) {
    if (length(x)) min(x, na.rm = TRUE) else 0L
  })
  
  maxs <- purrr::map_int(frames_mgf, function (x) {
    if (length(x)) max(x, na.rm = TRUE) else 0L
  })
  
  frames_theo <- as.integer(names(theopeps))
  
  # separates into intervals
  # (NAMEs are trivial node indexes from parallel processes)
  
  # DON'T: may cause uneven length between `mgf_frames` and `theopeps`
  # empties <- map_lgl(theopeps, is_empty)
  # theopeps[!empties]
  
  theopeps <- mapply(function (x, y) {
    theopeps[which(frames_theo >= (x - 1L) & frames_theo <= (y + 1L))]
  }, mins, maxs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  # (3) removes unused frames of `theopeps`
  # (NAMEs are trivial cluster node indexes...)
  theopeps <- mapply(subset_theoframes, mgf_frames, theopeps, 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)

  # (4) removes empties (zero overlap between mgf_frames and theopeps)
  oks <- purrr::map_lgl(mgf_frames, function (x) length(x) > 0L) |
    purrr::map_lgl(theopeps, function (x) length(x) > 0L)

  mgf_frames <- mgf_frames[oks]
  theopeps <- theopeps[oks]

  # (6) reverses the order (longer/heavier peptides towards the beginning)
  #     do the difficult ones first when paralleling with LB
  seqs <- rev(seq_along(theopeps))
  mgf_frames <- mgf_frames[seqs]
  theopeps <- theopeps[seqs]
  
  if (FALSE) {
    theopeps2 <- local({
      .path_bin2 <- gsub("ms1masses/", "ms2masses/", .path_bin)
      file <- file.path(.path_bin2, paste0("binned_ms2_", i, ".rds"))
      
      if (!file.exists(file))
        return(NULL)
      
      tps2 <- qs::qread(file)
      
      frames <- lapply(theopeps, names)
      frames <- lapply(frames, function (x) x[!is.na(x)])
      
      tps2 <- tps2[names(tps2) %in% unlist(frames)]
      
      # split by frames in accordance with theopeps
      tps2 <- lapply(frames, function (x) tps2[names(tps2) %in% x])
      
      # lapply(tps2, hash_frame_nums)
    })
  }

  invisible(list(mgf_frames = mgf_frames, theopeps = theopeps))
}


#' Subsets the frames of theoretical peptides.
#' 
#' @param mgf_frames MGFs in frames. Each frame contains one to multiple MGFs
#'   whose MS1 masses are in the same interval.
#' @param theopeps Binned theoretical peptides at a given combination of fixed
#'   and variable.
subset_theoframes <- function (mgf_frames = NULL, theopeps = NULL) 
{
  if (!(length(mgf_frames) && length(theopeps))) 
    return(NULL)

  frames <- as.integer(names(mgf_frames))
  breaks <- which(diff(frames) != 1L) + 1L
  grps   <- findInterval(frames, frames[breaks])
  frames <- split(frames, grps)
  
  frames <- lapply(frames, function (x) c(x[1] - 1, x, x[length(x)] + 1))
  frames <- unlist(frames, recursive = FALSE, use.names = FALSE)
  frames <- frames[!duplicated(frames)]
  
  theopeps[as.character(frames)]
}


#' Subsets theoretical peptides.
#'
#' Only entries containing the site of neuloss will be kept.
#'
#' @param pattern A regex of amino-acid residue(s).
#' @param theopeps Lists of theoretical peptides. A column of \code{pep_seq} is
#'   assumed.
subset_neuloss_peps <- function (pattern, theopeps) 
{
  rows <- lapply(theopeps, function (x) grepl(pattern, x$pep_seq))
  
  mapply(function (x, y) x[y, ], theopeps, rows, 
         SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


#' Finds MS2 N-terminal mass.
#'
#' When \code{length(ntmod) > 0}, \code{aa_mass["N-term"]} must be H given the
#' rule of no additive (terminal) modifications.
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
find_nterm_mass <- function (aa_masses) 
{
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  
  # (1) For the "if" part at "length(ntmod) > 0": 
  #   aa_mass["N-term"] must be hydrogen, provided the 
  #   rule of no additive (terminal) modifications. Thus, 
  #   aa_masses[names(ntmod)] + aa_masses["N-term"] - e is identical to
  #   aa_masses[names(ntmod)] + 1.00727647
  # 
  # (2) For the "else" part: 
  #   e.g. at fixed `TMT6plex (N-term)`
  #   aa_masses["N-term"] = 229 + H -> 230
  
  # hydrogen <- 1.007825
  # proton <- 1.00727647
  # electr <- 0.000549
  
  if (length(ntmod))
    aa_masses[names(ntmod)] + 1.00727647
  else
    aa_masses["N-term"] - 0.000549
}


#' Finds MS2 C-terminal mass.
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
find_cterm_mass <- function (aa_masses) 
{
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  # (1) For the "if" part at "length(ctmod) > 0": 
  #   aa_mass["C-term"] must be OH, provided the 
  #   rule of no additive (terminal) modifications. Thus, 
  #   aa_masses[names(ctmod)] + aa_masses["C-term"] + (H) + (H+) is identical to
  #   aa_masses[names(ctmod)] + 19

  if (length(ctmod))
    aa_masses[names(ctmod)] + 19.0178415
  else 
    aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
}


#' Right-joining of two data frames.
#'
#' Rows ordered by \code{y}, which is different to \link[dplyr]{right_join}.
#'
#' @param x The left data frame (to be proliferated by rows).
#' @param y The right data frame (the dominated one).
#' @param by The key.
#' @examples
#' \donttest{
#' df1 <- data.frame(A = c("a", "b", "c"), B = c(1, 1, 1))
#' df2 <- data.frame(A = c("a", "c", "d"), C = c(2, 2, "3"))
#'
#' x1 <- quick_rightjoin(df1, df2, by = "A")
#' x1 <- x1[, c("A", "B", "C")]
#' rownames(x1) <- seq_len(nrow(x1))
#'
#' x2 <- dplyr::right_join(df1, df2, by = "A")
#' x2 <- x2[, c("A", "B", "C")]
#'
#' stopifnot(identical(x1, x2))
#'
#' # row order may be different
#' x1 <- quick_rightjoin(df2, df1, by = "A")
#' x1 <- x1[, c("A", "B", "C")]
#' rownames(x1) <- seq_len(nrow(x1))
#'
#' x2 <- dplyr::right_join(df2, df1, by = "A")
#' x2 <- x2[, c("A", "B", "C")]
#'
#' # FALSE
#' identical(x1, x2)
#' }
quick_rightjoin <- function (x, y, by = NULL) 
{
  # the indexes of y in x;
  # NA rows in "x", if "y" not found in "x"

  rows <- match(y[[by]], x[[by]])
  x <- x[rows, ]
  x[[by]] <- NULL

  ## faster, but not for data.table
  # x <- x[, -which(names(x) == by), drop = FALSE]

  cbind2(x, y)
}


#' Left-joining of two data frames.
#'
#' @param x The left data frame.
#' @param y The right data frame.
#' @param by The key.
#' @examples
#' \donttest{
#' df1 <- data.frame(A = c("a", "b", "c"), B = c(1, 1, 1))
#' df2 <- data.frame(A = c("a", "c", "d"), C = c(2, 2, "3"))
#'
#' x1 <- quick_leftjoin(df1, df2, by = "A")
#' x1 <- x1[, c("A", "B", "C")]
#' rownames(x1) <- seq_len(nrow(x1))
#'
#' x2 <- dplyr::left_join(df1, df2, by = "A")
#'
#' stopifnot(identical(x1, x2))
#' }
quick_leftjoin <- function (x, y, by = NULL) 
{
  rows <- match(x[[by]], y[[by]])
  y <- y[rows, ]
  y[[by]] <- NULL

  ## faster, but not for data.table
  # y <- y[, -which(names(y) == by), drop = FALSE]

  cbind2(x, y)
}


#' Detects and suggests the number of CPU cores.
#'
#' @param max_n_cores The maximum number of cores for uses.
detect_cores <- function (max_n_cores = NULL) 
{
  n_cores <- parallel::detectCores()

  max_n_cores <- if (is.null(max_n_cores)) 
    n_cores
  else 
    min(max_n_cores, n_cores)

  max_n_cores <- if (n_cores > 128L) 
    min(max_n_cores, n_cores - 8L)
  else if (n_cores <= 128L && n_cores > 64L) 
    min(max_n_cores, n_cores - 4L)
  else if (n_cores <= 64L && n_cores > 32L) 
    min(max_n_cores, n_cores - 2L)
  else if (n_cores <= 32L && n_cores > 16L) 
    min(max_n_cores, n_cores - 1L)
  else 
    min(max_n_cores, n_cores)
}


#' Finds the amount of free system memory.
#' 
#' Outputs of free RAM in the unit of MB.
#' 
#' @param sys_ram The amount of system RAM; only for uses with Linux or Mas OS.
find_free_mem <- function (sys_ram = 32L) 
{
  nm_os <- Sys.info()['sysname']
  
  gc()
  
  if (nm_os == "Windows") {
    free_mem <- system('wmic OS get FreePhysicalMemory /Value', intern=TRUE)[3]
    free_mem <- gsub("^FreePhysicalMemory=(\\d+)\\r", "\\1", free_mem)
    free_mem <- as.numeric(free_mem)/1024
  } 
  else {
    # not yet tested for "Linux", "Darwin"
    # inaccurate e.g. if physical RAM is 32000 but in .RProfile 
    # `invisible(utils::memory.limit(64000))`
    # memory.limit() - memory.size(max = TRUE)
    
    warning("Cannot determine the amount of RAM with Linux or MAC OS.\n", 
            "To specify, use parameter \"sys_ram\".")
    free_mem <- sys_ram * .75
  }
  
  free_mem
}


#' Find the indexes of modifications.
#' 
#' @param file A full-path name of file where modifications are recorded.
find_mod_indexes <- function (file) 
{
  if (!file.exists(file)) 
    stop("File not found: ", file, call. = FALSE)
  
  mod_indexes <- readr::read_tsv(file, show_col_types = FALSE)
  inds <- mod_indexes$Abbr
  names(inds) <- mod_indexes$Desc
  
  inds
}


#' Finds if two sets are equal.
#' 
#' @param x A set.
#' @param y Another set.
is_equal_sets <- function(x, y) all(x %in% y) && all(y %in% x)


#' Purges decoy pep_seq(s) that are also found in target fasta.
#' 
#' Decoy sequences may be present in targets and thus removed.
#' 
#' @examples 
#' ## found in both forward and reverse fastas
#' # pep_seq  prot_acc
#' # RQEEELR NP_064522
#' # RQEEELR NP_997555
#' # RQEEELR NP_997553
#' # 
#' # pep_seq      prot_acc
#' # RQEEELR -NP_001073379
#' # RQEEELR -NP_001157031
#' # RQEEELR    -NP_796085
#' # RQEEELR    -NP_083410
#' 
#' @param target A target data frame with column \code{pep_seq}.
#' @param decoy A decoy data frame with column \code{pep_seq}.
purge_decoys <- function (target, decoy) 
{
  tpeps <- unique(target$pep_seq)
  dpeps <- unique(decoy$pep_seq)
  
  oks <- dpeps[! dpeps %in% tpeps]
  
  dplyr::filter(decoy, pep_seq %in% oks)
}


#' Expands grids.
#'
#' Outputs are vectors corresponding to rows in the the data.frame from the
#' original expand.grid.
#'
#' Some overhead in making row vectors but can avoid the expensive row
#' subsetting from a data.frame.
#'
#' @param use.names Logical; uses names or not.
#' @param ... Lists of data.
#' @examples
#' x <- list(`Oxidation (M)` = c(0.000000, 63.998285),
#'           `Carbamidomethyl (M)` = c(0.000000, 105.024835),
#'           `Oxidation (M)` = c(0.000000, 63.998285))
#'
#' expand_grid_rows(x)
#'
#' x <- list(`Bar (M)` = c(0, 3),
#'           `Foo (M)` = c(0, 5, 7),
#'           `Bar (M)` = c(0, 3))
#'
#' expand_grid_rows(x)
#'
#' x <- list(`Bar (M)` = c(0, 3))
#' expand_grid_rows(x)
expand_grid_rows <- function (..., use.names = TRUE) 
{
  args <- list(...)[[1]]
  nargs <- length(args)
  
  # if (!nargs) return(NULL)

  cargs <- vector("list", nargs)
  names(cargs) <- names(args)
  
  rep.fac <- 1L
  ds <- lengths(args)
  orep <- prod(ds)
  
  for (i in seq_len(nargs)) {
    x <- args[[i]]
    
    nx <- length(x)
    orep <- orep/nx
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)]
                                                
    cargs[[i]] <- x
    rep.fac <- rep.fac * nx
  }
  
  # vectors by rows
  ans <- vector("list", rep.fac)
  
  for (i in seq_len(rep.fac)) {
    x <- lapply(cargs, `[[`, i)
    x <- .Internal(unlist(x, recursive = FALSE, use.names = use.names))
    ans[[i]] <- x
  }

  ans
}


#' A simplified table utility.
#' 
#' A faster alternative to \code{table}.
#' 
#' @param vec A named vector.
#' @examples 
#' \donttest{
#' vec <- c("Carbamidomethyl (M)", "Carbamyl (M)", "Carbamidomethyl (M)")
#' 
#' microbenchmark::microbenchmark(count_elements(vec), table(vec))
#' }
count_elements <- function (vec) 
{
  vals <- unique(vec)
  len <- length(vals)
  
  out <- vector("integer", len)
  
  for (i in seq_len(len)) 
    out[i] <- sum(vec == vals[i])
  
  names(out) <- vals
  
  out
}


#' Split a named character vector to lists.
#' 
#' @param x A named character vector.
#' @examples 
#' \donttest{
#' x <- c(`Deamidated (N)` = "N", `Carbamidomethyl (S)` = "S")
#' 
#' identical(vec_to_list(x), split(x, x))
#' microbenchmark::microbenchmark(vec_to_list(x), split(x, x))
#' }
vec_to_list <- function (x) 
{
  len <- length(x)
  out <- vector("list", len)
  names(out) <- x
  
  for(i in seq_len(len)) 
    out[[i]] <- x[i]
  
  out
}


#' Split a vector by values
#' 
#' @param vec A vector.
#' 
#' @examples
#' \donttest{
#' ## M
#' vec <- c(`Carbamidomethyl (M)` = "M", 
#'          `Carbamyl (M)` = "M")
#' 
#' x <- split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark::microbenchmark(split_vec(vec), split(vec, vec))
#' 
#' ## N
#' vec <- c(`Deamidated (N)` = "N")
#' 
#' x <- split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark::microbenchmark(split_vec(vec), 
#'                                split(vec, vec))
#' 
#' ## M, N
#' vec <- c(`Carbamidomethyl (M)` = "M", 
#'          `Carbamyl (M)` = "M", 
#'          `Deamidated (N)` = "N")
#' 
#' x <- split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark::microbenchmark(split_vec(vec), split(vec, vec))
#' }
#' 
split_vec <- function (vec) 
{
  vals <- unique(vec)
  len <- length(vals)
  
  out <- vector("list", len)
  
  for (i in seq_len(len)) 
    out[[i]] <- vec[vec == vals[i]]
  
  names(out) <- vals
  
  out
}


#' Accumulates character strings.
#' 
#' @param x A vector of character strings.
#' @param f A function.
#' 
#' @examples 
#' \donttest{
#' x <- c("aa", "bb", "cc")
#' accumulate_char(x, paste0)
#' 
#' x <- c(a = 1, b = 2, c = 3)
#' accumulate_char(x, paste0)
#' 
#' x <- "-EDEIQDXI-"
#' accumulate_char(x, paste0)
#' 
#' }
#' 
#' @export
accumulate_char <- function(x, f) 
{
  len <- length(x)
  
  if (len == 1L) 
    return(x)
  
  out <- vector("character", len)
  
  out[1] <- x[1]
  
  for (i in seq(2, len)) 
    out[i] <- f(out[i-1], x[i])
  
  out
}


#' Populates the count matrix for combinations.
#' 
#' @param nb The number of balls.
#' @param ns The number of samplings (the number of columns). 
combi_mat <- function (nb = 5L, ns = 3L) 
{
  # stopifnot(ns >= 1L)
  
  m <- matrix(nrow = nb, ncol = ns)
  m[, 1] <- rep(1L, nb)
  
  if (ns == 1L) 
    return(m)
  
  for (i in seq(2, ns)) 
    m[, i] <- cumsum(m[, i-1])
  
  m
}


#' Makes a data frame of zero rows.
#' 
#' @param vec A vector of names.
make_zero_df <- function (vec) 
{
  df <- data.frame(matrix(ncol = length(vec), nrow = 0L))
  colnames(df) <- vec
  
  invisible(df)
}

#' Makes hash table.
#' 
#' Not yet used.
#' 
#' @param data A named vector or list.
#' @param r Numeric; bucket-to-key ratio.
hash_frame_nums <- function (data, r = 1.5) 
{
  vals <- names(data)
  
  if (is.null(vals))
    stop("Data need names.")
  
  n_bucks <- ceiling(r * length(vals))
  keys <- unlist(lapply(vals, digest::digest2int), recursive = FALSE, use.names = FALSE)
  coll_ids <- keys %% n_bucks + 1L
  
  ans <- vector("list", n_bucks)
  
  coll_data <- split(data, coll_ids) 
  uniq_ids <- as.integer(names(coll_data)) 
  
  for (i in seq_along(coll_data))
    ans[[uniq_ids[i]]] <- coll_data[[i]]

  invisible(ans)
}


#' Finds a value through a hash table.
#' 
#' Not yet used.
#' 
#' @param ht A hash table.
#' @param val A value.
#' @param n_bucks The number of buckets in \code{ht}.
#' @param offset An offset. Not yet used.
find_from_hash <- function (ht, val, n_bucks = length(ht), offset = 0L) 
{
  key <- digest::digest2int(val)
  coll_id <- key %% n_bucks + 1L
  x <- ht[[coll_id]]
  
  x[[val]]
}


#' Calculates the three-frame ppm.
#'
#' @param ppm A positive integer; the mass tolerance of MS1 species. The default
#'   is 20.
#' @param is_three_frame Logical; is a three-frame search or not. The value is
#'   always TRUE.
#' @param fct_ppm A factor to new ppm error. The value is always 0.5 for a
#'   three-frame search.
calc_threeframe_ppm <- function (ppm = 20L, is_three_frame = TRUE, fct_ppm = .5) 
{
  if (is_three_frame) as.integer(ceiling(ppm * fct_ppm)) else ppm
}


#' Checks the status of a prior execution of precursor mass calibration.
#' 
#' @inheritParams matchMS 
#' @return TRUE if without precursor mass calibration.
check_ms1calib <- function(out_path = NULL, calib_ms1mass = FALSE) 
{
  workflow_file <- file.path(out_path, "Calls", "workflow_info.rds")
  
  passed_ms1calib <- if (calib_ms1mass) {
    if (file.exists(workflow_file)) 
      qs::qread(workflow_file)[["passed_ms1calib"]]
    else 
      FALSE
  }
  else {
    TRUE
  }
  
  if (is.null(passed_ms1calib)) FALSE else passed_ms1calib
}


#' Saves the \code{ppm_ms1} before and after calibration.
#' 
#' @param ppm_ms1calib The mass error after calibration in ppm.
#' @inheritParams matchMS
save_ms1calib <- function (ppm_ms1, ppm_ms1calib, mgf_path)
{
  info_calib <- c(`ppm_ms1_bf` = ppm_ms1, `ppm_ms1_af` = ppm_ms1calib)
  qs::qsave(info_calib, file.path(mgf_path, "ppm_ms1calib.rds"))
}


#' Gets the MS1 charges.
#' 
#' @param charges A vector of \code{2+, 3+} etc.
get_ms1charges <- function (charges)
{
  charges <- lapply(charges, stringi::stri_reverse)
  charges <- .Internal(unlist(charges, recursive = FALSE, use.names = FALSE))
  
  as.integer(charges)
}


