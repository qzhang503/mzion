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
#' library(mzion)
#' 
#' mzion:::which_topx(c(1:5), 50)
#'
#' length(mzion:::which_topx(sample(5000, 500), 100))
#'
#' length(mzion:::which_topx(sample(100, 100, replace = TRUE), 100))
#' }
which_topx <- function(x, n = 50L, ...) 
{
  len <- length(x)
  p <- len - n
  
  if (p  <= 0L) 
    return(seq_along(x))
  
  xp <- sort(x, partial = p, ...)[p]
  
  .Internal(which(x > xp))
}


#' Finds the indexes of top-n entries without re-ordering.
#' 
#' Handles ties.
#' 
#' @param replace_na Logical; currently always TRUE; replaces NA or not.
#' @param vna The value of NA replacement.
#' @param exclude_na Logical; removes NA or not.
#' @inheritParams which_topx
#' @return The indexes of the top-n entries.
#' @examples 
#' library(mzion)
#' 
#' p <- 100
#' set.seed(1)
#' x <- sample(1:150, replace = TRUE)
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
#' ans <- mzion:::which_topx2(x, 2L, na.last = FALSE)
#' 
#' # OK
#' ans <- mzion:::which_topx2(x, 8L, na.last = FALSE)
#' 
#' # results can contain NA
#' x <- c(10, rep_len(NA_integer_, 3), 12, rep_len(NA_integer_, 3), 8, 10, 10)
#' ans <- mzion:::which_topx2(x, n = 7L)
#' x[ans]
#' 
#' # Bad
#' ans <- mzion:::which_topx2(x, 9L, na.last = FALSE)
#' # OK
#' ans <- mzion:::which_topx2(x, 9L, na.last = TRUE)
#'  
#' \donttest{
#' mzion:::which_topx2(5000, NA_integer_)
#' }
which_topx2 <- function(x, n = 50L, replace_na = TRUE, vna = 0, 
                        exclude_na = TRUE, ...) 
{
  if (is.na(n)) 
    return(NULL)
  
  len <- length(x)
  p <- len - n

  if (p  <= 0L) {
    if (exclude_na)
      return(.Internal(which(!is.na(x))))
    else
      return(seq_along(x))
  }

  if (replace_na)
    x[is.na(x)] <- vna
  
  # not yet work at replace_na = FALSE and x contains fewer non-NA than p
  xp <- sort(x, partial = p, ...)[p]
  inds <- .Internal(which(x > xp)) # NA drops by which()

  # in case of ties -> length(inds) < n
  # detrimental e.g. ms2_n = 500 and n = 100
  #   -> expect 100 `ms2_moverzs` guaranteed but may be only 99
  #
  # MGF `ms2_moverzs` is increasing
  # `inds2` goes first to ensure non-decreasing index for `ms2_moverzs`
  d <- n - length(inds)
  
  if (!d)
    return(inds)

  # must exist and length(ties) >= length(d) at replace_na = TRUE
  ties <- .Internal(which(x == xp))
  # better check the length(ties) >= length(d)
  for (i in seq_len(d)) {
    inds <- insVal(ties[[i]], inds)
  }
  
  inds
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
#' library(mzion)
#' sv <- c(4, 10)
#' 
#' sv <- mzion:::insVal(2, sv)
#' sv <- mzion:::insVal(6, sv)
#' sv <- mzion:::insVal(20, sv)
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


#' Sums elements across lists.
#'
#' Each list has the same length. NA values are removed.
#'
#' @param x A numeric value.
#' @param y A numeric value.
`%+%` <- function(x, y) mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))


#' Post frame advancing.
#'
#' Flattens mgfs within each frame.
#' 
#' @param res Results from frame-advanced searches.
#' @param mgf_frames Data of MGF frames.
post_frame_adv <- function (res, mgf_frames) 
{
  res  <- unlist(res, recursive = FALSE, use.names = FALSE)
  out <- do.call(rbind, mgf_frames)
  # stopifnot(nrow(out) == length(res))
  out <- dplyr::mutate(out, matches = res)
  out <- out[lengths(res, use.names = FALSE) > 0L, ]
}


#' Subsets the frames of theoretical peptides.
#' 
#' @param mgf_frames MGFs in frames. Each frame contains one to multiple MGFs
#'   whose MS1 masses are in the same interval.
#' @param theopeps Binned theoretical peptides at a given combination of fixed
#'   and variable.
subset_theoframes <- function (mgf_frames = NULL, theopeps = NULL) 
{
  if ((!length(mgf_frames)) || (!length(theopeps)))
    return(NULL)
  
  frames <- as.integer(names(mgf_frames))
  breaks <- which(diff(frames) != 1L) + 1L
  groups <- findInterval(frames, frames[breaks])
  frames <- split(frames, groups)
  
  frames <- lapply(frames, function (x) c(x[1] - 1, x, x[length(x)] + 1))
  frames <- unlist(frames, recursive = FALSE, use.names = FALSE)
  frames <- frames[!duplicated(frames)]

  # NA names and NULL contents if br_frames not in theopeps
  # (length determined by `br_frames`)
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
#' Not a general purpose utility. For simple circumstances that the number of
#' rows in the output is equal to that in \code{y}, and can be use as an
#' alternative to \link[dplyr]{right_join}.
#' 
#' Rows ordered by \code{y}, which is different to \link[dplyr]{right_join}. 
#'
#' @param x The left data frame (to be proliferated by rows).
#' @param y The right data frame (the dominated one).
#' @param by The key.
#' @examples
#' \donttest{
#' library(mzion)
#' library(dplyr)
#' 
#' df1 <- data.frame(A = c("a", "b", "c"), B = c(1, 1, 1))
#' df2 <- data.frame(A = c("a", "c", "d"), C = c(2, 2, "3"))
#'
#' x1 <- mzion:::quick_rightjoin(df1, df2, by = "A")
#' x1 <- x1[, c("A", "B", "C")]
#' rownames(x1) <- seq_len(nrow(x1))
#'
#' x2 <- dplyr::right_join(df1, df2, by = "A")
#' x2 <- x2[, c("A", "B", "C")]
#'
#' stopifnot(identical(x1, x2))
#'
#' # row order may be different
#' x1 <- mzion:::quick_rightjoin(df2, df1, by = "A")
#' x1 <- x1[, c("A", "B", "C")]
#' rownames(x1) <- seq_len(nrow(x1))
#'
#' x2 <- dplyr::right_join(df2, df1, by = "A")
#' x2 <- x2[, c("A", "B", "C")]
#'
#' # FALSE
#' !identical(x1, x2)
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
#' Not a general purpose utility. For simple circumstances that the number of
#' rows in the output is equal to that in \code{x}, and can be use as an
#' alternative to \link[dplyr]{left_join}.
#'
#' @param x The left data frame.
#' @param y The right data frame.
#' @param by The key.
#' @examples
#' \donttest{
#' library(mzion)
#' library(dplyr)
#'
#' df1 <- data.frame(A = c("a", "b", "c"), B = c(1, 1, 1))
#' df2 <- data.frame(A = c("a", "c", "d"), C = c(2, 2, "3"))
#'
#' x1 <- mzion:::quick_leftjoin(df1, df2, by = "A")
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
#' @param sys_ram The putative amount of system RAM, e.g., 32L.
find_free_mem <- function (sys_ram = NULL) 
{
  nm_os <- Sys.info()['sysname']

  if (nm_os == "Windows") {
    free_mem <- system('wmic OS get FreePhysicalMemory /Value', intern=TRUE)[3]
    free_mem <- gsub("^FreePhysicalMemory=(\\d+)\\r", "\\1", free_mem)
    free_mem <- as.numeric(free_mem)/1024
    
    if (!is.null(sys_ram)) 
      free_mem <- min(sys_ram, free_mem)
  } 
  else {
    # not yet tested for "Linux", "Darwin"
    # inaccurate e.g. if physical RAM is 32000 but in .RProfile 
    # `invisible(utils::memory.limit(64000))`
    # memory.limit() - memory.size(max = TRUE)
    
    warning("Cannot determine the amount of RAM with Linux or MAC OS.\n", 
            "To specify, use parameter \"sys_ram\".")
    free_mem <- 24L
  }
  
  free_mem
}


#' Find the indexes of modifications.
#' 
#' @param file A full-path name of file where modifications are recorded.
find_mod_indexes <- function (file) 
{
  if (!file.exists(file)) 
    stop("File not found: ", file)
  
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


#' Expands grids.
#'
#' Outputs are vectors corresponding to rows in the the data.frame from the
#' original expand.grid.
#'
#' Some overhead in making row vectors but can avoid the expensive row
#' subsetting from a data.frame.
#' 
#' @param nmax The maximum number of combinations allowed.
#' @param use.names Logical; uses names or not.
#' @param ... Lists of data.
#' @examples
#' library(mzion)
#' 
#' x <- list(`Oxidation (M)` = c(0.000000, 63.998285),
#'           `Carbamidomethyl (M)` = c(0.000000, 105.024835),
#'           `Oxidation (M)` = c(0.000000, 63.998285))
#'
#' mzion:::expand_grid_rows(x)
#'
#' x <- list(`Bar (M)` = c(0, 3),
#'           `Foo (M)` = c(0, 5, 7),
#'           `Bar (M)` = c(0, 3))
#'
#' mzion:::expand_grid_rows(x)
#'
#' x <- list(`Bar (M)` = c(0, 3))
#' mzion:::expand_grid_rows(x)
expand_grid_rows <- function (..., nmax = 3L, use.names = TRUE) 
{
  args <- list(...)[[1]]
  nargs <- length(args)
  tot <- nargs * nmax
  cargs <- vector("integer", )
  
  rep.fac <- 1L
  ds <- lengths(args)
  orep <- prod(ds)
  nmax <- min(nmax, orep)
  
  sta <- 1L
  end <- nmax
  
  for (i in seq_len(nargs)) {
    x <- args[[i]]
    
    nx <- length(x)
    orep <- orep/nx
    x <- x[rep_len(rep.int(seq_len(nx), rep.int(rep.fac, nx)), nmax)]
    
    cargs[sta:end] <- x
    sta <- sta + nmax
    end <- end + nmax
    rep.fac <- rep.fac * nx
  }
  
  # vectors by rows
  len <- length(x)
  ans <- vector("list", len)
  nms <- names(args)
  
  seqs <- 0:(nargs-1L)
  
  for (i in 1:len) {
    y <- cargs[seqs * nmax + i]
    names(y) <- nms
    ans[[i]] <- y
  }
  
  ans
}


#' Modified from expand.grid
#' 
#' Net yet used.
#' 
#' @param nmax The maximum number of combinations allowed.
#' @param ... Lists of data.
expand_grid <- function (..., nmax = 3L) 
{
  nargs <- length(args <- list(...))
  
  if (!nargs) 
    return(as.data.frame(list()))
  
  if (nargs == 1L && is.list(a1 <- args[[1L]])) 
    nargs <- length(args <- a1)
  if (nargs == 0L) 
    return(as.data.frame(list()))
  
  cargs <- vector("list", nargs)
  iArgs <- seq_len(nargs)
  nmc <- paste0("Var", iArgs)
  nm <- names(args)
  nm <- nmc
  names(cargs) <- nmc
  
  rep.fac <- 1L
  d <- lengths(args)
  orep <- prod(d)
  nmax <- min(nmax, orep)
  
  if (orep == 0L) {
    for (i in iArgs) cargs[[i]] <- args[[i]][FALSE]
  }
  else {
    for (i in iArgs) {
      x <- args[[i]]
      nx <- length(x)
      orep <- orep/nx
      
      x <- x[rep_len(rep.int(seq_len(nx), rep.int(rep.fac, nx)), nmax)]
      cargs[[i]] <- x
      rep.fac <- rep.fac * nx
    }
  }
  
  rn <- .set_row_names(length(x))
  structure(cargs, class = "data.frame", row.names = rn)
}


#' Modified from expand.grid
#' 
#' Inputs are in matrices
#' 
#' @param nmax The maximum number of combinations allowed.
#' @param ... Lists of data.
expand_gr <- function (..., nmax = 3L) 
{
  nargs <- length(args <- list(...))
  
  if (nargs == 1L && is.list(a1 <- args[[1L]])) 
    nargs <- length(args <- a1)
  
  nr <- .Internal(unlist(lapply(args, nrow), recursive = FALSE, use.names = FALSE))
  nc <- .Internal(unlist(lapply(args, ncol), recursive = FALSE, use.names = FALSE))
  cnc <- cumsum(nc)
  orep <- prod(nr)
  nmax <- min(nmax, orep)
  
  M <- matrix(nrow = nmax, ncol = sum(nc))
  rep.fac <- 1L
  c1 <- 1L
  
  for (i in seq_len(nargs)) {
    x <- args[[i]]
    nx <- nr[[i]]
    orep <- orep/nx
    c2 <- cnc[[i]]
    M[, c1:c2] <- x[rep_len(rep.int(seq_len(nx), rep.int(rep.fac, nx)), nmax), ]
    c1 <- c2 + 1L
    rep.fac <- rep.fac * nx
  }
  
  M
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
#' library(mzion)
#' 
#' x <- list(`Oxidation (M)` = c(0.000000, 63.998285),
#'           `Carbamidomethyl (M)` = c(0.000000, 105.024835),
#'           `Oxidation (M)` = c(0.000000, 63.998285))
#'
#' mzion:::expand_grid_rows0(x)
#'
#' x <- list(`Bar (M)` = c(0, 3),
#'           `Foo (M)` = c(0, 5, 7),
#'           `Bar (M)` = c(0, 3))
#'
#' mzion:::expand_grid_rows0(x)
#'
#' x <- list(`Bar (M)` = c(0, 3))
#' mzion:::expand_grid_rows0(x)
expand_grid_rows0 <- function (..., use.names = TRUE) 
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
#' \dontrun{
#' library(mzion)
#' library(microbenchmark)
#' 
#' vec <- c("Carbamidomethyl (M)", "Carbamyl (M)", "Carbamidomethyl (M)")
#' microbenchmark(mzion:::count_elements(vec), table(vec))
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
#' \dontrun{
#' library(mzion)
#' library(microbenchmark)
#' 
#' x <- c(`Deamidated (N)` = "N", `Carbamidomethyl (S)` = "S")
#' identical(mzion:::vec_to_list(x), split(x, x))
#' microbenchmark(mzion:::vec_to_list(x), split(x, x))
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


#' Splits a matrix
#' 
#' @param M A matrix.
#' @param by Split by matrix rows or columns.
split_matrix <- function (M, by = "row") 
{
  if (by == "row")
    len <- nrow(M)
  
  ans <- vector("list", len)
  
  for (i in seq_along(ans))
    ans[[i]] <- M[i, ]
  
  ans
}


#' Split a vector by values
#' 
#' @param vec A vector.
#' @param use_names Logical; uses names or not.
#' 
#' @examples
#' \dontrun{
#' ## M
#' library(mzion)
#' library(microbenchmark)
#' 
#' vec <- c(`Carbamidomethyl (M)` = "M", 
#'          `Carbamyl (M)` = "M")
#' 
#' x <- mzion:::split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark(mzion:::split_vec(vec), split(vec, vec))
#' 
#' ## N
#' vec <- c(`Deamidated (N)` = "N")
#' 
#' x <- mzion:::split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark(mzion:::split_vec(vec), split(vec, vec))
#' 
#' ## M, N
#' vec <- c(`Carbamidomethyl (M)` = "M", 
#'          `Carbamyl (M)` = "M", 
#'          `Deamidated (N)` = "N")
#' 
#' x <- mzion:::split_vec(vec)
#' y <- split(vec, vec)
#' 
#' stopifnot(identical(x, y))
#' 
#' microbenchmark(mzion:::split_vec(vec), split(vec, vec))
#' }
#' 
split_vec <- function (vec, use_names = TRUE) 
{
  vals <- unique(vec)
  len <- length(vals)
  
  out <- vector("list", len)
  
  for (i in seq_len(len)) 
    out[[i]] <- vec[vec == vals[i]]
  
  if (use_names)
    names(out) <- vals
  
  out
}


#' Folds a vector evenly
#'
#' Only for the special case that \code{length(vec)} is exactly n-times of the
#' \code{fold}.
#' 
#' @param vec A vector
#' @param fold The number of folds
#' @seealso \link{fold_vec2} for a more generalized solution.
fold_vec <- function (vec, fold = 5L) 
{
  # !!! only applicable when length(vec) is n-times of fold !!!

  ans <- vector("list", fold)
  # `r` must be integer
  r <- length(vec)/fold

  sta <- 1L
  end <- r
  
  for (i in 1:fold) {
    ans[[i]] <- vec[sta:end]
    sta <- sta + r
    end <- end + r
  }

  ans
}


#' Folds a vector evenly
#'
#' A more general form of \link{fold_vec} handling fractional
#' \code{length(vec)/fold}.
#'
#' @param vec A vector
#' @param fold The number of folds
#' @param remove_allna_subvec Logical; if TRUE, removes all NA sub vectors. This
#'   can occur with the rounding of lengths. E.g., at \code{length(vec) == 15}
#'   and \code{fold = 6}, the 6-th subvec will be all NA.
#' @examples
#' vec <- sort(rep(LETTERS[1:5], 1:5))
#' mzion:::fold_vec2(vec, 6L)
#' @seealso \link{find_group_breaks}
fold_vec2 <- function (vec, fold = 5L, remove_allna_subvec = TRUE) 
{
  len <- length(vec)
  r <- ceiling(len/fold)
  
  if (remove_allna_subvec) {
    dif <- r * fold -len
    
    if (dif > 0) {
      mod <- len %/% r
      fold <- if (dif %% r == 0L) mod else mod + 1L
    }
  }

  # faster than mapply
  ans <- vector("list", fold)
  sta <- 1L
  end <- r
  
  for (i in 1:fold) {
    ans[[i]] <- vec[sta:end]
    sta <- sta + r
    end <- end + r
  }

  alast <- ans[[fold]]
  ans[[fold]] <- alast[!is.na(alast)]
  
  ans
}


#' Separates a vector into chunks.
#' 
#' @param vec A vector.
#' @param fold The number of folds.
#' @examples
#' vec <- sort(rep(LETTERS[1:5], 1:5))
#' mzion:::sep_vec(vec, 6L)
#' @return A vector of fold indexes.
sep_vec <- function (vec, fold = 5L)
{
  len <- length(vec)
  seqs <- seq_along(vec)
  mod <- len / fold
  (seqs - 1L) %/% mod + 1L
}


#' Replicates a vector
#' 
#' Slower than rep
#' 
#' @param vec A vector
#' @param fold The number of folds
rep_vec <- function (vec, fold) 
{
  l <- length(vec)
  r <- l * fold
  ans <- vector("character", r)

  sta <- 1L
  end <- l
  
  for (i in 1:fold) {
    ans[sta:end] <- vec
    sta <- sta + l
    end <- end + l
  }
  
  ans
}


#' Accumulates character strings.
#' 
#' @param x A vector of character strings.
#' @param f A function.
#' 
#' @examples 
#' \donttest{
#' library(mzion)
#' 
#' x <- c("aa", "bb", "cc")
#' mzion:::accumulate_char(x, paste0)
#' 
#' x <- c(a = 1, b = 2, c = 3)
#' mzion:::accumulate_char(x, paste0)
#' 
#' x <- "-EDEIQDXI-"
#' mzion:::accumulate_char(x, paste0)
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


#' Gets the MS1 charges.
#' 
#' @param charges A vector of \code{2+, 3+} etc.
get_ms1charges <- function (charges)
{
  charges <- lapply(charges, stringi::stri_reverse)
  charges <- .Internal(unlist(charges, recursive = FALSE, use.names = FALSE))
  
  as.integer(charges)
}


#' Finds unique elements including the same name in a vector.
#' 
#' The same value at different name will be treated as different elements.
#' 
#' @param x A named vector.
#' @examples
#' \donttest{
#' library(mzion)
#' mzion:::finds_uniq_vec(c(a1 = 3, a2 = 3, a2 = 4, a2 = 3, b1 = 3, b2 = 1))
#' }
finds_uniq_vec <- function (x)
{
  v <- .Internal(paste0(list(x, names(x)), collapse = NULL, recycle0 = FALSE))
  oks <- !duplicated.default(v)
  x[oks]
}


#' Makes a data frame from lists.
#' 
#' Not yet used.
#' 
#' @param l A list of named vectors.
my_dataframe <- function (l)
{
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}


#' Alternative to \link[purrr]{flatten}
#' 
#' @param data A list of data.
#' @param use_names Logical; to use names or not.
flatten_list <- function (data, use_names = TRUE)
{
  len <- length(data)
  
  if (!len)
    return(data)
  
  lens <- .Internal(unlist(lapply(data, length), recursive = FALSE, 
                           use.names = FALSE))
  oks  <- which(lens > 0L)
  
  if (!length(oks))
    return(vector("list"))

  data <- data[oks]
  lens <- lens[oks]
  len  <- length(data)
  
  if (len == 1L)
    return(data[[1]])
  
  ends <- cumsum(lens)
  stas <- c(1, ends + 1L)

  ans <- vector("list", ends[len])
  
  for (i in 1:len) 
    ans[stas[[i]]:ends[[i]]] <- data[[i]]
  
  if (use_names) 
    names(ans) <- .Internal(unlist(lapply(data, names), recursive = FALSE, 
                                   use.names = FALSE))

  ans
}


#' Calculates the reversed MS2 from the forward
#' 
#' @param af An sequence of answer of the forward.
#' @param aas The sequence of amino acid residues.
calc_rev_ms2 <- function (af, aas) 
{
  l <- length(aas)
  l1 <- l - 1L
  l2 <- l - 2L
  sb <- af[2:l1]
  b <- c(af[1] + sb[l2] - sb[l2:1L], af[l1:l]) # 3.5 us
  nb <- c(aas[1], aas[l1:2], aas[l]) # 2.2 us
  names(b) <- nb
  
  ll <- l + l
  l3 <- ll - 1L
  sy <- af[(l+2L):l3]
  y <- c(af[l+1L] + sy[l2] - sy[l2:1L], af[l3:ll]) # 3.5 us
  names(y) <- nb[l:1L]
  
  c(b, y)
}


#' Binds data frames
#'
#' Assume identical column names across data frames. Also different to
#' \link[dplyr]{bind_rows} by row names.
#'
#' @param dfs A list of data frames.
bind_dfs <- function (dfs)
{
  nrs <- lapply(dfs, nrow)
  oks <- nrs > 0L
  dfs <- dfs[oks]
  nrs <- nrs[oks]
  
  lens <- unlist(nrs)
  tlen <- sum(lens)
  lend <- length(dfs)
  
  if (tlen == 1L || lend == 1L)
    return(dfs[[1]])
  
  cns <- colnames(dfs[[1]])
  out <- data.frame(matrix(ncol = length(cns), nrow = tlen))
  colnames(out) <- cns
  
  ends <- cumsum(lens)
  stas <- c(1, ends + 1L)
  
  for (i in seq_along(dfs))
    out[stas[[i]]:ends[[i]], ] <- dfs[[i]]

  out
}


#' Find the minimal number of cores.
#' 
#' @param x The number of files.
#' @param n The number of cores where \code{x > n > 1}.
find_min_ncores <- function (x = 25, n = 6)
{
  if (x <= 0)
    stop("`x` needs to be greater than zero.")
  
  if (n <= 1L)
    return(n)
  
  if (x <= n)
    return(x)
  
  len <- ceiling(x/n)
  len2 <- ceiling(x/(n2 <- n - 1L))
  # len2 >= len
  
  while(len2 == len) {
    n <- n2
    len <- len2
    n2 <- n - 1L # ok n2 == 0L
    len2 <- ceiling(x/n2)
  }
  
  n
}


