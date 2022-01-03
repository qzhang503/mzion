### Also in proteoQ

#' Prefix form of colnames(x)[c(2, 5, ...)] for use in pipes
#'
#' \code{names_pos<-} rename the columns at the indexes of \code{pos}.
#'
#' @param x A data frame.
#' @param pos Numeric.  The index of columns for name change.
#' @param value Characters.  The new column names.
#' @return The data frame with new names.
#'
#' @import dplyr
#' @importFrom magrittr %>% %T>% %$% %<>%
`names_pos<-` <- function(x, pos, value) 
{
  names(x)[pos] <- value
  x
}


#' Finds the columns of reporter-ion intensity.
#'
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
find_int_cols <- function (TMT_plex) 
{
  col_int <- if (TMT_plex == 18) 
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C",
      "I132N", "I132C", "I133N", "I133C", "I134N", "I134C", "I135N") 
  else if (TMT_plex == 16) 
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C",
      "I132N", "I132C", "I133N", "I133C", "I134N")
  else if (TMT_plex == 11) 
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C")
  else if (TMT_plex == 10) 
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131")
  else if(TMT_plex == 6) 
    c("I126", "I127", "I128", "I129", "I130", "I131")
  else 
    NULL
}


#' Re-order columns in a data frame
#'
#' \code{ins_cols_after} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param idx_bf A column index for the insertion of columns after.
#' @param idx_ins Column index(es) for the columns to be inserted after
#'   \code{idx_bf}.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>%
ins_cols_after <- function(df = NULL, idx_bf = ncol(df), idx_ins = NULL) 
{
  if (is.null(df)) 
    stop("`df` cannot be `NULL`.", call. = FALSE)
  if (is.null(idx_ins)) 
    return(df)
  if (idx_bf >= ncol(df)) 
    return(df)
  
  col_nms <- names(df)[idx_ins]
  
  dplyr::bind_cols(
    df[, 1:idx_bf, drop = FALSE],
    df[, idx_ins, drop = FALSE],
    dplyr::select(df[, (idx_bf + 1):ncol(df), drop = FALSE], -col_nms))
}


#' Pad columns to a placeholder data frame.
#'
#' @param df The original data frame.
#' @param df2 The data frame to be inserted.
#' @param idx The index of \code{df} column for \code{df2} to be inserted
#'   (after).
add_cols_at <- function(df, df2, idx) 
{
  stopifnot(idx >= 0L)
  
  bf <- if (idx == 0L) 
    NULL
  else 
    df[, seq_len(idx), drop = FALSE]
  
  af <- if ((idx + 1L) <= ncol(df)) 
    df[, (idx + 1L) : ncol(df), drop = FALSE]
  else 
    NULL
  
  dplyr::bind_cols(
    bf,
    df2,
    af)
}


#' Replace columns in the original PSM table.
#'
#' The column index(es) need to be continuous.
#'
#' @param df The original data frame.
#' @param df2 The data columns to replace those in \code{df}.
#' @param idxs The sequences of column indexes in \code{df}. Note that
#'   \code{idxs} need to be a continuous sequences.
replace_cols_at <- function(df, df2, idxs) 
{
  ncol <- ncol(df)
  stopifnot(all(idxs >= 1L), all(idxs <= ncol))
  
  idxs <- sort(idxs)
  stopifnot(all.equal(idxs - idxs[1] + 1L, seq_along(idxs)))
  
  bf <- if (idxs[1] >= 2L) 
    df[, 1:(idxs[1]-1L), drop = FALSE]
  else 
    NULL
  
  af <- if (idxs[length(idxs)] < ncol(df)) 
    df[, (idxs[length(idxs)]+1L):ncol(df), drop = FALSE]
  else 
    NULL
  
  dplyr::bind_cols(
    bf,
    df2,
    af)
}


#' Relocate column "to_move" immediately after column "col_before".
#'
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_before The anchor column to which the \code{to_move} will be moved
#'   after.
reloc_col_after <- function(df, to_move = "after_anchor", col_before = "anchor") 
{
  if (!(to_move %in% names(df) && col_before %in% names(df))) 
    return(df)
  
  if (to_move == col_before) 
    return(df)
  
  df2 <- dplyr::select(df, one_of(to_move))
  df <- dplyr::select(df, -one_of(to_move))
  
  idx <- which(names(df) == col_before)
  df <- add_cols_at(df, df2, idx)
  
  invisible(df)
}


#' Relocate column "to_move" immediately after the last column.
#'
#' @inheritParams reloc_col_after
reloc_col_after_last <- function (df, to_move = "after_anchor") 
{
  col_last <- names(df)[ncol(df)]
  reloc_col_after(df, to_move, col_last)
}


#' Relocate column "to_move" immediately after the first column.
#'
#' @inheritParams reloc_col_after
reloc_col_after_first <- function(df, to_move = "after_anchor") 
{
  col_first <- names(df)[1]
  reloc_col_after(df, to_move, col_first)
}


#' Relocate column "to_move" immediately before anchor column "col_after".
#'
#' The same as \code{reloc_col}.
#'
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_after The anchor column to which the \code{to_move} will be moved
#'   before.
reloc_col_before <- function(df, to_move = "before_anchor",
                             col_after = "anchor") 
{
  if (!(to_move %in% names(df) && col_after %in% names(df))) 
    return(df)
  
  df2 <- dplyr::select(df, one_of(to_move))
  df <- dplyr::select(df, -one_of(to_move))
  
  idx <- which(names(df) == col_after)
  
  df <- add_cols_at(df, df2, idx - 1)
  
  invisible(df)
}


#' Relocate column "to_move" immediately before the last column.
#'
#' @inheritParams reloc_col_after
reloc_col_before_last <- function(df, to_move = "after_anchor") 
{
  col_last <- names(df)[ncol(df)]
  reloc_col_before(df, to_move, col_last)
}


#' Relocate column "to_move" immediately before the first column.
#'
#' @inheritParams reloc_col_after
reloc_col_before_first <- function (df, to_move = "after_anchor") 
{
  col_first <- names(df)[1]
  reloc_col_before(df, to_move, col_first)
}


#' Helper: find the column name before \code{to_move}.
#'
#' To keep columns at the same order after descriptive summary.
#'
#' @inheritParams reloc_col_after
find_preceding_colnm <- function(df, to_move) 
{
  if (!to_move %in% names(df)) 
    stop("Column ", to_move, " not found.", call. = FALSE)

  ind_bf <- which(names(df) == to_move) - 1
  
  if (ind_bf == 0) 
    names(df)[1]
  else 
    names(df)[ind_bf]
}


#' Flatten lists recursively
#'
#' @param x Lists
recur_flatten <- function (x) 
{
  if (!inherits(x, "list")) 
    list(x)
  else 
    .Internal(unlist(c(lapply(x, recur_flatten)), recursive = FALSE, 
                     use.names = FALSE))
}


### End of also in proteoQ



###

#' Splits data into chunks by length.
#'
#' @param data Input data.
#' @param n_chunks The number of chunks.
#' @param type The type of data for splitting.
chunksplit <- function (data, n_chunks = 5L, type = c("list", "row") )
{
  type <- match.arg(type)

  stopifnot(type %in% c("list", "row"))
  
  if (n_chunks <= 1L) 
    return(data)
  
  if (type == "list") 
    len <- length(data)
  else if (type == "row") 
    len <- nrow(data)
  else 
    stop("Unknown type.", call. = TRUE)
  
  if (len == 0L) 
    return(data)
  
  labs <- levels(cut(1:len, n_chunks))
  
  x <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
             upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
  
  grps <- findInterval(1:len, x[, 1])
  
  split(data, grps)
}


#' Splits data into chunks with approximately equal sizes.
#'
#' @param nx Positive integer; an arbitrarily large number for data to be split
#'   into for estimating the cumulative sizes.
#' @inheritParams chunksplit
chunksplitLB <- function (data, n_chunks = 5L, nx = 100L, type = "list") 
{
  stopifnot(type %in% c("list", "row"))
  
  if (n_chunks <= 1L) 
    return(data)
  
  if (type == "list") 
    len <- length(data)
  else if (type == "row") 
    len <- nrow(data)
  else 
    stop("Unknown type.", call. = TRUE)
  
  if (len == 0L) 
    return(data)
  
  # The finer groups by 'nx'
  grps_nx <- local({
    labsx <- levels(cut(1:len, nx))
    
    xx <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1",
                                              labsx))),
                upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1",
                                                labsx))))
    
    findInterval(1:len, xx[, 1])
  })
  
  # The equated size for a chunk
  size_chunk <- local({
    size_nx <- data %>%
      split(., grps_nx) %>%
      lapply(object.size) %>%
      cumsum()
    
    size_nx[length(size_nx)]/n_chunks
  })
  
  #  Intervals
  grps <- local({
    size_data <- data %>%
      lapply(object.size) %>%
      cumsum()
    
    # the position indexes
    ps <- purrr::map_dbl(1:(n_chunks-1), function (x) {
      which(size_data < size_chunk * x) %>% `[`(length(.))
    })
    
    grps <- findInterval(1:len, ps)
  })
  
  split(data, grps)
}


#' Finds a file directory.
#'
#' With the option of creating a directory
#'
#' @param path A file path.
#' @param create Logical; if TRUE create the path if not yet existed.
find_dir <- function (path, create = FALSE) 
{
  stopifnot(length(path) == 1L)
  
  path <- gsub("\\\\", "/", path)
  
  p1 <- fs::path_expand_r(path)
  p2 <- fs::path_expand(path)
  
  if (fs::dir_exists(p1)) {
    path <- p1
  } else if (fs::dir_exists(p2)) {
    path <- p2
  } else {
    if (create) {
      dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)
      path <- p1
    } else {
      message(path, " not found.")
      path <- NULL
    }
  }
  
  path
}


#' Creates a file directory.
#'
#' @inheritParams find_dir
create_dir <- function (path) find_dir(path, create = TRUE)


#' Saves the arguments in a function call
#'
#' @param path A (system) file path.
#' @param fun The name of function being saved.
#' @param time The time stamp.
save_call2 <- function(path, fun, time = NULL) 
{
  stopifnot(length(path) == 1L, length(fun) == 1L)
  
  call_pars <- mget(names(formals(fun)), envir = parent.frame(), inherits = FALSE)
  call_pars[names(call_pars) == "..."] <- NULL
  
  if (is.null(time)) {
    p2 <- create_dir(path)
    save(call_pars, file = file.path(p2, paste0(fun, ".rda")))
  } 
  else {
    stopifnot(length(time) == 1L)
    p2 <- create_dir(file.path(path, fun))
    save(call_pars, file = file.path(p2, paste0(time, ".rda")))
  }
}


#' Finds the values of a list of arguments.
#'
#' @param args Arguments to be matched.
#' @inheritParams save_call2
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$%
#' @return The time stamp of a matched cache results.
find_callarg_vals <- function (time = NULL, path = NULL, fun = NULL,
                               args = NULL) 
{
  stopifnot(length(path) == 1L, length(fun) == 1L)
  
  if (is.null(time)) {
    file <- file.path(path, fun)
  } else {
    stopifnot(length(time) == 1L)
    file <- file.path(path, fun, time)
  }
  
  if (!file.exists(file)) 
    return(NULL)

  load(file = file)
  
  # --- back-compatibility
  call_pars <- local({
    new_args <- c(n_13c = 0L, min_ms1_charge = 2L, max_ms1_charge = 6L, 
                  min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                  min_ret_time = 0L, max_ret_time = .Machine$integer.max)
    
    for (i in seq_along(new_args)) {
      x <- new_args[i]
      nm <- names(x)
      val <- unname(x)
      
      if ((nm %in% args) && (! nm %in% names(call_pars)))
        call_pars[[nm]] <- val
    }
    
    call_pars
  })
  
  nots <- which(! args %in% names(call_pars))
  
  if (length(nots)) 
    stop("Arguments '", paste(args[nots], collapse = ", "),
         "' not found in the latest call to `", fun, "`.\n", 
         "This is probably due to new argument implementation in `matchMS`.\n", 
         "Delete `", file.path(path, fun), "`", " and try again.", 
         call. = FALSE)

  call_pars[args]
}


#' Finds the time stamp of a matched call from cached results.
#'
#' @param nms Names of arguments to be included in or excluded from matching.
#' @param type Logical; if TRUE, includes all arguments in \code{nms}. At FALSE,
#'   excludes all \code{nms}.
#' @inheritParams find_callarg_vals
#' @return An empty object if no matches.
match_calltime <- function (path = "~/proteoM/.MSearches/Cache/Calls",
                            fun = "calc_pepmasses2",
                            nms = c("parallel", "out_path"),
                            type = c(TRUE, FALSE)) 
{
  stopifnot(length(path) == 1L, length(fun) == 1L)
  
  if (length(type) > 1L) type <- TRUE
  
  # current
  args <- if (type) 
    mget(names(formals(fun)) %>% .[. %in% nms], 
         envir = parent.frame(), inherits = FALSE)
  else 
    mget(names(formals(fun)) %>% .[! . %in% nms], 
         envir = parent.frame(), inherits = FALSE)

  if (!length(args)) 
    stop("Arguments for matching is empty.", call. = FALSE)
  
  args <- lapply(args, sort)
  
  times <- list.files(path = file.path(path, fun),
                      pattern = "\\.rda$",
                      all.files = TRUE)
  
  # cached
  cached <- lapply(times, find_callarg_vals, path = path, fun = fun,
                   args = names(args))
  
  cached <- lapply(cached, function (x) lapply(x, sort))
  
  # matched
  oks <- lapply(cached, identical, args)
  oks <- unlist(oks, recursive = FALSE, use.names = FALSE)
  
  times[oks] %>% gsub("\\.rda$", "", .)
}


#' Deletes files under a directory.
#'
#' The directory will be kept.
#' 
#' @param path A file path.
#' @param ignores The file extensions to be ignored.
delete_files <- function (path, ignores = NULL, ...) 
{
  
  # Hadley Ch.19 Quasiquotation, note 97
  dots <- as.list(substitute(...()))
  recursive <- dots[["recursive"]]
  dots$recursive <- NULL
  
  if (is.null(recursive)) 
    recursive <- TRUE
  
  if (!is.logical(recursive)) 
    stop("Not logical `recursive`.")

  args <- c(list(path = file.path(path), recursive = recursive,
                 full.names = TRUE), dots)
  
  nms <- do.call(list.files, args)
  
  if (!is.null(ignores)) {
    nms <- local({
      dirs <- list.dirs(path, full.names = FALSE, recursive = recursive)
      dirs <- dirs[! dirs == ""]

      idxes_kept <- purrr::map_lgl(dirs, function (x) any(grepl(x, ignores)))
      
      nms_kept <- list.files(path = file.path(path, dirs[idxes_kept]),
                             recursive = TRUE, full.names = TRUE)
      
      nms[! nms %in% nms_kept]
    })
    
    nms <- local({
      exts <- gsub("^.*(\\.[^.]*)$", "\\1", nms)
      idxes_kept <- purrr::map_lgl(exts, function (x) any(grepl(x, ignores)))
      
      nms[!idxes_kept]
    })
  }
  
  if (length(nms) > 0L) 
    suppressMessages(file.remove(file.path(nms)))

  invisible(NULL)
}


