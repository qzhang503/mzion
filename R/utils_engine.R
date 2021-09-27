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
`names_pos<-` <- function(x, pos, value) {
  names(x)[pos] <- value
  x
}


#' Finds the columns of reporter-ion intensity.
#'
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
find_int_cols <- function (TMT_plex) {
  
  if (TMT_plex == 16) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C",
                 "I132N", "I132C", "I133N", "I133C", "I134N")
  } else if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }
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
ins_cols_after <- function(df = NULL, idx_bf = ncol(df), idx_ins = NULL) {
  
  if (is.null(df)) stop("`df` cannot be `NULL`.", call. = FALSE)
  if (is.null(idx_ins)) return(df)
  if (idx_bf >= ncol(df)) return(df)

  col_nms <- names(df)[idx_ins]

  dplyr::bind_cols(
    df[, 1:idx_bf, drop = FALSE],
    df[, idx_ins, drop = FALSE],
    df[, (idx_bf + 1):ncol(df), drop = FALSE] %>% dplyr::select(-col_nms),
  )
}


#' Pad columns to a placeholder data frame.
#'
#' @param df The original data frame.
#' @param df2 The data frame to be inserted.
#' @param idx The index of \code{df} column for \code{df2} to be inserted
#'   (after).
add_cols_at <- function(df, df2, idx) {

  stopifnot(idx >= 0L)

  if (idx == 0L) {
    bf <- NULL
  } else {
    bf <- df[, seq_len(idx), drop = FALSE]
  }

  if ((idx + 1L) <= ncol(df)) {
    af <- df[, (idx + 1L) : ncol(df), drop = FALSE]
  } else {
    af <- NULL
  }

  dplyr::bind_cols(
    bf,
    df2,
    af,
  )
}


#' Replace columns in the original PSM table.
#'
#' The column index(es) need to be continuous.
#'
#' @param df The original data frame.
#' @param df2 The data columns to replace those in \code{df}.
#' @param idxs The sequences of column indexes in \code{df}. Note that
#'   \code{idxs} need to be a continuous sequences.
replace_cols_at <- function(df, df2, idxs) {

  ncol <- ncol(df)
  stopifnot(all(idxs >= 1L), all(idxs <= ncol))

  idxs <- sort(idxs)
  stopifnot(all.equal(idxs - idxs[1] + 1L, seq_along(idxs)))

  if (idxs[1] >= 2L) {
    bf <- df[, 1:(idxs[1]-1L), drop = FALSE]
  } else {
    bf <- NULL
  }

  if (idxs[length(idxs)] < ncol(df)) {
    af <- df[, (idxs[length(idxs)]+1L):ncol(df), drop = FALSE]
  } else {
    af <- NULL
  }

  dplyr::bind_cols(
    bf,
    df2,
    af
  )
}


#' Relocate column "to_move" immediately after column "col_before".
#'
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_before The anchor column to which the \code{to_move} will be moved
#'   after.
reloc_col_after <- function(df, to_move = "after_anchor", col_before = "anchor") {

  if (!(to_move %in% names(df) && col_before %in% names(df))) return(df)

  if (to_move == col_before) return(df)

  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))

  idx <- which(names(df) == col_before)
  df <- add_cols_at(df, df2, idx)

  return(df)
}


#' Relocate column "to_move" immediately after the last column.
#'
#' @inheritParams reloc_col_after
reloc_col_after_last <- function (df, to_move = "after_anchor") {
  col_last <- names(df)[ncol(df)]
  reloc_col_after(df, to_move, col_last)
}


#' Relocate column "to_move" immediately after the first column.
#'
#' @inheritParams reloc_col_after
reloc_col_after_first <- function(df, to_move = "after_anchor") {
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
                             col_after = "anchor") {

  if (!(to_move %in% names(df) && col_after %in% names(df))) return(df)

  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))

  idx <- which(names(df) == col_after)

  df <- add_cols_at(df, df2, idx - 1)

  return(df)
}


#' Relocate column "to_move" immediately before the last column.
#'
#' @inheritParams reloc_col_after
reloc_col_before_last <- function(df, to_move = "after_anchor") {
  col_last <- names(df)[ncol(df)]
  reloc_col_before(df, to_move, col_last)
}


#' Relocate column "to_move" immediately before the first column.
#'
#' @inheritParams reloc_col_after
reloc_col_before_first <- function (df, to_move = "after_anchor") {
  col_first <- names(df)[1]
  reloc_col_before(df, to_move, col_first)
}


#' Helper: find the column name before \code{to_move}.
#'
#' To keep columns at the same order after descriptive summary.
#'
#' @inheritParams reloc_col_after
find_preceding_colnm <- function(df, to_move) {

  if (!to_move %in% names(df)) {
    stop("Column ", to_move, " not found.",
         call. = FALSE)
  }

  ind_bf <- which(names(df) == to_move) - 1

  if (ind_bf == 0) {
    names(df)[1]
  } else {
    names(df)[ind_bf]
  }
}


#' Flatten lists recursively
#'
#' @param x Lists
recur_flatten <- function (x) {
  
  if (!inherits(x, "list")) {
    list(x)
  } else {
    unlist(c(lapply(x, recur_flatten)), recursive = FALSE)
  }
}


### End of also in proteoQ




###

#' Finds a file directory.
#'
#' With the option of creating a direcotry.
#'
#' @param path A file path.
#' @param create Logical; if TRUE create the path if not yet existed.
find_dir <- function (path, create = FALSE) {

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
create_dir <- function (path) {
  find_dir(path, create = TRUE)
}


#' Saves the arguments in a function call
#'
#' @param path A (system) file path.
#' @param fun The name of function being saved.
#' @param time The time stamp.
#' @importFrom rlang caller_env
save_call2 <- function(path, fun, time = NULL) {

  stopifnot(length(path) == 1L, length(fun) == 1L)

  call_pars <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE)
  call_pars[names(call_pars) == "..."] <- NULL

  if (is.null(time)) {
    p2 <- create_dir(path)
    save(call_pars, file = file.path(p2, paste0(fun, ".rda")))
  } else {
    stopifnot(length(time) == 1L)
    p2 <- create_dir(file.path(path, fun))
    save(call_pars, file = file.path(p2, paste0(time, ".rda")))
  }
}


#' Finds the setting of a arguments.
#'
#' Not currently used.
#'
#' @param arg Argument to be matched.
#' @inheritParams save_call2
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$%
find_callarg_val <- function (time = ".2021-05-21_211227",
                              path = "~/proteoM/.MSearch/Cache/Calls",
                              fun = "calc_pepmasses2", arg = "fasta") {

  stopifnot(length(arg) == 1L)

  find_callarg_vals(time, path, fun, arg) %>%
    purrr::flatten()
}


#' Finds the values of a list of arguments.
#'
#' @param args Arguments to be matched.
#' @inheritParams save_call2
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$%
#' @return The time stamp of a matched cache results.
find_callarg_vals <- function (time = NULL, path = NULL, fun = NULL,
                               args = NULL) {

  stopifnot(length(path) == 1L, length(fun) == 1L)

  if (is.null(time)) {
    file <- file.path(path, fun)
  } else {
    stopifnot(length(time) == 1L)
    file <- file.path(path, fun, time)
  }

  if (!file.exists(file)) {
    # waring(file, " not found.", call. = FALSE)
    return(NULL)
  }

  load(file = file)

  nots <- which(! args %in% names(call_pars))

  if (length(nots) > 0L) {
    stop("Arguments '", paste(args[nots], collapse = ", "),
         "' not found in the latest call to ", fun, call. = FALSE)
  }

  call_pars[args]
}


#' Finds the time stamp of a matched call from cached results.
#'
#' @param nms Names of arguments to be included in or excluded from matching.
#' @param type Logical; if TRUE, includes all arguments in \code{nms}. At FALSE,
#'   excludes all \code{nms}.
#' @inheritParams find_callarg_vals
#' @importFrom rlang caller_env
#' @return An empty object if no matches.
match_calltime <- function (path = "~/proteoM/.MSearch/Cache/Calls",
                            fun = "calc_pepmasses2",
                            nms = c("parallel", "out_path"),
                            type = c(TRUE, FALSE)) {

  stopifnot(length(path) == 1L, length(fun) == 1L)

  if (length(type) > 1L) type <- TRUE

  # current
  args <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE) %>%
    { if (type) .[names(.) %in% nms] else .[! names(.) %in% nms] }

  if (!length(args)) stop("Arguments for matching is empty.", call. = FALSE)

  args <- args %>% map(sort)

  times <- list.files(path = file.path(path, fun),
                      pattern = "\\.rda$",
                      all.files = TRUE)

  # cached
  cached <- map(times, find_callarg_vals, path = path, fun = fun,
                args = names(args))

  cached <- cached %>% map(~ map(.x, sort))

  # matched
  oks <- map_lgl(cached, identical, args)

  times[oks] %>% gsub("\\.rda$", "", .)
}


#' Matches the value of a function argument.
#'
#' @param f The function where one of its argument value will be matched.
#' @param arg The argument where its value will be matched.
#' @param val The value to be matched.
#' @param default The default value of \code{arg}.
#' @examples
#' \donttest{
#' foo <- function (quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16")) {
#'  val <- rlang::enexpr(quant)
#'  val <- match_valexpr(f = !!match.call()[[1]], arg = quant, val = !!val)
#' }
#'
#' default <- foo()
#' custom <- foo(quant = tmt6)
#' }
match_valexpr <- function (f = NULL, arg = NULL, val = NULL, default = NULL) {

  f <- rlang::enexpr(f)
  arg <- rlang::enexpr(arg)

  stopifnot(length(f) == 1L, length(arg) == 1L)

  f <- rlang::as_string(f)
  arg <- rlang::as_string(arg)
  ok_vals <- as.character(formals(f)[[arg]])[-1]

  # may be plural, not yet to string
  val <- rlang::enexpr(val)

  if (length(val) == 0L) {
    stop("`val` cannot be empty.", call. = FALSE)
  }

  if (length(val) > 1) {
    if (is.null(default)) {
      if (is.call(val[1])) {
        val <- val[[2]]
      } else {
        val <- val[[1]]
      }
    } else {
      stopifnot(length(default) == 1L)
      val <- default
    }
  } else {
    val <- rlang::as_string(val)
  }

  if (! val %in% ok_vals) {
    stop("`", val, "` is not an option for `", f, "`.", call. = FALSE)
  }

  val
}


#' Deletes files under a directory.
#'
#' The directory will be kept.
#' @param path A file path.
#' @param ignores The file extensions to be ignored.
delete_files <- function (path, ignores = NULL, ...) {

  dots <- rlang::enexprs(...)
  recursive <- dots[["recursive"]]

  if (is.null(recursive)) recursive <- TRUE

  stopifnot(is.logical(recursive))

  nms <- list.files(path = file.path(path), recursive = recursive,
                    full.names = TRUE, ...)

  if (!is.null(ignores)) {
    nms <- local({
      dirs <- list.dirs(path, full.names = FALSE, recursive = recursive) %>%
        .[! . == ""]

      idxes_kept <- dirs %>% purrr::map_lgl(~ any(grepl(.x, ignores)))

      nms_kept <- list.files(path = file.path(path, dirs[idxes_kept]),
                             recursive = TRUE, full.names = TRUE)

      nms %>% .[! . %in% nms_kept]
    })

    nms <- local({
      exts <- nms %>% gsub("^.*(\\.[^.]*)$", "\\1", .)
      idxes_kept <- map_lgl(exts, ~ any(grepl(.x, ignores)))

      nms[!idxes_kept]
    })
  }

  if (length(nms) > 0L) {
    suppressMessages(file.remove(file.path(nms)))
  }

  invisible(NULL)
}


#' Sums elements across lists.
#'
#' Each list has the same length. NA values are removed.
#'
#' @param x A numeric value.
#' @param y A numeric value.
`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))


#' Post processing after ms2match.
#'
#' @param out An output from various ms2match(es).
#' @inheritParams ms2match_base
post_ms2match <- function (out, i, aa_masses, out_path) {

  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)

  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)

  if (grepl("^rev", i)) {
    is_decoy <- TRUE
  } else {
    is_decoy <- FALSE
  }

  out <- out %>%
    dplyr::mutate(pep_fmod = nm_fmods,
                  pep_vmod = nm_vmods,
                  pep_mod_group = as.character(i)) %>%
    { if (is_decoy) dplyr::mutate(., pep_isdecoy = TRUE) else
      dplyr::mutate(., pep_isdecoy = FALSE) }

  out <- out %>%
    reloc_col_after("raw_file", "scan_num") %>%
    reloc_col_after("pep_mod_group", "raw_file") %T>%
    saveRDS(file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")))
}


#' Post frame advancing.
#'
#' Flattens mgfs within each frame (the number of entries equals to the number
#' of mgfs).
#' 
#' @param res Results from frame-advanced searches.
#' @param mgf_frames Data of MGF frames.
post_frame_adv <- function (res, mgf_frames) {

  res <- res %>% unlist(recursive = FALSE)

  empties <- purrr::map_lgl(res, purrr::is_empty)

  res <- do.call(rbind, mgf_frames) %>%
    dplyr::mutate(matches = res)

  res[!empties, ]
}


#' Subsets the search space.
#' 
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' @inheritParams post_ms2match
purge_search_space <- function (i, aa_masses, mgf_path, n_cores, ppm_ms1 = 20L,
                                fmods_nl = NULL) {

  # loads freshly mgfs (as will be modified)
  mgf_frames <- readRDS(file.path(mgf_path, "mgf_queries.rds")) %>%
    dplyr::group_by(frame) %>%
    dplyr::group_split() %>%
    setNames(purrr::map_dbl(., ~ .x$frame[1]))

  mgf_frames <- local({
    labs <- levels(cut(1:length(mgf_frames), n_cores^2))

    x <- cbind(
      lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
      upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))

    grps <- findInterval(1:length(mgf_frames), x[, 1])

    split(mgf_frames, grps)
  })

  # parses aa_masses
  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
  msg_end <- if (grepl("^rev_", i)) " (decoy)." else "."

  message("Matching against: ",
          paste0(nm_fmods,
                 nm_vmods %>% { if (nchar(.) > 0) paste0(" + ", .) else . },
                 msg_end))

  # reads theoretical peptide data
  .path_fasta <- get(".path_fasta", envir = .GlobalEnv)
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv)

  theopeps <- readRDS(file.path(.path_fasta, "pepmasses", .time_stamp,
                                paste0("binned_theopeps_", i, ".rds")))

  theopeps <- lapply(theopeps, function (x) {
    x[, c("pep_seq", "mass")]
  })

  # purged by neuloss residues
  # (not used)
  if (!is.null(fmods_nl)) {
    sites <- names(fmods_nl)
    pattern <- paste(sites, collapse = "|")
    theopeps <- sub_neuloss_peps(pattern, theopeps)

    rm(list = c("sites", "pattern"))
  }

  # (1) for a given aa_masses_all[[i]], some mgf_frames[[i]]
  #     may not be found in theopeps[[i]]
  mgf_frames <- lapply(mgf_frames, function (x) {
    oks <- names(x) %in% names(theopeps)
    x <- x[oks]
    
    empties <- purrr::map_lgl(x, purrr::is_empty)
    x[!empties]
  })

  # (2) splits `theopeps` in accordance to `mgf_frames` with
  #     preceding and following frames: (o)|range of mgf_frames[[1]]|(o)
  theopeps <- local({
    frames <- lapply(mgf_frames, function (x) as.integer(names(x)))

    mins <- purrr::map_dbl(frames, ~ {
      if (!length(.x)) x <- 0 else x <- min(.x, na.rm = TRUE)
    })
    
    maxs <- purrr::map_dbl(frames, ~ {
      if (!length(.x)) x <- 0 else x <- max(.x, na.rm = TRUE)
    })
    
    nms <- as.integer(names(theopeps))

    # (NAMEs are trivial node indexes from parallel processes)
    theopeps <- mapply(function (x, y) {
      theopeps[which(nms >= (x - 1L) & nms <= (y + 1L))]
    }, mins, maxs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    # --- may cause uneven length between `mgf_frames` and `theopeps`
    # empties <- map_lgl(theopeps, is_empty)
    # theopeps[!empties]
  })

  # (3) removes unused frames of `theopeps`
  # (NAMEs are trivial cluster node indexes...)
  theopeps <- mapply(subset_theoframes, mgf_frames, theopeps, 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)

  # (4) removes empties (zero overlap between mgf_frames and theopeps)
  oks <- purrr::map_lgl(mgf_frames, function (x) length(x) > 0L) |
    purrr::map_lgl(theopeps, function (x) length(x) > 0L)

  mgf_frames <- mgf_frames[oks]
  theopeps <- theopeps[oks]

  return(list(mgf_frames = mgf_frames, theopeps = theopeps))
}


#' Subsets theoretical peptides.
#'
#' Only entries containing the site of neuloss will be kept.
#'
#' @param pattern A regex of amino-acid residue(s).
#' @param theopeps Lists of theoretical peptides. A column of \code{pep_seq} is
#'   assumed.
sub_neuloss_peps <- function (pattern, theopeps) {

  rows <- lapply(theopeps, function (x) grepl(pattern, x$pep_seq))
  purrr::map2(theopeps, rows, ~ .x[.y, ])
}


#' Finds N-terminal mass.
#'
#' Not yet used.
#'
#' @inheritParams hms2_base
find_nterm_mass <- function (aa_masses) {

  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)

  if (length(ntmod)) {
    ntmass <- aa_masses[names(ntmod)] - 0.000549
  } else {
    ntmass <- aa_masses["N-term"] - 0.000549 # - electron
  }

  ntmass
}


#' Finds C-terminal mass.
#'
#' Not yet used.
#'
#' @inheritParams hms2_base
find_cterm_mass <- function (aa_masses) {

  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)

  if (length(ctmod)) {
    ctmass <- aa_masses[names(ctmod)] + 2.01510147
  } else {
    ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
  }
  
  ctmass
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
quick_rightjoin <- function (x, y, by = NULL) {

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
quick_leftjoin <- function (x, y, by = NULL) {

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
detect_cores <- function (max_n_cores = NULL) {

  n_cores <- parallel::detectCores()

  if (is.null(max_n_cores)) {
    max_n_cores <- n_cores
  } else {
    max_n_cores <- min(max_n_cores, n_cores)
  }

  if (n_cores > 128L) {
    max_n_cores <- min(max_n_cores, n_cores - 8L)
  } else if (n_cores <= 128L && n_cores > 64L) {
    max_n_cores <- min(max_n_cores, n_cores - 4L)
  } else if (n_cores <= 64L && n_cores > 32L) {
    max_n_cores <- min(max_n_cores, n_cores - 2L)
  } else if (n_cores <= 32L && n_cores > 16L) {
    max_n_cores <- min(max_n_cores, n_cores - 1L)
  } else {
    max_n_cores <- min(max_n_cores, n_cores)
  }

  invisible(max_n_cores)
}


#' Finds the amount of free system memory.
#' 
#' In the unit of MB.
find_free_mem <- function () {
  
  nm_os <- Sys.info()['sysname']
  
  if (nm_os == "Windows") {
    free_mem <- 
      system('wmic OS get FreePhysicalMemory /Value', intern=TRUE)[3] %>% 
      gsub("^FreePhysicalMemory=(\\d+)\\r", "\\1", .) %>% 
      as.numeric() %>% 
      `/`(1024)
    # tot_ram <- free_mem + memory.size(max = TRUE)
  } else {
    # not yet tested for "Linux", "Darwin"
    # inaccurate e.g. if physical RAM is 32000 but in .RProfile 
    # `invisible(utils::memory.limit(64000))`
    free_mem <- memory.limit() - memory.size(max = TRUE)
  }
  
  free_mem
}


#' Find the indexes of modifications.
#' 
#' @param out_path A output path.
#' @param file The file name where modifications are recorded.
find_mod_indexes <- function (out_path, file = "mod_indexes.txt") {
  
  filepath <- file.path(out_path, file)
  
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  mod_indexes <- readr::read_tsv(filepath, show_col_types = FALSE)
  
  inds <- mod_indexes$Abbr
  names(inds) <- mod_indexes$Desc
  inds
}


#' Finds if two sets are equal.
#' 
#' @param x A set.
#' @param y Another set.
equal_sets <- function(x, y) all(x %in% y) && all(y %in% x)


#' Finds the \code{aa_masses} corresponding to the base modification.
#' 
#' @param out_path An output path.
find_base_aamasses_index <- function (out_path) {
  
  aa_masses_all <- readRDS(file.path(out_path, "temp", "aa_masses_all.rds"))
  
  base_mods <- readr::read_tsv(file.path(out_path, "mod_indexes.txt"), 
                               show_col_types = FALSE)
  base_mods <- base_mods[base_mods$Coerced, "Desc"]
  
  vmods_ps <- lapply(aa_masses_all, attr, "vmods_ps")
  vmods <- lapply(vmods_ps, names)
  
  idx <- lapply(vmods, equal_sets, base_mods)
  idx <- simplify2array(idx)
}



