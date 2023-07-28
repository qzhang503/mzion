#' Visualization of matched MS2 ions.
#'
#' @param in_name An input file name of PSMs, such as psmQ.txt, psmC.txt,
#'   psmT2.txt, psmT3.txt etc.
#' @param out_name An output file name for saving the the plot of MS2 ions.
#' @param out_path A file path where outputs of \link{matchMS} can be
#'   identified.
#' @param scan A scan number that can be found from the outputs. Positive
#'   integer, Thermo's Orbitrap data; character string, Bruker's timsTOF data.
#'
#'   The searches of scan number will proceed automatically through
#'   \code{psmQ.txt},  \code{psmT2.txt}, \code{psmT3.txt} and/or
#'   \code{psmC.txt}.
#' @param raw_file Character string; a RAW file name.
#' @param rank Positive integer; the rank of a match. The default is one for the
#'   best match.
#' @param is_decoy Logical; is the match from decoy. The default is FALSE.
#' @param width Plot width.
#' @param height Plot height.
#' @rawNamespace import(ggplot2, except = c("%+%"))
#' @examples
#' \dontrun{
#' library(mzion)
#' 
#' ans <- mapMS2ions(
#'     out_path = "a/mzion/output/folder", in_name = "psmQ.txt", 
#'     out_name = "bar.png", raw_file = "a-raw-file-name.raw", scan = 9933, 
#'     rank = 1L, is_decoy = FALSE, type_ms2ions = "by"
#'   )
#'
#' # Custom plots
#' library(ggplot2)
#'
#' mgf <- ans$mgf
#' duo <- ans$duo
#' duo2 <- ans$duo2
#'
#' p <- ggplot() +
#'   geom_segment(mgf, mapping = aes(x = ms2_moverz, y = ms2_int,
#'                                          xend = ms2_moverz, yend = 0),
#'                       color = "gray", size = .1) +
#'   geom_segment(duo, mapping = aes(x = ms2_moverz, y = ms2_int,
#'                                            xend = ms2_moverz, yend = 0,
#'                                            color = type),
#'                         size = 1, show.legend = FALSE) +
#'   geom_text(duo, mapping = aes(x = ms2_moverz, y = ms2_int,
#'                                         label = label, color = type),
#'                      size = 6, alpha = .5, hjust = 0, angle = 90, vjust = 0,
#'                      nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE,
#'                      show.legend = FALSE) +
#'   labs(x = "m/z", y = "Intensity")
#' 
#' if (nrow(duo2)) {
#'   p <- p + 
#'     ggplot2::geom_segment(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
#'                                               xend = ms2_moverz, yend = 0, 
#'                                               color = type), 
#'                           size = .5, show.legend = FALSE) + 
#'     ggplot2::geom_text(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
#'                                            label = label, color = type),
#'                        size = 4, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
#'                        nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
#'                        show.legend = FALSE) 
#' }
#' 
#' pep <- paste0(duo$site[1:(nrow(duo)/2L)], collapse = " ")
#' p + annotate("text", -Inf, Inf, label = pep, hjust = -.2, vjust = 2)
#' 
#' }
#' 
#' @export
mapMS2ions <- function (out_path = NULL, in_name = "psmQ.txt", 
                        out_name = "bar.png", raw_file = "foo.raw", 
                        scan = 1234, rank = 1L, is_decoy = FALSE, 
                        type_ms2ions = "by", width = 12.5, height = 6) 
{
  if (is.null(out_path) || is.na(out_path) || out_path == "") {
    warning("\"out_path\" cannot be empty.", call. = FALSE)
    return(NULL)
  }
  
  if (is.null(raw_file) || is.na(raw_file) || raw_file == "" ) {
    warning("\"raw_file\" cannot be empty.", call. = FALSE)
    return(NULL)
  }
  
  if (is.null(in_name) || is.na(in_name) || in_name == "" ) {
    warning("\"in_name\" is empty; assume `psmQ.txt`.", call. = FALSE)
    in_name <- "psmQ.txt"
  }
  
  if (is.null(out_name) || out_name == "") {
    warning("\"out_name\" is empty; use `bar.png`.", call. = FALSE)
    out_name <- "bar.png"
  }
  
  out_name <- check_ggname(out_name)
  
  # MGF
  mgf_path <- match_mgf_path(out_path)
  raw_id <- match_raw_id(raw_file, mgf_path)
  scan <- as.character(scan)
  mgf_ok <- find_mgf_query(mgf_path, raw_id, scan)
  
  if (is.null(mgf_ok) || !nrow(mgf_ok)) {
    warning("MGF query not found.")
    return(NULL)
  }

  mgf <- data.frame(ms2_moverz = mgf_ok$ms2_moverz[[1]], 
                    ms2_int = mgf_ok$ms2_int[[1]])
  mgf$iex <- seq_len(nrow(mgf))
  
  ## PSMs
  cols_excl <- c("pep_ms2_moverzs", "pep_ms2_ints", "pep_ms2_theos", 
                 "pep_ms2_theos2", "pep_ms2_exptints", "pep_ms2_exptints2", 
                 "pep_n_matches", "pep_n_matches2")
  
  if (file.exists(fi_psm <- file.path(out_path, in_name))) {
    gl_vals <- ls(all.names = TRUE, envir = .GlobalEnv)
    ok_psms <- any(gl_vals == ".psms")

    ok_file <- if (any(gl_vals == ".psm_file"))
      identical(get(".psm_file", envir = .GlobalEnv), fi_psm)
    else
      FALSE

    if (ok_psms && ok_file) {
      .psms <- get(".psms", envir = .GlobalEnv)
    }
    else {
      # some columns in psmQ.txt not in psmC.txt
      .psms <- suppressWarnings(
        readr::read_tsv(fi_psm, show_col_types = FALSE, 
                        col_types = get_mzion_coltypes()))
      .psms <- .psms[, -which(names(.psms) %in% cols_excl), drop = FALSE]
      assign(".psms", .psms, envir = .GlobalEnv)
      assign(".psm_file", file.path(out_path, in_name), envir = .GlobalEnv)
    }
    
    psm <- .psms |>
      dplyr::filter(pep_scan_num == scan, 
                    .data$raw_file == .env$raw_file, 
                    pep_rank == rank,
                    pep_isdecoy == is_decoy) 
    psm <- psm[, -grep("^prot_", names(psm)), drop = FALSE]
    # can be duplicated by prot_accs
    psm <- unique(psm)
  }
  else {
    warning("PSM file not found: ", fi_psm)
    return(NULL)
  }
  
  if (!(nrow <- nrow(psm))) {
    warning("PSM entry not found. Check the correctness of scan number etc.")
    return(NULL)
  }
  
  if (nrow > 1L) {
    warning("Multiple PSMs found and the the first match being used.")
    psm <- psm[1, , drop = FALSE]
  }
  
  cols_duo <- c("ms2_moverz", "theo", "ms2_int", "site", "type", "label")
  aas <- strsplit(psm$pep_seq, "")[[1]]
  naa <- length(aas)
  
  ## Primary
  duo <- local({
    ion_types <- unlist(strsplit(type_ms2ions, ""))
    
    if (length(ion_types) != 2L)
      stop("Not a two-character `type_ms2ions = ", type_ms2ions, "`.", 
           call. = FALSE)
    
    cols_pri <- c("pep_ms2_deltas", "pep_ms2_ideltas", "pep_ms2_iexs")

    if (!all(oks <- cols_pri %in% names(psm))) {
      warning("PSM columns not found: ", paste(cols_pri[!oks], collapse = ", "), 
              "\nPlease use the latest version of mzion.")
      return(NULL)
    }

    theoexpt <- lapply(psm[, cols_pri], function (x) strsplit(x, ";")[[1]]) |>
      dplyr::bind_cols() |>
      dplyr::mutate(pep_ms2_deltas = as.numeric(pep_ms2_deltas)/1E3, 
                    pep_ms2_ideltas = as.integer(pep_ms2_ideltas), 
                    pep_ms2_iexs = as.integer(pep_ms2_iexs)) |>
      dplyr::rename(ith = pep_ms2_ideltas, iex = pep_ms2_iexs)
    
    duo <- data.frame(site = c(aas, aas[naa:1L]))
    duo$ith <- seq_len(nrow(duo))
    duo$type <- rep(ion_types, each = naa)
    idxes <- paste0("(", seq_len(naa), ")")
    duo$label <- unlist(lapply(ion_types, function (x) paste0(x, idxes)))
    
    duo <- duo |>
      dplyr::left_join(theoexpt, by = "ith")  |> 
      dplyr::arrange(ith) |>
      dplyr::left_join(mgf, by = "iex") |> 
      dplyr::mutate(theo = ms2_moverz + pep_ms2_deltas) |> 
      dplyr::select(cols_duo)
  })
  
  ## Secondary
  duo2 <- local({
    cols_sec <- c("pep_ms2_deltas2", "pep_ms2_ideltas2", "pep_ms2_iexs2")
    
    if (!all(oks <- cols_sec %in% names(psm))) {
      warning("PSM columns not found: ", paste(cols_sec[!oks], collapse = ", "), 
              "\nPlease use the latest version of mzion.")
      return(NULL)
    }

    theoexpt2 <- lapply(psm[, cols_sec], function (x) strsplit(x, ";")[[1]]) |>
      dplyr::bind_cols() |>
      dplyr::mutate(pep_ms2_deltas2 = as.numeric(pep_ms2_deltas2)/1E3, 
                    pep_ms2_ideltas2 = as.integer(pep_ms2_ideltas2), 
                    pep_ms2_iexs2 = as.integer(pep_ms2_iexs2)) |>
      dplyr::rename(ith = pep_ms2_ideltas2, iex = pep_ms2_iexs2) 
    
    ion_types2 <- find_secion_types(type_ms2ions)
    nsec   <- length(ion_types2)/2L
    type2  <- rep(ion_types2, each = naa)
    idxes2 <- paste0("(", seq_len(naa), ")")
    labs2  <- unlist(lapply(ion_types2, function (x) paste0(x, idxes2)))
    
    duo2 <- data.frame(site = c(rep(aas, nsec), rep(aas[naa:1L], nsec)), 
                       ith = seq_len(nsec * 2L * naa), 
                       type = type2, label = labs2)
    
    duo2 <- duo2 |>
      dplyr::left_join(theoexpt2, by = "ith")  |>
      dplyr::arrange(ith) |>
      dplyr::left_join(mgf, by = "iex") |>
      dplyr::mutate(theo = ms2_moverz + pep_ms2_deltas2) |>
      dplyr::select(cols_duo)
  })
  
  duos <- list(duo = duo, duo2 = duo2, mgf = mgf[, c("ms2_moverz", "ms2_int")])
  plotMS2ions(duos, out_path = out_path, out_name = out_name, width = width, 
              height = height)
}


#' Plots matched MS2 ions
#' 
#' @param duos The duo of MGF data and matched data.
#' @param out_path An output path.
#' @param out_name An output file name for saving the the plot of MS2 ions.
#' @param width Plot width.
#' @param height Plot height.
plotMS2ions <- function (duos, out_path = "~", out_name = "bar.png", 
                         width = 12.5, height = 9)
{
  duo  <- duos$duo
  duo2 <- duos$duo2
  mgf  <- duos$mgf
  
  ## Visualizations
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("\n======================================================", 
            "\nPackage \"ggplot2\" required for visualization.",
            "\n======================================================",
            call. = FALSE)
  }
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_segment(mgf, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                             xend = ms2_moverz, yend = 0), 
                          color = "gray", size = .1) + 
    ggplot2::geom_segment(duo, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                             xend = ms2_moverz, yend = 0, 
                                             color = type), 
                          size = 1, show.legend = FALSE) + 
    ggplot2::geom_text(duo, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                          label = label, color = type),
                       size = 6, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
                       nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
                       show.legend = FALSE) + 
    ggplot2::labs(x = "m/z", y = "Intensity")
  
  if (nrow(duo2)) {
    p <- p + 
      ggplot2::geom_segment(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                                xend = ms2_moverz, yend = 0, 
                                                color = type), 
                            size = .5, show.legend = FALSE) + 
      ggplot2::geom_text(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                             label = label, color = type),
                         size = 4, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
                         nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
                         show.legend = FALSE) 
  }
  
  pep <- paste0(duo$site[1:(nrow(duo)/2L)], collapse = " ")
  
  if (nrow(duo2)) 
    ymax <- max(duo$ms2_int, duo2$ms2_int, na.rm = TRUE)
  else 
    ymax <- max(duo$ms2_int, na.rm = TRUE)
  
  p <- p + 
    annotate("text", -Inf, Inf, label = pep, hjust = -.2, vjust = 2) + 
    ylim(0, ymax * 1.05)

  ggplot2::ggsave(file.path(out_path, out_name), width = width, height = height,
                  dpi = 300)

  ans <- list(mgf = mgf, duo = duo, duo2 = duo2, p = p)
  rds <- paste0(gsub("\\.[^.]*$", "", out_name), ".rds")
  saveRDS(ans, file.path(out_path, rds))

  invisible(ans)
}


#' Matches the value of mgf_path.
#' 
#' @param out_path An output path.
match_mgf_path <- function (out_path) 
{
  rda <- file.path(out_path, "Calls", "matchMS.rda")
  
  if (!file.exists(rda))
    stop("Parameter file not found: ", rda, call. = FALSE)
  
  load(rda)
  
  call_pars$mgf_path
}


#' Matches the id of raw_file.
#' 
#' @param mgf_path An MGF path.
#' @inheritParams mapMS2ions
match_raw_id <- function (raw_file, mgf_path) 
{
  file <- file.path(mgf_path, "raw_indexes.rds")
  
  if (!file.exists(file))
    stop("File not found ", file)
  
  raw_lookup <- qs::qread(file)
  raw_id <- unname(raw_lookup[raw_file])
  
  if (is.na(raw_id)) {
    stop(raw_file, " not found in ", file, ".\n", 
         "Aside from the possibility of incorrect `raw_file`, ",
         "have the folder name been changed?", 
         call. = FALSE)
  }
  
  raw_id
}


#' Adds \code{raw_ids}.
#'
#' Currently only used with psmC.txt during matchMS. An inverse function of
#' \link{match_raw_id}.
#'
#' @param df A PSM table.
#' @inheritParams matchMS
add_raw_ids <- function (df, mgf_path) 
{
  if (!"raw_file" %in% names(df))
    stop("Column `raw_file` not found.")
  
  raw_files <- unique(df$raw_file)
  raw_ids <- unlist(lapply(raw_files, match_raw_id, mgf_path))
  raws_lookup <- data.frame(raw_file = raw_files, raw_id = raw_ids)
  
  dplyr::left_join(df, raws_lookup, by = "raw_file")
}


#' Finds the types of secondary ions.
#' 
#' The order of secondary ions was defined in \link{add_seions}.
#'
#' @inheritParams matchMS
find_secion_types <- function (type_ms2ions = "by") 
{
  switch(type_ms2ions, 
         by = c("b2", "b*", "b*2", "b0", "b02", "y2", "y*", "y*2", "y0", "y02"), 
         cz = c("a2", "a*", "a*2", "a0", "a02", "x2"), 
         ax = c("c2", "z2"), 
         stop("Unknown type.", call. = FALSE))
}


#' Finds matched MGF query.
#' 
#' @param mgf_path An MGF path.
#' @param raw_id The index of raw_file.
#' @param scan A scan number or identifier.
#' @param to_global Logical; if TRUE, assigned to the Global environment.
find_mgf_query <- function (mgf_path, raw_id, scan, to_global = TRUE) 
{
  ok <- any(ls(all.names = TRUE, envir = .GlobalEnv) == ".mgf_queries")
  
  if (ok) {
    .mgf_queries <- get(".mgf_queries", envir = .GlobalEnv)
  }
  else {
    files <- list.files(path = file.path(mgf_path), 
                        pattern = "mgf_queries[_]*[0-9]*\\.rds$")
    
    if (!length(files)) {
      warning("No parsed `mgf_queries.rds` under ", mgf_path, call. = FALSE)
      return(NULL)
    }

    .mgf_queries <- lapply(files, function (x) qs::qread(file.path(mgf_path, x)))
    .mgf_queries <- do.call(rbind, .mgf_queries)
    .mgf_queries <- dplyr::mutate(.mgf_queries, scan_num = as.character(scan_num))

    if (to_global)
      assign(".mgf_queries", .mgf_queries, envir = .GlobalEnv)
    
    rm(list = "files")
  }
  
  mgf <- .mgf_queries |>
    dplyr::filter(raw_file == raw_id, scan_num == scan)

  if (!(nrow <- nrow(mgf))) {
    warning("MGF entries not found.")
  }
  else if (nrow > 1L) {
    warning("Multiple `mgf_query` matches and used the first one.")
    mgf <- mgf[1, , drop = FALSE]
  }
  
  mgf
}


#' Gets column types.
#' 
#' @import readr
#' @export
get_mzion_coltypes <- function () 
{
  col_types_pq <- cols(
    prot_acc = col_character(), 
    prot_issig = col_logical(), 
    prot_isess = col_logical(),
    prot_tier = col_integer(), 
    prot_hit_num = col_integer(), 
    prot_family_member = col_integer(), 
    prot_es = col_number(), 
    prot_es_co = col_number(), 
    pep_seq = col_character(), 
    pep_n_ms2 = col_integer(), 
    pep_scan_title = col_character(), 
    pep_exp_mz = col_number(),
    pep_exp_mr = col_number(), 
    pep_exp_z = col_character(), 
    pep_calc_mr = col_number(), 
    pep_delta = col_number(),
    pep_tot_int = col_number(), 
    pep_ret_range = col_number(), 
    pep_scan_num = col_character(), # timsTOF
    pep_mod_group = col_integer(), 
    pep_fmod = col_character(),
    pep_vmod = col_character(),
    pep_isdecoy = col_logical(),
    pep_ivmod = col_character(),
    pep_len = col_integer(), 
    
    pep_ms2_moverzs = col_character(),
    pep_ms2_ints = col_character(), 
    pep_ms2_theos = col_character(),
    pep_ms2_theos2 = col_character(),
    pep_ms2_exptints = col_character(),
    pep_ms2_exptints2 = col_character(),
    
    pep_n_matches = col_integer(), 
    pep_n_matches2 = col_integer(),
    pep_ms2_deltas = col_character(),
    pep_ms2_ideltas = col_character(),
    pep_ms2_deltas2 = col_character(), 
    pep_ms2_ideltas2 = col_character(),
    pep_ms2_deltas_mean = col_double(),
    pep_ms2_deltas_sd = col_double(),
    
    pep_issig = col_logical(),
    pep_score = col_double(),
    # pep_score_co = col_double(),
    pep_rank = col_integer(), 
    pep_locprob = col_double(),
    pep_locdiff = col_double(),
    # pep_rank_nl = col_integer(), 
    pep_literal_unique = col_logical(),
    pep_razor_unique = col_logical(),
    raw_file = col_character(), 
  )
  
  nms <- names(col_types_pq$cols)
  
  stopifnot(length(nms) == length(unique(nms)))
  
  col_types_pq
}


#' Checks file names for ggsave
#' 
#' The same as proteoQ::gg_imgname.
#' 
#' @param filename Character string; An output file name.
check_ggname <- function(filename) 
{
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  exts <- c("png", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg") 
  
  if(!fn_suffix %in% exts) {
    warning("Unrecognized file extenstion: '", fn_suffix, 
            "'. Image will be saved as a '.png'.", 
            call. = FALSE)
    
    fn_suffix <- "png"
  }
  
  paste0(fn_prefix, ".", fn_suffix)
}
