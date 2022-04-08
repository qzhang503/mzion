#' Visualization of matched MS2 ions.
#'
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
#' @param filename An output file name for saving the the plot of MS2 ions.
#' @inheritParams matchMS
#' @rawNamespace import(ggplot2, except = c("%+%"))
#' @examples
#' \donttest{
#' ans <- mapMS2ions(out_path = "a/proteoM/output/folder",
#'                   scan = 9933,
#'                   raw_file = "a-raw-file-name.raw",
#'                   rank = 1L,
#'                   is_decoy = FALSE, 
#'                   filename = "bar.png")
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
#'                         size = .4, show.legend = FALSE) +
#'   geom_text(duo, mapping = aes(x = ms2_moverz, y = ms2_int,
#'                                         label = label, color = type),
#'                      size = 3, alpha = .5, hjust = 0, angle = 90, vjust = 0,
#'                      nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE,
#'                      show.legend = FALSE) +
#'   labs(x = "m/z", y = "Intensity")
#' 
#' if (nrow(duo2)) {
#'   p <- p + 
#'     ggplot2::geom_segment(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
#'                                               xend = ms2_moverz, yend = 0, 
#'                                               color = type), 
#'                           size = .2, show.legend = FALSE) + 
#'     ggplot2::geom_text(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
#'                                            label = label, color = type),
#'                        size = 2, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
#'                        nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
#'                        show.legend = FALSE) 
#' }
#' 
#' }
#' @export
mapMS2ions <- function (out_path = "~/proteoM/outs", scan = 1234, 
                        raw_file = "foo.raw", rank = 1L, is_decoy = FALSE, 
                        type_ms2ions = "by", filename = NULL) 
{
  if (is.null(filename)) 
    filename <- "bar.png"
  
  filename <- as.character(substitute(filename))
  raw_file <- as.character(substitute(raw_file))
  type_ms2ions <- as.character(substitute(type_ms2ions))

  mgf_path <- match_mgf_path(out_path)
  raw_id <- match_raw_id(raw_file, mgf_path)
  scan <- as.character(scan)

  fileQ <- list.files(path = file.path(out_path), pattern = "^psmQ.*\\.txt$")
  fileT2 <- list.files(path = file.path(out_path), pattern = "^psmT2.*\\.txt$")
  fileT3 <- list.files(path = file.path(out_path), pattern = "^psmT3.*\\.txt$")
  fileC <- list.files(path = file.path(out_path), pattern = "^psmC.*\\.txt$")
  
  lapply(c(fileQ, fileT2, fileT3, fileC), check_existed_psms)

  file_t1 <- file.path(out_path, fileQ)
  file_t2 <- file.path(out_path, fileT2)
  file_t3 <- file.path(out_path, fileT3)
  file_t0 <- file.path(out_path, fileC)

  ## PSMs
  psm <- find_psm_rows(file_t0 = file_t0, file_t1 = file_t1, file_t2 = file_t2, 
                       file_t3 = file_t3, scan = scan, raw_file = raw_file, 
                       rank = rank, is_decoy = is_decoy)

  ## Matches
  theoexpt_pair <- find_theoexpt_pair(psm = psm, out_path = out_path, 
                                      scan = scan, raw_id = raw_id, 
                                      is_decoy = is_decoy)
  theoexpt <- theoexpt_pair$theoexpt
  theoexpt2 <- theoexpt_pair$theoexpt2
  
  ## MGFs
  mgf_ok <- find_mgf_query(mgf_path, raw_id, scan)

  ## Trio: theo-expt-mgf
  th_ex_mgf <- list(theo = theoexpt$theo, 
                    expt = theoexpt$expt, 
                    theo2 = theoexpt2$theo, 
                    expt2 = theoexpt2$expt, 
                    ms2_moverz = mgf_ok$ms2_moverz[[1]], 
                    ms2_int = mgf_ok$ms2_int[[1]])

  duos <- combine_prisec_matches(th_ex_mgf, type_ms2ions)
  duo <- duos$duo
  duo2 <- duos$duo2
  mgf <- duos$mgf

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
                          size = .4, show.legend = FALSE) + 
    ggplot2::geom_text(duo, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                          label = label, color = type),
                       size = 3, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
                       nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
                       show.legend = FALSE) + 
    ggplot2::labs(x = "m/z", y = "Intensity")
  
  if (nrow(duo2)) {
    p <- p + 
      ggplot2::geom_segment(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                                xend = ms2_moverz, yend = 0, 
                                                color = type), 
                            size = .2, show.legend = FALSE) + 
      ggplot2::geom_text(duo2, mapping = aes(x = ms2_moverz, y = ms2_int, 
                                             label = label, color = type),
                         size = 2, alpha = .5, hjust = 0, angle = 90, vjust = 0, 
                         nudge_x = 0.05, nudge_y = 0.05, na.rm = TRUE, 
                         show.legend = FALSE) 
  }
  
  ggplot2::ggsave(file.path(out_path, filename), width = 8, height = 5)
  
  invisible(list(mgf = mgf, duo = duo, duo2 = duo2))
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
  
  raw_lookup <- readRDS(file)
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


#' Extracts the first row of matched PSMs.
#' 
#' @param file_t0 The filename of psmC results.
#' @param file_t1 The filename of tier-1 PSMs.
#' @param file_t2 The filename of tier-1 PSMs.
#' @param file_t3 The filename of tier-1 PSMs.
#' @param scan A scan number or identifier.
#' @inheritParams mapMS2ions
find_psm_rows <- function (file_t0, file_t1, file_t2, file_t3, scan, raw_file, 
                           rank = 1L, is_decoy = FALSE) 
{
  psm <- find_psm_rows1(file_t1 = file_t1, file_t2 = file_t2, file_t3 = file_t3, 
                        scan = scan, raw_file = raw_file, rank = rank, 
                        is_decoy = is_decoy)
  
  nrow <- nrow(psm)
  
  if (!nrow) {
    psm <- find_psm_rows2(file_t0 = file_t0, scan = scan, raw_file = raw_file, 
                          rank = rank, is_decoy = is_decoy)
    nrow <- nrow(psm)
  }
  
  if (!nrow)
    stop("PSM entry not found.", call. = FALSE)
  
  if (nrow > 1L) {
    warning("Multiple PSMs matched (e.g. at different neutral losses);", 
            " used the first match.")
    psm <- psm[1, ]
  }
  
  invisible(psm)
}


#' Extracts the first row of matched PSMs from tiers 1-3.
#' 
#' @inheritParams find_psm_rows
find_psm_rows1 <- function (file_t1, file_t2, file_t3, scan, raw_file, 
                            rank = 1L, is_decoy = FALSE) 
{
  ok <- any(ls(all.names = TRUE, envir = .GlobalEnv) == ".psms")
  
  if (ok) {
    .psms <- get(".psms", envir = .GlobalEnv)
  }
  else {
    psms_1 <- readr::read_tsv(file_t1, show_col_types = FALSE)
    
    if (length(file_t2)) 
      psms_2 <- readr::read_tsv(file_t2, show_col_types = FALSE)
    else 
      psms_2 <- NULL
    
    if (length(file_t3)) 
      psms_3 <- readr::read_tsv(file_t3, show_col_types = FALSE)
    else 
      psms_3 <- NULL
    
    # psms_1, _2 and _3 can have overlaps
    # (the same pep_seq to different prot_accs and thus different tiers)
    .psms <- dplyr::bind_rows(psms_1, psms_2, psms_3) %>% 
      dplyr::mutate(pep_scan_num = as.character(pep_scan_num))
    
    assign(".psms", .psms, envir = .GlobalEnv)
  }
  
  .psms %>% 
    dplyr::filter(pep_scan_num == scan, 
                  .data$raw_file == .env$raw_file, 
                  pep_rank == rank,
                  pep_isdecoy == is_decoy) %>% 
    # can be duplicated by prot_accs
    # no pep_start and pep_end in the data.frame -> no duplication per se
    dplyr::select(-grep("^prot_", names(.))) %>% 
    unique()
}


#' Extracts the first row of matched PSMs from psmC.
#' 
#' @inheritParams find_psm_rows
find_psm_rows2 <- function (file_t0, scan, raw_file, rank = 1L, 
                            is_decoy = FALSE) 
{
  ok <- any(ls(all.names = TRUE, envir = .GlobalEnv) == ".psmC")
  
  if (ok) {
    .psms <- get(".psmC", envir = .GlobalEnv)
  }
  else {
    .psms <- readr::read_tsv(file_t0, show_col_types = FALSE) %>% 
      dplyr::mutate(pep_scan_num = as.character(pep_scan_num))
    
    assign(".psmC", .psms, envir = .GlobalEnv)
  }
  
  .psms %>% 
    dplyr::filter(pep_scan_num == scan, 
                  .data$raw_file == .env$raw_file, 
                  pep_rank == rank,
                  pep_isdecoy == is_decoy) %>% 
    dplyr::select(-grep("^prot_", names(.))) %>% 
    unique()
}


#' Finds the pairs of theoretical and experimental values.
#' 
#' @param psm Matched PSM.
#' @param out_path An output path.
#' @param raw_id The index of raw_file.
#' @param scan A scan number or identifier.
#' @inheritParams mapMS2ions
find_theoexpt_pair <- function (psm, out_path, scan, raw_id, is_decoy = FALSE) 
{
  pep_seq <- psm$pep_seq
  pep_ivmod <- psm$pep_ivmod
  
  mod <- psm$pep_mod_group
  nm <- paste0(".ion_matches_", mod)
  nm2 <- paste0(".list_table_", mod)
  
  ok <- local({
    objs <- ls(all.names = TRUE, envir = .GlobalEnv)
    all(c(nm, nm2) %in% objs)
  })
  
  if (ok) {
    .ion_matches <- get(nm, envir = .GlobalEnv)
    .list_table <- get(nm2, envir = .GlobalEnv)
  }
  else {
    file <- file.path(out_path, "temp", paste0("ion_matches_", mod, ".rds"))
    file2 <- file.path(out_path, "temp", paste0("list_table_", mod, ".rds"))
    
    if (!file.exists(file))
      stop("Ion matches not found: ", file, call. = FALSE)
    
    if (!file.exists(file2))
      stop("Secondary ion matches not found: ", file, call. = FALSE)
    
    .ion_matches <- readRDS(file) %>% 
      dplyr::mutate(scan_num = as.character(scan_num))
    
    .list_table <- readRDS(file2) %>% 
      dplyr::mutate(scan_num = as.character(scan_num))
    
    assign(nm, .ion_matches, envir = .GlobalEnv)
    assign(nm2, .list_table, envir = .GlobalEnv)
  }
  
  ion_match <- .ion_matches %>% 
    dplyr::filter(scan_num == scan, 
                  raw_file == raw_id, 
                  pep_isdecoy == is_decoy)
  
  list_table <- .list_table %>% 
    dplyr::filter(scan_num == scan, 
                  raw_file == raw_id, )
  
  nrow <- nrow(ion_match)
  
  if (!nrow)
    stop("PSM not found.", call. = FALSE)
  
  if (nrow > 1L) {
    warning("Multiple PSMs found and the the first match being used.")
    ion_match <- ion_match[1, ]
  }
  
  # theo-expt pairs (unlist from list table)
  theoexpt <- ion_match$matches[[1]]
  
  # (1) matched by `pep_seq`
  theoexpt <- theoexpt[names(theoexpt) == pep_seq]
  
  if (length(theoexpt) > 1L) {
    warning("Multiple `pep_seq` matches and used the first one.")
    theoexpt <- theoexpt[1]
  }
  
  # unlist from the current list
  theoexpt <- theoexpt[[1]]
  
  # (2) matched by pep_ivmod
  # (pep_seq matched but can still have multiple pep_ivmod's)
  theoexpt <- theoexpt[names(theoexpt) == pep_ivmod]

  # (can have multiple NLs)
  if (length(theoexpt) > 1L) {
    warning("Multiple `pep_ivmod` matches and used the first one.")
    theoexpt <- theoexpt[1]
  }
  
  # unlist from the current list
  theoexpt <- theoexpt[[1]]
  
  ## list_table
  # (1) match by identical primary matches
  list_table <- local({
    # unlist
    pris <- lapply(list_table$pri_matches, `[[`, 1)
    
    # matched by the same ion ladder of `theo`s
    # (e.g. will exclude ladders at unmatched `pep_ivmod`)
    rows <- unlist(lapply(pris, function (x) identical(x$theo, theoexpt$theo)), 
                   use.names = FALSE)
    
    list_table <- list_table[rows, ]
    
    if (nrow(list_table) > 1L) {
      warning("More than one primary match; used the first one.")
      list_table <- list_table[1, ]
    }
    
    list_table
  })
  
  theoexpt2 <- list_table$sec_matches[[1]]
  
  # may still need to unlist
  while(is.list(theoexpt2) && length(theoexpt2) == 1L) {
    theoexpt2 <- theoexpt2[[1]]
  }
  
  names(theoexpt2$expt) <- names(theoexpt2$theo)
  
  list(theoexpt = theoexpt, 
       theoexpt2 = theoexpt2)
}


#' Finds matched MGF query.
#' 
#' @param mgf_path An MGF path.
#' @param raw_id The index of raw_file.
#' @param scan A scan number or identifier.
find_mgf_query <- function (mgf_path, raw_id, scan) 
{
  ok <- any(ls(all.names = TRUE, envir = .GlobalEnv) == ".mgf_queries")
  
  if (ok)
    .mgf_queries <- get(".mgf_queries", envir = .GlobalEnv)
  else {
    files <- list.files(path = file.path(mgf_path), 
                        pattern = "mgf_queries[_]*[0-9]*\\.rds$")
    
    if (!length(files))
      stop("No parsed `mgf_queries.rds` under ", mgf_path, call. = FALSE)
    
    .mgf_queries <- 
      lapply(files, function (x) readRDS(file.path(mgf_path, x))) %>% 
      do.call(rbind, .) %>% 
      dplyr::mutate(scan_num = as.character(scan_num))
    
    assign(".mgf_queries", .mgf_queries, envir = .GlobalEnv)
  }
  
  mgf <- .mgf_queries %>% 
    dplyr::filter(raw_file == raw_id, 
                  scan_num == scan)
  
  nrow <- nrow(mgf)
  
  if (!nrow) {
    stop("Corresponding MGF entries not found in ", 
         paste(files, collapse = ", "), 
         call. = FALSE)
  }
  else if (nrow > 1L) {
    warning("Multiple `mgf_query` matches and used the first one.")
    mgf <- mgf[1, ]
  }
  
  mgf
}


#' Combines primary and secondary matches.
#' 
#' @param th_ex_mgf The results of theoretical, experimental and MGF.
#' @inheritParams matchMS
combine_prisec_matches <- function (th_ex_mgf, type_ms2ions = "by") 
{
  mgf <- data.frame(ms2_moverz = th_ex_mgf$ms2_moverz, 
                    ms2_int = th_ex_mgf$ms2_int)
  
  # primary
  duo <- data.frame(ms2_moverz = th_ex_mgf$expt, theo = th_ex_mgf$theo) %>% 
    dplyr::left_join(mgf, by = "ms2_moverz") %>% 
    dplyr::mutate(site = names(th_ex_mgf$expt))
  
  ion_types <- unlist(strsplit(type_ms2ions, ""))
  n <- 2L
  
  if (!(length(ion_types) == n))
    stop("Not a two-character `type_ms2ions = ", type_ms2ions, "`.", 
         call. = FALSE)
  
  len <- nrow(duo)/n
  type <- rep(ion_types, each = len)
  idxes <- paste0("(", 1:len, ")")
  labs <- unlist(lapply(ion_types, function (x) paste0(x, idxes)))
  
  duo <- duo %>% dplyr::mutate(type = type, label = labs)

  # secondary
  if (!all(is.na(th_ex_mgf$expt2))) {
    duo2 <- data.frame(ms2_moverz = th_ex_mgf$expt2, theo = th_ex_mgf$theo2) %>% 
      dplyr::left_join(mgf, by = "ms2_moverz") %>% 
      dplyr::mutate(site = names(th_ex_mgf$expt2))
    
    ion_types2 <- find_secion_types(type_ms2ions)
    type2 <- rep(ion_types2, each = len)
    idxes2 <- idxes
    labs2 <- unlist(lapply(ion_types2, function (x) paste0(x, idxes2)))
    
    duo2 <- duo2 %>% dplyr::mutate(type = type2, label = labs2)
  }
  else {
    duo2 <- duo[0, ]
  }
  
  list(duo = duo, duo2 = duo2, mgf = mgf)
}


#' Checks the existing PSM files.
#' 
#' Required: psmQ, psmC; Optional: psmT2, psmT3
#' 
#' @param x A PSM file type.
check_existed_psms <- function (x) 
{
  type <- gsub("^(psm.*)\\.txt$", "\\1", x)
  len <- length(x)
  
  if (!len) {
    if (type %in% c("psmQ", "psmC")) {
      stop("No `", type, "` results.", call. = FALSE)
    }
  } 
  else if (len > 1L) {
    stop("No more than one `", type, "` results.", call. = FALSE)
  }
}
