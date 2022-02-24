#' Visualization of matched MS2 ions.
#'
#' Secondary ions of b0, y0 etc. were also used in ion matches though not
#' directly documented.
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
#'                   is_decoy = FALSE)
#'
#' # Custom plots
#' library(ggplot2)
#'
#' mgf <- ans$mgf
#' duo <- ans$duo
#'
#' ggplot() +
#' geom_segment(mgf, mapping = aes(x = ms2_moverz, y = ms2_int,
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
#' }
#' @export
mapMS2ions <- function (out_path = "~/proteoM/outs", scan = 1234, 
                        raw_file = "foo.raw", rank = 1L, is_decoy = FALSE, 
                        type_ms2ions = "by", filename = "bar.png") 
{
  ## MGF path
  mgf_path <- local({
    rda <- file.path(out_path, "Calls", "matchMS.rda")
    
    if (!file.exists(rda))
      stop("Parameter file not found: ", rda, call. = FALSE)
    
    load(rda)
    
    call_pars$mgf_path
  })
  
  ## RAW file look-up
  raw_id <- local({
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
  })
  
  fileQ <- list.files(path = file.path(out_path), pattern = "psmQ.*\\.txt$")
  fileT2 <- list.files(path = file.path(out_path), pattern = "psmT2.*\\.txt$")
  fileT3 <- list.files(path = file.path(out_path), pattern = "psmT3.*\\.txt$")
  fileC <- list.files(path = file.path(out_path), pattern = "psmC.*\\.txt$")
  
  # Required: psmQ, psmC; 
  # Optional: psmT2, psmT3
  lapply(c(fileQ, fileT2, fileT3, fileC), function (x) {
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
  })
  
  file_t1 <- file.path(out_path, fileQ)
  file_t2 <- file.path(out_path, fileT2)
  file_t3 <- file.path(out_path, fileT3)
  file_t0 <- file.path(out_path, fileC)
  
  col_types_pq <- cols(pep_scan_range = col_character()) # timsTOF
  scan <- as.character(scan)

  ## PSMs
  psm <- local({
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
  })
  
  nrow <- nrow(psm)
  
  if (!nrow) {
    psm <- local({
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
    })
    
    nrow <- nrow(psm)
  }
  
  if (!nrow)
    stop("PSM entry not found.", call. = FALSE)

  if (nrow > 1L) {
    warning("Multiple PSMs and use the first match.")
    psm <- psm[1, ]
  }
  
  pep_seq <- psm$pep_seq
  pep_ivmod <- psm$pep_ivmod
  
  ## Matches
  theoexpt <- local({
    mod <- psm$pep_mod_group
    nm_i <- paste0(".ion_matches_", mod)
    
    ok <- any(ls(all.names = TRUE, envir = .GlobalEnv) == nm_i)
    
    if (ok)
      .ion_matches <- get(nm_i, envir = .GlobalEnv)
    else {
      file <- file.path(out_path, "temp", paste0("ion_matches_", mod, ".rds"))
      
      if (!file.exists(file))
        stop("Ion matches not found: ", file, call. = FALSE)
      
      .ion_matches <- readRDS(file) %>% 
        dplyr::mutate(scan_num = as.character(scan_num))
      
      assign(nm_i, .ion_matches, envir = .GlobalEnv)
    }

    ion_match <- .ion_matches %>% 
      dplyr::filter(scan_num == scan, 
                    raw_file == raw_id, 
                    pep_isdecoy == is_decoy)
    
    nrow <- nrow(ion_match)
    
    if (!nrow)
      stop("PSM not found.", call. = FALSE)

    if (nrow > 1L) {
      warning("Multiple PSMs found and the the first match being used.")
      ion_match <- ion_match[1, ]
    }
    
    # theo-expt pairs
    theoexpt <- ion_match$matches[[1]]
    
    # by pep_seq
    theoexpt <- theoexpt[names(theoexpt) == pep_seq]
    
    if (length(theoexpt) > 1L) {
      warning("Multiple `pep_seq` matches and used the first one.")
      theoexpt <- theoexpt[1]
    }

    theoexpt <- theoexpt[[1]]
    
    # by pep_ivmod
    theoexpt <- theoexpt[names(theoexpt) == pep_ivmod]
    
    if (length(theoexpt) > 1L) {
      warning("Multiple `pep_ivmod` matches and used the first one.")
      theoexpt <- theoexpt[1]
    }
    
    theoexpt <- theoexpt[[1]]
  })
  
  # MGFs
  mgf <- local({
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
    
    if (!nrow)
      stop("Corresponding MGF entries not found in ", 
           paste(files, collapse = ", "), 
           call. = FALSE)
    else if (nrow > 1L) {
      warning("Multiple `mgf_query` matches and used the first one.")
      mgf <- mgf[1, ]
    }
    
    mgf
  })
  
  ms2_moverz <- mgf$ms2_moverz[[1]]
  ms2_int <- mgf$ms2_int[[1]]
  
  ans <- list(theo = theoexpt$theo, 
              expt = theoexpt$expt, 
              ms2_moverz = ms2_moverz, 
              ms2_int = ms2_int)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("\n======================================================", 
            "\nPackage \"ggplot2\" required for visualization.",
            "\n======================================================",
            call. = FALSE)
  }
  
  ans <- local({
    mgf <- data.frame(ms2_moverz = ans$ms2_moverz, ms2_int = ans$ms2_int)

    duo <- local({
      duo <- data.frame(ms2_moverz = ans$expt, theo = ans$theo) %>% 
        dplyr::left_join(mgf, by = "ms2_moverz") %>% 
        dplyr::mutate(site = names(ans$expt))
      
      ion_types <- unlist(strsplit(type_ms2ions, ""))
      n <- 2L
      
      if (!(length(ion_types) == n))
        stop("Not two characters `type_ms2ions = ", type_ms2ions, "`.", 
             call. = FALSE)
      
      len <- nrow(duo)/n
      type <- rep(ion_types, each = len)
      idxes <- paste0("(", 1:len, ")")
      labs <- c(paste0(ion_types[1], idxes), paste0(ion_types[2], idxes))
      
      duo %>% 
        dplyr::mutate(type = type, label = labs)
    })

    ggplot2::ggplot() + 
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

    ggplot2::ggsave(file.path(out_path, filename), width = 8, height = 5)
    
    list(mgf = mgf, duo = duo)
  })
  
  invisible(ans)
}


