#' Helper in loading MGFs.
#'
#' Currently, \code{index_ms2} is always FALSE.
#' 
#' @param min_mass A minimum mass of precursors for considerations.
#' @param max_mass A maximum mass of precursors for considerations.
#' @param min_ms2mass A minimum m/z of MS2 ions for considerations.
#' @param index_ms2 Logical; if TRUE, converts MS2 m/z values to indices.
#' @param is_ms1_three_frame Logical; is the searches by the three frames of
#'   preceeding, current and following.
#' @param is_ms2_three_frame Logical; is the searches by the three frames of
#'   preceeding, current and following.
#' @inheritParams matchMS
load_mgfs <- function (out_path, mgf_path, min_mass = 500L, max_mass = 6000L, 
                       min_ms2mass = 110L, topn_ms2ions = 100L, 
                       min_ms1_charge = 2L, max_ms1_charge = 6L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1 = 20L, ppm_ms2 = 25L, index_ms2 = FALSE, 
                       is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE) 
{
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    },
    add = TRUE
  )

  # ---
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  
  args_except <- NULL
  
  cache_pars <- find_callarg_vals(
    time = NULL, 
    path = file.path(out_path, "Calls"), 
    fun = paste0(fun, ".rda"), 
    args = names(formals(fun)) %>% 
      .[! . %in% args_except]
  ) 
  
  cache_pars <- cache_pars[sort(names(cache_pars))]
  
  call_pars <- mget(
    names(formals()) %>% .[! . %in% args_except], 
    envir = fun_env, 
    inherits = FALSE
  ) 
  
  call_pars <- call_pars[sort(names(call_pars))]

  ok_pars <- identical(call_pars, cache_pars)

  rds <- file.path(mgf_path, "mgf_queries.rds")
  
  if (ok_pars && file.exists(rds)) {
    message("Found cached MGFs: `", rds, "`.")
    .savecall <- FALSE
  } 
  else {
    message("Processing raw MGFs.")
    
    ppm_ms1_new <- if (is_ms1_three_frame) 
      as.integer(ceiling(ppm_ms1 * .5))
    else 
      ppm_ms1

    ppm_ms2_new <- if (is_ms2_three_frame) 
      as.integer(ceiling(ppm_ms2 * .5))
    else 
      ppm_ms2

    delete_files(
      out_path, 
      ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                  "\\.xls$", "\\.csv$", "\\.txt$", 
                  "^mgf$", "^mgfs$", "Calls")
    )

    readMGF(filepath = mgf_path,
            min_mass = min_mass,
            max_mass = max_mass, 
            min_ms2mass = min_ms2mass,
            topn_ms2ions = topn_ms2ions,
            ms1_charge_range = c(min_ms1_charge, max_ms1_charge), 
            ms1_scan_range = c(min_scan_num, max_scan_num), 
            ret_range = c(min_ret_time, max_ret_time),
            ppm_ms1 = ppm_ms1_new,
            ppm_ms2 = ppm_ms2_new,
            index_ms2 = index_ms2,
            out_path = rds)
    
    .savecall <- TRUE
  }

  invisible(NULL)
}


#' Reads MGF files in chunks.
#'
#' @param filepath The file path to a list of mgf files.
#' @param min_mass Numeric; the minimum mass of MS1 species. The value needs to
#'   match the one in  \link{binTheoSeqs}.
#' @param topn_ms2ions A non-negative integer; the top-n species for uses in
#'   MS2 ion searches. The default is to use the top-100 ions in an MS2 event.
#' @param ms1_charge_range The range of MS1 charge states.
#' @param ms1_scan_range The range of MS1 scan numbers.
#' @param ret_range The range of retention time in seconds.
#' @param out_path An output path.
#' @inheritParams load_mgfs
#' @inheritParams matchMS
#' @inheritParams frames_adv
#' @import stringi
#' @import readr
#' @import fs
#' @examples
#' \donttest{
#' mgf_queries <- proteoM:::readMGF()
#' }
readMGF <- function (filepath = "~/proteoM/mgf",
                     min_mass = 500L, max_mass = 6000L, min_ms2mass = 110L, 
                     topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                     ms1_scan_range = c(1L, .Machine$integer.max), 
                     ret_range = c(0, Inf), 
                     ppm_ms1 = 20L, ppm_ms2 = 25L, index_ms2 = FALSE,
                     out_path = file.path(filepath, "mgf_queries.rds")) 
{
  f <- function(x, pos) {
    nm <- file.path(filepath, "temp", paste0("chunk", "_", pos, ".mgf"))
    writeLines(x, nm)
  }

  # parsing rules
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (!length(filelist)) 
    stop("No '.mgf' files under ", filepath, call. = FALSE)

  pat_mgf <- find_mgf_type(file.path(filepath, filelist[[1]]))
  
  type_mgf <- pat_mgf$type
  n_bf_begin <- pat_mgf$n_bf_begin
  n_spacer <- pat_mgf$n_spacer
  n_hdr <- pat_mgf$n_hdr
  n_to_pepmass <- pat_mgf$n_to_pepmass
  n_to_title <- pat_mgf$n_to_title
  n_to_scan <- pat_mgf$n_to_scan
  n_to_rt <- pat_mgf$n_to_rt
  n_to_charge <- pat_mgf$n_to_charge
  sep_ms2s <- pat_mgf$sep_ms2s
  nfields_ms2s <- pat_mgf$nfields_ms2s
  sep_pepmass <- pat_mgf$sep_pepmass
  nfields_pepmass <- pat_mgf$nfields_pepmass
  raw_file <- pat_mgf$raw_file

  local({
    if (type_mgf == "msconv_thermo") {
      data_format <- "Thermo-RAW"
      mgf_format <- "MSconvert"
    }
    else if (type_mgf == "pd") {
      data_format <- "Thermo-RAW"
      mgf_format <- "Thermo-ProteomeDiscoverer"
    }
    else if (type_mgf == "msconv_pasef") {
      data_format <- "Bruker-D"
      mgf_format <- "MSconvert"
    }
    else if (type_mgf == "default_pasef") {
      data_format <- "Bruker-D"
      mgf_format <- "Bruker-DataAnalysis"
    }
    
    ans <- list(data_format = data_format, mgf_format = mgf_format)
    saveRDS(ans, file.path(filepath, "info_format.rds"))
  })
  
  rm(list = c("pat_mgf"))

  # chunks by mgf files
  out <- vector("list", length(filelist))

  for (i in seq_along(filelist)) {
    message("Loading '", filelist[i], "'.")

    # refresh at every `i`
    temp_dir <- local({
      dir <- find_dir(file.path(filepath, "temp"))

      if (!is.null(dir)) 
        fs::file_delete(dir)

      dir.create(file.path(filepath, "temp"), showWarnings = FALSE)
      find_dir(file.path(filepath, "temp"))
    })

    readr::read_lines_chunked(file = file.path(filepath, filelist[i]),
                              callback = SideEffectChunkCallback$new(f),
                              chunk_size = 1000000L)
    
    # for "default_pasef" format
    if (!is.null(raw_file)) {
      raw_file <- local({
        file <- file.path(filepath, "temp", "chunk_1.mgf")
        hdr <- readLines(file, 50L)
        pat <- "^COM="
        line_file <- hdr[grepl(pat, hdr)]
        gsub(pat, "", line_file)
      })
    }

    out[[i]] <- read_mgf_chunks(filepath = temp_dir,
                                topn_ms2ions = topn_ms2ions,
                                ms1_charge_range = ms1_charge_range, 
                                ms1_scan_range = ms1_scan_range, 
                                ret_range = ret_range,
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                index_ms2 = index_ms2,
                                type_mgf = type_mgf, 
                                n_bf_begin = n_bf_begin, 
                                n_spacer = n_spacer,
                                n_hdr = n_hdr,
                                n_to_pepmass = n_to_pepmass,
                                n_to_title = n_to_title,
                                n_to_scan = n_to_scan,
                                n_to_rt = n_to_rt,
                                n_to_charge = n_to_charge, 
                                sep_ms2s = sep_ms2s, 
                                nfields_ms2s = nfields_ms2s, 
                                sep_pepmass = sep_pepmass, 
                                nfields_pepmass = nfields_pepmass, 
                                raw_file = raw_file)

    local({
      dir2 <- file.path(filepath, gsub("\\.[^.]*$", "", filelist[i]))
      dir.create(dir2, showWarnings = FALSE)
      dir2 <- find_dir(dir2)

      if (fs::file_exists(dir2)) 
        fs::file_delete(dir2)

      fs::file_move(temp_dir, dir2)
    })
  }

  out <- out %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(ms1_mass) %>%
    dplyr::filter(ms1_mass >= min_mass, ms1_mass <= max_mass) %>% 
    dplyr::mutate(frame = find_ms1_interval(ms1_mass,
                                            from = min_mass,
                                            ppm = ppm_ms1)) 
  
  out <- local({
    raws_files <- out$raw_file
    raws <- raws_files[!duplicated.default(raws_files)]
    inds <- seq_along(raws)
    names(inds) <- raws
    saveRDS(inds, file.path(filepath, "raw_indexes.rds"))
    out$raw_file <- unname(inds[raws_files])
    
    scans <- out$scan_title
    inds <- seq_along(scans)
    names(inds) <- scans
    saveRDS(inds, file.path(filepath, "scan_indexes.rds"))
    out$scan_title <- unname(inds[scans])
    
    out
  })
  
  saveRDS(out, out_path)

  rm(list = c("out"))
  gc()

  invisible(NULL)
}


#' Reads mgfs in chunks.
#'
#' @param type_mgf The type of MGF format.
#' @param n_bf_begin The number of lines before \code{BEGIN IONS}. Zero for PD
#'   and MSConvert.
#' @param n_spacer The number of spacer lines between the preceding line END
#'   IONS and the following line BEGIN IONS. The value is 1 for Proteome
#'   Discoverer and 0 for MSConvert.
#' @param n_hdr The number of lines before MS2 data in an MGF. The value is +6
#'   for PD and +5 for MSConvert.
#' @param n_to_pepmass The number of lines from BEGIN to PEPMASS.
#' @param n_to_title The number of lines from BEGIN to TITLE. The value is the
#'   same between PD and MSConvert.
#' @param n_to_scan The number of lines from BEGIN to SCANS. The value is +5 for
#'   PD.
#' @param n_to_rt The number of lines from BEGIN to RTINSECONDS.
#' @param n_to_charge The number of lines from BEGIN to CHARGE.
#' @param sep_ms2s The separation character between MS2 m/z and intensity
#'   values.
#' @param nfields_ms2s The number of fields in MS2 entries. Mostely two and can
#'   be three for some Bruker MGFs.
#' @param sep_pepmass The separation character between MS1 m/z and intensity
#'   values.
#' @param nfields_pepmass The number of fields in \code{PEPMASS}.
#' @param raw_file The raw file name. Is NULL for PD and MSConvert.
#' @inheritParams readMGF
#' @inheritParams matchMS
read_mgf_chunks <- function (filepath = "~/proteoM/mgf",
                             topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                             ms1_scan_range = c(1L, .Machine$integer.max), 
                             ret_range = c(0, Inf), 
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                             type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                             n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L, 
                             n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L, 
                             n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                             sep_pepmass = " ", nfields_pepmass = 2L, 
                             raw_file = NULL) 
{
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (!length(filelist)) 
    stop("No mgf files under ", filepath, call. = FALSE)

  n_cores <- detect_cores(32L)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  
  parallel::clusterExport(
    cl, 
    c("stri_startswith_fixed", 
      "stri_endswith_fixed", 
      "stri_replace_first_fixed", 
      "stri_split_fixed"), 
    envir = environment(stringi::stri_startswith_fixed)
  )

  parallel::clusterExport(
    cl,
    c("proc_mgf_chunks", 
      "proc_mgfs", 
      "which_topx2", 
      "find_ms1_interval"), 
    envir = environment(proteoM:::proc_mgf_chunks)
  )

  out <- parallel::clusterApply(cl, file.path(filepath, filelist),
                                proc_mgf_chunks,
                                topn_ms2ions = topn_ms2ions,
                                ms1_charge_range = ms1_charge_range, 
                                ms1_scan_range = ms1_scan_range, 
                                ret_range = ret_range,
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                index_ms2 = index_ms2,
                                type_mgf = type_mgf,
                                n_bf_begin = n_bf_begin, 
                                n_spacer = n_spacer,
                                n_hdr = n_hdr,
                                n_to_pepmass = n_to_pepmass,
                                n_to_title = n_to_title,
                                n_to_scan = n_to_scan,
                                n_to_rt = n_to_rt,
                                n_to_charge = n_to_charge, 
                                sep_ms2s = sep_ms2s, 
                                nfields_ms2s = nfields_ms2s, 
                                sep_pepmass = sep_pepmass, 
                                nfields_pepmass = nfields_pepmass, 
                                raw_file = raw_file)
  
  parallel::stopCluster(cl)
  gc()

  out <- dplyr::bind_rows(out)

  # adds back broken mgf entries
  afs <- local({
    afs <- list.files(path = file.path(filepath), pattern = "^.*\\_af.mgf$")

    idxes <- afs %>%
      gsub("^chunk_(\\d+)_af\\.mgf", "\\1", .) %>%
      as.integer() %>%
      sort()

    paste0("chunk_", idxes, "_af.mgf") %>%
      .[-length(.)]
  })

  bfs <- local({
    bfs <- list.files(path = file.path(filepath), pattern = "^.*\\_bf.mgf$")

    idxes <- bfs %>%
      gsub("^chunk_(\\d+)_bf\\.mgf", "\\1", .) %>%
      as.integer() %>%
      sort()

    paste0("chunk_", idxes, "_bf.mgf") %>%
      .[-1]
  })

  # stopifnot(length(afs) == length(bfs))

  gaps <- purrr::map2(afs, bfs, function (x, y) {
    af <- stringi::stri_read_lines(file.path(filepath, x))
    bf <- stringi::stri_read_lines(file.path(filepath, y))
    append(af, bf)
  }) %>%
    unlist(use.names = FALSE) %T>%
    write(file.path(filepath, "gaps.mgf"))

  local({
    nms <- list.files(path = file.path(filepath), pattern = "^.*\\_[ab]f.mgf$")

    if (length(nms)) 
      suppressMessages(file.remove(file.path(filepath, nms)))
  })

  if (!is.null(gaps)) {
    out <- dplyr::bind_rows(
      out,
      proc_mgfs(gaps,
                topn_ms2ions = topn_ms2ions,
                ms1_charge_range = ms1_charge_range, 
                ms1_scan_range = ms1_scan_range, 
                ret_range = ret_range,
                ppm_ms2 = ppm_ms2,
                min_ms2mass = min_ms2mass,
                index_ms2 = index_ms2,
                type_mgf = type_mgf, 
                n_bf_begin = n_bf_begin, 
                n_spacer = n_spacer,
                n_hdr = n_hdr,
                n_to_pepmass = n_to_pepmass,
                n_to_title = n_to_title,
                n_to_scan = n_to_scan,
                n_to_rt = n_to_rt,
                n_to_charge = n_to_charge, 
                sep_ms2s = sep_ms2s, 
                nfields_ms2s = nfields_ms2s, 
                sep_pepmass = sep_pepmass, 
                nfields_pepmass = nfields_pepmass, 
                raw_file = raw_file)
    )
  }

  if (type_mgf == "default_pasef") {
    out <- out %>% 
      dplyr::mutate(scan_id = as.character(scan_num), 
                    scan_num = row_number())
  } 

  invisible(out)
}


#' Processes MGF entries in chunks.
#'
#' @param file A chunk of MGF (chunk_1.mgf etc.) with a prepending file path.
#' @inheritParams readMGF
#' @inheritParams read_mgf_chunks
proc_mgf_chunks <- function (file, topn_ms2ions = 100L, 
                             ms1_charge_range = c(2L, 6L), 
                             ms1_scan_range = c(1L, .Machine$integer.max), 
                             ret_range = c(0, Inf), 
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                             type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                             n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L, 
                             n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L, 
                             n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                             sep_pepmass = " ", nfields_pepmass = 2L, 
                             raw_file = NULL) 
{
  message("Parsing '", file, "'.")
  lines <- stringi::stri_read_lines(file)

  basename <- gsub("\\.[^.]*$", "", file)

  begins <- which(stringi::stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(lines, "END IONS"))

  af <- local({
    le <- ends[length(ends)]
    lb <- begins[length(begins)]

    af <- if (lb > le) 
      lines[(le + n_spacer + 1L):length(lines)]
    else 
      NULL

    write(af, file.path(paste0(basename, "_af.mgf")))

    af
  })

  bf <- local({
    le <- ends[1]
    lb <- begins[1]

    bf <- if (lb > le) 
      lines[1:(le + n_spacer)]
    else 
      NULL

    write(bf, file.path(paste0(basename, "_bf.mgf")))

    bf
  })

  if (!is.null(af)) 
    lines <- lines[1:(begins[length(begins)] - n_bf_begin- 1L)]

  if (!is.null(bf)) 
    lines <- lines[-c(1:(ends[1] + n_spacer))]

  out <- proc_mgfs(lines = lines,
                   topn_ms2ions = topn_ms2ions,
                   ms1_charge_range = ms1_charge_range, 
                   ms1_scan_range = ms1_scan_range, 
                   ret_range = ret_range,
                   ppm_ms2 = ppm_ms2,
                   min_ms2mass = min_ms2mass,
                   index_ms2 = index_ms2,
                   type_mgf = type_mgf, 
                   n_bf_begin = n_bf_begin,
                   n_spacer = n_spacer,
                   n_hdr = n_hdr,
                   n_to_pepmass = n_to_pepmass,
                   n_to_title = n_to_title,
                   n_to_scan = n_to_scan,
                   n_to_rt = n_to_rt,
                   n_to_charge = n_to_charge, 
                   sep_ms2s = sep_ms2s, 
                   nfields_ms2s = nfields_ms2s, 
                   sep_pepmass = sep_pepmass, 
                   nfields_pepmass = nfields_pepmass, 
                   raw_file = raw_file)
}


#' Helper in processing MGF entries in chunks.
#'
#' @inheritParams proc_mgf_chunks
proc_mgfs <- function (lines, topn_ms2ions = 100L, 
                       ms1_charge_range = c(2L, 6L), 
                       ms1_scan_range = c(1L, .Machine$integer.max), 
                       ret_range = c(0, Inf), 
                       ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                       type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                       n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L,
                       n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L,
                       n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                       sep_pepmass = " ", nfields_pepmass = 2L, raw_file = NULL) 
{
  options(digits = 9L)

  # MS2 ions
  begins <- which(stringi::stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(lines, "END IONS"))

  # (-1L: one line above "END IONS")
  ms2s <- mapply(function (x, y) lines[(x + n_hdr) : (y - 1L)], 
                 begins, ends, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ms2s <- lapply(ms2s, stringi::stri_split_fixed, pattern = sep_ms2s, 
                 n = nfields_ms2s, simplify = TRUE)

  ms2_moverzs <- lapply(ms2s, function (x) as.numeric(x[, 1]))

  # not to round the `ms2_ints` values here (to avoid ties and thus
  # shorter length after `which_topx` of a vector of integers;
  # more likely to have ties with integers than doubles)
  #
  # the latest `which_topx` handles ties and the length of the output
  #  equals `topn_ms2ions` even with ties;
  # except when the length of a vector is shorter than `topn_ms2ions`)

  ms2_ints <- lapply(ms2s, function (x) as.numeric(x[, 2]))

  lens <- lapply(ms2_moverzs, length)
  lens <- .Internal(unlist(lens, recursive = FALSE, use.names = FALSE))

  if (topn_ms2ions < Inf) {
    rows <- lapply(ms2_ints, which_topx2, topn_ms2ions)

    # OK to round `ms2_ints` now
    ms2_ints <- mapply(function (x, y) round(x[y], digits = 0L), 
                       ms2_ints, rows, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    ms2_moverzs <- mapply(function (x, y) round(x[y], digits = 5L), 
                          ms2_moverzs, rows, 
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    rm(list = c("rows", "ms2s"))
  }

  # MS1 ions
  ms1s <- stringi::stri_replace_first_fixed(lines[begins + n_to_pepmass], 
                                            "PEPMASS=", "")
  ms1s <- lapply(ms1s, stringi::stri_split_fixed, pattern = sep_pepmass, 
                 n = nfields_pepmass, simplify = TRUE)

  ms1_moverzs <- lapply(ms1s, function (x) round(as.numeric(x[, 1]), digits = 5L))
  ms1_moverzs <- .Internal(unlist(ms1_moverzs, recursive = FALSE, use.names = FALSE))

  # NA's if no MS1 intensities available
  ms1_ints <- lapply(ms1s, function (x) round(as.numeric(x[, 2]), digits = 0L))
  ms1_ints <- .Internal(unlist(ms1_ints, recursive = FALSE, use.names = FALSE))

  rm(list = c("ms1s"))
  gc()

  # Others
  scan_titles <- 
    stringi::stri_replace_first_fixed(lines[begins + n_to_title], "TITLE=", "")
  
  if (type_mgf %in% c("msconv_thermo", "msconv_pasef")) {
    raw_files <- 
      stringi::stri_replace_first_regex(scan_titles, "^.* File:\"([^\"]+)\".*", "$1")
    scan_nums <- as.integer(
      stringi::stri_replace_first_regex(scan_titles, 
                                        "^.*\\.(\\d+)\\.\\d+\\.\\d+ File:\".*", 
                                        "$1"))
  } 
  else if (type_mgf == "pd") {
    raw_files <- gsub("^.*File: \"([^\"]+)\".*", "\\1", scan_titles)
    raw_files <- gsub("\\\\", "/", raw_files)
    raw_files <- gsub("^.*/(.*)", "\\1", raw_files)
    scan_nums <- as.integer(gsub("^.* scans: \"([0-9]+)\"$", "\\1", scan_titles))
  } 
  else if (type_mgf == "default_pasef") {
    raw_files <- rep(raw_file, length(lens))
    scan_nums <- 
      stringi::stri_replace_first_fixed(lines[begins + n_to_scan], "RAWSCANS=", "")
  } 
  else {
    stop("Unknown MGF format.", call. = FALSE)
  }

  ret_times <- 
    stringi::stri_replace_first_fixed(lines[begins + n_to_rt], "RTINSECONDS=", "")
  ret_times <- as.numeric(ret_times)
  
  ms1_charges <- 
    stringi::stri_replace_first_fixed(lines[begins + n_to_charge], "CHARGE=", "")

  # MS1 neutral masses
  proton <- 1.00727647

  charges <- lapply(ms1_charges, stringi::stri_reverse)
  charges <- .Internal(unlist(charges, recursive = FALSE, use.names = FALSE))
  charges <- as.integer(charges)

  ms1_masses <- mapply(function (x, y) x * y - y * proton, 
                       ms1_moverzs, charges, 
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)
  ms1_masses <- round(ms1_masses, digits = 5L)

  # Subsetting
  # (no `scan_num` subsetting for PASEF)
  if (type_mgf == "default_pasef") 
    rows <- (charges >= ms1_charge_range[1] & charges <= ms1_charge_range[2] & 
               ret_times >= ret_range[1] & ret_times <= ret_range[2])
  else 
    rows <- (charges >= ms1_charge_range[1] & charges <= ms1_charge_range[2] & 
               scan_nums >= ms1_scan_range[1] & scan_nums <= ms1_scan_range[2] &
               ret_times >= ret_range[1] & ret_times <= ret_range[2])

  scan_titles <- scan_titles[rows]
  raw_files <- raw_files[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_masses <- ms1_masses[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- ret_times[rows]
  scan_nums <- scan_nums[rows]
  ms2_moverzs <- ms2_moverzs[rows]
  ms2_ints <- ms2_ints[rows]
  lens <- lens[rows]

  if (index_ms2) 
    ms2_moverzs <- lapply(ms2_moverzs,
                          find_ms1_interval,
                          from = min_ms2mass,
                          ppm = ppm_ms2)

  out <- tibble::tibble(scan_title = scan_titles,
                        raw_file = raw_files,
                        ms1_moverz = ms1_moverzs,
                        ms1_mass = ms1_masses,
                        ms1_int = ms1_ints,
                        ms1_charge = ms1_charges,
                        ret_time = ret_times,
                        scan_num = scan_nums,
                        ms2_moverz = ms2_moverzs,
                        ms2_int = ms2_ints,
                        ms2_n = lens)
}


#' Calculates the frame numbers for a list of experimental MS1 mass by
#' intervals.
#'
#' Needs correct \code{from}.
#' @param mass Numeric; a list of MS1 masses.
#' @inheritParams find_ms1_cutpoints
#' @examples
#' \donttest{
#' find_ms1_interval(c(500, 800.1))
#' }
#' @return Frame numbers.
#' @seealso find_ms1_cutpoints
find_ms1_interval <- function (mass = 1800.0, from = 350L, ppm = 20L) 
{
  d <- ppm/1e6
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+d))
}


#' Finds the type of MGF.
#'
#' @param file The path to an MGF file.
find_mgf_type <- function (file) 
{
  hdr <- readLines(file, 5000L)
  begins <- which(stringi::stri_startswith_fixed(hdr, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(hdr, "END IONS"))
  
  if (!length(begins)) 
    stop("Check corrupted files: the tag of `BEGIN IONS` not found in MGF.")
  
  # if (!length(ends))
  #   stop("The tag of `END IONS` not found in MGF.", call. = FALSE)

  
  ## MSConvert (Thermo)
  # <RunId>.<ScanNumber><ScanNumber><ChargeState> File:"<SourcePath>", NativID:"<Id>"
  # 
  # BEGIN IONS
  # TITLE=rawname.179.179.3 File:"rawname.raw", NativeID:"controllerType=0 controllerNumber=1 scan=179"
  # RTINSECONDS=63.4689
  # PEPMASS=482.224129434954 280125.927246099978
  # CHARGE=3+
  
  ## MSconvert (timsTOF)
  # BEGIN IONS
  # "TITLE=rawname.2.2.1 File:\"rawname.d\", 
  #   NativeID:\"merged=1 frame=2 scanStart=200 scanEnd=224\", IonMobility:\"1.3990913467070001\""
  # RTINSECONDS=2.950849933
  # PEPMASS=1221.991350361787
  # CHARGE=1+

  ## Proteome Discoverer
  # MASS=Monoisotopic
  # BEGIN IONS
  # TITLE=File: "Z:\Folder\rawname.raw"; SpectrumID: "1"; scans: "179"
  # PEPMASS=482.22421 110739.89844
  # CHARGE=3+
  # RTINSECONDS=63
  # SCANS=179
  
  ## RawConverter
  # BEGIN IONS
  # TITLE=Z:\Folder\rawname.raw
  # SCANS=179
  # RTINSECONDS=63.4689
  # CHARGE=3+
  # PEPMASS=482.2242
  
  # ###FS:    #m/z: 454.22925 #charge 2+
  # ###MS: 1
  # ###MSMS: 7, 9-11
  # ###Mobility: 0.7810
  # BEGIN IONS
  # TITLE=Cmpd 18, +MS2(454.2292), 27.7eV, 0.0min, 1/K0=0.781, #7-11
  # RTINSECONDS=1.04678
  # RAWSCANS=1,7,9-11
  # PEPMASS=454.22925	70295
  # CHARGE=2+
  #   221.10530	104	
  #   ...
  # END IONS
  # 
  # ###FS:    #m/z: 665.07643 #charge 1+
  # ###MS: 8
  # ###MSMS: 23
  # ###Mobility: 1.0572
  # BEGIN IONS
  # TITLE=Cmpd 49, +MS2(665.0764), 39.6eV, 0.1min, 1/K0=1.057, #23

  type <- local({
    ln_tit <- hdr[grepl("TITLE", hdr)][1]

    file_msconvert_pasef <- file_msconvert_thermo <- "File:\""
    file_pd <- "File: \""
    file_default_pasef <- "Cmpd "

    scan_msconv_thermo <- "scan=\\d+"
    scan_msconv_pasef <- "scanStart=\\d+"
    scan_pd <- "scans: \"\\d+\""
    scan_default_pasef <- NULL

    if (grepl(file_msconvert_thermo, ln_tit) && grepl(scan_msconv_thermo, ln_tit)) 
      "msconv_thermo"
    else if (grepl(file_msconvert_pasef, ln_tit) && grepl(scan_msconv_pasef, ln_tit))
      "msconv_pasef"
    else if (grepl(file_pd, ln_tit) && grepl(scan_pd, ln_tit)) 
      "pd"
    else if (grepl(file_default_pasef, ln_tit) && is.null(scan_default_pasef))
      "default_pasef"
    else 
      stop("Unkown format of MGFs.")
  })
  
  # n_bf_begin: the number of lines before `BEGIN IONS`
  # n_spacer: the number of white-space lines between two adjacent blocks
  # n_hdr: the number of hear lines from (including) `BEGIN IONS`
  # raw_file: if ".", an indicator to find RAW file names from the header
  
  if (type == "msconv_thermo") {
    n_bf_begin <- 0L
    n_spacer <- 0L
    n_hdr <- 5L
    n_to_pepmass <- 3L
    n_to_title <- 1L
    n_to_scan <- 0L
    n_to_rt <- 2L
    n_to_charge <- 4L
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " "
    nfields_pepmass <- 2L
    raw_file <- NULL
  } 
  else if (type == "msconv_pasef") {
    n_bf_begin <- 0L
    n_spacer <- 0L
    n_hdr <- 5L
    n_to_pepmass <- 3L
    n_to_title <- 1L
    n_to_scan <- 0L
    n_to_rt <- 2L
    n_to_charge <- 4L
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " " # as of 2021-12-27, missing MS1 intensity
    nfields_pepmass <- 2L # "NA" MS1 intensities after parsing
    raw_file <- NULL
  } 
  else if (type == "pd") {
    n_bf_begin <- 0L
    n_spacer <- 1L
    n_hdr <- 6L
    n_to_pepmass <- 2L
    n_to_title <- 1L
    n_to_scan <- 5L
    n_to_rt <- 4L
    n_to_charge <- 3L
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " "
    nfields_pepmass <- 2L
    raw_file <- NULL
  } 
  else if (type == "default_pasef") {
    n_bf_begin <- 4L
    n_spacer <- 1L
    n_hdr <- 6L
    n_to_pepmass <- 4L
    n_to_title <- 1L
    n_to_scan <- 3L
    n_to_rt <- 2L
    n_to_charge <- 5L
    
    sep_ms2s <- "\t"
    nfields_ms2s <- 3L
    sep_pepmass <- "\t"
    nfields_pepmass <- 2L
    raw_file <- "."
  }
  else {
    stop("Unkown format of MGFs.")
  }

  invisible(list(type = type, 
                 n_bf_begin = n_bf_begin, 
                 n_spacer = n_spacer,
                 n_hdr = n_hdr,
                 n_to_pepmass = n_to_pepmass,
                 n_to_title = n_to_title,
                 n_to_scan = n_to_scan,
                 n_to_rt = n_to_rt,
                 n_to_charge = n_to_charge, 
                 sep_ms2s = sep_ms2s, 
                 nfields_ms2s = nfields_ms2s, 
                 sep_pepmass = sep_pepmass, 
                 nfields_pepmass = nfields_pepmass, 
                 raw_file = raw_file))
}

