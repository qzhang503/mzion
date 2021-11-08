#' Helper in loading MGFs.
#'
#' Currently, \code{index_ms2} is always FALSE. Argument \code{ppm_ms2} is only
#' effective at \code{index_ms2 = TRUE}.
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
                       ppm_ms1 = 20L, ppm_ms2 = 25L, index_ms2 = FALSE, 
                       is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE) {
  
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    },
    add = TRUE
  )

  # ---
  fun <- as.character(match.call()[[1]])
  
  args_except <- NULL
  
  cache_pars <- find_callarg_vals(time = NULL, 
                                  path = file.path(out_path, "Calls"), 
                                  fun = paste0(fun, ".rda"), 
                                  args = names(formals(fun)) %>% 
                                    .[! . %in% args_except]) %>% 
    .[sort(names(.))]
  
  call_pars <- mget(names(formals()) %>% .[! . %in% args_except], 
                    envir = rlang::current_env(), 
                    inherits = FALSE) %>% 
    .[sort(names(.))]
  
  ok_pars <- identical(call_pars, cache_pars)

  rds <- file.path(mgf_path, "mgf_queries.rds")
  
  if (ok_pars && file.exists(rds)) {
    message("Found cached MGFs: `", rds, "`.")
    .savecall <- FALSE
  } else {
    message("Processing raw MGFs.")
    
    if (is_ms1_three_frame) {
      ppm_ms1_new <- as.integer(ceiling(ppm_ms1 * .5))
    } else {
      ppm_ms1_new <- ppm_ms1
    }
    
    if (is_ms2_three_frame) {
      ppm_ms2_new <- as.integer(ceiling(ppm_ms2 * .5))
    } else {
      ppm_ms2_new <- ppm_ms2
    }
    
    # dirs <- list.dirs(mgf_path, full.names = FALSE, recursive = TRUE) %>%
    #   .[! . == ""]
    # unlink(file.path(mgf_path, dirs), recursive = TRUE)

    delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                       "\\.xls$", "\\.csv$", "\\.txt$", 
                                       "^mgf$", "^mgfs$", "Calls"))

    readMGF(filepath = mgf_path,
            min_mass = min_mass,
            max_mass = max_mass, 
            min_ms2mass = min_ms2mass,
            topn_ms2ions = topn_ms2ions,
            ret_range = c(0, Inf),
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
readMGF <- function (filepath = "~/proteoM/mgfs",
                     min_mass = 500L, max_mass = 6000L, min_ms2mass = 110L, 
                     topn_ms2ions = 100L, ret_range = c(0, Inf), 
                     ppm_ms1 = 20L, ppm_ms2 = 25L, index_ms2 = FALSE,
                     out_path = file.path(filepath, "mgf_queries.rds")) {

  f <- function(x, pos) {
    nm <- file.path(filepath, "temp", paste0("chunk", "_", pos, ".mgf"))
    writeLines(x, nm)
  }

  # parsing rules
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (!length(filelist)) {
    stop("No '.mgf' files under ", filepath, call. = FALSE)
  }

  pat_mgf <- find_mgf_type(file.path(filepath, filelist[[1]]))

  pat_file = pat_mgf$pat_file
  pat_scan = pat_mgf$pat_scan
  n_spacer = pat_mgf$n_spacer
  n_hdr = pat_mgf$n_hdr
  n_to_pepmass = pat_mgf$n_to_pepmass
  n_to_title = pat_mgf$n_to_title
  n_to_scan = pat_mgf$n_to_scan
  n_to_rt = pat_mgf$n_to_rt
  n_to_charge = pat_mgf$n_to_charge

  rm(list = c("pat_mgf"))

  # chunks by mgf files
  out <- vector("list", length(filelist))

  for (i in seq_along(filelist)) {
    message("Loading '", filelist[i], "'.")

    # refresh at every `i`
    temp_dir <- local({
      temp_dir <- find_dir(file.path(filepath, "temp"))

      if (!is.null(temp_dir)) {
        fs::file_delete(temp_dir)
      }

      dir.create(file.path(filepath, "temp"), showWarnings = FALSE)
      temp_dir <- find_dir(file.path(filepath, "temp"))
    })

    readr::read_lines_chunked(file = file.path(filepath, filelist[i]),
                              callback = SideEffectChunkCallback$new(f),
                              chunk_size = 1000000L)

    out[[i]] <- read_mgf_chunks(filepath = temp_dir,
                                topn_ms2ions = topn_ms2ions,
                                ret_range = ret_range,
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                index_ms2 = index_ms2,
                                pat_file = pat_file,
                                pat_scan = pat_scan,
                                n_spacer = n_spacer,
                                n_hdr = n_hdr,
                                n_to_pepmass = n_to_pepmass,
                                n_to_title = n_to_title,
                                n_to_scan = n_to_scan,
                                n_to_rt = n_to_rt,
                                n_to_charge = n_to_charge)

    local({
      temp_dir2 <- file.path(filepath, gsub("\\.[^.]*$", "", filelist[i]))
      dir.create(temp_dir2, showWarnings = FALSE)
      temp_dir2 <- find_dir(temp_dir2)

      if (fs::file_exists(temp_dir2)) {
        fs::file_delete(temp_dir2)
      }

      fs::file_move(temp_dir, temp_dir2)
    })
  }

  out <- out %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(ms1_mass) %>%
    dplyr::filter(ms1_mass >= min_mass, ms1_mass <= max_mass) %>% 
    dplyr::mutate(frame = find_ms1_interval(ms1_mass,
                                            from = min_mass,
                                            ppm = ppm_ms1)) %T>%
    saveRDS(out_path)

  rm(list = c("out"))
  gc()

  invisible(NULL)
}


#' Reads mgfs in chunks.
#'
#' @param pat_file The pattern for parsing the \code{File} field in MGF.
#' @param pat_scan  The pattern for parsing the \code{scan} field in MGF.
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
#' @inheritParams readMGF
#' @inheritParams matchMS
read_mgf_chunks <- function (filepath = "~/proteoM/mgfs",
                             topn_ms2ions = 100L, ret_range = c(0L, Inf),
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                             pat_file, pat_scan, n_spacer, n_hdr, n_to_pepmass,
                             n_to_title, n_to_scan, n_to_rt, n_to_charge) {

  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (!length(filelist)) {
    stop("No mgf files under ", filepath, call. = FALSE)
  }

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
    c("proc_mgf_chunks_i", 
      "proc_mgf_chunks", 
      "proc_mgfs", 
      "which_topx2", 
      "find_ms1_interval"), 
    envir = environment(proteoM:::proc_mgf_chunks)
  )
  
  out <- parallel::clusterApply(cl, file.path(filepath, filelist),
                                proc_mgf_chunks_i,
                                topn_ms2ions = topn_ms2ions,
                                ret_range = ret_range,
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                index_ms2 = index_ms2,
                                pat_file = pat_file,
                                pat_scan = pat_scan,
                                n_spacer = n_spacer,
                                n_hdr = n_hdr,
                                n_to_pepmass = n_to_pepmass,
                                n_to_title = n_to_title,
                                n_to_scan = n_to_scan,
                                n_to_rt = n_to_rt,
                                n_to_charge = n_to_charge)

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

  stopifnot(length(afs) == length(bfs))

  gaps <- purrr::map2(afs, bfs, function (x, y) {
    af <- stri_read_lines(file.path(filepath, x))
    bf <- stri_read_lines(file.path(filepath, y))
    append(af, bf)
  }) %>%
    unlist(use.names = FALSE) %T>%
    write(file.path(filepath, "gaps.mgf"))

  local({
    nms <- list.files(path = file.path(filepath), pattern = "^.*\\_[ab]f.mgf$")

    if (length(nms)) {
      suppressMessages(file.remove(file.path(filepath, nms)))
    }
  })

  if (!is.null(gaps)) {
    out <- dplyr::bind_rows(
      out,
      proc_mgfs(gaps,
                topn_ms2ions = topn_ms2ions,
                ret_range = ret_range,
                ppm_ms2 = ppm_ms2,
                min_ms2mass = min_ms2mass,
                index_ms2 = index_ms2,
                pat_file = pat_file,
                pat_scan = pat_scan,
                n_spacer = n_spacer,
                n_hdr = n_hdr,
                n_to_pepmass = n_to_pepmass,
                n_to_title = n_to_title,
                n_to_scan = n_to_scan,
                n_to_rt = n_to_rt,
                n_to_charge = n_to_charge)
    )
  }

  invisible(out)
}


#' Helper of \link{proc_mgf_chunks}.
#'
#' @param file A mgf chunk (chunk_1.mgf etc.) with prepending file path.
#' @inheritParams read_mgf_chunks
proc_mgf_chunks_i <- function (file, topn_ms2ions = 100L, ret_range = c(0L, Inf),
                               ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                               pat_file, pat_scan, n_spacer, n_hdr, n_to_pepmass,
                               n_to_title, n_to_scan, n_to_rt, n_to_charge) {

  message("Parsing '", file, "'.")

  x <- proc_mgf_chunks(lines = stringi::stri_read_lines(file), 
                       topn_ms2ions = topn_ms2ions,
                       ret_range = ret_range,
                       ppm_ms2 = ppm_ms2,
                       min_ms2mass = min_ms2mass,
                       index_ms2 = index_ms2,
                       filepath = file,
                       pat_file = pat_file,
                       pat_scan = pat_scan,
                       n_spacer = n_spacer,
                       n_hdr = n_hdr,
                       n_to_pepmass = n_to_pepmass,
                       n_to_title = n_to_title,
                       n_to_scan = n_to_scan,
                       n_to_rt = n_to_rt,
                       n_to_charge = n_to_charge)
}


#' Processes MGF entries in chunks.
#'
#' @param lines MGF lines.
#' @inheritParams readMGF
#' @inheritParams read_mgf_chunks
proc_mgf_chunks <- function (lines, topn_ms2ions = 100L, ret_range = c(0L, Inf),
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                             filepath = file.path("~/proteoM/mgfs/temp"),
                             pat_file, pat_scan, n_spacer, n_hdr, n_to_pepmass,
                             n_to_title, n_to_scan, n_to_rt, n_to_charge) {

  basename <- gsub("\\.[^.]*$", "", filepath)

  begins <- which(stringi::stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(lines, "END IONS"))

  af <- local({
    le <- ends[length(ends)]
    lb <- begins[length(begins)]

    if (lb > le) {
      af <- lines[(le + n_spacer + 1L):length(lines)]
    } else {
      af <- NULL
    }

    write(af, file.path(paste0(basename, "_af.mgf")))

    af
  })

  bf <- local({
    le <- ends[1]
    lb <- begins[1]

    if (lb > le) {
      bf <- lines[1:(le + n_spacer)]
    } else {
      bf <- NULL
    }

    write(bf, file.path(paste0(basename, "_bf.mgf")))

    bf
  })

  if (!is.null(af)) {
    lines <- lines[1:(begins[length(begins)] - 1L)]
  }

  if (!is.null(bf)) {
    lines <- lines[-c(1:(ends[1] + n_spacer))]
  }

  out <- proc_mgfs(lines,
                   topn_ms2ions = topn_ms2ions,
                   ret_range = ret_range,
                   ppm_ms2 = ppm_ms2,
                   min_ms2mass = min_ms2mass,
                   index_ms2 = index_ms2,
                   pat_file = pat_file,
                   pat_scan = pat_scan,
                   n_spacer = n_spacer,
                   n_hdr = n_hdr,
                   n_to_pepmass = n_to_pepmass,
                   n_to_title = n_to_title,
                   n_to_scan = n_to_scan,
                   n_to_rt = n_to_rt,
                   n_to_charge = n_to_charge)
}


#' Helper in processing MGF entries in chunks.
#'
#' @inheritParams proc_mgf_chunks
proc_mgfs <- function (lines, topn_ms2ions = 100L, ret_range = c(0L, Inf),
                       ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE,
                       pat_file, pat_scan,
                       n_spacer = 1L, n_hdr = 6L, n_to_pepmass = 2L,
                       n_to_title = 1L, n_to_scan = 5L, n_to_rt = 4L,
                       n_to_charge = 3L) {

  options(digits = 9L)

  # MS2 ions
  begins <- which(stringi::stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(lines, "END IONS"))

  # (-1L: one line above "END IONS")
  ms2s <- mapply(function (x, y) lines[(x + n_hdr) : (y - 1L)], 
                 begins, ends, 
                 SIMPLIFY = FALSE, 
                 USE.NAMES = FALSE)
  ms2s <- lapply(ms2s, stringi::stri_split_fixed, " ", n = 2, simplify = TRUE)

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
  ms1s <- lapply(ms1s, stringi::stri_split_fixed, " ", n = 2, simplify = TRUE)

  ms1_moverzs <- lapply(ms1s, function (x) round(as.numeric(x[, 1]), digits = 5L))
  ms1_moverzs <- .Internal(unlist(ms1_moverzs, recursive = FALSE, use.names = FALSE))

  ms1_ints <- lapply(ms1s, function (x) round(as.numeric(x[, 2]), digits = 0L))
  ms1_ints <- .Internal(unlist(ms1_ints, recursive = FALSE, use.names = FALSE))

  rm(list = c("ms1s"))
  gc()

  # Others
  scan_titles <- stringi::stri_replace_first_fixed(lines[begins + n_to_title], 
                                                   "TITLE=", "")

  if (pat_file == "^.* File:\"([^\"]+)\".*") {
    raw_files <- gsub(pat_file, "\\1", scan_titles)
    raw_files <- gsub("\\\\", "/", raw_files)
  } else if (pat_file == "^.*File: \"([^\"]+)\".*") {
    raw_files <- gsub(pat_file, "\\1", scan_titles)
    raw_files <- gsub("\\\\", "/", raw_files)
    raw_files <- gsub("^.*/(.*)", "\\1", raw_files)
  } else {
    stop("Unknown MGF format.", call. = FALSE)
  }

  scan_nums <- as.integer(gsub(pat_scan, "\\1", scan_titles))

  ret_times <- stringi::stri_replace_first_fixed(lines[begins + n_to_rt], "RTINSECONDS=", "")
  ret_times <- as.numeric(ret_times)
  
  ms1_charges <- stringi::stri_replace_first_fixed(lines[begins + n_to_charge], "CHARGE=", "")

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
  rows <- (ret_times >= ret_range[1] & ret_times <= ret_range[2])

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

  if (index_ms2) {
    ms2_moverzs <- lapply(ms2_moverzs,
                          find_ms1_interval,
                          from = min_ms2mass,
                          ppm = ppm_ms2)
  }

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
find_ms1_interval <- function (mass = 1800.0, from = 350L, ppm = 20L) {
  d <- ppm/1e6
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+d))
}


#' Add the MS file name to timsTOF mgf files.
#'
#' @param path The file path to mgf files.
#' @param n Integer; the first \code{n} lines for quick checks of the MS file
#'   name.
#' @export
proc_mgf_timstof <- function (path, n = 1000L) {
  filelist <- list.files(path = file.path(path),
                         pattern = "^.*\\.mgf$")

  if (!length(filelist)) {
    stop("No mgf files under ", path,
         call. = FALSE)
  }

  purrr::walk(filelist, ~ {
    lines <- readLines(file.path(path, .x))

    hdr <- lines[1:n]

    fn <- hdr %>%
      .[grepl("^COM\\=.*\\.d$", .)] %>%
      gsub("^COM\\=(.*)$", "\\1", .) %>%
      paste0("File:~", ., "~,")

    stopifnot(length(fn) == 1L)

    lntit <- hdr %>%
      .[grepl("^TITLE", .)] %>%
      `[`(1)

    if (grepl("File:~", lntit)) {
      warning("MS file names already in ", .x,
              call. = FALSE)
    } else {
      rows <- grepl("^TITLE", lines)

      lines[rows] <- lines[rows]  %>%
        gsub("(Cmpd \\d+,)", paste("\\1", fn, sep =" "), .)

      writeLines(lines, file.path(path, .x))
    }
  })
}


#' Finds the type of MGF.
#'
#' @param file The path to an MGF file.
find_mgf_type <- function (file) {

  hdr <- readLines(file, 1000L)

  ln_tit <- hdr %>%
    .[grepl("^TITLE", .)] %>%
    `[`(1)

  file_msconvert <- "File:\""
  file_pd <- "File: \""

  scan_msonvert <- "scan=\\d+"
  scan_pd <- "scans: \"\\d+\""

  if (grepl(file_msconvert, ln_tit) && grepl(scan_msonvert, ln_tit)) {
    type <- "msconvert"
  } else if (grepl(file_pd, ln_tit) && grepl(scan_pd, ln_tit)) {
    type <- "pd"
  }

  if (type == "msconvert") {
    pat_file2 <- "^.* File:\"([^\"]+)\".*"
    pat_scan2 <- "^.* scan=([0-9]+)\"$"

    n_spacer = 0L
    n_hdr = 5L
    n_to_pepmass = 3L
    n_to_title = 1L
    n_to_scan = 0L
    n_to_rt = 2L
    n_to_charge = 4L
  } else if (type == "pd") {
    pat_file2 <- "^.*File: \"([^\"]+)\".*"
    pat_scan2 <- "^.* scans: \"([0-9]+)\"$"

    n_spacer = 1L
    n_hdr = 6L
    n_to_pepmass = 2L
    n_to_title = 1L
    n_to_scan = 5L
    n_to_rt = 4L
    n_to_charge = 3L
  } else {
    stop("The `File` lines in MGF needs to be in the format of ",
         "File:\"my.raw\" or File: \"my.raw\".\n",
         "The `scan` lines in MGF needs to be in the format of ",
         "scan=123 or scans: \"123\".")
  }

  invisible(list(pat_file = pat_file2,
                 pat_scan = pat_scan2,
                 n_spacer = n_spacer,
                 n_hdr = n_hdr,
                 n_to_pepmass = n_to_pepmass,
                 n_to_title = n_to_title,
                 n_to_scan = n_to_scan,
                 n_to_rt = n_to_rt,
                 n_to_charge = n_to_charge))
}


