#' Helper in loading MGFs.
#'
#' @param min_mass A minimum mass of precursors for considerations.
#' @param max_mass A maximum mass of precursors for considerations.
#' @param min_ms2mass A minimum m/z of MS2 ions for considerations.
#' @param is_ms1_three_frame Logical; is the searches by the three frames of
#'   preceding, current and following.
#' @param is_ms2_three_frame Logical; is the searches by the three frames of
#'   preceding, current and following.
#' @param mgf_cutmzs Cut points of MS1 m-over-z values in peak picking.
#' @param mgf_cutpercs The counts of MS2 features in each region of
#'   \code{mgf_cutmzs}.
#' @inheritParams matchMS
load_mgfs <- function (out_path, mgf_path, min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, topn_ms2ions = 100L, 
                       min_ms1_charge = 2L, max_ms1_charge = 6L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1 = 20L, ppm_ms2 = 20L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       enzyme = "trypsin_p", quant = "none", digits = 4L) 
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
  this_call <- match.call()
  fun <- as.character(this_call[[1]])
  fun <- fun[length(fun)] # may be called as mzion:::load_mgfs
  fun_env <- environment()
  
  args_except <- c("out_path")
  args <- names(formals(fun))
  args_must <- args[! args %in% args_except]
  
  cache_pars <- find_callarg_vals(
    time = NULL, 
    path = file.path(out_path, "Calls"), 
    fun = paste0(fun, ".rda"), 
    args = args_must, 
    new_args = unlist(formals(load_mgfs)[c("enzyme")])
  ) 
  
  cache_pars <- cache_pars[sort(names(cache_pars))]
  call_pars <- mget(args_must, envir = fun_env, inherits = FALSE)
  call_pars <- call_pars[sort(names(call_pars))]
  ok_pars <- identical(call_pars, cache_pars)
  
  # suboptimal for handling matchMS_noenzyme()
  if ((!ok_pars) && enzyme == "noenzyme") 
    ok_pars <- TRUE
  
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
                  "\\.xls$", "\\.csv$", "\\.txt$", "\\.tsv$", 
                  "^mgf$", "^mgfs$", "Calls", 
                  
                  # in case calling from proteoQ with MSGF workflows
                  "fraction_scheme.rda", "label_scheme.rda", 
                  "label_scheme_full.rda"))

    fi_mgf <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mgf$")
    fi_mzml <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mzML$")
    len_mgf <- length(fi_mgf)
    len_mzml <- length(fi_mzml)
    
    if (len_mgf && len_mzml)
      stop("Peak lists need to be in either MGF or mzML, but not both.")
    
    filelist <- if (len_mgf) fi_mgf else fi_mzml
    
    if (len_mgf) {
      readMGF(filepath = mgf_path,
              filelist = filelist, 
              min_mass = min_mass,
              max_mass = max_mass, 
              min_ms2mass = min_ms2mass,
              max_ms2mass = max_ms2mass, 
              topn_ms2ions = topn_ms2ions,
              quant = quant, 
              ms1_charge_range = c(min_ms1_charge, max_ms1_charge), 
              ms1_scan_range = c(min_scan_num, max_scan_num), 
              ret_range = c(min_ret_time, max_ret_time),
              ppm_ms1 = ppm_ms1_new,
              ppm_ms2 = ppm_ms2_new,
              tmt_reporter_lower = tmt_reporter_lower, 
              tmt_reporter_upper = tmt_reporter_upper, 
              exclude_reporter_region = exclude_reporter_region, 
              index_mgf_ms2 = index_mgf_ms2, 
              mgf_cutmzs = mgf_cutmzs, 
              mgf_cutpercs = mgf_cutpercs, 
              out_path = rds, 
              digits = digits)
    }
    else if (len_mzml) {
      warning("Please uncheck \"Use zlib compression\" with mzML from MSConvert.", 
              call. = FALSE)
      
      readmzML(filepath = mgf_path,
               filelist = filelist, 
               min_mass = min_mass,
               max_mass = max_mass, 
               min_ms2mass = min_ms2mass,
               max_ms2mass = max_ms2mass, 
               topn_ms2ions = topn_ms2ions,
               quant = quant, 
               ms1_charge_range = c(min_ms1_charge, max_ms1_charge), 
               ms1_scan_range = c(min_scan_num, max_scan_num), 
               ret_range = c(min_ret_time, max_ret_time),
               ppm_ms1 = ppm_ms1_new,
               ppm_ms2 = ppm_ms2_new,
               tmt_reporter_lower = tmt_reporter_lower, 
               tmt_reporter_upper = tmt_reporter_upper, 
               exclude_reporter_region = exclude_reporter_region, 
               index_mgf_ms2 = index_mgf_ms2, 
               mgf_cutmzs = mgf_cutmzs, 
               mgf_cutpercs = mgf_cutpercs, 
               out_path = rds, 
               digits = digits)
    }
    else {
      stop("No files of peak lists found.")
    }
    
    .savecall <- TRUE
  }

  invisible(NULL)
}


#' Reads MGF files in chunks.
#'
#' @param filepath The file path to a list of MGF or mzML files.
#' @param filelist A list of MGF or mzML files.
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
readMGF <- function (filepath = NULL, filelist = NULL, 
                     min_mass = 200L, max_mass = 4500L, 
                     min_ms2mass = 115L, max_ms2mass = 4500L, 
                     topn_ms2ions = 100L, quant = "none", 
                     ms1_charge_range = c(2L, 6L), 
                     ms1_scan_range = c(1L, .Machine$integer.max), 
                     ret_range = c(0, Inf), ppm_ms1 = 10L, ppm_ms2 = 10L, 
                     tmt_reporter_lower = 126.1, 
                     tmt_reporter_upper = 135.2, 
                     exclude_reporter_region = FALSE, 
                     index_mgf_ms2 = FALSE, 
                     mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                     out_path = file.path(filepath, "mgf_queries.rds"), 
                     digits = 4L) 
{
  ## Parsing rules
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
    qs::qsave(ans, file.path(filepath, "info_format.rds"), preset = "fast")
  })

  ## Reads MGF into chunks
  # separate parallel process: 
  # (1) one large MGF file and parallel chunks
  # (2) parallel five MGF files and parallel chunks in each
  len <- length(filelist)
  n_cores <- min(len, detect_cores(32L))

  if (n_cores == 1L)
    raw_files <- readlineMGFs(1, filelist, filepath, raw_file)
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    raw_files <- parallel::clusterMap(cl, readlineMGFs, 
                                      1:len, filelist, 
                                      MoreArgs = list(filepath = filepath, 
                                                      raw_file = raw_file), 
                                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    parallel::stopCluster(cl)
  }

  ## Reads from chunks
  out <- vector("list", len)
  
  for (i in seq_along(filelist)) {
    file <- filelist[i]
    temp_dir <- file.path(filepath, paste0("temp_", i))
    
    message("Loading '", file, "'.")
    
    out[[i]] <- read_mgf_chunks(filepath = temp_dir,
                                topn_ms2ions = topn_ms2ions,
                                ms1_charge_range = ms1_charge_range, 
                                ms1_scan_range = ms1_scan_range, 
                                ret_range = ret_range,
                                min_mass = min_mass, 
                                max_mass = max_mass, 
                                ppm_ms1 = ppm_ms1, 
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                max_ms2mass = max_ms2mass, 
                                mgf_cutmzs = mgf_cutmzs, 
                                mgf_cutpercs = mgf_cutpercs, 
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
                                raw_file = raw_files[[i]], 
                                tmt_reporter_lower = tmt_reporter_lower, 
                                tmt_reporter_upper = tmt_reporter_upper, 
                                exclude_reporter_region = exclude_reporter_region, 
                                index_mgf_ms2 = index_mgf_ms2, 
                                quant = quant, 
                                digits = digits)
    
    local({
      dir2 <- file.path(filepath, gsub("\\.[^.]*$", "", file))
      dir.create(dir2, showWarnings = FALSE)
      dir2 <- find_dir(dir2)
      
      if (fs::file_exists(dir2)) 
        fs::file_delete(dir2)
      
      fs::file_move(temp_dir, dir2)
    })
  }
  
  ## Clean up
  out <- dplyr::bind_rows(out)
  
  post_readmgf(out, min_mass = min_mass, max_mass = max_mass, ppm_ms1 = ppm_ms1, 
               filepath = filepath, out_path = out_path)
}


#' Post-processing of MGF or mzML
#' 
#' Calculates mass \code{frame}s etc.
#' 
#' @param df A data frame of processed peak lists.
#' @inheritParams readMGF
post_readmgf <- function (df, min_mass = 200L, max_mass = 4500L, ppm_ms1 = 10L, 
                          filepath, out_path) 
{
  df <- dplyr::arrange(df, ms1_mass)
  # df <- dplyr::filter(df, ms1_mass >= min_mass, ms1_mass <= max_mass)
  df <- dplyr::mutate(df, frame = find_ms1_interval(ms1_mass, from = min_mass, ppm = ppm_ms1))

  raws_files <- df$raw_file
  raws <- raws_files[!duplicated.default(raws_files)]
  inds <- seq_along(raws)
  names(inds) <- raws
  qs::qsave(inds, file.path(filepath, "raw_indexes.rds"), preset = "fast")
  df$raw_file <- unname(inds[raws_files])
  
  scans <- df$scan_title
  inds2 <- seq_along(scans)
  names(inds2) <- scans
  qs::qsave(inds2, file.path(filepath, "scan_indexes.rds"), preset = "fast")
  df$scan_title <- unname(inds2[scans])
  
  qs::qsave(df, out_path, preset = "fast")
  
  invisible(NULL)
}


#' Helper of \link{readMGF}.
#'
#' @param i An index of the i-th file.
#' @param file An MGF file name.
#' @inheritParams readMGF
#' @inheritParams read_mgf_chunks
#' @return Updated raw_file (for Bruker's timsTOF). Otherwise, raw_file remains
#'   NULL.
readlineMGFs <- function (i, file, filepath, raw_file) 
{
  f <- function(x, pos) {
    nm <- file.path(filepath, temp_i, paste0("chunk", "_", pos, ".mgf"))
    writeLines(x, nm)
  }
  
  message("Loading '", file, "'.")
  
  temp_i <- paste0("temp_", i)
  
  # refresh at every `i`
  temp_dir <- local({
    path <- file.path(filepath, temp_i)
    ok <- find_dir(path)
    
    if (!is.null(ok)) 
      fs::file_delete(ok)
    
    dir.create(path, showWarnings = FALSE)
    find_dir(path)
  })
  
  readr::read_lines_chunked(file = file.path(filepath, file),
                            callback = SideEffectChunkCallback$new(f),
                            chunk_size = 1000000L)
  
  # for "default_pasef" format
  if (!is.null(raw_file)) {
    raw_file <- local({
      file <- file.path(temp_dir, "chunk_1.mgf")
      hdr <- readLines(file, 50L)
      pat <- "^COM="
      line_file <- hdr[grepl(pat, hdr)]
      gsub(pat, "", line_file)
    })
  }
  
  invisible(raw_file)
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
read_mgf_chunks <- function (filepath = "~/mzion/mgf/temp_1",
                             topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                             ms1_scan_range = c(1L, .Machine$integer.max), 
                             ret_range = c(0, Inf), min_mass = 200L, 
                             max_mass = 4500L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                             min_ms2mass = 115L, max_ms2mass = 4500L, 
                             mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                             type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                             n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L, 
                             n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L, 
                             n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                             sep_pepmass = " ", nfields_pepmass = 2L, 
                             raw_file = NULL, 
                             tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                             exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                             quant = "none", digits = 4L) 
{
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")
  len <- length(filelist)

  if (!len) 
    stop("No mgf files under ", filepath, call. = FALSE)

  n_cores <- min(detect_cores(32L), len)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

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
      "index_mz", 
      "integerize_ms2ints", 
      "find_ms1_interval"), 
    envir = environment(mzion:::proc_mgf_chunks)
  )

  out <- parallel::clusterApply(cl, file.path(filepath, filelist),
                                proc_mgf_chunks,
                                topn_ms2ions = topn_ms2ions,
                                ms1_charge_range = ms1_charge_range, 
                                ms1_scan_range = ms1_scan_range, 
                                ret_range = ret_range,
                                min_mass = min_mass, 
                                max_mass = max_mass, 
                                ppm_ms1 = ppm_ms1,
                                ppm_ms2 = ppm_ms2,
                                min_ms2mass = min_ms2mass,
                                max_ms2mass = max_ms2mass, 
                                mgf_cutmzs = mgf_cutmzs, 
                                mgf_cutpercs = mgf_cutpercs, 
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
                                raw_file = raw_file, 
                                tmt_reporter_lower = tmt_reporter_lower, 
                                tmt_reporter_upper = tmt_reporter_upper, 
                                exclude_reporter_region = exclude_reporter_region, 
                                index_mgf_ms2 = index_mgf_ms2, 
                                quant = quant, 
                                digits = digits)
  
  parallel::stopCluster(cl)
  # gc()

  out <- dplyr::bind_rows(out)

  # adds back broken mgf entries
  afs <- local({
    afs <- list.files(path = file.path(filepath), pattern = "^.*\\_af.mgf$")
    idxes <- sort(as.integer(gsub("^chunk_(\\d+)_af\\.mgf", "\\1", afs)))
    afs <- paste0("chunk_", idxes, "_af.mgf")
    afs <- afs[-length(afs)]
  })

  bfs <- local({
    bfs <- list.files(path = file.path(filepath), pattern = "^.*\\_bf.mgf$")
    idxes <- sort(as.integer(gsub("^chunk_(\\d+)_bf\\.mgf", "\\1", bfs)))
    bfs <- paste0("chunk_", idxes, "_bf.mgf")
    bfs <- bfs[-1]
  })

  # stopifnot(length(afs) == length(bfs))

  gaps <- purrr::map2(afs, bfs, function (x, y) {
    af <- stringi::stri_read_lines(file.path(filepath, x))
    bf <- stringi::stri_read_lines(file.path(filepath, y))
    ab <- append(af, bf)
    
    # perfect case of no gaps: two lines of "" and ""
    if (length(ab) > 2L) ab else NULL
  })
  
  gaps <- unlist(gaps, use.names = FALSE)
  write(gaps, file.path(filepath, "gaps.mgf"))

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
                min_mass = min_mass, 
                max_mass = max_mass, 
                ppm_ms2 = ppm_ms2,
                min_ms2mass = min_ms2mass,
                max_ms2mass = max_ms2mass, 
                mgf_cutmzs = mgf_cutmzs, 
                mgf_cutpercs = mgf_cutpercs, 
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
                raw_file = raw_file, 
                tmt_reporter_lower = tmt_reporter_lower, 
                tmt_reporter_upper = tmt_reporter_upper, 
                exclude_reporter_region = exclude_reporter_region, 
                index_mgf_ms2 = index_mgf_ms2, 
                quant = quant, 
                digits = digits)
    )
  }

  if (type_mgf == "default_pasef") {
    out <- dplyr::mutate(out, scan_id = as.character(scan_num), 
                         scan_num = as.character(row_number()))
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
                             ret_range = c(0, Inf), min_mass = 200L, 
                             max_mass = 4500L, ppm_ms1 = 10L, ppm_ms2 = 10L, 
                             min_ms2mass = 115L, max_ms2mass = 4500L, 
                             mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                             type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                             n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L, 
                             n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L, 
                             n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                             sep_pepmass = " ", nfields_pepmass = 2L, 
                             raw_file = NULL, 
                             tmt_reporter_lower = 126.1, 
                             tmt_reporter_upper = 135.2, 
                             exclude_reporter_region = FALSE, 
                             index_mgf_ms2 = FALSE, 
                             quant = "none", digits = 4L) 
{
  message("Parsing '", file, "'.")
  lines <- stringi::stri_read_lines(file)
  basename <- gsub("\\.[^.]*$", "", file)
  
  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends   <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))

  af <- local({
    le <- ends[length(ends)]
    lb <- begins[length(begins)]
    af <- if (lb > le) lines[(le + n_spacer + 1L):length(lines)] else NULL
    write(af, file.path(paste0(basename, "_af.mgf")))

    af
  })

  bf <- local({
    le <- ends[1]
    lb <- begins[1]
    bf <- if (lb > le) lines[1:(le + n_spacer)] else NULL
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
                   min_mass = min_mass, 
                   max_mass = max_mass, 
                   ppm_ms1 = ppm_ms1,
                   ppm_ms2 = ppm_ms2,
                   min_ms2mass = min_ms2mass,
                   max_ms2mass = max_ms2mass, 
                   mgf_cutmzs = mgf_cutmzs, 
                   mgf_cutpercs = mgf_cutpercs, 
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
                   raw_file = raw_file, 
                   tmt_reporter_lower = tmt_reporter_lower, 
                   tmt_reporter_upper = tmt_reporter_upper, 
                   exclude_reporter_region = exclude_reporter_region, 
                   index_mgf_ms2 = index_mgf_ms2, 
                   quant = quant, 
                   digits = digits)
}


#' Helper in processing MGF entries in chunks.
#'
#' @param lines Lines of MGF.
#' @inheritParams proc_mgf_chunks
proc_mgfs <- function (lines, topn_ms2ions = 100L, 
                       ms1_charge_range = c(2L, 6L), 
                       ms1_scan_range = c(1L, .Machine$integer.max), 
                       ret_range = c(0, Inf), min_mass = 200L, max_mass = 4500L, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                       n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L,
                       n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L,
                       n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                       sep_pepmass = " ", nfields_pepmass = 2L, raw_file = NULL, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       quant = "none", digits = 4L) 
{
  options(digits = 9L)

  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))

  ## MS1 
  # (1) m-over-z and intensity
  ms1s <- stringi::stri_replace_first_fixed(lines[begins + n_to_pepmass], 
                                            "PEPMASS=", "")
  ms1s <- lapply(ms1s, stringi::stri_split_fixed, pattern = sep_pepmass, 
                 n = nfields_pepmass, simplify = TRUE)

  ms1_moverzs <- lapply(ms1s, function (x) as.numeric(x[, 1]))
  ms1_moverzs <- .Internal(unlist(ms1_moverzs, recursive = FALSE, use.names = FALSE))
  # not as.integer; intensity may be > .Machine$integer.max (2147483647)
  ms1_ints <- lapply(ms1s, function (x) round(as.numeric(x[, 2]), digits = 1L))
  ms1_ints <- .Internal(unlist(ms1_ints, recursive = FALSE, use.names = FALSE))

  rm(list = c("ms1s"))
  gc()
  
  # (2) retention time
  ret_times <- stringi::stri_replace_first_fixed(lines[begins + n_to_rt], "RTINSECONDS=", "")
  ret_times <- as.numeric(ret_times)

  # (3) MS1 charges and masses
  ms1_charges <- stringi::stri_replace_first_fixed(lines[begins + n_to_charge], "CHARGE=", "")
  charges <- lapply(ms1_charges, stringi::stri_reverse)
  charges <- .Internal(unlist(charges, recursive = FALSE, use.names = FALSE))
  # timsTOF no CHARGE line -> NAs introduced by coercion
  charges <- as.integer(charges)
  
  ms1_masses <- mapply(function (x, y) x * y - y * 1.00727647, 
                       ms1_moverzs, charges, 
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)
  ms1_masses <- round(ms1_masses, digits = digits)
  ms1_moverzs <- round(ms1_moverzs, digits = digits)

  rows <- (charges >= ms1_charge_range[1] & charges <= ms1_charge_range[2] & 
             ret_times >= ret_range[1] & ret_times <= ret_range[2] & 
             ms1_masses >= min_mass & ms1_masses <= max_mass & 
             # timsTOF: no MS1 masses
             !is.na(ms1_masses))
  
  # timsTOF data may have undetermined charge states
  na_rows <- .Internal(which(is.na(rows)))
  if (length(na_rows)) rows[na_rows] <- FALSE

  begins <- begins[rows]
  ends <- ends[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- round(ret_times[rows], digits = 2L)
  ms1_masses <- ms1_masses[rows]
  charges <- charges[rows]

  # Others
  scan_titles <- 
    stringi::stri_replace_first_fixed(lines[begins + n_to_title], "TITLE=", "")
  
  if (type_mgf %in% c("msconv_thermo", "msconv_pasef")) {
    raw_files <- 
      stringi::stri_replace_first_regex(scan_titles, "^.* File:\"([^\"]+)\".*", "$1")
    scan_nums <- 
      stringi::stri_replace_first_regex(scan_titles, 
                                        "^.*\\.(\\d+)\\.\\d+\\.\\d+ File:\".*", 
                                        "$1")
  } 
  else if (type_mgf == "pd") {
    raw_files <- gsub("^.*File: \"([^\"]+)\".*", "\\1", scan_titles)
    raw_files <- gsub("\\\\", "/", raw_files)
    raw_files <- gsub("^.*/(.*)", "\\1", raw_files)
    scan_nums <- gsub("^.* scans: \"([0-9]+)\"$", "\\1", scan_titles)
  } 
  else if (type_mgf == "default_pasef") {
    # one raw_file one .d file guaranteed
    raw_files <- rep(raw_file, length(begins))
    scan_nums <- stringi::stri_replace_first_fixed(lines[begins + n_to_scan], "RAWSCANS=", "")
  } 
  else {
    stop("Unknown MGF format.", call. = FALSE)
  }

  ## MS2
  # (-1L: one line above "END IONS")
  ms2s <- mapply(function (x, y) lines[(x + n_hdr) : (y - 1L)], 
                 begins, ends, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ms2s <- lapply(ms2s, stringi::stri_split_fixed, pattern = sep_ms2s, 
                 n = nfields_ms2s, simplify = TRUE)

  ms2_moverzs <- lapply(ms2s, function (x) as.numeric(x[, 1]))
  # not as.integer; intensity may be > .Machine$integer.max
  ms2_ints <- lapply(ms2s, function (x) as.numeric(x[, 2]))
  rm(list = c("ms2s"))
  
  # extract the TMT region of MS2 moverz and intensity
  # (also convert reporter-ion intensities to integers)
  restmt <- extract_mgf_rptrs(ms2_moverzs = ms2_moverzs, 
                              ms2_ints = ms2_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  
  ms2_moverzs <- restmt[["ms2_moverzs"]]
  ms2_ints <- restmt[["ms2_ints"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  rm(list = "restmt")

  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = ms2_moverzs, 
                          ms2_ints = ms2_ints, 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  
  ms2_moverzs <- mz_n_int[["ms2_moverzs"]]
  ms2_ints <- mz_n_int[["ms2_ints"]]
  lens <- mz_n_int[["lens"]]
  rm(list = "mz_n_int")
  
  if (index_mgf_ms2) {
    ms2_moverzs <- lapply(ms2_moverzs, index_mz, min_ms2mass, ppm_ms2/1E6)
    # dups <- lapply(ms2_moverzs, duplicated.default)
    # ms2_moverzs <- mapply(function (x, y) x[!y], ms2_moverzs, dups, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    # ms2_ints <- mapply(function (x, y) x[!y], ms2_ints, dups, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    # rm(list = "dups")
  }
  else
    ms2_moverzs <- lapply(ms2_moverzs, round, digits = digits)

  tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_files,
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_int = ms1_ints,
    ms1_charge = ms1_charges,
    ret_time = ret_times,
    scan_num = scan_nums,
    ms2_moverz = ms2_moverzs,
    ms2_int = ms2_ints,
    ms2_n = lens, 
    # charge = charges, 
    rptr_moverz = rptr_moverzs, 
    rptr_int = rptr_ints, 
    )
}


#' Subsets MGFs by top-n
#' 
#' \code{lens} after filtered by \code{min_ms2mass} but before subset by 
#' \code{topn_ms2ions} to reflect noise levels.
#' 
#' @param ms2_moverzs Lists of MS2 moverz values.
#' @param ms2_ints Lists of MS2 intensities
#' @inheritParams load_mgfs
sub_mgftopn <- function (ms2_moverzs, ms2_ints, topn_ms2ions = 100L, 
                         mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                         min_ms2mass = 115L, max_ms2mass = 4500L) 
{
  options(digits = 9L)
  
  ## subsets by min_ms2mass
  oks <- lapply(ms2_moverzs, function (x) x >= min_ms2mass)
  
  ms2_moverzs <- mapply(function (x, y) x[y], ms2_moverzs, oks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ms2_ints <- mapply(function (x, y) x[y], ms2_ints, oks, 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  rm(list = c("oks"))
  
  ## subsets by topn
  lens <- lapply(ms2_moverzs, length)
  lens <- .Internal(unlist(lens, recursive = FALSE, use.names = FALSE))
  
  if (topn_ms2ions < Inf) {
    is_long <- lens > topn_ms2ions
    
    if (length(mgf_cutmzs)) {
      m_long <- ms2_moverzs[is_long]
      i_long <- ms2_ints[is_long]
      
      for (i in seq_along(m_long)) {
        x <- m_long[[i]]
        y <- i_long[[i]]
        
        # (`<` not `<=`)
        ok_ms2 <- x < max_ms2mass
        x <- x[ok_ms2]
        y <- y[ok_ms2]
        
        idxes <- findInterval(x, mgf_cutmzs)
        xs <- split(x, idxes)
        ys <- split(y, idxes)
        
        # some zones may have no entries
        ok_idxes <- as.integer(names(xs)) + 1L
        ok_percs <- mgf_cutpercs[ok_idxes]
        
        # e.g. the last interval from ms2masses >= max_ms2mass -> c(100, 10, NA)
        # but no need to remove NA since already `ok_ms2`;
        # 
        # `which_topx2` also guard against NA
        # 
        # ok_percs <- ok_percs[!is.na(ok_percs)]
        
        ## (little benefit with padding) 
        if (FALSE) {
          percs_no_one <- ok_percs[-1]
          ys_no_one <- ys[-1]
          rows_no_one <- mapply(which_topx2,ys_no_one, percs_no_one, 
                                SIMPLIFY = FALSE, USE.NAMES = FALSE)
          
          cts_no_one <- lapply(rows_no_one, length)
          cts_no_one <- .Internal(unlist(cts_no_one, recursive = FALSE, use.names = FALSE))
          cts_delta <- sum(percs_no_one) - sum(cts_no_one)
          
          cts_one <- ok_percs[1] + cts_delta
          rows_one <- which_topx2(ys[[1]], cts_one)
          
          rows <- c(list(rows_one), rows_no_one)
        }
        ##
        
        rows <- mapply(which_topx2, ys, ok_percs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        ans_x <- mapply(function (x, y) x[y], xs, rows, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        ans_y <- mapply(function (x, y) x[y], ys, rows, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        
        ans_x <- .Internal(unlist(ans_x, recursive = FALSE, use.names = FALSE))
        ans_y <- .Internal(unlist(ans_y, recursive = FALSE, use.names = FALSE))
        
        ms2_moverzs[is_long][[i]] <- ans_x
        ms2_ints[is_long][[i]] <- ans_y
        
        if (i %% 5000L == 0) gc()
      }
      
      rm(list = c("m_long", "i_long"))
    }
    else {
      rows <- lapply(ms2_ints[is_long], which_topx2, topn_ms2ions)
      
      ms2_ints[is_long] <- mapply(function (x, y) x[y], 
                                  ms2_ints[is_long], rows, 
                                  SIMPLIFY = FALSE, USE.NAMES = FALSE)
      
      ms2_moverzs[is_long] <- mapply(function (x, y) x[y], 
                                     ms2_moverzs[is_long], rows, 
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
      
      rm(list = c("rows"))
    }
    
    rm(list = c("is_long"))
  }
  
  # also handles MS2 intensity max-outs, which usually don't happen
  ms2_ints <- integerize_ms2ints(ms2_ints)

  list(ms2_moverzs = ms2_moverzs, ms2_ints = ms2_ints, lens = lens)
}


#' Integerizes the MS2 intensities.
#'
#' Also guards against intensity integers above machine maximum (32-bit).
#'
#' @param ms2_ints Lists of MS2 intensities
#' @param max_intdbl Maximum 32-bit integer as double:
#'   \code{as.double(.Machine$integer.max)}.
integerize_ms2ints <- function(ms2_ints, max_intdbl = 2147483647.0) 
{
  lapply(ms2_ints, function (x) {
    x[x > max_intdbl] <- max_intdbl
    as.integer(x)
  })
}


#' Extracts reporter-ion data from MGF.
#' 
#' Also purges MS2 m-over-z and intensity when applicable.
#' 
#' @param ms2_moverzs Lists of MS2 m-over-z values.
#' @param ms2_ints Lists of MS2 intensity values. 
#' @inheritParams matchMS
extract_mgf_rptrs <- function (ms2_moverzs, ms2_ints, quant = "none", 
                               tmt_reporter_lower = 126.1, 
                               tmt_reporter_upper = 135.2, 
                               exclude_reporter_region = FALSE) 
{
  if (grepl("^tmt.*\\d+", quant)) {
    ok_rptrs <- lapply(ms2_moverzs, function (x) x > tmt_reporter_lower & 
                         x < tmt_reporter_upper)
    
    no_rptrs <- lapply(ok_rptrs, `!`)
    
    rptr_moverzs <- mapply(function (x, y) x[y], ms2_moverzs, ok_rptrs, 
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    rptr_ints <- mapply(function (x, y) x[y], ms2_ints, ok_rptrs, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    rptr_ints <- integerize_ms2ints(rptr_ints)
    
    if (exclude_reporter_region) {
      ms2_moverzs <- mapply(function (x, y) x[y], ms2_moverzs, no_rptrs, 
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
      
      ms2_ints <- mapply(function (x, y) x[y], ms2_ints, no_rptrs, 
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
  }
  else {
    rptr_moverzs <- NA_real_
    rptr_ints <- NA_integer_
  }
  
  list(ms2_moverzs = ms2_moverzs, 
       ms2_ints = ms2_ints,  
       rptr_moverzs = rptr_moverzs, 
       rptr_ints = rptr_ints)
}


#' Calculates the frame numbers for a list of experimental MS1 mass by
#' intervals.
#'
#' Needs correct \code{from}.
#' @param mass Numeric; a list of MS1 masses.
#' @inheritParams find_ms1_cutpoints
#' @examples
#' \donttest{
#' library(mzion)
#' mzion:::find_ms1_interval(c(500, 800.1))
#' }
#' @return Frame numbers.
#' @seealso find_ms1_cutpoints
find_ms1_interval <- function (mass = 1800.0, from = 115L, ppm = 10L) 
{
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+ppm/1e6))
}


#' Converts ms2_moverz to integers.
#'
#' @param from Numeric; the starting MS1 mass.
#' @param x Numeric; MS2 mass.
#' @param d Numeric; \eqn{ppm * 10E-6}.
#' @seealso find_ms1_interval
index_mz <- function (x, from = 115L, d = 1E-5) ceiling(log(x/from)/log(1+d))


#' Finds the type of MGF.
#'
#' @param file The path to an MGF file.
find_mgf_type <- function (file) 
{
  hdr <- readLines(file, 5000L)
  begins <- which(stringi::stri_startswith_fixed(hdr, "BEGIN IONS"))
  ends <- which(stringi::stri_endswith_fixed(hdr, "END IONS"))
  
  len_h <- length(hdr)
  len_b <- length(begins)
  
  if (!len_b) 
    stop("Check corrupted files: the tag of `BEGIN IONS` not found in MGF.")
  else if (len_b == 1L) {
    b2 <- len_h + 1L
    hdr[b2] <- "BEGIN IONS"
    begins <- c(begins, b2)
    rm(list = "b2")
  }
    
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
  
  # between ends[1] and begins[2]
  n_end <- 1L
  n_begin <- 1L
  
  end <- ends[1]
  begin2 <- begins[2]
  lines_bf <- hdr[end:begin2]
  
  n_spacer <- sum(unlist(lapply(lines_bf, function (x) x == "")))
  n_bf_begin <- begin2 - end - n_spacer - n_end
  
  # between begins[1] and ends[1]
  begin <- begins[1]
  lines <- hdr[begin:end]
  
  n_hdr <- grep("^[0-9]+", lines)[1] - 1L # 5
  n_to_pepmass <- grep("^PEPMASS", lines)[1] - n_begin
  n_to_title <- grep("^TITLE", lines)[1] - n_begin
  n_to_rt <- grep("^RTINSECONDS", lines)[1] - n_begin
  n_to_charge <- grep("^CHARGE", lines)[1] - n_begin

  if (type == "msconv_thermo") {
    # n_bf_begin <- 0L
    # n_spacer <- 0L
    # n_hdr <- 5L
    # n_to_pepmass <- 3L
    # n_to_title <- 1L
    # n_to_rt <- 2L
    # n_to_charge <- 4L
    n_to_scan <- 0L
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " "
    nfields_pepmass <- 2L
    raw_file <- NULL
  } 
  else if (type == "msconv_pasef") {
    warning("Suggest MGFs of manufacturer's default")
    
    n_bf_begin <- 0L
    n_spacer <- 0L
    n_hdr <- 5L
    n_to_pepmass <- 3L
    n_to_title <- 1L
    n_to_rt <- 2L
    n_to_charge <- 4L
    n_to_scan <- 0L
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " " # as of 2021-12-27, missing MS1 intensity
    nfields_pepmass <- 2L # "NA" MS1 intensities after parsing
    raw_file <- NULL
  } 
  else if (type == "pd") {
    # n_bf_begin <- 0L
    # n_spacer <- 1L
    # n_hdr <- 6L
    # n_to_pepmass <- 2L
    # n_to_title <- 1L
    # n_to_rt <- 4L
    # n_to_charge <- 3L
    # n_to_scan <- 5L
    
    n_to_scan <- grep("^SCANS", lines)[1] - n_begin
    if (is.na(n_to_scan)) stop("`SCANS` not found.")
    
    sep_ms2s <- " "
    nfields_ms2s <- 2L
    sep_pepmass <- " "
    nfields_pepmass <- 2L
    raw_file <- NULL
  } 
  else if (type == "default_pasef") {
    # n_bf_begin <- 4L
    # n_spacer <- 1L
    # n_hdr <- 6L
    # n_to_pepmass <- 4L
    # n_to_title <- 1L
    # n_to_scan <- 3L
    # n_to_rt <- 2L
    # n_to_charge <- 5L
    
    n_to_scan <- grep("^RAWSCANS", lines)[1] - n_begin
    if (is.na(n_to_scan)) stop("`RAWSCANS` not found.")
    
    sep_ms2s <- "\t"
    nfields_ms2s <- 3L
    sep_pepmass <- "\t"
    nfields_pepmass <- 2L
    raw_file <- "."
  }
  else {
    stop("Unkown format of MGFs.")
  }

  if (is.na(n_hdr)) stop("MS2 data not found.")
  if (is.na(n_to_pepmass)) stop("`PEPMASS` not found.")
  if (is.na(n_to_title)) stop("`TITLE` not found.")
  if (is.na(n_to_rt)) stop("`RTINSECONDS` not found.")
  if (is.na(n_to_charge)) stop("`CHARGE` not found.")
  
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


#' Reads mzML files.
#'
#' @inheritParams readMGF
readmzML <- function (filepath = NULL, filelist = NULL, 
                      min_mass = 200L, max_mass = 4500L, 
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      topn_ms2ions = 100L, quant = "none", 
                      ms1_charge_range = c(2L, 6L), 
                      ms1_scan_range = c(1L, .Machine$integer.max), 
                      ret_range = c(0, Inf), ppm_ms1 = 10L, ppm_ms2 = 10L, 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, 
                      index_mgf_ms2 = FALSE, 
                      mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                      out_path, digits = 4L)
{
  len <- length(filelist)
  out <- vector("list", len)
  
  # parallel here
  files <- file.path(filepath, filelist)
  sizes <- max(unlist(lapply(files, file.size)))/1024^3
  # n_cores <- min(detect_cores(16L), floor((max_ram <- 32)/(sizes * 8)), len)
  n_cores <- min(detect_cores(32L), floor((find_free_mem()/1024)/(sizes * 8)), len)
  n_cores <- max(1L, n_cores)
  
  if (n_cores == 1L) {
    for (i in 1:len) {
      out[[i]] <- proc_mzml(files[[i]], 
                            topn_ms2ions = topn_ms2ions, 
                            ms1_charge_range = ms1_charge_range, 
                            ret_range = ret_range, 
                            min_mass = min_mass, 
                            max_mass = max_mass, 
                            ppm_ms1 = ppm_ms1, 
                            ppm_ms2 = ppm_ms2, 
                            min_ms2mass = min_ms2mass, 
                            max_ms2mass = max_ms2mass, 
                            mgf_cutmzs = mgf_cutmzs, 
                            mgf_cutpercs = mgf_cutpercs, 
                            tmt_reporter_lower = tmt_reporter_lower, 
                            tmt_reporter_upper = tmt_reporter_upper, 
                            exclude_reporter_region = exclude_reporter_region, 
                            index_mgf_ms2 = index_mgf_ms2, 
                            quant = quant, 
                            digits = digits)
    }
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    out <- parallel::clusterApply(cl, files, proc_mzml, 
                                  topn_ms2ions = topn_ms2ions, 
                                  ms1_charge_range = ms1_charge_range, 
                                  ret_range = ret_range, 
                                  min_mass = min_mass, 
                                  max_mass = max_mass, 
                                  ppm_ms1 = ppm_ms1, 
                                  ppm_ms2 = ppm_ms2, 
                                  min_ms2mass = min_ms2mass, 
                                  max_ms2mass = max_ms2mass, 
                                  mgf_cutmzs = mgf_cutmzs, 
                                  mgf_cutpercs = mgf_cutpercs, 
                                  tmt_reporter_lower = tmt_reporter_lower, 
                                  tmt_reporter_upper = tmt_reporter_upper, 
                                  exclude_reporter_region = exclude_reporter_region, 
                                  index_mgf_ms2 = index_mgf_ms2, 
                                  quant = quant, 
                                  digits = digits)
    parallel::stopCluster(cl)
  }
  
  out <- dplyr::bind_rows(out)
  
  ## frame number
  post_readmgf(out, min_mass = min_mass, max_mass = max_mass, ppm_ms1 = ppm_ms1, 
               filepath = filepath, out_path = out_path)
}


#' Helper of \link{readmzML}
#' 
#' No scan range subsetting with PASEF timsTOF.
#' 
#' @param file A file name to mzML with a prepending path.
#' @inheritParams readmzML
proc_mzml <- function (file, topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                       ret_range = c(0, Inf), min_mass = 200L, max_mass = 4500L, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       quant = "none", digits = 4L) 
{
  message("Loading ", file)
  
  df <- read_mzml(file, tmt_reporter_lower = tmt_reporter_lower, 
                  tmt_reporter_upper = tmt_reporter_upper, 
                  exclude_reporter_region = exclude_reporter_region, 
                  index_mgf_ms2 = index_mgf_ms2, ppm_ms1 = ppm_ms1, 
                  ppm_ms2 = ppm_ms2, min_ms2mass = min_ms2mass, 
                  max_ms2mass = max_ms2mass, quant = quant, digits = digits)
  df <- df[with(df, !is.na(ms1_mass)), ]
  .Internal(gc(verbose = FALSE, reset = FALSE, full = TRUE))
  
  if (!"charge" %in% names(df)) {
    ms1_charges <- df[["ms1_charges"]]
    charges <- lapply(ms1_charges, stringi::stri_reverse)
    charges <- .Internal(unlist(charges, recursive = FALSE, use.names = FALSE))
    charges <- as.integer(charges)
    df <- dplyr::mutate(df, charge = charges)
    rm(list = c("charges", "ms1_charges"))
  }
  
  df <- dplyr::filter(df, 
                      charge >= ms1_charge_range[1], charge <= ms1_charge_range[2], 
                      ret_time >= ret_range[1], ret_time <= ret_range[2], 
                      ms1_mass >= min_mass, ms1_mass <= max_mass, )
  
  df[["charge"]] <- NULL
  
  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = df[["ms2_moverz"]], 
                          ms2_ints = df[["ms2_int"]], 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  
  # can be integers if "index_mgf_ms2 = TRUE"
  df[["ms2_moverz"]] <- mz_n_int[["ms2_moverzs"]]
  df[["ms2_int"]] <- mz_n_int[["ms2_ints"]]
  df[["ms2_n"]] <- mz_n_int[["lens"]]

  invisible(df)
}


#' Reads mzML from MSConvert
#' 
#' zlib compression need to be disabled when creating mzML from MSConvert.
#' 
#' @param xml_file A file name of mzML.
#' @inheritParams matchMS
read_mzml <- function (xml_file, tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, min_ms2mass = 115L, 
                       max_ms2mass = 4500L, quant = "none", digits = 4L)
{
  ## spectrum
  xml_root <- xml2::read_xml(xml_file)
  mzML <- xml2::xml_child(xml_root)
  idx_run <- which(xml2::xml_name(xml2::xml_children(mzML)) == "run")
  run <- xml2::xml_children(mzML)[[idx_run]]
  idx_specs <- which(xml2::xml_name(xml2::xml_children(run)) == "spectrumList")
  spec <- xml2::xml_children(xml2::xml_children(run)[[idx_specs]])
  rm(list = c("mzML", "idx_run", "run", "idx_specs"))
  
  len <- length(spec)
  scan_nums <- ms1_charges <- raw_files <- scan_titles <- character(len)
  ret_times <- ms1_ints <- ms1_masses <- ms1_moverzs <- numeric(len)
  ms2_moverzs <- ms2_ints <- vector("list", len)
  charges <- ms2_ns <- frames <- integer(len)
  
  for (i in seq_along(spec)) {
    x <- spec[[i]]
    scan_nums[i] <- gsub(".* scan=(.*)$", "\\1", xml2::xml_attr(x, "id"))
    xc <- xml2::xml_children(x)
    idx_precursor <- grep("precursorList", xc)
    rm(list = c("x"))
    
    if (length(idx_precursor)) {
      nms <- xml2::xml_attr(xc, "name")
      idx_title <- .Internal(which(nms == "spectrum title"))
      idx_scanList <- grep("scanList", xc) # 11
      idx_bin <- grep("binaryDataArrayList", xc)
      
      ## title
      title <- xml2::xml_attr(xc[[idx_title]], "value")
      scan_titles[i] <- title
      raw_files[i] <- gsub("(.*)\\.[0-9]+\\.[0-9]+\\.[0-9]+ File:.*", "\\1", title)
      
      ## retention
      scanList <- xml2::xml_children(xc[[idx_scanList]])
      idx_rt <- grep("scan", scanList) # 2
      scanList_scan <- xml2::xml_children(scanList[[idx_rt]])
      idx_scan_start <- 
        .Internal(which(xml2::xml_attr(scanList_scan, "name") == "scan start time"))
      ret_times[i] <- xml2::xml_attr(scanList_scan[[idx_scan_start]], "value")
      rm(list = c("nms", "title", "scanList_scan", "scanList"))
      
      ## precursorList
      precursorList <- xml2::xml_children(xc[[idx_precursor]])
      
      # (assume one precursor, not yet chimeric)
      precursor <- precursorList[[1]]
      precursorc <- xml2::xml_children(precursor)
      idx_selectedIonList <- grep("selectedIonList", precursorc)
      
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      ms1_moverzs[i] <- xml2::xml_attr(selectedIonc[[1]], "value")
      ms1_charges[i] <- xml2::xml_attr(selectedIonc[[2]], "value")
      
      # may be no precursor intensity
      ms1_ints[i] <- if (length(selectedIonc) > 2L)
        xml2::xml_attr(selectedIonc[[3]], "value")
      else
        numeric(1)
      
      rm(list = c("precursor", "precursorc", "idx_selectedIonList", 
                  "selectedIon", "selectedIonc"))
      
      ## binaryDataArrayList
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin]]))
      ms2s <- xml2::xml_contents(binData)
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(ms2s[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(ms2s[[2]]))
      ms2_n <- length(r1)/8L
      ms2_moverzs[[i]] <- readBin(r1, "double", n = ms2_n, size = 8L)
      ms2_ints[[i]] <- readBin(r2, "double", n = ms2_n, size = 8L)
      ms2_ns[i] <- ms2_n
      rm(list = c("r1", "r2", "ms2s", "binData", "ms2_n"))
    }
  }
  
  rows <- !(is.na(ms1_charges) | ms1_charges == "")
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
  ms2_ns <- ms2_ns[rows]
  
  charges <- as.integer(ms1_charges)
  ms1_charges <- paste0(ms1_charges, "+") # assume always "+" for now
  ms1_moverzs <- round(as.numeric(ms1_moverzs), digits = digits)
  ms1_masses <- round(ms1_moverzs * charges - charges * 1.00727647, digits = digits)
  # ms1_ints not "as.integer": may be > .Machine$integer.max (2147483647)
  ms1_ints <- round(as.numeric(ms1_ints), digits = 1L)
  ret_times <- round(as.numeric(ret_times) * 60, digits = 2L)
  
  # extract the TMT region of MS2 moverz and intensity
  # (also convert reporter-ion intensities to integers)
  restmt <- extract_mgf_rptrs(ms2_moverzs = ms2_moverzs, 
                              ms2_ints = ms2_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  
  ms2_moverzs <- restmt[["ms2_moverzs"]]
  ms2_ints <- restmt[["ms2_ints"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  rm(list = "restmt")

  ms2_moverzs <- if (index_mgf_ms2) 
    lapply(ms2_moverzs, index_mz, min_ms2mass, ppm_ms2/1E6)
  else
    lapply(ms2_moverzs, round, digits = digits)
  
  tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_files,
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_int = ms1_ints,
    ms1_charge = ms1_charges,
    ret_time = ret_times,
    scan_num = scan_nums,
    ms2_moverz = ms2_moverzs,
    ms2_int = ms2_ints,
    
    # before subset by min_ms2mass
    ms2_n = ms2_ns, 
    
    # temporarily kept
    charge = charges,
    
    rptr_moverz = rptr_moverzs, 
    rptr_int = rptr_ints, )
}


