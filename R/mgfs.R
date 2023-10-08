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
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, topn_ms2ions = 150L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1 = 20L, ppm_ms2 = 20L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       enzyme = "trypsin_p", 
                       is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       quant = "none", digits = 4L) 
{
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    }, add = TRUE)

  # ---
  this_call <- match.call()
  fun <- as.character(this_call[[1]])
  fun <- fun[length(fun)] # may be called as mzion:::load_mgfs
  fun_env <- environment()
  
  args_except <- c("out_path")
  args <- names(formals(fun))
  args_must <- args[!args %in% args_except]
  
  cache_pars <- find_callarg_vals(
    time = NULL, 
    path = file.path(out_path, "Calls"), 
    fun = paste0(fun, ".rda"), 
    args = args_must, 
    new_args = unlist(formals(load_mgfs)[c("enzyme")]))

  cache_pars <- cache_pars[sort(names(cache_pars))]
  call_pars <- mget(args_must, envir = fun_env, inherits = FALSE)
  call_pars <- call_pars[sort(names(call_pars))]
  ok_pars <- identical(call_pars, cache_pars)
  
  # suboptimal for handling matchMS_noenzyme()
  if ((!ok_pars) && enzyme == "noenzyme") 
    ok_pars <- TRUE
  
  # checks processed mgfs
  raws_indexes <- file.path(mgf_path, "raw_indexes.rds")
  
  if (file.exists(raws_indexes)) {
    raws <- qs::qread(raws_indexes)
    ques <- list.files(mgf_path, pattern = "^mgf_queries_\\d+\\.rds$")
    ok_mgfs <- if (length(raws) == length(ques)) TRUE else FALSE
    rm(list = c("raws", "ques", "raws_indexes"))
  }
  else {
    ok_mgfs <- FALSE
    rm(list = c("raws_indexes"))
  }

  if (ok_pars && ok_mgfs) {
    message("Found cached MGFs.")
    .savecall <- FALSE
    return(NULL)
  }
  
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
    ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.mzML$", "\\.xlsx$", 
                "\\.xls$", "\\.csv$", "\\.txt$", "\\.tsv$", "\\.pars$", 
                "^mgf$", "^mgfs$", "Calls", 
                # in case of reprocessing after proteoQ
                "fraction_scheme.rda", "label_scheme.rda", 
                "label_scheme_full.rda"))
  
  fi_mgf   <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mgf$")
  fi_mzml  <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mzML$")
  len_mgf  <- length(fi_mgf)
  len_mzml <- length(fi_mzml)
  
  if (len_mgf && len_mzml)
    stop("Peak lists need to be in either MGF or mzML, but not both.")
  
  filelist <- if (len_mgf) fi_mgf else fi_mzml
  
  if (len_mgf) {
    readMGF(filepath = mgf_path,
            out_path = out_path, 
            filelist = filelist, 
            min_mass = min_mass,
            max_mass = max_mass, 
            min_ms2mass = min_ms2mass,
            max_ms2mass = max_ms2mass, 
            topn_ms2ions = topn_ms2ions,
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
            is_mdda = FALSE, 
            use_defpeaks = use_defpeaks, 
            deisotope_ms2 = deisotope_ms2, 
            max_ms2_charge = max_ms2_charge, 
            maxn_dia_precurs = maxn_dia_precurs, 
            maxn_mdda_precurs = maxn_mdda_precurs, 
            n_mdda_flanks = n_mdda_flanks, 
            ppm_ms1_deisotope = ppm_ms1_deisotope, 
            ppm_ms2_deisotope = ppm_ms2_deisotope, 
            quant = quant, 
            digits = digits)
  }
  else if (len_mzml) {
    readmzML(filepath = mgf_path,
             out_path = out_path, 
             filelist = filelist, 
             min_mass = min_mass,
             max_mass = max_mass, 
             min_ms2mass = min_ms2mass,
             max_ms2mass = max_ms2mass, 
             topn_ms2ions = topn_ms2ions,
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
             
             deisotope_ms2 = deisotope_ms2, 
             max_ms2_charge = max_ms2_charge, 
             maxn_dia_precurs = maxn_dia_precurs, 
             
             is_mdda = is_mdda, 
             use_defpeaks = use_defpeaks, 
             maxn_mdda_precurs = maxn_mdda_precurs, 
             n_mdda_flanks = n_mdda_flanks, 
             ppm_ms1_deisotope = ppm_ms1_deisotope, 
             ppm_ms2_deisotope = ppm_ms2_deisotope, 
             
             quant = quant, 
             digits = digits)
  }
  else {
    stop("No peak lists of mzML or MGF formats found.")
  }
  
  .savecall <- TRUE

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
                     topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                     ms1_scan_range = c(1L, .Machine$integer.max), 
                     ret_range = c(0, Inf), ppm_ms1 = 10L, ppm_ms2 = 10L, 
                     tmt_reporter_lower = 126.1, 
                     tmt_reporter_upper = 135.2, 
                     exclude_reporter_region = FALSE, 
                     index_mgf_ms2 = FALSE, 
                     mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                     out_path = file.path(filepath, "mgf_queries.rds"), 
                     is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                     use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                     maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                     ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                     quant = "none", 
                     digits = 4L) 
{
  if (is_mdda) {
    warning("No multi-precursor DDA with MGF. Use mzML to enable the feature.")
    is_mdda <- FALSE
  }

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
  
  if (type_mgf == "default_pasef")
    mprepBrukerMGF(filepath)
  
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
    
    qs::qsave(list(data_format = data_format, mgf_format = mgf_format), 
              file.path(filepath, "info_format.rds"), preset = "fast")
  })

  ## Reads MGF into chunks
  # separate parallel process: 
  # (1) one large MGF file and parallel chunks
  # (2) parallel five MGF files and parallel chunks in each
  n_cores <- min(len <- length(filelist), detect_cores(32L))

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
                                is_mdda = is_mdda, 
                                use_defpeaks = use_defpeaks, 
                                deisotope_ms2 = deisotope_ms2, 
                                max_ms2_charge = max_ms2_charge, 
                                maxn_dia_precurs = maxn_dia_precurs, 
                                maxn_mdda_precurs = maxn_mdda_precurs, 
                                n_mdda_flanks = n_mdda_flanks, 
                                ppm_ms1_deisotope = ppm_ms1_deisotope, 
                                ppm_ms2_deisotope = ppm_ms2_deisotope, 
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
               filepath = filepath)
}


#' Post-processing of MGF or mzML
#' 
#' Calculates mass \code{frame}s etc.
#' 
#' @param df A data frame of processed peak lists.
#' @inheritParams readMGF
post_readmgf <- function (df, min_mass = 200L, max_mass = 4500L, ppm_ms1 = 10L, 
                          filepath) 
{
  if (is.atomic(df[1, "ms1_charge", drop = TRUE])) {
    df <- dplyr::arrange(df, ms1_mass)
    # df <- dplyr::filter(df, ms1_mass >= min_mass, ms1_mass <= max_mass)
  }

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
  
  df  <- split(df, df$raw_file)
  nms <- names(df)
  
  for (i in seq_along(df)) {
    qs::qsave(df[[i]], file.path(filepath, paste0("mgf_queries_", nms[i], ".rds")), 
              preset = "fast")
  }

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
                             deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                             maxn_dia_precurs = 300L, 
                             is_mdda = FALSE, use_defpeaks = FALSE, 
                             maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                             ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                             quant = "none", digits = 4L) 
{
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (!(len <- length(filelist))) 
    stop("No mgf files under ", filepath, call. = FALSE)
  
  if (len == 1L) {
    out <- proc_mgf_chunks(
      file.path(filepath, filelist),
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
      is_mdda = is_mdda, 
      deisotope_ms2 = deisotope_ms2, 
      max_ms2_charge = max_ms2_charge, 
      use_defpeaks = use_defpeaks, 
      maxn_dia_precurs = maxn_dia_precurs, 
      maxn_mdda_precurs = maxn_mdda_precurs, 
      n_mdda_flanks = n_mdda_flanks, 
      ppm_ms1_deisotope = ppm_ms1_deisotope, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      quant = quant, 
      digits = digits
    )
  }
  else {
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
      envir = environment(mzion::matchMS)
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
                                  is_mdda = is_mdda, 
                                  deisotope_ms2 = deisotope_ms2, 
                                  max_ms2_charge = max_ms2_charge, 
                                  use_defpeaks = use_defpeaks, 
                                  maxn_dia_precurs = maxn_dia_precurs, 
                                  maxn_mdda_precurs = maxn_mdda_precurs, 
                                  n_mdda_flanks = n_mdda_flanks, 
                                  ppm_ms1_deisotope = ppm_ms1_deisotope, 
                                  ppm_ms2_deisotope = ppm_ms2_deisotope, 
                                  quant = quant, 
                                  digits = digits)
    
    parallel::stopCluster(cl)
    
    out <- dplyr::bind_rows(out)
  }

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
                is_mdda = is_mdda, 
                deisotope_ms2 = deisotope_ms2, 
                max_ms2_charge = max_ms2_charge, 
                use_defpeaks = use_defpeaks, 
                maxn_dia_precurs = maxn_dia_precurs, 
                maxn_mdda_precurs = maxn_mdda_precurs, 
                n_mdda_flanks = n_mdda_flanks, 
                ppm_ms1_deisotope = ppm_ms1_deisotope, 
                ppm_ms2_deisotope = ppm_ms2_deisotope, 
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
                             index_mgf_ms2 = FALSE, is_mdda = FALSE, 
                             deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                             use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                             maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                             ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
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
                   is_mdda = is_mdda, 
                   deisotope_ms2 = deisotope_ms2, 
                   max_ms2_charge = max_ms2_charge, 
                   use_defpeaks = use_defpeaks, 
                   maxn_dia_precurs = maxn_dia_precurs, 
                   maxn_mdda_precurs = maxn_mdda_precurs, 
                   n_mdda_flanks = n_mdda_flanks, 
                   ppm_ms1_deisotope = ppm_ms1_deisotope, 
                   ppm_ms2_deisotope = ppm_ms2_deisotope, 
                   quant = quant, 
                   digits = digits)
}


#' Helper in processing MGF entries in chunks.
#'
#' @param lines Lines of MGF.
#' @inheritParams proc_mgf_chunks
proc_mgfs <- function (lines, topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
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
                       is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       quant = "none", digits = 4L) 
{
  options(digits = 9L)

  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends   <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))

  ## MS1 
  # (1) m-over-z and intensity
  ms1s <- stringi::stri_replace_first_fixed(lines[begins + n_to_pepmass], 
                                            "PEPMASS=", "")
  ms1s <- lapply(ms1s, stringi::stri_split_fixed, pattern = sep_pepmass, 
                 n = nfields_pepmass, simplify = TRUE)

  ms1_moverzs <- lapply(ms1s, function (x) as.numeric(x[, 1]))
  ms1_moverzs <- .Internal(unlist(ms1_moverzs, recursive = FALSE, use.names = FALSE))
  # not as.integer; intensity may be > .Machine$integer.max (2147483647)
  ms1_ints <- lapply(ms1s, function (x) as.numeric(x[, 2]))
  ms1_ints <- .Internal(unlist(ms1_ints, recursive = FALSE, use.names = FALSE))
  rm(list = c("ms1s"))

  # (2) retention time
  ret_times <- stringi::stri_replace_first_fixed(lines[begins + n_to_rt], "RTINSECONDS=", "")
  ret_times <- as.numeric(ret_times)

  # (3) MS1 charges and masses
  ms1_charges <- stringi::stri_replace_first_fixed(lines[begins + n_to_charge], "CHARGE=", "")
  ms1_charges <- gsub("[\\+]$", "", ms1_charges)
  # if ((polarity <- gsub(".*([\\+-])$", "\\1", ms1_charges[[1]])) == "+")
  #   ms1_charges <- gsub("[\\+]$", "", ms1_charges)
  # else if (polarity == "-")
  #   ms1_charges <- gsub("-$", "", ms1_charges)
  ms1_charges <- as.integer(ms1_charges)

  # timsTOF no CHARGE line -> NAs introduced by coercion

  ms1_masses <- mapply(function (x, y) x * y - y * 1.00727647, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)

  rows <- (ms1_charges >= ms1_charge_range[1] & ms1_charges <= ms1_charge_range[2] & 
             ret_times >= ret_range[1] & ret_times <= ret_range[2] & 
             ms1_masses >= min_mass & ms1_masses <= max_mass & 
             # timsTOF: no MS1 masses
             !is.na(ms1_masses))
  
  # timsTOF data may have undetermined charge states
  if (length(na_rows <- .Internal(which(is.na(rows))))) 
    rows[na_rows] <- FALSE

  begins <- begins[rows]
  ends <- ends[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- round(ret_times[rows], digits = 2L)
  ms1_masses <- ms1_masses[rows]

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
    raw_files <- rep_len(raw_file, length.out = length(begins))
    scan_nums <- stringi::stri_replace_first_fixed(lines[begins + n_to_scan], "RAWSCANS=", "")
  } 
  else {
    stop("Unknown MGF format.")
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
  
  ms2_charges <- vector("list", length(ms2_moverzs))
  
  is_tmt <- if (grepl("^tmt.*\\d+", quant)) TRUE else FALSE
  
  if (deisotope_ms2) {
    mics <- mapply(deisotope, ms2_moverzs, ms2_ints, 
                   MoreArgs = list(
                     exclude_reporter_region = is_tmt, 
                     tmt_reporter_lower = tmt_reporter_lower, 
                     tmt_reporter_upper = tmt_reporter_upper, 
                     ppm = ppm_ms2_deisotope, 
                     ms_lev = 2L, maxn_feats = topn_ms2ions, 
                     max_charge = max_ms2_charge, n_fwd = 10L, 
                     offset_upr = 30L, offset_lwr = 30L, order_mz = TRUE
                   ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ms2_moverzs <- lapply(mics, `[[`, "masses")
    ms2_ints <- lapply(mics, `[[`, "intensities")
    ms2_charges <- lapply(mics, `[[`, "charges")
  }

  # extract the TMT region of MS2 moverz and intensity
  # (also convert reporter-ion intensities to integers)
  restmt <- extract_mgf_rptrs(ms2_moverzs, 
                              ms2_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  
  ms2_moverzs <- restmt[["xvals"]]
  ms2_ints <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  rm(list = "restmt")

  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = ms2_moverzs, 
                          ms2_ints = ms2_ints, 
                          ms2_charges = ms2_charges, 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  
  ms2_moverzs <- mz_n_int[["ms2_moverzs"]]
  ms2_ints <- mz_n_int[["ms2_ints"]]
  ms2_charges <- mz_n_int[["ms2_charges"]]
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

  df <- tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_files,
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_int = ms1_ints,
    ms1_charge = ms1_charges,
    ret_time = ret_times,
    scan_num = scan_nums,
    ms2_moverzs= ms2_moverzs,
    ms2_ints = ms2_ints,
    ms2_charges = ms2_charges, 
    ms2_n = lens, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints, )
}


#' Subsets MGFs by top-n entries.
#' 
#' \code{lens} after filtered by \code{min_ms2mass} but before subset by 
#' \code{topn_ms2ions} to reflect noise levels.
#' 
#' @param ms2_moverzs Lists of MS2 moverz values.
#' @param ms2_ints Lists of MS2 intensities.
#' @param ms2_charges Lists of MS2 charges.
#' @inheritParams load_mgfs
sub_mgftopn <- function (ms2_moverzs = NULL, ms2_ints = NULL, ms2_charges = NULL, 
                         topn_ms2ions = 100L, mgf_cutmzs = numeric(), 
                         mgf_cutpercs = numeric(), min_ms2mass = 115L, 
                         max_ms2mass = 4500L) 
{
  options(digits = 9L)
  
  ## subsets by min_ms2mass
  oks <- lapply(ms2_moverzs, function (x) x >= min_ms2mass)
  
  ms2_moverzs <- mapply(function (x, y) x[y], ms2_moverzs, oks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ms2_ints <- mapply(function (x, y) x[y], ms2_ints, oks, 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ms2_charges <- mapply(function (x, y) x[y], ms2_charges, oks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)

  ## subsets by topn
  lens <- lengths(ms2_moverzs)
  
  if (topn_ms2ions < Inf) {
    is_long <- lens > topn_ms2ions
    
    if (length(mgf_cutmzs)) {
      m_long <- ms2_moverzs[is_long]
      i_long <- ms2_ints[is_long]
      z_long <- ms2_charges[is_long]
      
      for (i in seq_along(m_long)) {
        x <- m_long[[i]]
        y <- i_long[[i]]
        z <- z_long[[i]]
        
        # (`<` not `<=`)
        ok_ms2 <- x < max_ms2mass
        x <- x[ok_ms2]
        y <- y[ok_ms2]
        z <- z[ok_ms2]
        
        idxes <- findInterval(x, mgf_cutmzs)
        xs <- split(x, idxes)
        ys <- split(y, idxes)
        zs <- split(z, idxes)
        
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
        
        rows  <- mapply(which_topx2, ys, ok_percs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        ans_x <- mapply(function (x, y) x[y], xs, rows, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        ans_y <- mapply(function (x, y) x[y], ys, rows, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        ans_z <- mapply(function (x, y) x[y], zs, rows, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        
        ans_x <- .Internal(unlist(ans_x, recursive = FALSE, use.names = FALSE))
        ans_y <- .Internal(unlist(ans_y, recursive = FALSE, use.names = FALSE))
        ans_z <- .Internal(unlist(ans_z, recursive = FALSE, use.names = FALSE))
        
        ms2_moverzs[is_long][[i]] <- ans_x
        ms2_ints[is_long][[i]] <- ans_y
        ms2_charges[is_long][[i]] <- ans_z
        
        # if (i %% 5000L == 0L) gc()
      }
    }
    else {
      rows <- lapply(ms2_ints[is_long], which_topx2, topn_ms2ions)
      
      ms2_ints[is_long] <- mapply(function (x, y) x[y], 
                                  ms2_ints[is_long], rows, 
                                  SIMPLIFY = FALSE, USE.NAMES = FALSE)
      ms2_moverzs[is_long] <- mapply(function (x, y) x[y], 
                                     ms2_moverzs[is_long], rows, 
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
      ms2_charges[is_long] <- mapply(function (x, y) x[y], 
                                     ms2_charges[is_long], rows, 
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
  }
  
  # also handles MS2 intensity max-outs, which usually don't happen
  ms2_ints <- integerize_ms2ints(ms2_ints)

  list(ms2_moverzs = ms2_moverzs, ms2_ints = ms2_ints, 
       ms2_charges = ms2_charges, lens = lens)
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
#' @param xvals Lists of MS2 m-over-z values.
#' @param yvals Lists of MS2 intensity values. 
#' @inheritParams matchMS
extract_mgf_rptrs <- function (xvals, yvals, quant = "none", 
                               tmt_reporter_lower = 126.1, 
                               tmt_reporter_upper = 135.2, 
                               exclude_reporter_region = FALSE) 
{
  if (grepl("^tmt.*\\d+", quant)) {
    ok_rptrs <- lapply(xvals, function (x) x > tmt_reporter_lower & 
                         x < tmt_reporter_upper)
    rptr_moverzs <- mapply(function (x, y) x[y], xvals, ok_rptrs, 
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
    rptr_ints <- mapply(function (x, y) x[y], yvals, ok_rptrs, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    rptr_ints <- integerize_ms2ints(rptr_ints)
    
    if (exclude_reporter_region) {
      no_rptrs <- lapply(ok_rptrs, `!`)
      xvals <- mapply(function (x, y) x[y], xvals, no_rptrs, 
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
      yvals <- mapply(function (x, y) x[y], yvals, no_rptrs, 
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
  }
  else {
    rptr_moverzs <- NA_real_
    rptr_ints <- NA_integer_
  }
  
  list(xvals = xvals, 
       yvals = yvals,  
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
find_ms1_interval <- function (mass = 1800.0, from = 200L, ppm = 10L) 
{
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+ppm/1e6))
}


#' Converts ms2_moverzs to integers.
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
readmzML <- function (filepath = NULL, filelist = NULL, out_path = NULL, 
                      min_mass = 200L, max_mass = 4500L, 
                      min_ms2mass = 115L, max_ms2mass = 4500L, 
                      topn_ms2ions = 100L, ms1_charge_range = c(2L, 6L), 
                      ms1_scan_range = c(1L, .Machine$integer.max), 
                      ret_range = c(0, Inf), ppm_ms1 = 10L, ppm_ms2 = 10L, 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, 
                      index_mgf_ms2 = FALSE, 
                      mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                      is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                      use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                      maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                      ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                      quant = "none", digits = 4L)
{
  out <- vector("list", len <- length(filelist))
  
  files <- file.path(filepath, filelist)
  sizes <- max(unlist(lapply(files, file.size)))/1024^3
  n_cores <- min(detect_cores(32L), floor((find_free_mem()/1024)/(sizes * 8)) + 1L, len)
  n_cores <- max(1L, n_cores)
  
  if (n_cores == 1L) {
    for (i in 1:len)
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
                            is_mdda = is_mdda, 
                            deisotope_ms2 = deisotope_ms2, 
                            max_ms2_charge = max_ms2_charge, 
                            use_defpeaks = use_defpeaks, 
                            maxn_dia_precurs = maxn_dia_precurs, 
                            maxn_mdda_precurs = maxn_mdda_precurs, 
                            n_mdda_flanks = n_mdda_flanks, 
                            ppm_ms1_deisotope = ppm_ms1_deisotope, 
                            ppm_ms2_deisotope = ppm_ms2_deisotope, 
                            quant = quant, 
                            digits = digits)
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
                                  is_mdda = is_mdda, 
                                  deisotope_ms2 = deisotope_ms2, 
                                  max_ms2_charge = max_ms2_charge, 
                                  use_defpeaks = use_defpeaks, 
                                  maxn_dia_precurs = maxn_dia_precurs, 
                                  maxn_mdda_precurs = maxn_mdda_precurs, 
                                  n_mdda_flanks = n_mdda_flanks, 
                                  ppm_ms1_deisotope = ppm_ms1_deisotope, 
                                  ppm_ms2_deisotope = ppm_ms2_deisotope, 
                                  quant = quant, 
                                  digits = digits)
    parallel::stopCluster(cl)
  }
  
  out <- dplyr::bind_rows(out)

  post_readmgf(out, min_mass = min_mass, max_mass = max_mass, ppm_ms1 = ppm_ms1, 
               filepath = filepath)
}


#' Helper of \link{readmzML}
#' 
#' No scan range subsetting with PASEF timsTOF.
#' 
#' @param file A file name to mzML with a prepending path.
#' @inheritParams readmzML
proc_mzml <- function (file, topn_ms2ions = 100L, ms1_charge_range = c(2L, 4L), 
                       ret_range = c(0, Inf), min_mass = 200L, max_mass = 4500L, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       quant = "none", digits = 4L) 
{
  min_ms1_charge <- ms1_charge_range[1]
  max_ms1_charge <- ms1_charge_range[2]
  
  df <- read_mzml(file, topn_ms2ions = topn_ms2ions, 
                  tmt_reporter_lower = tmt_reporter_lower, 
                  tmt_reporter_upper = tmt_reporter_upper, 
                  exclude_reporter_region = exclude_reporter_region, 
                  index_mgf_ms2 = index_mgf_ms2, ppm_ms1 = ppm_ms1, 
                  ppm_ms2 = ppm_ms2, min_ms2mass = min_ms2mass, 
                  max_ms2mass = max_ms2mass, max_ms1_charge = max_ms1_charge, 
                  is_mdda = is_mdda, 
                  deisotope_ms2 = deisotope_ms2, 
                  max_ms2_charge = max_ms2_charge, 
                  use_defpeaks = use_defpeaks, 
                  maxn_dia_precurs = maxn_dia_precurs, 
                  maxn_mdda_precurs = maxn_mdda_precurs, 
                  n_mdda_flanks = n_mdda_flanks, 
                  ppm_ms1_deisotope = ppm_ms1_deisotope, 
                  ppm_ms2_deisotope = ppm_ms2_deisotope, 
                  quant = quant, digits = digits)
  
  if (is_mdda && maxn_mdda_precurs == 1L) {
    # empties <- lapply(df$ms1_moverz, is.null)
    # empties <- unlist(empties, recursive = FALSE, use.names = FALSE)
    # df <- df[!empties, ]
    df$ms1_moverz <- unlist(df$ms1_moverz, recursive = TRUE, use.names = FALSE)
    df$ms1_mass <- unlist(df$ms1_mass, recursive = TRUE, use.names = FALSE)
    df$ms1_charge <- unlist(df$ms1_charge, recursive = TRUE, use.names = FALSE)
    df$ms1_int <- unlist(df$ms1_int, recursive = TRUE, use.names = FALSE)
  }
  
  if (is.atomic(df[1, "ms1_charge", drop = TRUE])) {
    df <- df[with(df, !is.na(ms1_mass)), ]
    
    df <- dplyr::filter(df, 
                        ms1_charge >= min_ms1_charge, 
                        ms1_charge <= max_ms1_charge, 
                        ret_time >= ret_range[1], ret_time <= ret_range[2], 
                        ms1_mass >= min_mass, ms1_mass <= max_mass, )
  }

  # subsets by top-n and min_ms2mass
  # (also convert non reporter-ion MS2 intensities to integers)
  mz_n_int <- sub_mgftopn(ms2_moverzs = df[["ms2_moverzs"]], 
                          ms2_ints = df[["ms2_ints"]], 
                          ms2_charges = df[["ms2_charges"]], 
                          topn_ms2ions = topn_ms2ions, 
                          mgf_cutmzs = mgf_cutmzs, 
                          mgf_cutpercs = mgf_cutpercs, 
                          min_ms2mass = min_ms2mass, 
                          max_ms2mass = max_ms2mass)
  
  # can be integers if "index_mgf_ms2 = TRUE"
  df[["ms2_moverzs"]] <- mz_n_int[["ms2_moverzs"]]
  df[["ms2_ints"]] <- mz_n_int[["ms2_ints"]]
  df[["ms2_charges"]] <- mz_n_int[["ms2_charges"]]
  df[["ms2_n"]] <- mz_n_int[["lens"]]

  invisible(df)
}


#' Reads mzML from MSConvert
#'
#' Only for Thermo's RAW format.
#'
#' @param xml_file A file name of mzML.
#' @inheritParams matchMS
read_mzml <- function (xml_file, topn_ms2ions = 100L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, min_ms2mass = 115L, 
                       max_ms2mass = 4500L, max_ms1_charge = 4L, 
                       is_mdda = FALSE, deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                       maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                       quant = "none", digits = 4L)
{
  ## spectrum
  xml_root <- xml2::read_xml(xml_file)
  mzML <- xml2::xml_child(xml_root)
  
  raw_file <- local({
    idx_file <- which(xml2::xml_name(xml2::xml_children(mzML)) == "fileDescription")
    file_des <- xml2::xml_children(mzML)[[idx_file]]
    idx_srcl <- which(xml2::xml_name(xml2::xml_children(file_des)) == "sourceFileList")
    info_raw <- xml2::xml_children(file_des)[[idx_srcl]]
    idx_srcf <- which(xml2::xml_name(xml2::xml_children(info_raw)) == "sourceFile")
    info_fi  <- xml2::xml_children(info_raw)[[idx_srcf]]
    raw_file <- xml2::xml_attr(info_fi, "name")
  })
  
  idx_run <- which(xml2::xml_name(xml2::xml_children(mzML)) == "run")
  run <- xml2::xml_children(mzML)[[idx_run]]
  idx_specs <- which(xml2::xml_name(xml2::xml_children(run)) == "spectrumList")
  spec <- xml2::xml_children(xml2::xml_children(run)[[idx_specs]])
  rm(list = c("mzML", "idx_run", "run", "idx_specs"))
  
  # - mzML
  #   - ...
  #   - fileDescription
  #     - ...
  #     - sourceFileList
  #       - sourceFile
  
  # - run
  #  - spectrumList
  #   - spectrum # x
  #     - ...
  #     - <cvParam ... name="ms level" value="2"/>
  #     - <cvParam ... name="total ion current" value="2578.8672"/>
  #     - <cvParam ... name="spectrum title" value="23aug2017_hela_serum_timecourse_4mz_narrow_1.1.1. 
  #                    File:"FILENAME.raw", ...originalScan=1 demux=0 scan=1/>
  #     - scanList # xc
  #       - ...
  #       - <scan ... originalScan=3 demux=1 scan=4> # originalScan=3 demux=1 scan=4 only only for MS2
  #         - <cvParam ... name="scan start time" value="0.012728725" .../>
  #         - <cvParam ... name="ion injection time" value="54.999999701977" .../>
  #     - precursorList # xc
  #       - <precursor ... originalScan=3 demux=1 scan=4"> # redundant?
  #         - <cvParam ... name="isolation window target m/z" value="401.432342529297" .../>
  #         - selectedIonList
  #           - selectedIon
  #             -  <cvParam ... name="selected ion m/z" value="401.432342529297" .../>
  #             -  <cvParam ... name="charge state" value="3" .../> # No charge state if DIA
  #             -  <cvParam ... name="peak intensity" value="0" .../> # 0 if DIA
  #     - binaryDataArrayList # xc
  #       - binaryDataArray
  #         - binary
  #       - binaryDataArray
  #         - binary
  #  - chromatogramList
  #   - chromatogram
  #    - binaryDataArrayList
  #     - binaryDataArray # name="time array"
  #     - binaryDataArray # name="intensity array"
  #     - binaryDataArray # name="non-standard data array"
  
  ## the first scan
  x <- spec[[1]]
  ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))[[1]]
  ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))
  id_nms <- lapply(ids, `[[`, 1)
  idx_demux <- which(id_nms == "demux")
  is_demux <- if (length(idx_demux)) TRUE else FALSE
  idx_sc <- which(id_nms == "scan")
  if (!length(idx_osc <- which(id_nms == "originalScan"))) idx_osc <- idx_sc
  
  xc <- xml2::xml_children(x)
  xcp_attrs <- xml2::xml_attrs(xc[which(xml2::xml_name(xc) == "cvParam")])
  xcp_names <- lapply(xcp_attrs, `[[`, "name")
  xcp_vals <- lapply(xcp_attrs, `[[`, "value")
  idx_mslev <- which(xcp_names == "ms level")
  idx_title <- which(xcp_names == "spectrum title")
  rm(list = c("x", "ids", "xc", "xcp_attrs", "xcp_names", "xcp_vals", "id_nms"))
  
  # MS2 indexes (no MS1 with DDA -> MS2 goes first)
  for (i in 1:(len <- length(spec))) {
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (xml2::xml_attr(xc[[idx_mslev]], "value") == "2") {
      idx_precursor_2 <- grep("precursorList", xc)
      idx_scanList_2 <- grep("scanList", xc)
      idx_bin_2 <- grep("binaryDataArrayList", xc)
      
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      idx_rt_2 <- which(xml2::xml_name(scanList) == "scan")
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      scanList_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      idx_scan_start_2 <- which(scanList_ret_attrs == "scan start time")
      ms2_reso <- scanList_ret[[which(scanList_ret_attrs == "mass resolving power")]]
      ms2_reso <- as.integer(xml2::xml_attr(ms2_reso, "value"))
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # assume one precursor
      precursorc <- xml2::xml_children(precursor)
      idx_selectedIonList <- grep("selectedIonList", precursorc)
      
      idx_isolationWindow <- grep("isolationWindow", precursorc)
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_nms <- lapply(isolationWindowc, function (x) xml2::xml_attr(x, "name"))
      idx_ctrmz <- which(iso_nms == "isolation window target m/z")
      idx_lwrmz <- which(iso_nms == "isolation window lower offset")
      idx_uprmz <- which(iso_nms == "isolation window upper offset")
      
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      selion_nms <- lapply(selectedIonc, function (x) xml2::xml_attr(x, "name"))
      idx_moverz <- which(selion_nms == "selected ion m/z")
      idx_ms1int <- which(selion_nms == "peak intensity") # DIA: zero intensity
      idx_charge <- which(selion_nms == "charge state") # DIA: no "charge state"
      is_dia <- if (length(idx_charge)) FALSE else TRUE
      
      rm(list = c("scanList", "scanList_ret", "scanList_ret_attrs", 
                  "precursorList", "precursor", "precursorc", 
                  "selectedIon", "selectedIonc", "selion_nms"))
      break
    }
  }
  
  # MS1 indexes
  for (i in seq_along(spec)) {
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (as.integer(xml2::xml_attr(xc[[idx_mslev]], "value")) == 1) {
      idx_scanList_1 <- grep("scanList", xc)
      idx_bin_1 <- grep("binaryDataArrayList", xc)
      
      local({
        binData <- xml2::xml_children(xc[[idx_bin_1]])[[1]]
        binDataC <- xml2::xml_children(binData)
        binDataC <- binDataC[which(xml2::xml_name(binDataC) == "cvParam")]
        
        lapply(binDataC, function (v) {
          attrs <- xml2::xml_attrs(v)
          oks <- grepl("compression", attrs)
          
          if (any(oks))
            if (attrs[oks][[1]] == "zlib compression")
              stop("Please uncheck zlib compression in  mzML generations.")
        })
      })
      
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      idx_rt_1 <- which(xml2::xml_name(scanList) == "scan")
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      
      scan_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      # later remove the assumption of a fixed resolution with Astral
      # ms1_reso <- scanList_ret[[which(scan_ret_attrs == "mass resolving power")]]
      # ms1_reso <- as.integer(xml2::xml_attr(ms1_reso, "value"))
      idx_scan_start_1 <- which(scan_ret_attrs == "scan start time")
      
      rm(list = c("scanList", "scanList_ret", "scan_ret_attrs"))
      break
    }
  }
  
  rm(list = c("x", "xc", "i"))
  
  # msx_: both ms1 and ms2, differentiated by ms_lev
  # ms1_: ms1 only
  # ms0_: by other peak-pickings, e.g., MSConvert
  
  if (is_mdda) {
    df <- proc_mdda(spec, raw_file = raw_file, idx_sc = idx_sc, idx_osc = idx_osc, 
                    idx_mslev = idx_mslev, idx_title = idx_title, 
                    idx_scanList_2 = idx_scanList_2, idx_rt_2 = idx_rt_2, 
                    idx_scan_start_2 = idx_scan_start_2, 
                    idx_precursor_2 = idx_precursor_2, 
                    idx_isolationWindow = idx_isolationWindow, 
                    idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz,
                    idx_uprmz = idx_uprmz, 
                    idx_selectedIonList = idx_selectedIonList, 
                    idx_moverz = idx_moverz, idx_charge = idx_charge, 
                    idx_ms1int = idx_ms1int, idx_bin_2 = idx_bin_2, 
                    maxn_precurs = maxn_mdda_precurs, 
                    max_ms1_charge = max_ms1_charge, 
                    topn_ms2ions = topn_ms2ions, quant = quant, 
                    tmt_reporter_lower = tmt_reporter_lower, 
                    tmt_reporter_upper = tmt_reporter_upper, 
                    exclude_reporter_region = exclude_reporter_region, 
                    idx_scanList_1 = idx_scanList_1, 
                    idx_rt_1 = idx_rt_1, idx_scan_start_1 = idx_scan_start_1, 
                    idx_bin_1 = idx_bin_1, deisotope_ms2 = deisotope_ms2, 
                    n_mdda_flanks = n_mdda_flanks, 
                    ppm_ms1_deisotope = ppm_ms1_deisotope, 
                    ppm_ms2_deisotope = ppm_ms2_deisotope, 
                    use_defpeaks = use_defpeaks)
  }
  else if (is_dia) {
    df <- proc_dia(spec, raw_file = raw_file, is_demux = is_demux, 
                   idx_sc = idx_sc, idx_osc = idx_osc, idx_mslev = idx_mslev, 
                   idx_title = idx_title, idx_scanList_1 = idx_scanList_1, 
                   idx_scanList_2 = idx_scanList_2, idx_rt_1 = idx_rt_1, 
                   idx_rt_2 = idx_rt_2, idx_scan_start_1 = idx_scan_start_1, 
                   idx_scan_start_2 = idx_scan_start_2, 
                   idx_precursor_2 = idx_precursor_2, 
                   idx_isolationWindow = idx_isolationWindow, 
                   idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz, 
                   idx_uprmz = idx_uprmz, idx_selectedIonList = idx_selectedIonList, 
                   idx_demux = idx_demux, idx_bin_1 = idx_bin_1, 
                   idx_bin_2 = idx_bin_2, 
                   maxn_dia_precurs = maxn_dia_precurs, 
                   max_ms1_charge = max_ms1_charge, 
                   ppm_ms1_deisotope = ppm_ms1_deisotope, 
                   ppm_ms2_deisotope = ppm_ms2_deisotope, 
                   deisotope_ms2 = deisotope_ms2, 
                   max_ms2_charge = max_ms2_charge, 
                   topn_ms2ions = topn_ms2ions, quant = quant, 
                   tmt_reporter_lower = tmt_reporter_lower, 
                   tmt_reporter_upper = tmt_reporter_upper, 
                   exclude_reporter_region = exclude_reporter_region)
    
    
  }
  else {
    df <- proc_dda(spec, raw_file = raw_file, idx_sc = idx_sc, idx_osc = idx_osc, 
                   idx_mslev = idx_mslev, idx_title = idx_title, 
                   idx_scanList_2 = idx_scanList_2, idx_rt_2 = idx_rt_2, 
                   idx_scan_start_2 = idx_scan_start_2, 
                   idx_precursor_2 = idx_precursor_2, 
                   idx_isolationWindow = idx_isolationWindow, 
                   idx_selectedIonList = idx_selectedIonList, 
                   idx_moverz = idx_moverz, idx_charge = idx_charge, 
                   idx_ms1int = idx_ms1int, idx_bin_2 = idx_bin_2, 
                   deisotope_ms2 = deisotope_ms2, max_ms2_charge = max_ms2_charge, 
                   ppm_ms2_deisotope = ppm_ms2_deisotope, 
                   topn_ms2ions = topn_ms2ions, quant = quant, 
                   tmt_reporter_lower = tmt_reporter_lower, 
                   tmt_reporter_upper = tmt_reporter_upper, 
                   exclude_reporter_region = exclude_reporter_region, 
                   index_mgf_ms2 = index_mgf_ms2)
  }
}


#' Proc mzML data for mDDA workflows.
#' 
#' @param spec Spectrum list.
#' @param idx_sc Index of scan numbers.
#' @param idx_osc Index of original scan numbers.
#' @param idx_mslev Index of MS levels.
#' @param idx_title Index of scan titles.
#' @param idx_scanList_2 Index of MS2 scanList.
#' @param idx_rt_2 Index of MS2 retention times.
#' @param idx_scan_start_2 Index of MS2 scan start.
#' @param idx_precursor_2 Index of MS2 precursor.
#' @param idx_isolationWindow Index of isolationWindow.
#' @param idx_ctrmz Index of center m-over-z.
#' @param idx_lwrmz Index of lower m-over-z.
#' @param idx_uprmz Index of upper m-over-z.
#' @param idx_selectedIonList Index of selectedIonList.
#' @param idx_moverz Index of m-over-z sequences.
#' @param idx_charge Index of charge states.
#' @param idx_ms1int Index of MS1 intensities.
#' @param idx_bin_2 Index of MS2 binary data.
#' @param maxn_precurs Maximum number of chimeric precursors.
#' @param max_ms1_charge Maximum precursor charge states.
#' @param topn_ms2ions Top-n MS2 ions for consideration.
#' @param idx_scanList_1 Index of MS1 scanList.
#' @param idx_rt_1 Index of precursor retention times.
#' @param idx_scan_start_1 Index of MS1 scan starts.
#' @param idx_bin_1 Index of MS1 binary data.
#' @param deisotope_ms2 Logical; to deisotope MS2 features or not.
#' @param raw_file The RAW file name of an mzML.
#' @inheritParams matchMS
proc_mdda <- function (spec, raw_file, idx_sc = 3L, idx_osc = 3L, idx_mslev = 2L, 
                       idx_title = 10L, idx_scanList_2 = 11L, idx_rt_2 = 2L, 
                       idx_scan_start_2 = 1L, idx_precursor_2 = 12L, 
                       idx_isolationWindow = 1L, idx_ctrmz = 1L, idx_lwrmz = 2L,
                       idx_uprmz = 3L, idx_selectedIonList = 2L, 
                       idx_moverz = 1L, idx_charge = 2L, idx_ms1int = 3L, 
                       idx_bin_2 = 13L, maxn_precurs = 5L, max_ms1_charge = 4L, 
                       topn_ms2ions = 150L, quant = "none", 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       idx_scanList_1 = 11L, idx_rt_1 = 2L, 
                       idx_scan_start_1 = 1L, idx_bin_1 = 12L, 
                       deisotope_ms2 = TRUE, ppm_ms1_deisotope = 10L, 
                       ppm_ms2_deisotope = 10L, n_mdda_flanks = 6L, 
                       use_defpeaks = FALSE)
{
  len <- length(spec)
  
  ret_times <- orig_scans <- scan_nums <- scan_titles <- 
    iso_ctr <- iso_lwr <- iso_upr <- character(len)
  
  # msx_: both ms1 and ms2, differentiated by ms_lev
  # ms1_: ms1 only
  # ms0_: by other peak-pickings, e.g., MSConvert
  
  ms_levs <- msx_ns <- integer(len)
  msx_moverzs <- msx_ints <- ms2_charges <- vector("list", len)
  ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", len)
  ms0_moverzs <- ms0_ints <- ms0_charges <- character(len)
  
  for (i in 1:len) {
    x <- spec[[i]]
    ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))[[1]]
    ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))
    scan_nums[[i]]  <- ids[[idx_sc]][[2]]
    orig_scans[[i]] <- ids[[idx_osc]][[2]]

    xc <- xml2::xml_children(x)
    ms_levs[[i]] <- ms_lev <- as.integer(xml2::xml_attr(xc[[idx_mslev]], "value"))
    scan_titles[[i]] <- xml2::xml_attr(xc[[idx_title]], "value")
    
    is_tmt <- if (grepl("^tmt.*\\d+", quant)) TRUE else FALSE
    
    if (ms_lev == 2L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_2]], "value")
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # (assume one precursor, not yet chimeric)
      precursorc <- xml2::xml_children(precursor)
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_ctr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_ctrmz]], "value")
      iso_lwr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value")
      iso_upr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value")
      
      selectedIon  <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      
      ms0_moverzs[[i]] <- xml2::xml_attr(selectedIonc[[idx_moverz]], "value")
      ms0_charges[[i]] <- xml2::xml_attr(selectedIonc[[idx_charge]], "value")
      # no precursor intensity: 01CPTAC3_Benchmarking_P_BI_20170523_BL_f12.mzML; scan=54837
      if (length(selectedIonc) >= idx_ms1int)
        ms0_ints[[i]] <- xml2::xml_attr(selectedIonc[[idx_ms1int]], "value")
      else
        ms0_ints[[i]] <- character(1)
      
      ## binaryDataArrayList
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_2]]))
      msData  <- xml2::xml_contents(binData)
      
      if (length(msData) == 2L) {
        r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
        r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
        msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
        msx_xs <- readBin(r1, "double", n = msx_n, size = 8L)
        msx_ys <- readBin(r2, "double", n = msx_n, size = 8L)
        
        if (deisotope_ms2) {
          mic <- deisotope(msx_xs, msx_ys, exclude_reporter_region = is_tmt, 
                           tmt_reporter_lower = tmt_reporter_lower, 
                           tmt_reporter_upper = tmt_reporter_upper, 
                           ppm = ppm_ms2_deisotope, ms_lev = ms_lev, 
                           maxn_feats = topn_ms2ions, max_charge = 3L, 
                           n_fwd = 10L, offset_upr = 30L, offset_lwr = 30L, 
                           order_mz = TRUE)
          msx_moverzs[[i]] <- mic[["masses"]]
          msx_ints[[i]] <- mic[["intensities"]]
          ms2_charges[[i]] <- mic[["charges"]]
        }
        else {
          msx_moverzs[[i]] <- msx_xs
          msx_ints[[i]] <- msx_ys
        }
      }
      else {
        msx_ns[[i]] <- 0L
        msx_moverzs[[i]] <- NULL
        msx_ints[[i]] <- NULL
      }
    }
    else if (ms_lev == 1L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_1]], "value")
      
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_1]]))
      msData <- xml2::xml_contents(binData)
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_moverzs[[i]] <- msx_xs <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ints[[i]] <- msx_ys <- readBin(r2, "double", n = msx_n, size = 8L)
    }
  }
  
  restmt <- extract_mgf_rptrs(msx_moverzs, 
                              msx_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  msx_moverzs <- restmt[["xvals"]]
  msx_ints <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  
  df <- tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_file,
    ms_level = ms_levs, 
    ms1_moverzs = as.numeric(ms0_moverzs), 
    ms1_masses = list(NULL),
    ms1_charges = as.integer(ms0_charges), 
    ms1_ints = as.numeric(ms0_ints), 
    ret_time = as.numeric(ret_times) * 60, 
    scan_num = as.integer(scan_nums), 
    orig_scan = orig_scans,
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    ms2_charges = ms2_charges, 
    msx_ns = msx_ns, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints, 
    iso_ctr = as.numeric(iso_ctr), 
    iso_lwr = as.numeric(iso_lwr), 
    iso_upr = as.numeric(iso_upr), 
  )
  
  if (use_defpeaks) {
    df$ms1_moverzs <- as.list(df$ms1_moverzs)
    df$ms1_charges <- as.list(df$ms1_charges)
    df$ms1_ints <- as.list(df$ms1_ints)
  } 
  else {
    df$ms1_ints <- df$ms1_charges <- df$ms1_moverzs <- list(NULL)
  }

  idxes_ms1 <- which(df$ms_level == 1L)
  diff_ms1 <- c(0, diff(idxes_ms1))
  oks <- which(diff_ms1 > 1) # non-consecutive MS1s
  ms1_stas <- idxes_ms1[oks - 1L]
  ms2_stas <- ms1_stas + 1L
  ms2_ends <- idxes_ms1[oks] - 1L
  len <- length(ms1_stas)
  rm(list = c("oks", "diff_ms1", "idxes_ms1"))
  
  cols <- match(c("ms1_moverzs", "ms1_masses", "ms1_charges", "ms1_ints"), 
                names(df))

  for (i in 1:len) {
    stacr <- ms1_stas[i]
    stas1 <- ms1_stas[max(1L, i - n_mdda_flanks):min(len, i + n_mdda_flanks)]
    stas2 <- ms2_stas[i]
    ends2 <- ms2_ends[i]
    df[stas2:ends2, cols] <- 
      find_mdda_mms1s(df1 = df[stas1, ], df2 = df[stas2:ends2, ], 
                      stas1 = stas1, stas2 = stas2, ends2 = ends2, 
                      ppm = ppm_ms1_deisotope, maxn_precurs = maxn_precurs, 
                      max_ms1_charge = max_ms1_charge, n_fwd = 20L)
  }
  rm(list = c("stacr", "stas1", "stas2", "ends2", "ms1_stas", "ms2_stas", 
              "ms2_ends", "len"))
  
  df$iso_lwr <- df$iso_upr <- df$iso_ctr <- NULL
  df <- df[with(df, ms_level != 1L), ]
  bads <- unlist(lapply(df$ms1_moverzs, is.null))
  df <- df[!bads, ]

  df <- dplyr::rename(
    df, 
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_charge = ms1_charges,
    ms1_int = ms1_ints,
    ms2_moverzs = msx_moverzs,
    ms2_ints = msx_ints,
    ms2_n = msx_ns,
  )
}


#' Proc mzML data for DIA workflows.
#' 
#' @param is_demux Is demux DIA or not.
#' @param idx_demux Index of demux.
#' @inheritParams proc_mdda
#' @inheritParams matchMS
proc_dia <- function (spec, raw_file, is_demux = FALSE, idx_sc = 5L, idx_osc = 3L, 
                      idx_mslev = 2L, idx_title = 10L, idx_scanList_1 = 11L, 
                      idx_scanList_2 = 11L, idx_rt_1 = 2L, idx_rt_2 = 2L, 
                      idx_scan_start_1 = 1L, idx_scan_start_2 = 1L, 
                      idx_precursor_2 = 12L, idx_isolationWindow = 1L, 
                      idx_ctrmz = 1L, idx_lwrmz = 2L, idx_uprmz = 3L, 
                      idx_selectedIonList = 2L, idx_demux = 4L, 
                      idx_bin_1 = 12L, idx_bin_2 = 13L, 
                      maxn_dia_precurs = 300L, max_ms1_charge = 4L, 
                      ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                      deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                      topn_ms2ions = 100L, quant = "none", 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE)
{
  len <- length(spec)
  
  ret_times <- orig_scans <- scan_nums <- scan_titles <- 
    iso_ctr <- iso_lwr <- iso_upr <- character(len)
  
  ms_levs <- msx_ns <- integer(len)
  msx_moverzs <- msx_ints <- ms2_charges <- vector("list", len)
  ms1_moverzs <- ms1_ints <- ms1_charges <- vector("list", len)
  demux <- if (is_demux) rep("0", len) else NULL
  
  is_tmt <- if (grepl("^tmt.*\\d+", quant)) TRUE else FALSE
  
  for (i in 1:len) {
    x <- spec[[i]]
    ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))[[1]]
    ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))
    scan_nums[[i]] <- ids[[idx_sc]][[2]]
    orig_scans[[i]] <- ids[[idx_osc]][[2]]
    
    xc <- xml2::xml_children(x)
    ms_levs[[i]] <- ms_lev <- as.integer(xml2::xml_attr(xc[[idx_mslev]], "value"))
    scan_titles[[i]] <- xml2::xml_attr(xc[[idx_title]], "value")
    
    if (ms_lev == 2L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_2]], "value")
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # (assume one precursor, not yet chimeric)
      precursorc <- xml2::xml_children(precursor)
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_ctr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_ctrmz]], "value")
      iso_lwr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value")
      iso_upr[[i]] <- xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value")
      
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      
      # TBA: DIA ms1 moverzs, charges and intensity
      if (is_demux)
        demux[[i]] <- ids[[idx_demux]][[2]]
      
      ## binaryDataArrayList
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_2]]))
      msData <- xml2::xml_contents(binData)
      
      if (length(msData) == 2L) {
        r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
        r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
        msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
        msx_xs <- readBin(r1, "double", n = msx_n, size = 8L)
        msx_ys <- readBin(r2, "double", n = msx_n, size = 8L)
        
        if (deisotope_ms2) {
          mic <- deisotope(msx_xs, msx_ys, exclude_reporter_region = is_tmt, 
                           tmt_reporter_lower = tmt_reporter_lower, 
                           tmt_reporter_upper = tmt_reporter_upper, 
                           ppm = 10L, ms_lev = ms_lev, maxn_feats = topn_ms2ions, 
                           max_charge = max_ms2_charge, n_fwd = 10L, 
                           offset_upr = 30L, offset_lwr = 30L, order_mz = TRUE)
          msx_moverzs[[i]] <- mic[["masses"]]
          msx_ints[[i]] <- mic[["intensities"]]
          ms2_charges[[i]] <- mic[["charges"]]
        }
        else {
          msx_moverzs[[i]] <- msx_xs
          msx_ints[[i]] <- msx_ys
        }
      } 
      else {
        msx_ns[[i]] <- 0L
        # msx_moverzs[[i]] <- NA_real_
        # msx_ints[[i]] <- NA_real_
        msx_moverzs[[i]] <- NULL
        msx_ints[[i]] <- NULL
      }
    }
    else if (ms_lev == 1L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_1]], "value")
      
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_1]]))
      msData <- xml2::xml_contents(binData)
      r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
      r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
      msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
      msx_xs <- readBin(r1, "double", n = msx_n, size = 8L)
      msx_ys <- readBin(r2, "double", n = msx_n, size = 8L)
      
      mic <- deisotope(msx_xs, msx_ys, exclude_reporter_region = FALSE, 
                       ppm = ppm_ms1_deisotope, ms_lev = ms_lev, 
                       maxn_feats = maxn_dia_precurs, max_charge = max_ms1_charge, 
                       n_fwd = 20L, offset_upr = 30L, offset_lwr = 30L, 
                       order_mz = TRUE)
      msx_moverzs[[i]] <- mic[["masses"]]
      msx_ints[[i]] <- mic[["intensities"]]
      ms1_charges[[i]] <- mic[["charges"]] # MS2: NULL; MS1: integer vectors
    }
  }
  
  df <- tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_file,
    ms_level = ms_levs, 
    ret_time = as.numeric(ret_times) * 60, 
    scan_num = as.integer(scan_nums), 
    orig_scan = orig_scans,
    # either MS1 or MS2
    msx_moverzs = msx_moverzs, 
    msx_ints = msx_ints, 
    ms2_charges = ms2_charges, 
    msx_ns = msx_ns, 
    # MS1 placeholder
    ms1_moverzs = ms1_moverzs, 
    ms1_ints = ms1_ints, 
    # MS1: vector or empty vector for MS1; MS2: NULL
    ms1_charges = ms1_charges, 
    iso_ctr = as.numeric(iso_ctr), 
    iso_lwr = as.numeric(iso_lwr), 
    iso_upr = as.numeric(iso_upr), 
    demux = demux,
  )
  
  # correspondence of row indexes between MS1 and MS2 scans
  idxes_ms1 <- which(df$ms_level == 1L)
  diff_ms1 <- c(0L, diff(idxes_ms1))
  oks <- which(diff_ms1 > 1L)
  ms1_stas <- idxes_ms1[oks - 1L]
  ms2_stas <- ms1_stas + 1L
  ms2_ends <- idxes_ms1[oks] - 1L
  rm(list = c("oks", "diff_ms1", "idxes_ms1"))
  
  for (i in seq_along(ms1_stas)) {
    idx <- ms1_stas[i]
    df1 <- df[idx, ]
    mzs <- df1$msx_moverzs[[1]]
    lenmz <- length(mzs)
    
    if (lenmz == 0L || (lenmz == 1L && is.na(mzs)))
      next
    
    ys  <- df1$msx_ints[[1]]
    css <- df1$ms1_charges[[1]]
    
    sta <- ms2_stas[i]
    end <- ms2_ends[i]
    df2 <- df[sta:end, ]
    maxrow <- nrow(df2)
    
    half <- df2$iso_upr[[1]]
    edges <- c(df2$iso_ctr[[1]] - half, df2$iso_ctr + half)
    cuts <- findInterval(mzs, edges) # 3.9 us
    mz_cuts <- split(mzs, cuts)
    ys_cuts <- split(ys, cuts)
    css_cuts <- split(css, cuts)
    
    # excludes precursors outside the isolation window
    rows <- as.integer(names(mz_cuts))
    oks <- rows > 0L & rows <= maxrow
    rows <- rows[oks]
    
    df[sta:end, ]$ms1_moverzs[rows] <- mz_cuts[oks]
    df[sta:end, ]$ms1_ints[rows] <- ys_cuts[oks]
    df[sta:end, ]$ms1_charges[rows] <- css_cuts[oks]
  }
  
  df <- df[with(df, ms_level != 1L), ]
  bads <- unlist(lapply(df$ms1_moverzs, is.null))
  df <- df[!bads, ]
  df$iso_lwr <- df$iso_upr <- df$iso_ctr <- NULL
  
  # use singular names for consistency with DDA
  df$ms1_mass <- mapply(function (x, y) (x - 1.00727647) * y, 
                        df$ms1_moverzs, df$ms1_charges, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  restmt <- extract_mgf_rptrs(df$msx_moverzs, 
                              df$msx_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  df$msx_moverzs <- restmt[["xvals"]]
  df$msx_ints <- restmt[["yvals"]]
  df$rptr_moverzs <- restmt[["rptr_moverzs"]]
  df$rptr_ints <- restmt[["rptr_ints"]]
  
  df <- dplyr::rename(df, ms2_moverzs = msx_moverzs, ms2_ints = msx_ints, 
                      ms2_n = msx_ns, ms1_moverz = ms1_moverzs, 
                      ms1_int = ms1_ints, ms1_charge = ms1_charges, )
}


#' Proc mzML data for DDA workflows.
#' 
#' @inheritParams proc_mdda
#' @inheritParams matchMS
proc_dda <- function (spec, raw_file, idx_sc = 3L, idx_osc = 3L, 
                      idx_mslev = 2L, idx_title = 10L, idx_scanList_2 = 11L, 
                      idx_rt_2 = 2L, idx_scan_start_2 = 1L, 
                      idx_precursor_2 = 12L, idx_isolationWindow = 1L, 
                      idx_selectedIonList = 2L, idx_moverz = 1L, idx_charge = 2L, 
                      idx_ms1int = 3L, idx_bin_2 = 13L, deisotope_ms2 = TRUE, 
                      max_ms2_charge = 3L, ppm_ms2_deisotope = 10L, 
                      topn_ms2ions = 100L, quant = "none", 
                      tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                      exclude_reporter_region = FALSE, index_mgf_ms2 = FALSE)
{
  len <- length(spec)
  ret_times <- orig_scans <- scan_nums <- scan_titles <- character(len)
  ms_levs <- msx_ns <- integer(len)
  msx_moverzs <- msx_ints <- ms2_charges <- vector("list", len)
  ms1_moverzs <- ms1_ints <- ms1_charges <- character(len)
  
  is_tmt <- if (grepl("^tmt.*\\d+", quant)) TRUE else FALSE
  
  for (i in 1:len) {
    x <- spec[[i]]
    ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))[[1]]
    ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                              perl = FALSE, useBytes = FALSE))
    scan_nums[[i]] <- ids[[idx_sc]][[2]]
    orig_scans[[i]] <- ids[[idx_osc]][[2]]
    
    xc <- xml2::xml_children(x)
    ms_levs[[i]] <- ms_lev <- as.integer(xml2::xml_attr(xc[[idx_mslev]], "value"))
    scan_titles[[i]] <- xml2::xml_attr(xc[[idx_title]], "value")
    
    if (ms_lev == 2L) {
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      ret_times[[i]] <- xml2::xml_attr(scanList_ret[[idx_scan_start_2]], "value")
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # (assume one precursor, not yet chimeric)
      precursorc <- xml2::xml_children(precursor)
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      ms1_moverzs[[i]] <- xml2::xml_attr(selectedIonc[[idx_moverz]], "value")
      ms1_charges[[i]] <- xml2::xml_attr(selectedIonc[[idx_charge]], "value")
      # no precursor intensity: 01CPTAC3_Benchmarking_P_BI_20170523_BL_f12.mzML; scan=54837
      if (length(selectedIonc) >= idx_ms1int)
        ms1_ints[[i]] <- xml2::xml_attr(selectedIonc[[idx_ms1int]], "value")
      else
        ms1_ints[[i]] <- character(1)
      
      ## binaryDataArrayList
      binData <- xml2::xml_children(xml2::xml_children(xc[[idx_bin_2]]))
      msData <- xml2::xml_contents(binData)
      
      if (length(msData) == 2L) {
        r1 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[1]]))
        r2 <- .Call(base64enc:::B64_decode, xml2::xml_text(msData[[2]]))
        msx_ns[[i]] <- msx_n <- as.integer(length(r1)/8L)
        msx_xs <- readBin(r1, "double", n = msx_n, size = 8L)
        msx_ys <- readBin(r2, "double", n = msx_n, size = 8L)
        
        if (deisotope_ms2) {
          mic <- deisotope(msx_xs, msx_ys, exclude_reporter_region = is_tmt, 
                           tmt_reporter_lower = tmt_reporter_lower, 
                           tmt_reporter_upper = tmt_reporter_upper, 
                           ppm = ppm_ms2_deisotope, ms_lev = ms_lev, 
                           maxn_feats = topn_ms2ions, max_charge = max_ms2_charge, 
                           n_fwd = 10L, offset_upr = 30L, offset_lwr = 30L, 
                           order_mz = TRUE)
          msx_moverzs[[i]] <- mic[["masses"]]
          msx_ints[[i]] <- mic[["intensities"]]
          ms2_charges[[i]] <- mic[["charges"]]
        }
        else {
          msx_moverzs[[i]] <- msx_xs
          msx_ints[[i]] <- msx_ys
        }
      } 
      else {
        msx_ns[[i]] <- 0L
        msx_moverzs[[i]] <- NULL
        msx_ints[[i]] <- NULL
      }
    }
  }
  
  rows <- ms1_charges != ""
  scan_titles <- scan_titles[rows]
  ms1_moverzs <- as.numeric(ms1_moverzs[rows])
  # not "as.integer": may be > .Machine$integer.max (2147483647)
  ms1_ints <- as.numeric(ms1_ints[rows])
  ms1_charges <- as.integer(ms1_charges[rows])
  ret_times <- as.numeric(ret_times[rows]) * 60
  scan_nums <- as.integer(scan_nums[rows])
  orig_scans <- orig_scans[rows]
  ms_levs <- as.integer(ms_levs[rows])
  msx_moverzs <- msx_moverzs[rows]
  msx_ints <- msx_ints[rows]
  ms2_charges <- ms2_charges[rows]
  msx_ns <- msx_ns[rows]
  ms1_masses <- (ms1_moverzs - 1.00727647) * ms1_charges
  
  # extract the TMT region of MS2 moverz and intensity
  # (also convert reporter-ion intensities to integers)
  restmt <- extract_mgf_rptrs(msx_moverzs, 
                              msx_ints, 
                              quant = quant, 
                              tmt_reporter_lower = tmt_reporter_lower, 
                              tmt_reporter_upper = tmt_reporter_upper, 
                              exclude_reporter_region = exclude_reporter_region)
  msx_moverzs <- restmt[["xvals"]]
  msx_ints <- restmt[["yvals"]]
  rptr_moverzs <- restmt[["rptr_moverzs"]]
  rptr_ints <- restmt[["rptr_ints"]]
  
  if (index_mgf_ms2) 
    msx_moverzs <- lapply(msx_moverzs, index_mz, min_ms2mass, ppm_ms2/1E6)
  
  df <- tibble::tibble(
    scan_title = scan_titles,
    raw_file = raw_file,
    ms_level = ms_levs, 
    ms1_moverz = ms1_moverzs,
    ms1_mass = ms1_masses,
    ms1_charge = ms1_charges,
    ms1_int = ms1_ints,
    ret_time = ret_times,
    scan_num = scan_nums,
    orig_scan = orig_scans,
    ms2_moverzs = msx_moverzs,
    ms2_ints = msx_ints,
    ms2_charges = ms2_charges, 
    # before subset by min_ms2mass
    ms2_n = msx_ns, 
    rptr_moverzs = rptr_moverzs, 
    rptr_ints = rptr_ints)
}


#' Finds the bracketing positions of MS1 and MS2 levels along the scan indexes.
#'
#' @param vals A vector of MS levels. Assume that \code{vals} contains only
#'   levels 1L or 2L and always start with 1L.
#' @param out_type The type of outputs.
#'
#' @examples
#' # trailing MS1s
#' vals <- c(rep(1L, 1), rep(2L, 4), rep(1L, 2), rep(2L, 3), rep(1L, 3))
#' # find_mslev_brackets(vals)
#'
#' # multiple MS1s at the beginning
#' vals <- c(rep(1L, 2), rep(2L, 4), rep(1L, 2), rep(2L, 3), rep(1L, 3), rep(2L, 4))
#' # find_mslev_brackets(vals)
#'
#' # One MS1 at the beginning
#' vals <- c(rep(1L, 1), rep(2L, 4), rep(1L, 2), rep(2L, 3), rep(1L, 3), rep(2L, 4))
find_mslev_brackets <- function (vals, out_type = c("sim_vec", "sim_list", "full_list"))
{
  out_type <- match.arg(out_type)
  
  len <- length(vals)
  oks2 <- vals == 2L
  
  if (!oks2[len]) {
    ids2 <- which(oks2)
    vals <- vals[-c((ids2[length(ids2)] + 1L):len)]
  }
  
  ds <- c(-1L, diff(vals))
  ms2_stas <- which(ds == 1L)
  ms1_ends <- ms2_stas - 1L 
  ms1_stas <- which(ds == -1L)
  ms2_ends <- c(ms1_stas[2:length(ms1_stas)] - 1L, length(vals))
  
  if (out_type == "sim_vec")
    list(ms1_ends = ms1_ends, ms2_stas = ms2_stas, ms2_ends = ms2_ends)
  else {
    ms1_scans <- if (out_type == "sim_list")
      as.list(ms1_ends)
    else
      mapply(`:`, ms1_stas, ms1_ends)
    
    ms2_scans <- mapply(`:`, ms2_stas, ms2_ends)
    list(ms1_scans = ms1_scans, ms2_scans = ms2_scans)
  }
}


#' Finds chimeric precursors.
#' 
#' Not used; by single MS1 scan.
#' 
#' @param df1 A one-row data frame of precursor data.
#' @param df2 The data frame of MS2 corresponding to df1.
#' @param stas1 The indexes of df1 starting lines.
#' @param stas2 The indexes of df2 starting lines.
#' @param ends2 The indexes of df2 ending lines.
#' @param ppm Mass error tolerance.
#' @param maxn_precurs The maximum number of chimeric precursors.
#' @param max_ms1_charge The maximum charge state of precursors.
find_mdda_ms1s <- function (df1, df2, stas1, stas2, ends2, ppm = 6L, 
                            maxn_precurs = 5L, max_ms1_charge = 4L)
{
  if (!nrow(df2))
    return(NULL)
  
  x <- df1$msx_moverzs[[1]]
  y <- df1$msx_ints[[1]]
  
  # (1) deisotope +/- 2.01 window
  moks <- lapply(df2$iso_ctr, function (m) x > m - 2.01 & x < m + 2.01)
  xs <- lapply(moks, function (m) x[m])
  ys <- lapply(moks, function (m) y[m])
  
  mics <- mapply(
    deisotope,
    xs, ys,
    MoreArgs = list(
      ppm = ppm, ms_lev = 1L, maxn_feats = maxn_precurs, 
      max_charge = max_ms1_charge, offset_upr = 30L, 
      offset_lwr = 30L, order_mz = TRUE, bound = FALSE
    ), 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  masses <- lapply(mics, `[[`, "masses")
  charges <- lapply(mics, `[[`, "charges")
  intensities <- lapply(mics, `[[`, "intensities")
  
  # (2) subset by isolation window
  moks <- mapply(function (x, m, w) x > m - w & x < m + w, 
                 x = masses, m = df2$iso_ctr, w = df2$iso_lwr, 
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
  masses <- mapply(function (x, i) x[i], x = masses, i = moks, 
                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
  charges <- mapply(function (x, i) x[i], x = charges, i = moks, 
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  intensities <- mapply(function (x, i) x[i], x = intensities, i = moks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # impurities <- lapply(intensities, function (x) x/sum(x))
  rm(list = "moks")
  
  # update df2
  df2$ms1_moverzs <- masses
  df2$ms1_charges <- charges
  df2$ms1_ints <- intensities
  
  df2[, c("ms1_moverzs", "ms1_charges", "ms1_ints")]
}


#' Finds the positions of logical gates.
#' 
#' @param vec A logical vector.
#' 
#' @examples
#' library(mzion)
#' 
#' # starts at low and ends at low
#' vec <- c(rep(0L, 2), rep(1L, 4L), rep(0L, 3L), rep(1L, 3L), rep(0L, 1L))
#' # find_gatepos(vec)
#' 
#' # low-high
#' vec <- c(rep(0L, 2), rep(1L, 4L), rep(0L, 3L), rep(1L, 3L), rep(0L, 1L), rep(1L, 2L))
#' 
#' # high-low
#' vec <- c(rep(1L, 4L), rep(0L, 3L), rep(1L, 3L), rep(0L, 1L), rep(1L, 2L), rep(0L, 2))
#' 
#' # high-high
#' vec <- c(rep(1L, 4L), rep(0L, 3L), rep(1L, 3L), rep(0L, 1L), rep(1L, 3L))
#' 
#' # all-zeros
#' # find_gatepos(diff(c(6L, 12L, 18L, 164L), 1L) == 1L)
#' 
#' # all-ones
#' # find_gatepos(rep(1L, 5L))
#' 
#' # up-width == 1L
#' # find_gatepos(diff(c(143L, 159L, 310L, 311L, 316L), 1L) == 1L)
#' 
#' # single TRUE
#' # find_gatepos(diff(c(310L, 311L), 1L) == 1L)
find_gatepos <- function (vec)
{
  # stopifnot(all(vec %in% c(0L, 1L)))
  lenv <- length(vec)
  
  # if (lenv == 1L) {
  #   if (vec)
  #     return(list(up_starts = 1L, up_ends = 1L))
  #   else
  #     return(NULL)
  # }
  
  if (all(vec))
    return(list(up_starts = 1L, up_ends = lenv))
  
  if (!any(vec)) 
    return(NULL)

  ds <- diff(vec)
  
  # if (!any(ds != 0L)) # +/-1L
  #   return(NULL)

  ups <- which(ds == 1L) + 1L
  dns <- which(ds == -1L)
  lenu <- length(ups)
  lend <- length(dns)
  
  if (lenu == lend) {
    if (ups[[1]] <= dns[[1]])
      list(up_starts = ups, up_ends = dns)
    else
      list(up_starts = c(1L, ups), up_ends = c(dns, lenv))
  } 
  else if (lenu > lend)
    list(up_starts = ups, up_ends = c(dns, lenv))
  else
    list(up_starts = c(1L, ups), up_ends = dns)
}


#' Collapse MS1 intensities.
#' 
#' Allow adjacent values in unv and later collapse adjacent columns/values.
#' 
#' @param xs Vectors of m-over-z values.
#' @param ys Vectors of intensity values.
#' @param lwr The lower mass limit.
#' @param step The step size for mass bins.
#' 
#' @examples
#' # Twos adjacent bins of xs: 392.1796, 392.1845
#' # xs longer than the universe
#' xs <- list(c(391.1883,391.2848,391.6627,391.6995,392.1646,392.1796,392.1845,
#'            392.2030,392.2887,392.6641,392.7812,393.0833,393.2975))
#' ys <- list(c(12827.41,337002.19,617819.69,18045.10,205851.53,15194.98,11318.61,
#'              12970.02,118604.48,75726.89,11676.51,23723.18,55749.93))
#' # collapse_mms1ints(xs, ys, 389.6529)
#' 
#' xs <- list(c(400.6596,401.7076,402.1813,402.1944,402.1969,402.2094,402.5438,402.7112,403.1812,404.1777), 
#'            c(400.6599,401.7075,402.1954,402.1975,402.7112,403.1822,404.2777))
#' ys <- list(c(24003.98,53431.96,110619.26,10988.55,12291.00,140045.06,67601.16,11413.04,21651.61,16686.06), 
#'            c(10000.1,40000.1,20000.1,50000.1,2500.2,5000.1,30000.1))
#' # collapse_mms1ints(xs, ys, 400.1994)
#' 
#' xs <- ys <- vector("list", 13L)
#' xs[[7]] <- 954.607849; xs[[8]] <- 954.630249; xs[[10]] <- 954.622925
#' ys[[7]] <- 15706.2627; ys[[8]] <- 19803.5879; ys[[10]] <- 31178.9648
#' # collapse_mms1ints(xs, ys, 951.089731)
collapse_mms1ints <- function (xs, ys, lwr, step = 1e-5)
{
  # all xs can be NULL
  oks <- lengths(xs) > 0L

  if (!any(oks))
    return(list(xmeans = NULL, ysum = NULL))
  
  xs <- xs[oks]
  ys <- ys[oks]
  ixs <- lapply(xs, index_mz, lwr, step)
  
  # removes duplicated ixs
  for (i in seq_along(ixs)) {
    ix <- ixs[[i]]
    x <- xs[[i]]
    y <- ys[[i]]
    oks <- !duplicated(ix)
    ixs[[i]] <- ix[oks]
    xs[[i]]  <- x[oks]
    ys[[i]]  <- y[oks]
  }

  unv <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv <- sort(unique(unv))
  x0 <- y0 <- rep(NA_real_, length(unv))
  
  # maps ixs onto unv (presence or absence)
  # also note one-to-one correspondence between ixs and xs
  xps <- lapply(ixs, function (x) unv %in% x)
  
  xout <- mapply(function (i, v) {
    x0[i] <- v
    x0
  }, xps, xs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  yout <- mapply(function (i, v) {
    y0[i] <- v
    y0
  }, xps, ys, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  xout <- do.call(rbind, xout)
  yout <- do.call(rbind, yout)
  
  # collapses adjacent (+/-1 in index) columns in unv
  if (any(adjs <- diff(unv, 1L) == 1L)) {
    ps <- find_gatepos(adjs)
    ps1 <- ps[[1]]
    ps2 <- ps[[2]] <- ps[[2]] + 1L
    ps12 <- mapply(function (x, y) x:y, ps1, ps2, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    exs <- mapply(function (x, y) (x + 1L):y, ps1, ps2, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    exs <- .Internal(unlist(exs, recursive = FALSE, use.names = FALSE))
    
    ysum <- colSums(yout, na.rm = TRUE)
    xmeans <- colSums(xout * yout, na.rm = TRUE)/ysum
    # xmeans <- colMeans(xout, na.rm = TRUE)
    ms <- lapply(ps12, function (cols) sum(xmeans[cols])/length(cols))
    xmeans[ps1] <- .Internal(unlist(ms, recursive = FALSE, use.names = FALSE))
    xmeans <- xmeans[-exs]
    
    # ysum <- colSums(yout, na.rm = TRUE)
    us <- lapply(ps12, function (cols) sum(ysum[cols], na.rm = TRUE))
    ysum[ps1] <- .Internal(unlist(us, recursive = FALSE, use.names = FALSE))
    ysum <- ysum[-exs]
  }
  else {
    # xmeans <- colMeans(xout, na.rm = TRUE)
    ysum <- colSums(yout, na.rm = TRUE)
    xmeans <- colSums(xout * yout, na.rm = TRUE)/ysum
  }
  
  list(xmeans = xmeans, ysum = ysum)
}


#' Finds mDDA precursors.
#'
#' Averages of multiple MS1 scans.
#'
#' @param df1 A sub data frame of MS1.
#' @param df2 A sub data frame of MS2.
#' @param stas1 Indexes of MS1 starting lines (in the full data frame).
#' @param stas2 Indexes of MS2 starting lines.
#' @param ends2 Indexes of MS2 ending lines.
#' @param ppm Mass error tolerance.
#' @param maxn_precurs Maximum number of precursors for consideration.
#' @param max_ms1_charge Maximum charge state of precursors for consideration.
#' @param width The width of an MS1 window. A wide window is used for containing
#'   isotope envelops.
#' @param n_fwd Forward looking up to \code{n_fwd} mass entries.
#' @param step The bin size in converting numeric m-over-z values to integers.
find_mdda_mms1s <- function (df1, df2, stas1, stas2, ends2, ppm = 10L, 
                             maxn_precurs = 5L, max_ms1_charge = 4L, 
                             n_fwd = 20L, width = 2.01, step = ppm/1e6)
{
  # for all (6+1+6) MS1 frames subset by one MS2 iso-window
  ansx1 <- ansy1 <- vector("list", len1 <- length(stas1))
  
  # for all MS2s from averaged (6+1+6 -> 1) MS1s 
  ansx2 <- ansy2 <- vector("list", len2 <- length(df2$iso_ctr))
  
  # go through MS2
  for (i in 1:len2) {
    m2  <- df2$iso_ctr[[i]]
    lwr <- m2 - width
    upr <- m2 + width
    
    # go through MS1s
    for (j in 1:len1) {
      x1s <- df1$msx_moverzs[[j]]
      y1s <- df1$msx_ints[[j]]
      oks <- x1s > lwr & x1s < upr
      ansx1[[j]] <- x1s[oks]
      ansy1[[j]] <- y1s[oks]
    }
    
    # collapse MS1s
    ans <- collapse_mms1ints(ansx1, ansy1, lwr, step)
    ansx <- ans[["xmeans"]]
    ansy <- ans[["ysum"]]
    
    if (!is.null(ansx)) {
      ansx2[[i]] <- ansx
      ansy2[[i]] <- ansy
    }
  }
  
  # ansx2[[i]] can be NULL (no precursor found in the isolation window)
  mics <- mapply(
    deisotope,
    ansx2, ansy2, df2$iso_ctr, 
    MoreArgs = list(
      exclude_reporter_region = FALSE, 
      ppm = ppm, ms_lev = 1L, maxn_feats = maxn_precurs, 
      max_charge = max_ms1_charge, n_fwd = n_fwd, offset_upr = 30L, 
      offset_lwr = 30L, order_mz = TRUE, bound = FALSE
    ), 
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  masses <- lapply(mics, `[[`, "masses")
  charges <- lapply(mics, `[[`, "charges")
  intensities <- lapply(mics, `[[`, "intensities")
  
  # (2) subset by isolation window
  moks <- mapply(function (x, m, w) {
    if (any(oks <- x > m - w & x < m + w)) oks else rep(TRUE, length(x))
    # if ((len <- length(x)) == 1L && is.na(x))
    #   FALSE
    # else {
    #   if (any(oks <- x > m - w & x < m + w)) oks else rep(TRUE, len)
    # }
  }, x = masses, m = df2$iso_ctr, w = df2$iso_lwr, 
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  # masses may be NULL and the corresponding moks is logical(0)
  masses <- mapply(function (x, i) x[i], x = masses, i = moks, 
                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
  charges <- mapply(function (x, i) x[i], x = charges, i = moks, 
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  intensities <- mapply(function (x, i) x[i], x = intensities, i = moks, 
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # impurities <- lapply(intensities, function (x) x/sum(x))
  rm(list = "moks")
  
  # (3) update df2
  rows <- lengths(masses) > 0L
  df2$ms1_moverzs[rows] <- masses[rows]
  df2$ms1_charges[rows] <- charges[rows]
  df2$ms1_ints[rows] <- intensities[rows]
  
  df2$ms1_masses <- mapply(function (x, y) (x - 1.00727647) * y, 
                           df2$ms1_moverzs, df2$ms1_charges, 
                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  df2[, c("ms1_moverzs", "ms1_masses", "ms1_charges", "ms1_ints")]
}


#' Preprocessing Bruker's MGF
#' 
#' Some entries may have no CHARGE line
#'
#' @param file A file name
#' @param begin_offset The number of lines before a BEGIN line.
#' @param charge_offset The number lines after a BEGIN line to a following
#'   CHARGE line.
prepBrukerMGF <- function (file = NULL, begin_offset = 5L, charge_offset = 5L)
{
  if (is.null(file))
    stop("`file` cannot be NULL.")
  
  message("Processing: ", file)
  
  lines <- readLines(file)
  hd <- lines[1:100]
  
  ### Some MGFs may have additional "SCANS=" lines
  pat_sc <- "SCANS"
  sc <- .Internal(which(stringi::stri_startswith_fixed(hd, pat_sc)))
  
  if (length(sc)) {
    lscs  <- .Internal(which(stringi::stri_startswith_fixed(lines, pat_sc)))
    lines <- lines[-lscs]
    rm(list = "lscs")
  }
  rm(list = c("pat_sc", "sc"))

  ###PrecursorID: 
  pat_pr <- "###PrecursorID"
  pr <- .Internal(which(stringi::stri_startswith_fixed(hd, pat_pr)))
  
  if (length(pr)) {
    lprs  <- .Internal(which(stringi::stri_startswith_fixed(lines, pat_pr)))
    lines <- lines[-lprs]
    rm(list = "lprs")
  }
  rm(list = c("pat_pr", "pr"))

  ###CCS: 
  pat_ccs <- "###CCS"
  ccs <- .Internal(which(stringi::stri_startswith_fixed(hd, pat_ccs)))
  
  if (length(ccs)) {
    ccs  <- .Internal(which(stringi::stri_startswith_fixed(lines, pat_ccs)))
    lines <- lines[-ccs]
    rm(list = "ccs")
  }
  rm(list = c("pat_ccs", "ccs"))

  ## Processing
  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends   <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))
  hdrs   <- 1:(begins[1]-begin_offset-1L)
  
  zls <- lines[begins+5L]
  # oks <- grepl("^CHARGE", zls) & (zls != "CHARGE=1+")
  oks <- grepl("^CHARGE", zls) & (zls != "CHARGE=1+") & grepl("^###Mobility", lines[begins-1L])
  rm(list = "zls")
  
  b_oks <- begins[oks] - begin_offset
  e_oks <- ends[oks]
  
  ranges <- mapply(function (x, y) x:y, b_oks, e_oks, SIMPLIFY = TRUE)
  ranges <- do.call(`c`, ranges)
  ranges <- c(hdrs, ranges)
  
  writeLines(lines[ranges], file)
}


#' Parallel \link{prepBrukerMGF}
#' 
#' @param filepath A file path to MGF.
#' @param n_cores The number of CPU cores.
#' @export
mprepBrukerMGF <- function (filepath, n_cores = 32L) 
{
  message("Preparing Bruker's MGFs.")
  
  files <- list.files(filepath, pattern = "\\.mgf$", full.names = TRUE, 
                      recursive = TRUE)
  
  len <- length(files)
  
  if (!len)
    stop("No MGF files found.")
  
  if (n_cores > 1L)
    n_cores <- min(parallel::detectCores(), n_cores, len)
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterApply(cl, files, prepBrukerMGF)
    parallel::stopCluster(cl)
  }
  else 
    lapply(files, prepBrukerMGF)
  
  message("Done preparing Bruker's MGFs.")
}


