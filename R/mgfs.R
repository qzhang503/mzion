#' Helper in loading MGFs.
#' 
#' @param topn_dia_ms2ions The top-n DIA-MS2 ions for deisotoping.
#' @param is_ms1_three_frame Logical; is the searches by the three frames of
#'   preceding, current and following.
#' @param is_ms2_three_frame Logical; is the searches by the three frames of
#'   preceding, current and following.
#' @param mgf_cutmzs Cut points of MS1 m-over-z values in peak picking.
#' @param mgf_cutpercs The counts of MS2 features in each region of
#'   \code{mgf_cutmzs}.
#' @inheritParams matchMS
load_mgfs <- function (out_path = NULL, mgf_path = NULL, topn_ms2ions = 150L, 
                       maxn_dia_precurs = 1000L, # max MS1 deisotoping features
                       topn_dia_ms2ions = 500L, # max MS2 deisotoping features
                       delayed_diams2_tracing = FALSE, 
                       n_dia_ms2bins = 3L, n_dia_scans = 4L, 
                       min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       ppm_ms1 = 20L, ppm_ms2 = 20L, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       is_ms1_three_frame = TRUE, is_ms2_three_frame = TRUE, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       enzyme = "trypsin_p", 
                       deisotope_ms2 = TRUE, grad_isotope = 1.6, fct_iso2 = 3.0,
                       max_ms2_charge = 3L, use_defpeaks = FALSE, 
                       maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                       ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
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
  if ((!ok_pars) && isTRUE(enzyme == "noenzyme")) 
    ok_pars <- TRUE
  
  scns <- list.files(mgf_path, pattern = "^scan_map_.*\\.rds$")
  ques <- list.files(mgf_path, pattern = "^mgf_queries_.*\\.rds$")
  n_scns <- length(scns)
  n_ques <- length(ques)
  
  if (n_scns) {
    ok_mgfs <- if (n_scns == n_ques) TRUE else FALSE
  }
  else {
    # backward compatible
    ok_mgfs <- local({
      raws_indexes <- file.path(mgf_path, "raw_indexes.rds")
      
      if (file.exists(raws_indexes)) {
        raws <- qs::qread(raws_indexes)
        ques <- list.files(mgf_path, pattern = "^mgf_queries_\\d+\\.rds$")
        if (length(raws) == length(ques)) TRUE else FALSE
      }
      else {
        FALSE
      }
    })
  }

  rm(list = c("scns", "ques", "n_scns", "n_ques"))
  
  if (ok_pars && ok_mgfs) {
    message("Found cached MGFs.")
    .savecall <- FALSE
    return(NULL)
  }
  
  message("Processing raw MGFs.")
  
  ppm_ms1_bin <- calc_threeframe_ppm(ppm_ms1)
  ppm_ms2_bin <- calc_threeframe_ppm(ppm_ms2)

  delete_files(
    out_path, 
    ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.(mzML|mzml)$", 
                "\\.xlsx$", "\\.xls$", "\\.csv$", "\\.txt$", "\\.pars$", 
                "^mgf$", "^mgfs$", "^mzML$", "^mzMLs$", 
                "Calls", "^PSM$", "^Peptide$", "^Protein$", 
                "fraction_scheme.rda", "label_scheme.rda", 
                "label_scheme_full.rda"))

  fi_mgf   <- list.files(path = mgf_path, pattern = "^.*\\.(mgf|MGF)$")
  fi_mzml  <- list.files(path = mgf_path, pattern = "^.*\\.(mzML|mzml)$")
  len_mgf  <- length(fi_mgf)
  len_mzml <- length(fi_mzml)
  
  if (len_mgf && len_mzml)
    stop("Peak lists need to be in either MGF or mzML, but not both.")
  
  filelist <- if (len_mgf) fi_mgf else fi_mzml

  if (len_mgf) {
    readMGF(
      filepath = mgf_path,
      filelist = filelist, 
      out_path = out_path, 
      topn_ms2ions = topn_ms2ions,
      min_mass = min_mass,
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass,
      max_ms2mass = max_ms2mass, 
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge,
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      ppm_ms1 = ppm_ms1_bin, # change arg name from ppm_ms1 to ppm_ms1_bin
      ppm_ms2 = ppm_ms2_bin, # change arg name from ppm_ms2 to ppm_ms2_bin
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region, 
      mgf_cutmzs = mgf_cutmzs, 
      mgf_cutpercs = mgf_cutpercs, 
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
    readmzML(
      filelist = filelist, 
      mgf_path = mgf_path, 
      topn_ms2ions = topn_ms2ions, 
      topn_dia_ms2ions = topn_dia_ms2ions, 
      maxn_dia_precurs = maxn_dia_precurs, 
      n_dia_ms2bins = n_dia_ms2bins, 
      n_dia_scans = n_dia_scans, 
      delayed_diams2_tracing = delayed_diams2_tracing, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass, 
      max_ms2mass = max_ms2mass, 
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge, 
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      ppm_ms1 = ppm_ms1_bin, # change arg name from ppm_ms1 to ppm_ms1_bin
      ppm_ms2 = ppm_ms2_bin, # change arg name from ppm_ms2 to ppm_ms2_bin
      tmt_reporter_lower = tmt_reporter_lower, 
      tmt_reporter_upper = tmt_reporter_upper, 
      exclude_reporter_region = exclude_reporter_region, 
      is_ms1_three_frame = is_ms1_three_frame, 
      is_ms2_three_frame = is_ms2_three_frame, 
      mgf_cutmzs = mgf_cutmzs, 
      mgf_cutpercs = mgf_cutpercs, 
      enzyme = enzyme, 
      deisotope_ms2 = deisotope_ms2, 
      grad_isotope = grad_isotope, 
      fct_iso2 = fct_iso2,
      max_ms2_charge = max_ms2_charge, 
      use_defpeaks = use_defpeaks, 
      maxn_mdda_precurs = maxn_mdda_precurs, 
      n_mdda_flanks = n_mdda_flanks, 
      ppm_ms1_deisotope = ppm_ms1_deisotope, 
      ppm_ms2_deisotope = ppm_ms2_deisotope, 
      quant = quant, 
      digits = digits)
  }
  else {
    stop("No peak lists at an mzML or MGF format found.")
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
#' @param out_path An output path.
#' @inheritParams load_mgfs
#' @inheritParams matchMS
#' @inheritParams frames_adv
#' @import stringi
#' @import readr
#' @import fs
readMGF <- function (filepath = NULL, filelist = NULL, out_path = NULL, 
                     topn_ms2ions = 150L, min_mass = 200L, max_mass = 4500L, 
                     min_ms2mass = 115L, max_ms2mass = 4500L, 
                     min_ms1_charge = 2L, max_ms1_charge = 4L,
                     min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                     min_ret_time = 0, max_ret_time = Inf, 
                     ppm_ms1 = 10L, ppm_ms2 = 10L, 
                     tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                     exclude_reporter_region = FALSE, 
                     mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                     deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                     use_defpeaks = FALSE, maxn_dia_precurs = 1000L, 
                     maxn_mdda_precurs = 1L, n_mdda_flanks = 6L, 
                     ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                     quant = "none", digits = 4L) 
{
  if (maxn_mdda_precurs >= 1L) {
    warning("No multi-precursor DDA with MGF. Use mzML to enable the feature.")
    maxn_mdda_precurs <- 0L
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
  warning("An mzML or MGF with multiple RAWs not supported since v1.3.4. ", 
          "Each peaklist file need contain exactly one RAW file.")
  
  out <- vector("list", len)
  
  for (i in seq_along(filelist)) {
    file <- filelist[i]
    temp_dir <- file.path(filepath, paste0("temp_", i))
    
    message("Loading '", file, "'.")
    
    out[[i]] <- read_mgf_chunks(
      filepath = filepath, 
      temp_dir = temp_dir, 
      raw_id = i, 
      topn_ms2ions = topn_ms2ions,
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge, 
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass,
      max_ms2mass = max_ms2mass, 
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2,
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
  
  raws <- unlist(out, recursive = FALSE, use.names = TRUE)
  qs::qsave(raws, file.path(filepath, "raw_indexes.rds"), preset = "fast")

  invisible(NULL)
}


#' Post-processing of MGF or mzML
#' 
#' Calculates mass \code{frame}s etc.
#' 
#' @param df A data frame of processed peak lists.
#' @param raw_id An ID to replace the original RAW file name.
#' @param mgf_path A path to mzML or MGF files.
post_readmgf <- function (df, raw_id, mgf_path) 
{
  raw <- unique(df$raw_file)
  
  if (length(raw) > 1L)
    stop("An mzML or MGF with multiple RAWs not supported since v1.3.4. ", 
         "Each peaklist file need contain exactly one RAW file.")
  
  raw_map <- raw_id
  names(raw_map) <- raw
  df$raw_file <- raw_id
  
  scans <- df$scan_title
  scans_map <- df$scan_title <- seq_along(scans)
  names(scans_map) <- scans
  
  qs::qsave(df, file.path(mgf_path, paste0("mgf_queries_", raw, ".rds")), 
            preset = "fast")
  qs::qsave(scans_map, file.path(mgf_path, paste0("scan_map_", raw, ".rds")), 
            preset = "fast")

  invisible(raw_map)
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
    if (!is.null(ok)) fs::file_delete(ok)
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


#' Reads MGFs in chunks.
#'
#' @param temp_dir A temporary path of MGFs.
#' @param raw_id An ID to RAW file name.
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
read_mgf_chunks <- function (filepath, temp_dir, raw_id = 1L, topn_ms2ions = 150L, 
                             min_ms1_charge = 2L, max_ms1_charge = 4L,
                             min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                             min_ret_time = 0, max_ret_time = Inf, 
                             min_mass = 200L, max_mass = 4500L, 
                             min_ms2mass = 115L, max_ms2mass = 4500L, 
                             ppm_ms1 = 10L, ppm_ms2 = 10L, 
                             mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                             type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                             n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L, 
                             n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L, 
                             n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                             sep_pepmass = " ", nfields_pepmass = 2L, 
                             raw_file = NULL, 
                             tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                             exclude_reporter_region = FALSE, 
                             deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                             maxn_dia_precurs = 1000L, 
                             use_defpeaks = FALSE, 
                             maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                             ppm_ms1_deisotope = 10L, ppm_ms2_deisotope = 10L, 
                             quant = "none", digits = 4L) 
{
  filelist <- list.files(path = temp_dir, pattern = "^.*\\.mgf$")
  
  if (!(len <- length(filelist))) 
    stop("No mgf files under ", temp_dir)
  
  if (len == 1L) {
    out <- proc_mgf_chunks(
      file.path(temp_dir, filelist),
      topn_ms2ions = topn_ms2ions,
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge,
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass,
      max_ms2mass = max_ms2mass, 
      ppm_ms1 = ppm_ms1,
      ppm_ms2 = ppm_ms2,
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
        "integerize_ms2ints"), 
      envir = environment(mzion::matchMS)
    )
    
    out <- parallel::clusterApply(
      cl, file.path(temp_dir, filelist),
      proc_mgf_chunks,
      topn_ms2ions = topn_ms2ions,
      min_ms1_charge = min_ms1_charge, 
      max_ms1_charge = max_ms1_charge,
      min_scan_num = min_scan_num, 
      max_scan_num = max_scan_num, 
      min_ret_time = min_ret_time, 
      max_ret_time = max_ret_time, 
      min_mass = min_mass, 
      max_mass = max_mass, 
      min_ms2mass = min_ms2mass,
      max_ms2mass = max_ms2mass, 
      ppm_ms1 = ppm_ms1,
      ppm_ms2 = ppm_ms2,
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
    afs <- list.files(path = temp_dir, pattern = "^.*\\_af.mgf$")
    idxes <- sort(as.integer(gsub("^chunk_(\\d+)_af\\.mgf", "\\1", afs)))
    afs <- paste0("chunk_", idxes, "_af.mgf")
    afs <- afs[-length(afs)]
  })
  
  bfs <- local({
    bfs <- list.files(path = temp_dir, pattern = "^.*\\_bf.mgf$")
    idxes <- sort(as.integer(gsub("^chunk_(\\d+)_bf\\.mgf", "\\1", bfs)))
    bfs <- paste0("chunk_", idxes, "_bf.mgf")
    bfs <- bfs[-1]
  })
  
  # stopifnot(length(afs) == length(bfs))
  
  gaps <- purrr::map2(afs, bfs, function (x, y) {
    af <- stringi::stri_read_lines(file.path(temp_dir, x))
    bf <- stringi::stri_read_lines(file.path(temp_dir, y))
    ab <- append(af, bf)
    # perfect case of no gaps: two lines of "" and ""
    if (length(ab) > 2L) ab else NULL
  })
  
  gaps <- unlist(gaps, use.names = FALSE)
  write(gaps, file.path(temp_dir, "gaps.mgf"))
  
  local({
    nms <- list.files(path = file.path(temp_dir), pattern = "^.*\\_[ab]f.mgf$")
    if (length(nms)) suppressMessages(file.remove(file.path(temp_dir, nms)))
  })
  
  if (!is.null(gaps)) {
    out <- dplyr::bind_rows(
      out,
      proc_mgfs(lines = gaps,
                topn_ms2ions = topn_ms2ions,
                min_ms1_charge = min_ms1_charge, 
                max_ms1_charge = max_ms1_charge,
                min_scan_num = min_scan_num, 
                max_scan_num = max_scan_num, 
                min_ret_time = min_ret_time, 
                max_ret_time = max_ret_time, 
                min_mass = min_mass, 
                max_mass = max_mass, 
                min_ms2mass = min_ms2mass,
                max_ms2mass = max_ms2mass, 
                ppm_ms1 = ppm_ms1,
                ppm_ms2 = ppm_ms2,
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
  
  post_readmgf(out, raw_id = raw_id, mgf_path = filepath) 
}


#' Processes MGF entries in chunks.
#'
#' @param file A chunk of MGF (chunk_1.mgf etc.) with a prepending file path.
#' @inheritParams readMGF
#' @inheritParams read_mgf_chunks
proc_mgf_chunks <- function (file, topn_ms2ions = 150L, 
                             min_ms1_charge = 2L, max_ms1_charge = 4L, 
                             min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                             min_ret_time = 0, max_ret_time = Inf, 
                             min_mass = 200L, max_mass = 4500L, 
                             min_ms2mass = 115L, max_ms2mass = 4500L, 
                             ppm_ms1 = 10L, ppm_ms2 = 10L, 
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
                             deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                             use_defpeaks = FALSE, maxn_dia_precurs = 1000L, 
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
                   min_ms1_charge = min_ms1_charge, 
                   max_ms1_charge = max_ms1_charge,
                   min_scan_num = min_scan_num, 
                   max_scan_num = max_scan_num, 
                   min_ret_time = min_ret_time, 
                   max_ret_time = max_ret_time, 
                   min_mass = min_mass, 
                   max_mass = max_mass, 
                   min_ms2mass = min_ms2mass,
                   max_ms2mass = max_ms2mass, 
                   ppm_ms1 = ppm_ms1,
                   ppm_ms2 = ppm_ms2,
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
proc_mgfs <- function (lines, topn_ms2ions = 150L, 
                       min_ms1_charge = 2L, max_ms1_charge = 4L, 
                       min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                       min_ret_time = 0, max_ret_time = Inf, 
                       min_mass = 200L, max_mass = 4500L, 
                       min_ms2mass = 115L, max_ms2mass = 4500L, 
                       ppm_ms1 = 10L, ppm_ms2 = 10L, 
                       mgf_cutmzs = numeric(), mgf_cutpercs = numeric(), 
                       type_mgf = "msconv_thermo", n_bf_begin = 0L, 
                       n_spacer = 0L, n_hdr = 5L, n_to_pepmass = 3L,
                       n_to_title = 1L, n_to_scan = 0L, n_to_rt = 2L,
                       n_to_charge = 4L, sep_ms2s = " ", nfields_ms2s = 2L, 
                       sep_pepmass = " ", nfields_pepmass = 2L, raw_file = NULL, 
                       tmt_reporter_lower = 126.1, tmt_reporter_upper = 135.2, 
                       exclude_reporter_region = FALSE, 
                       deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                       use_defpeaks = FALSE, maxn_dia_precurs = 1000L, 
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
  if (FALSE) {
    if ((polarity <- gsub(".*([\\+-])$", "\\1", ms1_charges[[1]])) == "+")
      ms1_charges <- gsub("[\\+]$", "", ms1_charges)
    else if (polarity == "-")
      ms1_charges <- gsub("-$", "", ms1_charges)
  }
  ms1_charges <- as.integer(ms1_charges)

  # timsTOF no CHARGE line -> NAs introduced by coercion
  ms1_masses <- mapply(function (x, y) x * y - y * 1.00727647, 
                       ms1_moverzs, ms1_charges, 
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)

  ans_rts <- reset_rettimes(ret_times = ret_times, min_ret_time = min_ret_time, 
                            max_ret_time = max_ret_time)
  min_ret_time <- ans_rts$min_ret_time
  max_ret_time <- ans_rts$max_ret_time

  rows <- (ms1_charges >= min_ms1_charge & ms1_charges <= max_ms1_charge & 
             ret_times >= min_ret_time & ret_times <= max_ret_time & 
             ms1_masses >= min_mass & ms1_masses <= max_mass & 
             # timsTOF: no MS1 masses
             !is.na(ms1_masses))
  
  # timsTOF may have undetermined charge states
  if (length(na_rows <- .Internal(which(is.na(rows))))) 
    rows[na_rows] <- FALSE

  begins <- begins[rows]
  ends <- ends[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- round(ret_times[rows], digits = 2L)
  ms1_masses <- ms1_masses[rows]
  rm(list = "rows")

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
  
  ## Low priority: no data subsetting by scan_nums; use retention time instead

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
    mics <- mapply(find_ms1stat, moverzs = ms2_moverzs, msxints = ms2_ints, 
                   MoreArgs = list(
                     center = 0, 
                     exclude_reporter_region = is_tmt, 
                     tmt_reporter_lower = tmt_reporter_lower, 
                     tmt_reporter_upper = tmt_reporter_upper, 
                     ppm = ppm_ms2_deisotope, 
                     ms_lev = 2L, maxn_feats = topn_ms2ions, 
                     max_charge = max_ms2_charge, n_fwd = 10L, 
                     offset_upr = 30L, offset_lwr = 30L, order_mz = FALSE
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


#' Resets the range of retention times by experimental observations.
#' 
#' Needed since \code{max_ret_time} can be negative.
#' 
#' @param ret_times A vector of experimental retention times.
#' @inheritParams matchMS
reset_rettimes <- function (ret_times, min_ret_time = 0, max_ret_time = Inf) 
{
  min_rt <- min(ret_times, na.rm = TRUE)
  max_rt <- max(ret_times, na.rm = TRUE)
  
  if (min_ret_time < min_rt)
    min_ret_time <- min_rt 
  
  if (max_ret_time > 0) {
    if (max_ret_time > max_rt)
      max_ret_time <- max_rt
  }
  else {
    max_ret_time <- max_rt + max_ret_time
    
    if (max_ret_time < min_ret_time) {
      warning("Invalid \"max_ret_time\". Choose a less negative value." )
      max_ret_time <- max_rt
    }
  }
  
  list(min_ret_time = min_ret_time, max_ret_time = max_ret_time)
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
                         topn_ms2ions = 150L, mgf_cutmzs = numeric(), 
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
  if (isTRUE(grepl("^tmt.*\\d+", quant))) {
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


#' Converts ms2_moverzs to integers.
#'
#' Note that as.integer is needed. When an indexed x equals 0 but is double, it
#' can cause %fin% to fail.
#'
#' @param from Numeric; the starting MS1 mass.
#' @param x Numeric; MS2 mass.
#' @param d Numeric; \eqn{ppm * 10E-6}.
#' @examples
#' xs <- ceiling(log(c(114.0916, 114.9999)/115)/log(1+8E-6))
#' match(xs, xs); fastmatch::fmatch(xs, xs)
#' 
#' ys <- as.integer(xs)
#' match(ys, ys); fastmatch::fmatch(ys, ys)
index_mz <- function (x, from = 115L, d = 1E-5) 
  as.integer(ceiling(log(x/from)/log(1+d)))


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
  #   stop("The tag of `END IONS` not found in MGF.")

  
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

    if (isTRUE(grepl(file_msconvert_thermo, ln_tit)) && 
        isTRUE(grepl(scan_msconv_thermo, ln_tit))) 
      "msconv_thermo"
    else if (isTRUE(grepl(file_msconvert_pasef, ln_tit)) && 
             isTRUE(grepl(scan_msconv_pasef, ln_tit)))
      "msconv_pasef"
    else if (isTRUE(grepl(file_pd, ln_tit)) && 
             isTRUE(grepl(scan_pd, ln_tit))) 
      "pd"
    else if (isTRUE(grepl(file_default_pasef, ln_tit)) && 
             isTRUE(is.null(scan_default_pasef)))
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

  # the first entry in the original MGF may have no CHARGE line
  if (is.na(n_to_charge) && type == "default_pasef")
    n_to_charge <- 5L
  
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


