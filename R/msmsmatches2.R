#' Matches theoretical peptides (parallel by mgf chunks).
#'
#' All files under `out_path` are removed if incur \code{calc_pepmasses} in the
#' upstream.
#'
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @param .path_bin The file path to binned precursor masses.
#' @param reframe_mgfs Logical; if TRUE, recalculates the frame indexes of MGFs.
#' @param first_search Logical; is the first search (for MGF mass calibration)
#'   or not.
#' @param .savecall Logical; if TRUE, saves the current call.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @inheritParams frames_adv
#' @inheritParams add_var_masses
#' @import parallel
ms2match <- function (mgf_path, aa_masses_all, out_path, .path_bin, 
                      mod_indexes, type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                      maxn_sites_per_vmod = 3L, maxn_fnl_per_seq = 64L, 
                      maxn_vnl_per_seq = 64L, maxn_vmods_sitescombi_per_pep = 64L, 
                      minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 20L, 
                      min_mass = 200L, max_mass = 4500L, min_ms2mass = 115L, 
                      quant = "none", ppm_reporters = 10L, 
                      by_modules = TRUE, reframe_mgfs = FALSE, 
                      n_13c = NULL, ms1_notches = 0, 

                      # dummies
                      fasta, acc_type, acc_pattern,
                      topn_ms2ions, fixedmods, varmods, 
                      enzyme, 
                      maxn_fasta_seqs, maxn_vmods_setscombi, 
                      min_len, max_len, max_miss, 
                      
                      index_mgf_ms2 = FALSE, first_search = FALSE, 
                      .savecall = TRUE) 
{
  options(digits = 9L)
  
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    }, add = TRUE
  )
  
  # Check cached 
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  fml_nms <- names(formals(fun))
  faa <- file.path(out_path, "aa_masses_all.rds")
  
  # (OK as `argument` not for users)
  # min_mass and max_mass only for calib_ms1mass, not to be changed by users
  # args_except <- c("quant", "min_mass", "max_mass", "by_modules")
  args_except <- c("by_modules", "first_search")
  fml_incl    <- fml_nms[!fml_nms %in% args_except]
  cache_pars  <- find_callarg_vals(time = NULL, 
                                   path = file.path(out_path, "Calls"), 
                                   fun = paste0(fun, ".rda"), 
                                   args = fml_incl) 
  cache_pars  <- cache_pars[sort(names(cache_pars))]
  call_pars   <- mget(fml_incl, envir = fun_env, inherits = FALSE)
  call_pars   <- call_pars[sort(names(call_pars))]
  
  if (identical(cache_pars, call_pars)) {
    fions   <- list.files(path = file.path(out_path, "temp"), 
                          pattern = "ion_matches_[0-9]+\\.rds$")

    if (length(fions)) {
      message("Found ", length(fions), " cached ion matches.")

      if (!file.exists(faa))
        qs::qsave(aa_masses_all, faa)
      
      .savecall <- FALSE
      
      return(NULL)
    }
  }
  
  rm(list = c("args_except", "cache_pars", "call_pars"))
  
  delete_files(
    out_path, 
    ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                "\\.xls$", "\\.csv$", "\\.txt$", 
                "^mgf$", "^mgfs$", "Calls"))

  # pairs expts and theos
  files_a  <-  list.files(mgf_path, pattern = "^expttheo_", full.names = TRUE)
  files_b  <-  list.files(mgf_path, pattern = "^mgftheo_",  full.names = TRUE)
  nfiles_a <- length(files_a)
  nfiles_b <- length(files_b)

  if (nfiles_a && nfiles_b) {
    warning("Both cached `expttheo_` and `mgftheo_` under", mgf_path, ".\n", 
            "Deleting ", paste(files_a, collapse = "\n"))
    file.remove(files_a)
    nfiles_a <- 0L
  }

  # For three-frame searches
  # (matches of secondary ions may use `outer` products and no adjustments)
  ppm_ms1_bin <- calc_threeframe_ppm(ppm_ms1)
  ppm_ms2_bin <- calc_threeframe_ppm(ppm_ms2)
  ms1_offsets <- find_ms1_offsets(n_13c = n_13c, ms1_notches = ms1_notches)
  
  pair_mgftheos(mgf_path = mgf_path, n_modules = length(aa_masses_all), 
                ms1_offsets = ms1_offsets, by_modules = by_modules, 
                min_mass = min_mass, max_mass = max_mass, 
                ppm_ms1_bin = ppm_ms1_bin, .path_bin = .path_bin, 
                reframe_mgfs = reframe_mgfs, first_search = first_search)

  rm(list = c("files_a", "files_b", "nfiles_a", "nfiles_b"))

  # MS2 generation functions
  types  <- unlist(lapply(aa_masses_all, attr, "type", exact = TRUE))
  
  funs_ms2 <- lapply(types, function (x) {
    if (x %in% c("amods- tmod- vnl- fnl-", "amods- tmod+ vnl- fnl-"))
      "gen_ms2ions_base"
    else if (x %in% c("amods- tmod- vnl- fnl+", "amods- tmod+ vnl- fnl+"))
      "gen_ms2ions_a0_vnl0_fnl1"
    else if (x %in% c("amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"))
      "gen_ms2ions_a1_vnl0_fnl0"
    else if (x %in% c("amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"))
      "gen_ms2ions_a1_vnl1_fnl0"
    else if (x %in% c("amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"))
      "gen_ms2ions_a1_vnl0_fnl1"
    else
      stop("Unknown modification type.")
  })
  
  # other attributes
  ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
                         maxn_vmods_per_pep = maxn_vmods_per_pep,
                         maxn_sites_per_vmod = maxn_sites_per_vmod)
  ms2vmods_all <- lapply(ms1vmods_all, lapply, make_ms2vmods)
  
  # Searches
  df0 <- tibble::tibble(scan_title = integer(), ms1_moverz = numeric(), 
                        ms1_mass = numeric(), ms1_int = numeric(), 
                        ms1_charge = character(), ret_time = numeric(), 
                        scan_num = character(), raw_file = integer(), 
                        ms2_moverz = list(list()), 
                        ms2_int = list(list()), 
                        ms2_n = integer(), frame  = numeric(), 
                        matches = list(list()), 
                        pep_isdecoy = logical())

  hms2match(aa_masses_all = aa_masses_all, 
            funs_ms2 = funs_ms2, 
            ms1vmods_all = ms1vmods_all, 
            ms2vmods_all = ms2vmods_all, 
            mod_indexes = mod_indexes, 
            mgf_path = mgf_path, 
            out_path = out_path, 
            type_ms2ions = type_ms2ions, 
            maxn_vmods_per_pep = maxn_vmods_per_pep, 
            maxn_sites_per_vmod = maxn_sites_per_vmod, 
            maxn_fnl_per_seq = maxn_fnl_per_seq, 
            maxn_vnl_per_seq = maxn_vnl_per_seq, 
            maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
            minn_ms2 = minn_ms2, 
            ppm_ms1 = ppm_ms1_bin, 
            ppm_ms2 = ppm_ms2_bin, 
            min_ms2mass = min_ms2mass, 
            index_mgf_ms2 = index_mgf_ms2, 
            by_modules = by_modules, 
            df0 = df0)

  qs::qsave(aa_masses_all, faa)
  
  invisible(NULL)
}


#' Helper of \link{reverse_seqs}.
#' 
#' Reverses \code{pep_seq} in a frame.
#' 
#' @param pep_frame A frame of data.
#' 
#' @examples 
#' ## pep_frame
#' #                                  pep_seq      mass  frame    prot_acc
#' # 2148391     CSCNNGEMCDRFQGCLCSPGWQGLQCER 3694.4923 100001   TIE2_HUMAN
#' # 2148392 EGSARASEQPENAESPDNEDGDCEETTEEAGR 3694.5248 100001  TXLNB_HUMAN
reverse_peps_in_frame <- function (pep_frame) 
{
  nms <- names(pep_frame)
  
  if ("pep_seq" %in% nms) 
    pep_frame[["pep_seq"]]  <- reverse_seqs(pep_frame[["pep_seq"]])

  if ("prot_acc" %in% nms) 
    pep_frame[["prot_acc"]] <- paste0("-", pep_frame[["prot_acc"]])

  pep_frame
}


#' Reverses peptide sequences.
#' 
#' @param seqs Lists of peptide sequences.
#' @examples 
#' \donttest{
#' seqs = c(paste0(LETTERS[1:10], collapse = ""), paste0(letters[1:10], collapse = ""))
#' }
reverse_seqs <- function (seqs) 
{
  fis  <- stringi::stri_sub(seqs, 1L, 1L, use_matrix = FALSE)
  las  <- stringi::stri_sub(seqs, -1L, -1L, use_matrix = FALSE)
  lens <- stringi::stri_length(seqs)
  revs <- stringi::stri_reverse(seqs)
  
  substring(revs, 1)    <- fis
  substring(revs, lens) <- las
  
  revs
}


#' MGF precursor mass calibration.
#'
#' \code{ppm_ms1} only for the calculation of frame indexes of precursors.
#'
#' @param aa_masses_all List(1); The first list of all amino-acid look-ups.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @param .path_bin The file path to binned precursor masses.
#' @param reframe_mgfs Logical; if TRUE, recalculates the frame indexes of MGFs
#' @param knots The number of knots for spline fits.
#' @inheritParams matchMS
calib_mgf <- function (mgf_path = NULL, aa_masses_all = NULL,out_path = NULL, 
                       .path_bin, mod_indexes = NULL, type_ms2ions = "by", 
                       maxn_vmods_per_pep = 5L,maxn_sites_per_vmod = 3L, 
                       maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                       maxn_vmods_sitescombi_per_pep = 64L, minn_ms2 = 6L, 
                       ppm_ms1 = 20L, reframe_mgfs = TRUE, 
                       ppm_ms2 = 20L, min_mass = 200L, 
                       max_mass = 4500L, min_ms2mass = 115L, quant = "none", 
                       ppm_reporters = 10L, index_mgf_ms2 = FALSE, 
                       by_modules = TRUE, fasta = NULL, acc_type = NULL, 
                       acc_pattern = NULL, topn_ms2ions = 100L, 
                       fixedmods = NULL, varmods = NULL, enzyme = "trypsin_p", 
                       maxn_fasta_seqs = 200000L, maxn_vmods_setscombi = 512L,
                       min_len = 7L, max_len = 40L, max_miss = 2L, knots = 50L)
{
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) save_call2(path = file.path(out_path, "Calls"), fun = fun)
    }, add = TRUE
  )
  
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  args <- names(formals(fun))
  args_except <- NULL
  args_must <- if (length(args_except)) args[!args %in% args_except] else args

  cache_pars <- find_callarg_vals(
    time = NULL, 
    path = file.path(out_path, "Calls"), 
    fun = paste0(fun, ".rda"), 
    args = args_must)
  
  cache_pars <- cache_pars[sort(names(cache_pars))]
  call_pars  <- mget(args_must, envir = fun_env, inherits = FALSE)
  call_pars  <- call_pars[sort(names(call_pars))]
  ok_pars    <- identical(call_pars, cache_pars)

  if (ok_pars) {
    message("Mass calibration performed previously. ", 
            "Delete `", paste0(fun, ".rda"), "` to recalibrate.")
    .savecall <- FALSE
    return(NULL)
  }
  
  # may need to delete mgf_queries_[...].rds when changing, e.g., from
  #   ppm_ms1 = 20 to 10; or save a copy of the original mgf_queries
  
  ## the first search
  tempdir <- file.path(out_path, "temp")
  pat_th <- if (by_modules) "^expttheo_\\d+.*\\.rds$" else "^mgftheo_\\d+.*\\.rds$"
  pat_im <- "^ion_matches_\\d+.*\\.rds$"
  fs_th <- list.files(mgf_path, pattern = pat_th, full.names = TRUE)
  fs_im <- list.files(tempdir,  pattern = pat_im, full.names = TRUE)
  file.remove(fs_th, recursive = TRUE)
  file.remove(fs_im, recursive = TRUE)
  
  if (!dir.exists(tempdir))
    create_dir(tempdir)

  fi_aa <- file.path(out_path, "aa_masses_all.rds")
  fi_mi <- file.path(out_path, "mod_indexes.txt")

  if (!file.exists(fi_aa))
    stop("Amino-acid look-ups not found: ", fi_aa)
  if (!file.exists(fi_mi))
    stop("Amino-acid look-ups not found: ", fi_mi)

  fi_aa2 <- file.path(out_path, "Calls", "aa_masses_all.rds")
  fi_mi2 <- file.path(out_path, "Calls", "mod_indexes.txt")
  file.rename(fi_aa, fi_aa2)
  file.rename(fi_mi, fi_mi2)

  ms2match(mgf_path = mgf_path,
           aa_masses_all = aa_masses_all,
           out_path = out_path,
           .path_bin = .path_bin, 
           mod_indexes = mod_indexes, 
           type_ms2ions = type_ms2ions,
           maxn_vmods_per_pep = maxn_vmods_per_pep,
           maxn_sites_per_vmod = maxn_sites_per_vmod,
           maxn_fnl_per_seq = maxn_fnl_per_seq, 
           maxn_vnl_per_seq = maxn_vnl_per_seq, 
           maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep,
           minn_ms2 = minn_ms2,
           ppm_ms1 = ppm_ms1,
           ppm_ms2 = ppm_ms2,
           min_mass = min_mass, 
           max_mass = max_mass, 
           min_ms2mass = min_ms2mass,
           quant = quant,
           ppm_reporters = ppm_reporters,
           reframe_mgfs = reframe_mgfs, 
           index_mgf_ms2 = index_mgf_ms2, 
           by_modules = by_modules, 
           fasta = fasta,
           acc_type = acc_type,
           acc_pattern = acc_pattern,
           topn_ms2ions = topn_ms2ions,
           fixedmods = fixedmods,
           varmods = varmods,
           enzyme = enzyme,
           maxn_fasta_seqs = maxn_fasta_seqs,
           maxn_vmods_setscombi = maxn_vmods_setscombi,
           min_len = min_len,
           max_len = max_len,
           max_miss = max_miss,
           first_search = TRUE, 
           .savecall = FALSE)
  
  file.rename(fi_aa2, fi_aa)
  file.rename(fi_mi2, fi_mi)
  
  ## mass calibration
  fs_mgf <- list.files(mgf_path, "^mgf_queries.*\\.rds$")
  fi_ion <- file.path(out_path, "temp", "ion_matches_1.rds")
  
  if (!length(fs_mgf))
    stop("No `mgf_queries` files found for calibrations.")
  if (!file.exists(fi_ion))
    stop("No `ion_matches` files found for calibrations.")
  
  df  <- qs::qread(fi_ion)
  
  if (!"raw_file" %in% names(df))
    stop("Column not found in search results: `raw_file`")

  dfs <- split(df, df[["raw_file"]])
  ord <- sort(as.integer(gsub("^mgf_queries_(\\d+)\\.rds", "\\1", fs_mgf)))
  dfs <- dfs[ord]
  fs_mgf <- fs_mgf[ord]
  rm(list = c("df", "ord"))
  
  len <- length(dfs)
  n_cores <- min(len, detect_cores(32L))
  
  if (len <= 2L) {
    mapply(calib_ms1, fs_mgf, dfs, 
           MoreArgs = list(
             mgf_path = mgf_path, out_path = out_path, ppm_ms1 = ppm_ms1, 
             min_mass = min_mass, max_mass = max_mass, knots = knots), 
           SIMPLIFY = FALSE, 
           USE.NAMES = FALSE)
  }
  else {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterExport(cl, "calib_ms1", envir = environment(mzion::matchMS))
    
    parallel::clusterMap(
      cl, calib_ms1, fs_mgf, dfs, 
      MoreArgs = list(
        mgf_path = mgf_path, out_path = out_path, ppm_ms1 = ppm_ms1, 
        min_mass = min_mass, max_mass = max_mass, knots = knots), 
      SIMPLIFY = FALSE, 
      USE.NAMES = FALSE)
    
    parallel::stopCluster(cl)
  }
  
  message("Completed precursor mass calibration.\n")

  fs_th <- list.files(mgf_path, pattern = pat_th, full.names = TRUE)
  fs_im <- list.files(tempdir,  pattern = pat_im, full.names = TRUE)
  file.remove(fs_th, recursive = TRUE)
  file.remove(fs_im, recursive = TRUE)
  
  .savecall <- TRUE

  invisible(NULL)
}


#' Calibrates precursor masses (by individual RAW_Files)
#' 
#' @param filename An MGF file name
#' @param df A data frame of \code{ion_matches_1.rds}
#' @inheritParams calib_mgf
calib_ms1 <- function (filename, df = NULL, mgf_path = NULL, out_path = NULL, 
                       ppm_ms1 = 20L, min_mass = 200L, max_mass = 4500L, 
                       knots = 50L)
{
  mgfs <- qs::qread(file.path(mgf_path, filename))
  
  # subsets by minn_ms2 and ms1_int
  if (FALSE) {
    ms1int_co <- quantile(df$ms1_int, probs = .25, na.rm = TRUE)
    df <- df[with(df, ms1_int >= ms1int_co), ]
    
    minn_ok <- lapply(df[["matches"]], function (x) {
      m <- x[[1]][[1]][["m"]]
      (!is.null(m)) && m >= 8L
    })
    minn_ok  <- .Internal(unlist(minn_ok, recursive = FALSE, use.names = FALSE))
    
    df <- df[minn_ok, ]
  }
  
  # x[[1]]: no nested structures, as the first search is always 
  # against the all-fixed modifications
  theo_ms1 <- lapply(df[["matches"]], function (x) {
    attr(x[[1]], "theo_ms1", exact = TRUE)
  })
  theo_ms1  <- .Internal(unlist(theo_ms1, recursive = FALSE, use.names = FALSE))
  
  diff_ms1  <- (df[["pep_exp_mr"]] - theo_ms1)/theo_ms1 * 1E6
  ret_time  <- df[["pep_ret_range"]]
  ppm_ms1_bin  <- calc_threeframe_ppm(ppm_ms1)
  
  fit_ns <- tryCatch(
    lm(diff_ms1 ~ splines::ns(ret_time, knots)),
    error = function(e) NA)
  
  fit_bs <- tryCatch(
    lm(diff_ms1 ~ splines::bs(ret_time, knots)),
    error = function(e) NA)
  
  res_ns <- if (class(fit_ns) == "lm") sum(resid(fit_ns)^2, na.rm = TRUE) else Inf
  res_bs <- if (class(fit_bs) == "lm") sum(resid(fit_bs)^2, na.rm = TRUE) else Inf
  fit    <- if (res_ns <= res_bs) fit_ns else fit_bs

  # (keeps the original df$ms1_mass -> can later infer mass deltas)
  # charges <- get_ms1charges(df[["ms1_charge"]])
  # df[["ms1_moverz"]] <- (df[["ms1_mass"]] + 1.00727647 * charges)/charges
  
  ## Update mgf
  rt <- mgfs[["ret_time"]]
  oks_le <- rt >= min(ret_time, na.rm = TRUE)
  oks_gr <- rt <= max(ret_time, na.rm = TRUE)
  oks <- oks_le & oks_gr
  ms1err <- predict.lm(fit, newdata = data.frame(ret_time = rt[oks])) / 1E6
  ms1oks <- mgfs[["ms1_mass"]][oks]
  mgfs[["ms1_mass"]][oks] <- ms1oks - ms1oks * ms1err
  
  # beyond the boundary of RT
  if (FALSE) {
    rts_ok <- rt[oks]
    err_le <- ms1err[which.min(rts_ok)]
    err_gr <- ms1err[which.max(rts_ok)]
    
    ms1_le <- mgfs[["ms1_mass"]][!oks_le]
    ms1_gr <- mgfs[["ms1_mass"]][!oks_gr]
    mgfs[["ms1_mass"]][!oks_le] <- ms1_le - ms1_le * err_le
    mgfs[["ms1_mass"]][!oks_gr] <- ms1_gr - ms1_gr * err_gr
  }
  
  ## update MGF
  mgfs <- mgfs %>%
    dplyr::arrange(ms1_mass) %>% 
    dplyr::filter(ms1_mass >= min_mass, ms1_mass <= max_mass)

  # charges <- get_ms1charges(mgfs[["ms1_charge"]])
  # mgfs[["ms1_moverz"]] <- (mgfs[["ms1_mass"]] + 1.00727647 * charges)/charges
  
  .savecall <- TRUE

  qs::qsave(mgfs, file.path(mgf_path, filename), preset = "fast")
}


#' Finds offsets in precursor masses.
#' 
#' @inheritParams matchMS
find_ms1_offsets <- function (n_13c = 0L, ms1_notches = 0) 
{
  if (length(n_13c))
    n_13c <- n_13c[n_13c != 0L]
  
  offsets_13c <- if (length(n_13c)) n_13c * 1.00335483 else NULL
  ms1_offsets <- unique(c(offsets_13c, ms1_notches))
  round(ms1_offsets, digits = 4L)
}


