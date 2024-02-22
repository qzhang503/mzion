#' Matches theoretical peptides (parallel by mgf chunks).
#'
#' All files under `out_path` are removed if incur \code{calc_pepmasses} in the
#' upstream.
#'
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @param ms1_offsets The MS1 off-sets.
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
                      maxn_sites_per_vmod = 3L, maxn_fnl_per_seq = 3L, 
                      maxn_vnl_per_seq = 3L, maxn_vmods_sitescombi_per_pep = 64L, 
                      minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 20L, 
                      min_mass = 200L, max_mass = 4500L, min_ms2mass = 115L, 
                      quant = "none", ppm_reporters = 10L, 
                      by_modules = TRUE, reframe_mgfs = FALSE, ms1_offsets = 0, 
                      ms1_neulosses = NULL, maxn_neulosses_fnl = 1L, 
                      maxn_neulosses_vnl = 1L, deisotope_ms2 = TRUE, 

                      # dummies
                      fasta, acc_type, acc_pattern, topn_ms2ions, fixedmods, 
                      varmods, enzyme, maxn_fasta_seqs, maxn_vmods_setscombi, 
                      min_len, max_len, max_miss, 
                      first_search = FALSE, .savecall = TRUE) 
                      
{
  options(digits = 9L)
  
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    }, add = TRUE)

  # Check cached 
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  fml_nms <- names(formals(fun))
  faa <- file.path(out_path, "aa_masses_all.rds")
  
  # (OK as `argument` not for users)
  # min_mass and max_mass only for calib_ms1mass, not to be changed by users
  # args_except <- c("ms1_offsets")
  args_except <- c("by_modules")
  fml_incl    <- fml_nms[!fml_nms %in% args_except]
  cache_pars  <- find_callarg_vals(time = NULL, 
                                   path = file.path(out_path, "Calls"), 
                                   fun = paste0(fun, ".rda"), 
                                   args = fml_incl) 
  cache_pars  <- cache_pars[sort(names(cache_pars))]
  call_pars   <- mget(fml_incl, envir = fun_env, inherits = FALSE)
  call_pars   <- call_pars[sort(names(call_pars))]
  
  # backward compatibility of old cached parameters
  if (".path_bin" %in% names(cache_pars) && ".path_bin" %in% names(call_pars)) {
    if (!(is.null(cache_pars$.path_bin) || 
          "fs_path" %in% class(cache_pars$.path_bin))) {
      cache_pars$.path_bin <- fs::fs_path(cache_pars$.path_bin)
    }
    
    if (!(is.null(call_pars$.path_bin) || 
          "fs_path" %in% class(call_pars$.path_bin))) {
      call_pars$.path_bin <- fs::fs_path(call_pars$.path_bin)
    }
  }
  
  if (identical(cache_pars, call_pars)) {
    fions   <- list.files(path = file.path(out_path, "temp"), 
                          pattern = "ion_matches_[0-9]+\\.rds$")

    if (n_fi <- length(fions)) {
      message("Found ", n_fi, " cached ion matches.")

      if (!file.exists(faa))
        qs::qsave(aa_masses_all, faa)
      
      .savecall <- FALSE
      
      return(NULL)
    }
  }
  
  rm(list = c("args_except", "cache_pars", "call_pars"))
  
  delete_files(
    out_path, 
    ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.(mzML|mzml)$", 
                "\\.xlsx$", "\\.xls$", "\\.csv$", "\\.txt$", "\\.pars$", 
                "^mgf$", "^mgfs$", "^mzML$", "^mzMLs$", 
                "Calls", "^PSM$", "^Peptide$", "^Protein$", 
                "fraction_scheme.rda", "label_scheme.rda", 
                "label_scheme_full.rda"), 
    paths_excluded = mgf_path)

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
  
  pair_mgftheos(mgf_path = mgf_path, n_modules = length(aa_masses_all), 
                ms1_offsets = comb_ms1_offsets(ms1_offsets = ms1_offsets, 
                                               ms1_neulosses = ms1_neulosses), 
                quant = quant, by_modules = by_modules, 
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
  
  df0 <- tibble::tibble(
    scan_title = integer(), raw_file = integer(), 
    pep_mod_group = integer(), pep_exp_mz = numeric(), 
    pep_exp_mr = numeric(), pep_tot_int = numeric(), 
    pep_exp_z = numeric(), pep_ret_range = numeric(), 
    pep_scan_num = character(), 
    pep_ms2_moverzs = list(list()), 
    pep_ms2_ints = list(list()), 
    pep_n_ms2 = integer(), 
    rptr_moverz = list(list()), 
    rptr_int = list(list()), 
    pep_ms1_offset = numeric(), 
    matches = list(list()), 
    pep_fmod = character(), 
    pep_vmod = character())

  hms2match(
    aa_masses_all = aa_masses_all, 
    funs_ms2 = funs_ms2, 
    ms1vmods_all = ms1vmods_all, 
    ms2vmods_all = ms2vmods_all, 
    ms1_neulosses = ms1_neulosses, 
    maxn_neulosses_fnl = maxn_neulosses_fnl, 
    maxn_neulosses_vnl = maxn_neulosses_vnl, 
    deisotope_ms2 = deisotope_ms2, 
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
#' @inheritParams matchMS
calib_mgf <- function (mgf_path, aa_masses_all, out_path, .path_bin, 
                       mod_indexes = NULL, type_ms2ions = "by", 
                       maxn_vmods_per_pep = 5L,maxn_sites_per_vmod = 3L, 
                       maxn_fnl_per_seq = 3L, maxn_vnl_per_seq = 3L, 
                       maxn_vmods_sitescombi_per_pep = 64L, minn_ms2 = 6L, 
                       ppm_ms1 = 20L, reframe_mgfs = TRUE, 
                       ppm_ms2 = 20L, min_mass = 200L, 
                       max_mass = 4500L, min_ms2mass = 115L, quant = "none", 
                       ppm_reporters = 10L, 
                       by_modules = TRUE, fasta = NULL, acc_type = NULL, 
                       acc_pattern = NULL, topn_ms2ions = 150L, 
                       fixedmods = NULL, varmods = NULL, enzyme = "trypsin_p", 
                       maxn_fasta_seqs = 200000L, maxn_vmods_setscombi = 512L,
                       min_len = 7L, max_len = 40L, max_miss = 2L)
{
  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) save_call2(path = file.path(out_path, "Calls"), fun = fun)
    }, add = TRUE)

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
  
  ## the first search
  tempdir <- file.path(out_path, "temp")
  pat_th <- if (by_modules) "^(theo|expt)_\\d+.*\\.rds$" else "^mgftheo_\\d+.*\\.rds$"
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
           maxn_vmods_per_pep = 1L,
           maxn_sites_per_vmod = 1L,
           maxn_fnl_per_seq = 1L, 
           maxn_vnl_per_seq = 1L, 
           ms1_offsets = 0, 
           ms1_neulosses = NULL, 
           maxn_neulosses_fnl = 1L, 
           maxn_neulosses_vnl = 1L, 
           maxn_vmods_sitescombi_per_pep = 1L,
           minn_ms2 = minn_ms2,
           ppm_ms1 = ppm_ms1,
           ppm_ms2 = ppm_ms2,
           min_mass = min_mass, 
           max_mass = max_mass, 
           min_ms2mass = min_ms2mass,
           quant = "none",
           ppm_reporters = ppm_reporters,
           reframe_mgfs = reframe_mgfs, 
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
  
  calc_pepscores(topn_ms2ions = topn_ms2ions,
                 type_ms2ions = type_ms2ions,
                 target_fdr = .01,
                 min_len = min_len,
                 max_len = max_len,
                 ppm_ms2 = ppm_ms2,
                 soft_secions = FALSE, 
                 out_path = out_path,
                 min_ms2mass = min_ms2mass,
                 tally_ms2ints = TRUE, 
                 
                 # dummies
                 mgf_path = mgf_path,
                 maxn_vmods_per_pep = 1L,
                 maxn_sites_per_vmod = 1L,
                 maxn_vmods_sitescombi_per_pep = 1L,
                 minn_ms2 = minn_ms2,
                 ppm_ms1 = ppm_ms1,
                 quant = quant,
                 ppm_reporters = ppm_reporters,
                 fasta = fasta,
                 acc_type = acc_type,
                 acc_pattern = acc_pattern,
                 fixedmods = fixedmods,
                 varmods = varmods,
                 enzyme = enzyme,
                 maxn_fasta_seqs = maxn_fasta_seqs,
                 maxn_vmods_setscombi = maxn_vmods_setscombi,
                 add_ms2theos = FALSE, 
                 add_ms2theos2 = FALSE, 
                 add_ms2moverzs = FALSE, 
                 add_ms2ints = FALSE,
                 by_modules = by_modules, 
                 digits = 4L)
  
  file.rename(fi_aa2, fi_aa)
  file.rename(fi_mi2, fi_mi)
  
  ## mass calibration
  fs_mgf <- list.files(mgf_path, "^mgf_queries_.*\\.rds$")
  fi_ion <- file.path(out_path, "temp", "prescores_1_1.rds")

  if (!length(fs_mgf))
    stop("No `mgf_queries` files found for calibrations.")
  if (!file.exists(fi_ion))
    stop("No `ion_matches` files found for calibrations.")
  
  df <- qs::qread(fi_ion)
  df <- df[with(df, pep_prob <= .01), ]
  
  if (!"raw_file" %in% names(df))
    stop("Column not found in search results: `raw_file`")

  dfs <- split(df, df[["raw_file"]])
  raws <- qs::qread(file.path(mgf_path, "raw_indexes.rds"))
  ord <- match(names(dfs), as.character(raws))
  names(dfs) <- names(raws)[ord]

  ord_mgf <- match(gsub("^mgf_queries_(.*)\\.rds$", "\\1", fs_mgf), names(dfs))
  fs_mgf <- fs_mgf[ord_mgf]

  len <- length(dfs)
  n_cores <- min(len, detect_cores(32L))
  
  if (len <= 2L) {
    mapply(calib_ms1, fs_mgf, dfs, 
           MoreArgs = list(
             mgf_path = mgf_path, out_path = out_path, ppm_ms1 = ppm_ms1, 
             min_mass = min_mass, max_mass = max_mass), 
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
        min_mass = min_mass, max_mass = max_mass), 
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
#' @param range The range of spine knots
#' @inheritParams calib_mgf
calib_ms1 <- function (filename, df = NULL, mgf_path = NULL, out_path = NULL, 
                       ppm_ms1 = 20L, min_mass = 200L, max_mass = 4500L, 
                       range = 3:6)
{
  n_row <- nrow(df)
  mgfs  <- qs::qread(file.path(mgf_path, filename))
  
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
  
  diff_ms1 <- (df[["pep_exp_mr"]] - df[["theo_ms1"]])/df[["theo_ms1"]] * 1E6
  mdiff <- median(diff_ms1, na.rm = TRUE)/1E6

  if (n_row <= 100L || mdiff <= 1e-6) {
    mgfs[["ms1_mass"]] <- mgfs[["ms1_mass"]] - mdiff
    post_calib(mgfs, min_mass, max_mass, mgf_path, filename)
    .savecall <- TRUE
    return(NULL)
  }
  else {
    cvs <- lapply(range, cv_ms1err, k = 10, df = df) 
    cvs <- unlist(cvs, recursive = FALSE, use.names = FALSE)
    # stopifnot(length(cvs) == length(range))
    
    if (all(is.na(cvs))) {
      mgfs[["ms1_mass"]] <- mgfs[["ms1_mass"]] - mdiff
      post_calib(mgfs, min_mass, max_mass, mgf_path, filename)
      .savecall <- TRUE
      return(NULL)
    }
    
    cvs[is.na(cvs)] <- Inf
    knots <- range[which.min(cvs)]
    
    if (length(knots) > 1L)
      knots <- knots[[1]]
  }
  
  ret_time <- df[["pep_ret_range"]]
  fit_ns <- lm(diff_ms1 ~ splines::ns(ret_time, knots))
  fit_bs <- lm(diff_ms1 ~ splines::bs(ret_time, knots))

  if (all(is.na(fit_ns))) {
    if (all(is.na(fit_bs))) {
      mgfs[["ms1_mass"]] <- mgfs[["ms1_mass"]] - mdiff
      post_calib(mgfs, min_mass, max_mass, mgf_path, filename)
      .savecall <- TRUE
      return(NULL)
    }
    else
      fit_ns <- fit_bs
  }
  else if (all(is.na(fit_bs)))
    fit_bs <- fit_ns
  
  bad_ns <- anyNA(coef(fit_ns))
  bad_bs <- anyNA(coef(fit_bs))
  
  if (bad_ns) {
    if (bad_bs) {
      mgfs[["ms1_mass"]] <- mgfs[["ms1_mass"]] - mdiff
      post_calib(mgfs, min_mass, max_mass, mgf_path, filename)
      .savecall <- TRUE
      return(NULL)
    }
    else
      fit_ns <- fit_bs
  }
  else if (bad_bs)
    fit_bs <- fit_ns

  res_ns <- if (class(fit_ns) == "lm") mean(resid(fit_ns)^2, na.rm = TRUE) else Inf
  res_bs <- if (class(fit_bs) == "lm") mean(resid(fit_bs)^2, na.rm = TRUE) else Inf
  fit    <- if (res_ns <= res_bs) fit_ns else fit_bs
  
  # (keeps the original df$ms1_mass -> can later infer mass deltas)
  # charges <- get_ms1charges(df[["ms1_charge"]])
  # df[["ms1_moverz"]] <- (df[["ms1_mass"]] + 1.00727647 * charges)/charges
  
  ## Update mgf
  rt <- mgfs[["ret_time"]]
  min_rt <- min(ret_time, na.rm = TRUE)
  max_rt <- max(ret_time, na.rm = TRUE)
  
  oks <- if ((min_rt2 <- min_rt * 1.1) < (max_rt2 <- max_rt / 1.2)) {
    rt >= min_rt2 & rt <= max_rt2
  }
  else {
    rt >= min_rt & rt <= max_rt
  }

  ms1err <- predict.lm(fit, newdata = data.frame(ret_time = rt[oks])) / 1E6
  mgfs[["ms1_mass"]][oks] <- mgfs[["ms1_mass"]][oks] * (1 - ms1err)

  uls <- rt[!oks]
  mgfs[["ms1_mass"]][uls] <- mgfs[["ms1_mass"]][uls] - mdiff
  
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
  # charges <- get_ms1charges(mgfs[["ms1_charge"]])
  # mgfs[["ms1_moverz"]] <- (mgfs[["ms1_mass"]] + 1.00727647 * charges)/charges
  post_calib(mgfs, min_mass, max_mass, mgf_path, filename)
  .savecall <- TRUE

  invisible(NULL)
}


#' Cross-validation of mass error at a given number of knots
#' 
#' @param df A data frame of search results.
#' @param m The number of knots for fitting.
#' @param k The fold of cross-valications
cv_ms1err <- function(m = 3L, k = 5L, df)
{
  if (!is.data.frame(df))
    stop("Input is a not data frame.")
  
  if ((nr <- nrow(df)) < 50L)
    return(NA_real_)
  
  cv_errs <- vector("list", k)
  folds <- create_folds(seq_len(nr), k = k)
  
  for (i in seq_len(k)) {
    fdi <- folds[[i]]
    tei <- df[fdi, ]
    tri <- df[-fdi, ]
    
    ret_time <- tei[["pep_ret_range"]]
    deltas <- (tei[["pep_exp_mr"]] - tei[["theo_ms1"]])/tei[["theo_ms1"]] * 1E6
    fit <- lm(deltas ~ splines::ns(ret_time, m))
    
    if (all(is.na(fit)) || anyNA(coef(fit)))
      cv_errs[i] <- list(NULL)
    
    tri_rt <- tri[["pep_ret_range"]]
    min_rt <- min(tri_rt, na.rm = TRUE)
    max_rt <- max(tri_rt, na.rm = TRUE)
    
    oks <- if ((min_rt2 <- min_rt * 1.1) < (max_rt2 <- max_rt / 1.2)) {
      tri_rt >= min_rt2 & tri_rt <= max_rt2
    }
    else {
      tri_rt >= min_rt & tri_rt <= max_rt
    }
    
    tri_rt <- tri_rt[oks]
    prd <- predict.lm(fit, newdata = data.frame(ret_time = tri_rt)) / 1E6
    cv_errs[[i]] <- mean(prd^2, na.rm = TRUE)
  }
  
  # return NA at all NUll
  mean(unlist(cv_errs, recursive = FALSE, use.names = FALSE), na.rm = TRUE)
}


#' Post MS1 calibrations.
#' 
#' @param mgfs MGF data.
#' @param filename An MGF file name
#' @inheritParams matchMS
post_calib <- function (mgfs, min_mass, max_mass, mgf_path, filename)
{
  mgfs <- mgfs |>
    dplyr::arrange(ms1_mass) |> 
    dplyr::filter(ms1_mass >= min_mass, ms1_mass <= max_mass)
  
  qs::qsave(mgfs, file.path(mgf_path, filename), preset = "fast")
}


#' Finds off-sets in precursor masses.
#'
#' @inheritParams matchMS
#' @return A vector of mass off-sets.
find_ms1_offsets <- function (n_13c = 0L, ms1_notches = 0) 
{
  if (dups <- anyDuplicated(n_13c)) {
    warning("At least one duplicated value in `n_13c`: ", n_13c[dups])
    n_13c <- unique(n_13c)
  }

  if (dups <- anyDuplicated(ms1_notches)) {
    warning("At least one duplicated values in `ms1_offsets`: ", 
            ms1_notches[dups])
    ms1_notches <- unique(ms1_notches)
  }

  offsets_13c <- if (length(n_13c)) n_13c * 1.00335483 else NULL
  ms1_offsets <- unique(c(0, offsets_13c, ms1_notches))
  ms1_offsets <- round(ms1_offsets, digits = 4L)
}


#' Combines off-sets in precursor masses (notches and neutral losses).
#'
#' The return contains no information of Unimod titles and positions since it is
#' only used for pairing with experimental MGF data.
#'
#' @param ms1_offsets Precursor mass off-sets (notches, not neutral losses).
#' @inheritParams matchMS
#' @return A vector of mass off-sets.
comb_ms1_offsets <- function (ms1_offsets = 0, ms1_neulosses = NULL)
{
  if (dups <- anyDuplicated(ms1_neulosses))
    stop("At least one duplicated values in `ms1_neulosses`: ", 
         ms1_neulosses[dups])

  if ((!(nnl <- length(ms1_neulosses))) || (nnl == 1L && nnl == 0)) 
    return(ms1_offsets)
  
  nls <- extract_umods(ms1_neulosses)
  nls <- unique(unlist(lapply(nls, `[[`, "nl")))
  nls <- -nls[nls != 0]
  
  round(unique(c(ms1_offsets, nls)), digits = 4L)
}


