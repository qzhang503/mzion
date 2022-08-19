#' Matches theoretical peptides (parallel by mgf chunks).
#'
#' All files under `out_path` are removed if incur \code{calc_pepmasses} in the
#' upstream.
#'
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @param use_first_rev Logical; if TRUE, uses the first set of \code{aa_masses}
#'   for reversed searches; otherwise, uses the set with the maximum output.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @inheritParams frames_adv
#' @inheritParams add_var_masses
#' @import parallel
ms2match <- function (mgf_path, aa_masses_all, out_path, 
                      mod_indexes, type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                      maxn_sites_per_vmod = 3L, maxn_vmods_sitescombi_per_pep = 64L, 
                      minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                      min_ms2mass = 115L, quant = "none", ppm_reporters = 10L, 
                      use_first_rev = FALSE, 

                      # dummies
                      fasta, acc_type, acc_pattern,
                      topn_ms2ions, fixedmods, varmods, 
                      enzyme, 
                      maxn_fasta_seqs, maxn_vmods_setscombi, 
                      min_len, max_len, max_miss, 

                      digits) 
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

  # (OK as `use_first_rev` is not for users)
  args_except <- c("use_first_rev", "quant")
  fml_incl <- fml_nms[!fml_nms %in% args_except]
  
  cache_pars <- find_callarg_vals(time = NULL, 
                                  path = file.path(out_path, "Calls"), 
                                  fun = paste0(fun, ".rda"), 
                                  args = fml_incl) 
  cache_pars <- cache_pars[sort(names(cache_pars))]
  
  call_pars <- mget(fml_incl, envir = fun_env, inherits = FALSE)
  call_pars <- call_pars[sort(names(call_pars))]
  
  if (identical(cache_pars, call_pars)) {
    len <- length(aa_masses_all)
    
    fions <- list.files(path = file.path(out_path, "temp"), 
                        pattern = "ion_matches_[0-9]+\\.rds$")
    
    ok_ions <- (length(fions) == len)
    
    if (grepl("^tmt[0-9]+$", quant)) {
      ftmt <- list.files(path = file.path(out_path, "temp"), 
                          pattern = "reporters_[0-9]+\\.rds$")
      
      ok_tmt <- (length(ftmt) == len)
    } 
    else {
      ok_tmt <- TRUE
    }
    
    if (ok_ions && ok_tmt) {
      message("Found cached ion matches.")
      .savecall <- FALSE
      
      return(NULL)
    }
  }
  
  rm(list = c("args_except", "cache_pars", "call_pars"))
  
  delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                     "\\.xls$", "\\.csv$", "\\.txt$", 
                                     "^mgf$", "^mgfs$", "Calls"))
  
  out0 <- tibble::tibble(scan_title = integer(), ms1_moverz = numeric(), 
                         ms1_mass = numeric(), ms1_int = numeric(), 
                         ms1_charge = character(), ret_time = numeric(), 
                         scan_num = character(), raw_file = integer(), 
                         pep_mod_group = character(), 
                         ms2_moverz = list(list()), 
                         ms2_int = list(list()), 
                         ms2_n = integer(), frame  = numeric(), 
                         matches = list(list()), 
                         pep_fmod = character(), pep_vmod = character(), 
                         pep_isdecoy = logical())
  
  # For three-frame searches
  # (matches of secondary ions using `outer` and no adjustments)
  is_ms2_three_frame <- is_ms1_three_frame <- TRUE
  
  ppm_ms1_new <- if (is_ms1_three_frame) 
    as.integer(ceiling(ppm_ms1 * .5))
  else 
    ppm_ms1
  
  ppm_ms2_new <- if (is_ms2_three_frame) 
    as.integer(ceiling(ppm_ms2 * .5))
  else 
    ppm_ms2
  
  ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
                         maxn_vmods_per_pep = maxn_vmods_per_pep,
                         maxn_sites_per_vmod = maxn_sites_per_vmod)
  
  ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
  
  
  message("\n===  MS2 ion searches started at ", Sys.time(), ". ===\n")
  
  ## Targets 
  obj_sizes <- numeric(length(aa_masses_all))
  types <- purrr::map_chr(aa_masses_all, attr, "type", exact = TRUE)
  
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl-", 
                             "amods- tmod+ vnl- fnl-"))
  
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      ms1vmods <- ms1vmods_all[[i]]
      ms2vmods <- ms2vmods_all[[i]]
      
      # need ntmod and ctmod for `amod+_...` for excluding 
      #   additive terminal mod and anywhere mod
      # 
      # aa_masses["N-term"] not necessary H; 
      #   e.g., `Acetyl + hydrogen` in case of FIXED Protein N-term 
      #   or fixed N-term `TMT + hydrogen`
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      ntmass <- find_nterm_mass(aa_masses)
      ctmass <- find_cterm_mass(aa_masses)
      
      out <- ms2match_base(
        i = i, 
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1_new, 
        ppm_ms2 = ppm_ms2_new, 
        min_ms2mass = min_ms2mass, 
        df0 = out0, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          qs::qsave(file.path(out_path, "temp", paste0("reporters_", i, ".rds")), 
                    preset = "fast")
      }
      
      rm(list = c("out", "aa_masses", "ms1vmods", "ms2vmods", 
                  "ntmod", "ntmass", "ctmod", "ctmass"))
      gc()
    }
  }
  
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  #        (mutual exclusive btw. (1, 2) and (5, 6)
  #         "ANY" fmod has neuloss -> 5, 6;
  #         "ALL" fmods have no neuloss -> 1, 2)
  
  inds <- which(types %in% c("amods- tmod- vnl- fnl+", 
                             "amods- tmod+ vnl- fnl+"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      ms1vmods <- ms1vmods_all[[i]]
      ms2vmods <- ms2vmods_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      ntmass <- find_nterm_mass(aa_masses)
      ctmass <- find_cterm_mass(aa_masses)
      
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a0_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1_new, 
        ppm_ms2 = ppm_ms2_new, 
        min_ms2mass = min_ms2mass, 
        df0 = out0, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          qs::qsave(file.path(out_path, "temp", paste0("reporters_", i, ".rds")), 
                    preset = "fast")
      }
      
      rm(list = c("out", "aa_masses", "ms1vmods", "ms2vmods", 
                  "ntmod", "ntmass", "ctmod", "ctmass"))
      gc()
    }
  }
  
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #        (ALL amods are vnl-)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl-", 
                             "amods+ tmod+ vnl- fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      ms1vmods <- ms1vmods_all[[i]]
      ms2vmods <- ms2vmods_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      ntmass <- find_nterm_mass(aa_masses)
      ctmass <- find_cterm_mass(aa_masses)
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      
      out <- ms2match_a1_vnl0_fnl0(
        i = i, 
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1_new, 
        ppm_ms2 = ppm_ms2_new, 
        min_ms2mass = min_ms2mass, 
        df0 = out0, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          qs::qsave(file.path(out_path, "temp", paste0("reporters_", i, ".rds")), 
                    preset = "fast")
      }
      
      rm(list = c("out", "aa_masses", "ms1vmods", "ms2vmods", 
                  "ntmod", "ntmass", "ctmod", "ctmass"))
      gc()
    }
  }
  
  
  # (9, 10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #         (ANY amod is vnl+)
  
  inds <- which(types %in% c("amods+ tmod- vnl+ fnl-", 
                             "amods+ tmod+ vnl+ fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      ms1vmods <- ms1vmods_all[[i]]
      ms2vmods <- ms2vmods_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      ntmass <- find_nterm_mass(aa_masses)
      ctmass <- find_cterm_mass(aa_masses)
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl1_fnl0(
        i = i, 
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        vmods_nl = vmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1_new, 
        ppm_ms2 = ppm_ms2_new, 
        min_ms2mass = min_ms2mass, 
        df0 = out0, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          qs::qsave(file.path(out_path, "temp", paste0("reporters_", i, ".rds")), 
                    preset = "fast")
      }
      
      rm(list = c("out", "aa_masses", "ms1vmods", "ms2vmods", 
                  "ntmod", "ntmass", "ctmod", "ctmass"))
      gc()
    }
  }
  
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  #          (mutual exclusive btw. (11, 12) and (7, 8);
  #           logicial ANY versus ALL)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl+", 
                             "amods+ tmod+ vnl- fnl+"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      ms1vmods <- ms1vmods_all[[i]]
      ms2vmods <- ms2vmods_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      ntmass <- find_nterm_mass(aa_masses)
      ctmass <- find_cterm_mass(aa_masses)
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
        ms1vmods = ms1vmods, 
        ms2vmods = ms2vmods, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1_new, 
        ppm_ms2 = ppm_ms2_new, 
        min_ms2mass = min_ms2mass, 
        df0 = out0, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          qs::qsave(file.path(out_path, "temp", paste0("reporters_", i, ".rds")), 
                    preset = "fast")
      }
      
      rm(list = c("out", "aa_masses", "ms1vmods", "ms2vmods", 
                  "ntmod", "ntmass", "ctmod", "ctmass"))
      gc()
    }
  }
  
  ## Decoys
  # (1) makes binned_theopeps_rev_[i_max].rds
  i_max <- if (use_first_rev) 1L else which.max(obj_sizes)
  i_max2 <- paste0("rev_", i_max)

  .path_bin <- get(".path_bin", envir = .GlobalEnv, inherits = FALSE)
  bin_file <- file.path(.path_bin, paste0("binned_theopeps_", i_max, ".rds"))
  bin_file2 <- file.path(.path_bin, paste0("binned_theopeps_", i_max2, ".rds"))
  
  if (!file.exists(bin_file2)) {
    rev_peps <- 
      qs::qread(bin_file) %>% 
      lapply(reverse_peps_in_frame) %T>% 
      qs::qsave(bin_file2, preset = "fast")

    rm(list = c("rev_peps"))
  }
  
  rm(list = c("bin_file", "bin_file2"))
  
  # (2) makes MS2 ions 
  aa_masses <- aa_masses_all[[i_max]]
  ms1vmods <- ms1vmods_all[[i_max]]
  ms2vmods <- ms2vmods_all[[i_max]]

  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  ntmass <- find_nterm_mass(aa_masses)
  ctmass <- find_cterm_mass(aa_masses)

  amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
  
  if (length(amods)) { # (7, 8)
    out <- ms2match_a1_vnl0_fnl0(
      i = i_max2, 
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      mod_indexes = mod_indexes, 
      mgf_path = mgf_path, 
      out_path = out_path, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      minn_ms2 = minn_ms2, 
      ppm_ms1 = ppm_ms1_new, 
      ppm_ms2 = ppm_ms2_new, 
      min_ms2mass = min_ms2mass, 
      df0 = out0, 
      digits = digits)
  } 
  else { # (1, 2)
    out <- ms2match_base(
      i = i_max2, 
      aa_masses = aa_masses, 
      ms1vmods = ms1vmods, 
      ms2vmods = ms2vmods, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      mod_indexes = mod_indexes, 
      mgf_path = mgf_path, 
      out_path = out_path, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      minn_ms2 = minn_ms2, 
      ppm_ms1 = ppm_ms1_new, 
      ppm_ms2 = ppm_ms2_new, 
      min_ms2mass = min_ms2mass, 
      df0 = out0, 
      digits = digits)
  }
  
  # if (is.null(out)) out <- out0
  qs::qsave(out, file.path(out_path, "temp", paste0("ion_matches_", i_max2, ".rds")), 
            preset = "fast") 
  
  .savecall <- TRUE
  
  message("\n===  MS2 ion searches completed at ", Sys.time(), ". ===\n")

  invisible(NULL)
}


#' Helper of \link{reverse_seqs}.
#' 
#' Reverses `pep_seq` in a frame.
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
    pep_frame[["pep_seq"]] <- reverse_seqs(pep_frame[["pep_seq"]])

  if ("prot_acc" %in% nms) 
    pep_frame[["prot_acc"]] <- paste0("-", pep_frame[["prot_acc"]])

  pep_frame
}


#' Reverses peptide sequences.
#' 
#' @param seqs Lists of peptide sequences.
reverse_seqs <- function (seqs) 
{
  fis <- stringi::stri_sub(seqs, 1, 1)
  las <- stringi::stri_sub(seqs, -1, -1)
  lens <- stringi::stri_length(seqs)
  
  revs <- stringi::stri_reverse(seqs)
  substring(revs, 1) <- fis
  substring(revs, lens) <- las
  
  revs
}

