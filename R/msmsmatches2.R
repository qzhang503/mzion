#' Matches theoretical peptides (parallel by mgf chunks).
#'
#' All files under `out_path` are removed if incur \code{calc_pepmasses} in the
#' upstream.
#' 
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @inheritParams matchMS
#' @inheritParams load_mgfs
#' @inheritParams hms2_base
#' @inheritParams add_fixvar_masses
#' @import parallel
ms2match <- function (mgf_path, aa_masses_all, out_path, 
                      mod_indexes, type_ms2ions, maxn_vmods_per_pep, 
                      maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                      minn_ms2, ppm_ms1, ppm_ms2, min_ms2mass, 
                      quant, ppm_reporters, 

                      # dummies
                      fasta, acc_type, acc_pattern,
                      topn_ms2ions, fixedmods, varmods, 
                      include_insource_nl, enzyme, 
                      maxn_fasta_seqs, maxn_vmods_setscombi, 
                      min_len, max_len, max_miss, 

                      digits) {
  
  options(digits = 9L)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = fun)
      }
    }, 
    add = TRUE
  )
  
  # Check cached 
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
  
  if (identical(cache_pars, call_pars)) {
    len <- length(aa_masses_all)
    
    fions <- list.files(path = file.path(out_path, "temp"), 
                        pattern = "ion_matches_[0-9]+\\.rds$")
    
    ok_ions <- (length(fions) == len)
    
    if (grepl("^tmt[0-9]+$", quant)) {
      ftmt <- list.files(path = file.path(out_path, "temp"), 
                          pattern = "reporters_[0-9]+\\.rds$")
      
      ok_tmt <- (length(ftmt) == len)
    } else {
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
  
  # For three-frame searches
  # (matches of secondary ions using `outer` and no adjustments)
  is_ms2_three_frame <- is_ms1_three_frame <- TRUE

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
      
      # need ntmod and ctmod for `amod+_...` for 
      #   excluding additive terminal mod and anywhere mod
      
      # aa_masses["N-term"] not necessary H; 
      #   e.g., `Acetyl + hydrogen` in case of FIXED Protein N-term 
      #   or fixed N-term `TMT + hydrogen`
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod)) {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647 # + proton
      } else {
        ntmass <- aa_masses["N-term"] - 0.000549 # - electron
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod)) {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      } else {
        ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
      }
      
      # (`map` against groups of frames)
      out <- ms2match_base(
        i = i, 
        aa_masses = aa_masses, 
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
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
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
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod)) {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      } else {
        ntmass <- aa_masses["N-term"] - 0.000549
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod)) {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      } else {
        ctmass <- aa_masses["C-term"] + 2.01510147
      }
      
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a0_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
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
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
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
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod)) {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      } else {
        ntmass <- aa_masses["N-term"] - 0.000549
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod)) {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      } else {
        ctmass <- aa_masses["C-term"] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      
      out <- ms2match_a1_vnl0_fnl0(
        i = i, 
        aa_masses = aa_masses, 
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
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
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
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod)) {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      } else {
        ntmass <- aa_masses["N-term"] - 0.000549
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod)) {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      } else {
        ctmass <- aa_masses["C-term"] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl1_fnl0(
        i = i, 
        aa_masses = aa_masses, 
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
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
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
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod)) {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      } else {
        ntmass <- aa_masses["N-term"] - 0.000549
      }

      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod)) {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      } else {
        ctmass <- aa_masses["C-term"] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
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
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  ## Decoys
  # (1) makes binned_theopeps_rev_[i_max].rds
  i_max <- which.max(obj_sizes)
  i_max2 <- paste0("rev_", i_max)
  
  .path_bin <- get(".path_bin", envir = .GlobalEnv, inherits = FALSE)

  bin_file <- file.path(.path_bin, paste0("binned_theopeps_", i_max, ".rds"))
  bin_file2 <- file.path(.path_bin, paste0("binned_theopeps_", i_max2, ".rds"))
  
  if (!file.exists(bin_file2)) {
    rev_peps <- readRDS(bin_file) %>% 
      lapply(reverse_peps_in_frame) %T>% 
      saveRDS(bin_file2)
    
    rm(list = c("rev_peps"))
  }
  
  rm(list = c("bin_file", "bin_file2"))
  
  # (2) makes MS2 ions 
  aa_masses <- aa_masses_all[[i_max]]
  
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  if (length(ntmod)) {
    ntmass <- aa_masses[names(ntmod)] + 1.00727647
  } else {
    ntmass <- aa_masses["N-term"] - 0.000549
  }
  
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  if (length(ctmod)) {
    ctmass <- aa_masses[names(ctmod)] + 2.01510147
  } else {
    ctmass <- aa_masses["C-term"] + 2.01510147
  }
  
  amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
  
  if (length(amods)) { # (7, 8)
    rev <- ms2match_a1_vnl0_fnl0(
      i = i_max2, 
      aa_masses = aa_masses, 
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
      digits = digits)
  } else { # (1, 2)
    rev <- ms2match_base(
      i = i_max2, 
      aa_masses = aa_masses, 
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
      digits = digits)
  }
  
  saveRDS(rev, file.path(out_path, "temp", paste0("ion_matches_", i_max2, ".rds"))) 
  
  rm(list = c("rev"))
  gc()
  
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
reverse_peps_in_frame <- function (pep_frame) {
  
  nms <- names(pep_frame)
  
  if ("pep_seq" %in% nms) {
    pep_frame[["pep_seq"]] <- reverse_seqs(pep_frame[["pep_seq"]])
  }
  
  if ("prot_acc" %in% nms) {
    pep_frame[["prot_acc"]] <- paste0("-", pep_frame[["prot_acc"]])
  }
  
  
  
  pep_frame
}


#' Reverses peptide sequences.
#' 
#' @param seqs Lists of peptide sequences.
reverse_seqs <- function (seqs) {
  
  fis <- stringi::stri_sub(seqs, 1, 1)
  las <- stringi::stri_sub(seqs, -1, -1)
  lens <- stringi::stri_length(seqs)
  
  revs <- stringi::stri_reverse(seqs)
  substring(revs, 1) <- fis
  substring(revs, lens) <- las
  
  revs
}

