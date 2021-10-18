#' Batch processing of MS2 ions
#' 
#' @inheritParams matchMS
#' @export
batch_ms2ions <- function (fasta = c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
                                     "~/proteoM/dbs/fasta/crap/crap.fasta"),
                           fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                         "Carbamidomethyl (C)"),
                           varmods = c("Acetyl (Protein N-term)",
                                       "Oxidation (M)", "Deamidated (N)",
                                       "Gln->pyro-Glu (N-term = Q)"),
                           include_insource_nl = FALSE,
                           exclude_phospho_nl = TRUE, 
                           type_ms2ions = "by", 
                           maxn_vmods_setscombi = 64L,
                           maxn_vmods_per_pep = 5L, 
                           maxn_sites_per_vmod = 3L, 
                           maxn_vmods_sitescombi_per_pep = 32L, 
                           .path_cache = NULL,
                           .path_fasta = NULL,
                           digits = 5L) {
  
  options(digits = 9L)
  
  # ---
  if (is.null(.path_cache)) {
    .path_cache <- "~/proteoM/.MSearches/Cache/Calls/"
  }
  
  .path_cache <- find_dir(.path_cache, create = FALSE)
  
  if (is.null(.path_cache)) {
    stop("Cache path not existed: ", .path_cache, call. = FALSE)
  }
  
  # ---
  if (is.null(.path_fasta)) {
    .path_fasta <- file.path(gsub("(.*)\\.[^\\.]*$", "\\1", fasta[1]))
  }
  
  .path_fasta <- find_dir(.path_fasta, create = FALSE)
  
  if (is.null(.path_fasta)) {
    stop("Fasta path not existed: ", .path_fasta, call. = FALSE)
  }
  
  # ---
  file_aa <- file.path(.path_fasta, "aa_masses_all.rds")
  
  if (!file.exists(file_aa)) {
    aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                                   varmods = varmods,
                                   maxn_vmods_setscombi = maxn_vmods_setscombi,
                                   add_varmasses = FALSE,
                                   add_nlmasses = FALSE, 
                                   exclude_phospho_nl = exclude_phospho_nl, 
                                   out_path = .path_fasta) %T>%
      saveRDS(file_aa)
  } else {
    aa_masses_all <- readRDS(file_aa)
  }
  
  types <- purrr::map_chr(aa_masses_all, attr, "type", exact = TRUE)
  mod_indexes <- find_mod_indexes(.path_fasta)
  
  # ---
  path_ms1masses <- create_dir(file.path(.path_fasta, "ms1masses"))
  ms1_times <- dir(path_ms1masses, all.files = TRUE, no.. = TRUE)
  
  if (!length(ms1_times)) {
    stop("Time stamps of MS1 masses not found.", call. = FALSE)
  }
  
  lapply(ms1_times, hbatch_ms2ions, 
         types = types, 
         mod_indexes = mod_indexes, 
         type_ms2ions = type_ms2ions, 
         maxn_vmods_per_pep = maxn_vmods_per_pep, 
         maxn_sites_per_vmod = maxn_sites_per_vmod, 
         maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
         .path_cache = .path_cache, 
         .path_fasta = .path_fasta, 
         digits = digits)
  
  invisible(NULL)
}


#' Helper of \link{batch_ms2ions}
#'
#' For a given time stamp.
#'
#' length(pepmasses_i.rds) == length(ms2masses_i.rds). Note tht ms2masses_i can
#' be nested due to the permutation and NLs of AA residues under the same
#' pep_seq and ms1mass.
#'
#' @param ms1_time A cached MS1 time (directory).
#' @param types The types of modification, e.g., "amods- tmod+ vnl- fnl-" etc.
#' @inheritParams batch_ms2ions
#' @inheritParams ms2match
hbatch_ms2ions <- function (ms1_time = NULL, types = NULL, mod_indexes = NULL, 
                            type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                            maxn_sites_per_vmod = 3L, 
                            maxn_vmods_sitescombi_per_pep = 32L, 
                            .path_cache = NULL, .path_fasta = NULL, 
                            digits = 5L) {
  
  path_ms2time <- create_dir(file.path(.path_fasta, "ms2masses", ms1_time))
  path_bin2 <- create_dir(file.path(path_ms2time, "bin_ms2masses"))
  
  new_bins <- check_ms2frames(.path_fasta, ms1_time, .path_cache, path_bin2)
  
  if (!length(new_bins)) {
    return(NULL)
  }
  
  # --- MS2 Not yet binned
  path_ms1masses <- file.path(.path_fasta, "ms1masses")
  path_ms1time <- file.path(path_ms1masses, ms1_time)
  
  file_ms1time <- local({
    files <- list.files(path = path_ms1time, 
                        pattern = paste0("^pepmasses_", "\\d+\\.rds$"))
    
    idxes <- gsub("^pepmasses_(\\d+)\\.rds$", "\\1", files)
    idxes <- as.integer(idxes)
    ords <- order(idxes)
    files[ords]
  })

  ## Targets
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl-", 
                             "amods- tmod+ vnl- fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(path_ms1time, paste0("pepmasses_", i, ".rds")))
      ms1s <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
      ms2i <- file.path(path_ms2time, paste0("ms2masses_", i, ".rds"))
      
      if (file.exists(ms2i)) {
        ms2s <- readRDS(ms2i)
      } else {
        ntmass <- find_nterm_mass(aa_masses)
        ctmass <- find_cterm_mass(aa_masses)
        ntmod = NULL
        ctmod = NULL
        amods = NULL
        vmods_nl = NULL
        fmods_nl = NULL
        
        n_cores <- detect_cores(16L)
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("gen_ms2ions_base", 
            "ms2ions_by_type", 
            "byions", "czions", "axions"), 
          envir = environment(proteoM:::gen_ms2ions_base)
        )
        
        ms2s <- parallel::clusterApply(
          cl, chunksplit(ms1s, n_cores, "list"), 
          mgen_ms2ions, 
          aa_masses = aa_masses, 
          ntmod = ntmod, 
          ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, 
          vmods_nl = vmods_nl, 
          fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits, 
          FUN = gen_ms2ions_base
        )
        
        parallel::stopCluster(cl)
        
        ms2s <- purrr::flatten(ms2s)
        saveRDS(ms2s, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      }

      ms2s <- make_ms2frames(ms1_time, ms1s, ms2s, .path_cache, path_bin2, i)
      rm(list = c("aa_masses", "ms1s", "ms2s"))
      gc()
    }
  }
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl+", 
                             "amods- tmod+ vnl- fnl+"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(path_ms1time, paste0("pepmasses_", i, ".rds")))
      ms1s <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
      ms2i <- file.path(path_ms2time, paste0("ms2masses_", i, ".rds"))
      
      if (file.exists(ms2i)) {
        ms2s <- readRDS(ms2i)
      } else {
        ntmass <- find_nterm_mass(aa_masses)
        ctmass <- find_cterm_mass(aa_masses)
        ntmod = NULL
        ctmod = NULL
        amods = NULL
        vmods_nl = NULL
        fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
        
        n_cores <- detect_cores(16L)
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("gen_ms2ions_a0_vnl0_fnl1", 
            "gen_ms2ions_base", 
            "ms2ions_by_type", 
            "byions", "czions", "axions"), 
          envir = environment(proteoM:::gen_ms2ions_a0_vnl0_fnl1)
        )
        
        ms2s <- parallel::clusterApply(
          cl, chunksplit(ms1s, n_cores, "list"), 
          mgen_ms2ions, 
          aa_masses = aa_masses, 
          ntmod = ntmod, 
          ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, 
          vmods_nl = vmods_nl, 
          fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits, 
          FUN = gen_ms2ions_a0_vnl0_fnl1
        )
        
        parallel::stopCluster(cl)
        
        ms2s <- purrr::flatten(ms2s)
        saveRDS(ms2s, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      }

      ms2s <- make_ms2frames(ms1_time, ms1s, ms2s, .path_cache, path_bin2, i)
      rm(list = c("aa_masses", "ms1s", "ms2s"))
      gc()
    }
  }
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #        (ALL amods are vnl-)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl-", 
                             "amods+ tmod+ vnl- fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(path_ms1time, paste0("pepmasses_", i, ".rds")))
      ms1s <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
      ms2i <- file.path(path_ms2time, paste0("ms2masses_", i, ".rds"))
      
      if (file.exists(ms2i)) {
        ms2s <- readRDS(ms2i)
      } else {
        ntmass <- find_nterm_mass(aa_masses)
        ctmass <- find_cterm_mass(aa_masses)
        ntmod <- attr(aa_masses, "ntmod", exact = TRUE) # for checking MS1 masses
        ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
        amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
        vmods_nl = NULL
        fmods_nl <- NULL
        
        n_cores <- detect_cores(16L)
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("gen_ms2ions_a1_vnl0_fnl0", 
            "combi_mvmods2", 
            "combi_vmods2", 
            "find_intercombi_p2", 
            "check_ms1_mass_vmods2", 
            "calc_ms2ions_a1_vnl0_fnl0", 
            "ms2ions_by_type", 
            "byions", "czions", "axions", 
            "add_hexcodes"), 
          envir = environment(proteoM:::gen_ms2ions_a1_vnl0_fnl0)
        )
        
        ms2s <- parallel::clusterApply(
          cl, chunksplit(ms1s, n_cores, "list"), 
          mgen_ms2ions, 
          aa_masses = aa_masses, 
          ntmod = ntmod, 
          ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, 
          vmods_nl = vmods_nl, 
          fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits, 
          FUN = gen_ms2ions_a1_vnl0_fnl0
        )
        
        parallel::stopCluster(cl)
        
        ms2s <- purrr::flatten(ms2s)
        saveRDS(ms2s, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      }

      ms2s <- make_ms2frames(ms1_time, ms1s, ms2s, .path_cache, path_bin2, i)
      rm(list = c("aa_masses", "ms1s", "ms2s"))
      gc()
    }
  }
  
  # (9, 10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #         (ANY amod is vnl+)
  
  inds <- which(types %in% c("amods+ tmod- vnl+ fnl-", 
                             "amods+ tmod+ vnl+ fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(path_ms1time, paste0("pepmasses_", i, ".rds")))
      ms1s <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
      ms2i <- file.path(path_ms2time, paste0("ms2masses_", i, ".rds"))
      
      if (file.exists(ms2i)) {
        ms2s <- readRDS(ms2i)
      } else {
        ntmass <- find_nterm_mass(aa_masses)
        ctmass <- find_cterm_mass(aa_masses)
        ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
        ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
        amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
        vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
        fmods_nl <- NULL
        
        n_cores <- detect_cores(16L)
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("gen_ms2ions_a1_vnl1_fnl0", 
            "combi_mvmods2", 
            "combi_vmods2", 
            "find_intercombi_p2", 
            "check_ms1_mass_vmods2", 
            "calc_ms2ions_a1_vnl1_fnl0", 
            "ms2ions_by_type", 
            "byions", "czions", "axions", 
            "add_hexcodes_vnl2"), 
          envir = environment(proteoM:::gen_ms2ions_a1_vnl1_fnl0)
        )
        
        ms2s <- parallel::clusterApply(
          cl, chunksplit(ms1s, n_cores, "list"), 
          mgen_ms2ions, 
          aa_masses = aa_masses, 
          ntmod = ntmod, 
          ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, 
          vmods_nl = vmods_nl, 
          fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits, 
          FUN = gen_ms2ions_a1_vnl1_fnl0
        )
        
        parallel::stopCluster(cl)
        
        ms2s <- purrr::flatten(ms2s)
        saveRDS(ms2s, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      }
      
      ms2s <- make_ms2frames(ms1_time, ms1s, ms2s, .path_cache, path_bin2, i)
      rm(list = c("aa_masses", "ms1s", "ms2s"))
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
      aa_masses <- readRDS(file.path(path_ms1time, paste0("pepmasses_", i, ".rds")))
      ms1s <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
      ms2i <- file.path(path_ms2time, paste0("ms2masses_", i, ".rds"))
      
      if (file.exists(ms2i)) {
        ms2s <- readRDS(ms2i)
      } else {
        ntmass <- find_nterm_mass(aa_masses)
        ctmass <- find_cterm_mass(aa_masses)
        ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
        ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
        amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
        vmods_nl <- NULL
        fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
        
        n_cores <- detect_cores(16L)
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("gen_ms2ions_a1_vnl0_fnl1", 
            "combi_mvmods2", 
            "combi_vmods2", 
            "find_intercombi_p2", 
            "check_ms1_mass_vmods2", 
            "calc_ms2ions_a1_vnl0_fnl1", 
            "ms2ions_by_type", 
            "byions", "czions", "axions", 
            "add_hexcodes_fnl2"), 
          envir = environment(proteoM:::gen_ms2ions_a1_vnl0_fnl1)
        )
        
        ms2s <- parallel::clusterApply(
          cl, chunksplit(ms1s, n_cores, "list"), 
          mgen_ms2ions, 
          aa_masses = aa_masses, 
          ntmod = ntmod, 
          ctmod = ctmod, 
          ntmass = ntmass, 
          ctmass = ctmass, 
          amods = amods, 
          vmods_nl = vmods_nl, 
          fmods_nl = fmods_nl, 
          mod_indexes = mod_indexes, 
          type_ms2ions = type_ms2ions, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits, 
          FUN = gen_ms2ions_a1_vnl0_fnl1
        )
        
        parallel::stopCluster(cl)
        
        ms2s <- purrr::flatten(ms2s)
        saveRDS(ms2s, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      }

      ms2s <- make_ms2frames(ms1_time, ms1s, ms2s, .path_cache, path_bin2, i)
      rm(list = c("aa_masses", "ms1s", "ms2s"))
      gc()
    }
  }
  
  ## Decoys...
  
  invisible(NULL)
}





#' Multiple generations of MS2 ions.
#'
#' @param ms1s Lists of named values. Peptide sequences in names and precursor
#'   masses in values.
#' @param FUN A function pointer to, e.g., \link{gen_ms2ions_base}.
#' @inheritParams gen_ms2ions_base
mgen_ms2ions <- function (ms1s = NULL, aa_masses = NULL, 
                          ntmod = NULL, ctmod = NULL, 
                          ntmass = NULL, ctmass = NULL, 
                          amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                          mod_indexes = NULL, type_ms2ions = "by", 
                          maxn_vmods_per_pep = 5L, 
                          maxn_sites_per_vmod = 3L, 
                          maxn_vmods_sitescombi_per_pep = 32L, 
                          digits = 4L, 
                          FUN) {
  
  seqs <- names(ms1s)
  masses <- unname(ms1s)
  
  ms2s <- mapply(
    FUN, 
    aa_seq = seqs, 
    ms1_mass = masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      vmods_nl = vmods_nl, 
      fmods_nl = fmods_nl,
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  names(ms2s) <- seqs

  invisible(ms2s)
}


#' Bin MS2 ion series by precursor masses.
#'
#' The same MS2 results but different frame numbers according to the paramter
#' set of min_mass, max_mass and ppm_ms1.
#' 
#' @param ms1_time A ms1_time stamp.
#' @param ms2s Lists of MS2 ion series.
#' @param path_bin2 A file path to binned MS2 results. 
#' @param i An index of aa_masses
#' @inheritParams hbatch_ms2ions
#' @inheritParams mgen_ms2ions
make_ms2frames <- function (ms1_time, ms1s, ms2s, .path_cache, path_bin2, i) {

  path_time_cache <- file.path(.path_cache, "calc_pepmasses2", ms1_time, 
                               "bin_ms1masses")
  
  files_ms1bins <- list.files(path = path_time_cache, 
                              pattern = "\\.[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]+.*\\.rda$", 
                              all.files = TRUE, no.. = TRUE)
  
  if (!length(files_ms1bins)) {
    stop("Parameter files not found under ", path_time_cache, call. = FALSE)
  }
  
  pars_ms1bins <- lapply(files_ms1bins, function (x) {
    load(file.path(path_time_cache, x))
    
    list(min_mass = call_pars$min_mass, 
         max_mass = call_pars$max_mass, 
         ppm_ms1 = call_pars$ppm_ms1, 
         .time_bin = gsub("\\.rda", "", x))
  })
  
  len <- length(pars_ms1bins)
  
  if (len == 1L) {
    make_ms2frames_bypars(pars_ms1bins, ms1s, ms2s, path_bin2, i)
  } else {
    n_cores <- detect_cores(len)
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterApply(
      cl, pars_ms1bins, 
      make_ms2frames_bypars, 
      ms1s, ms2s, path_bin2, i
    )
    
    parallel::stopCluster(cl)
  }
}


#' Bins MS2 ion series at a set of parameters.
#' 
#' The set of parameters include min_mass, max_mass and ppm_ms1. 
#' 
#' @param pars A set of cached parameters from \code{.Msearches}.
#' @param is_ms1_three_frame Logical; is the searches by the three frames of
#'   preceeding, current and following.
#' @inheritParams make_ms2frames
make_ms2frames_bypars <- function (pars, ms1s, ms2s, path_bin2, i, 
                                   is_ms1_three_frame = TRUE) {
  
  min_mass <- pars$min_mass
  max_mass <- pars$max_mass
  ppm_ms1 <- pars$ppm_ms1
  .time_bin <- pars$.time_bin
  
  out_path <- create_dir(file.path(path_bin2, .time_bin))
  
  if (is_ms1_three_frame) {
    ppm_ms1_new <- as.integer(ceiling(ppm_ms1 * .5))
  } else {
    ppm_ms1_new <- ppm_ms1
  }

  ps <- find_ms1_cutpoints(min_mass, max_mass, ppm_ms1_new)
  frames <- findInterval(ms1s, ps)

  ans <- split(ms2s, frames)
  saveRDS(ans, file.path(out_path, paste0("binned_ms2_", i, ".rds")))
  
  invisible(NULL)
}


#' Checks for pre-existed MS2 bins.
#' 
#' @inheritParams hbatch_ms2ions
#' @inheritParams make_ms2frames
check_ms2frames <- function (.path_fasta, ms1_time, .path_cache, path_bin2) {
  
  # finds `.times_bin` under a given `ms1_time`
  .times_bin <- local({
    path_cache <- file.path(.path_cache, "calc_pepmasses2", ms1_time, "bin_ms1masses")
    
    rdas <- list.files(path = path_cache, 
                       pattern = "\\.[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]+.*\\.rda$", 
                       all.files = TRUE, no.. = TRUE)
    gsub("\\.rda", "", rdas)
  })
  
  
  if (!length(.times_bin)) {
    stop("Parameter files not found under ", path_cache, call. = FALSE)
  }
  
  # compare files
  .path_bin1 <- file.path(.path_fasta, "ms1masses", ms1_time, "bin_ms1masses")
  
  nots <- lapply(.times_bin, function (x) {
    ms1files <- list.files(path = file.path(.path_bin1, x), 
                           pattern = "^binned_theopeps_\\d+\\.rds$")
    
    ms2files <- list.files(path = file.path(path_bin2, x), 
                           pattern = "^binned_ms2_\\d+\\.rds$")
    
    (length(ms2files) != length(ms1files))
  })
  
  nots <- unlist(nots, recursive = FALSE, use.names = FALSE)
  newbins <- .times_bin[nots]
}

