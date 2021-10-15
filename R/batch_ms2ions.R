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
                           .path_fasta = NULL,
                           digits = 5L) {
  
  options(digits = 9L)
  
  if (is.null(.path_fasta)) {
    .path_fasta <- file.path(gsub("(.*)\\.[^\\.]*$", "\\1", fasta[1]))
  }

  path_ms1masses <- create_dir(file.path(.path_fasta, "ms1masses"))

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
  time_stamps <- dir(path_ms1masses, all.files = TRUE, no.. = TRUE)
  
  # ---
  # time_stamps = time_stamps[2]
  # ---
  
  lapply(time_stamps, hbatch_ms2ions, 
         types = types, 
         mod_indexes = mod_indexes, 
         type_ms2ions = type_ms2ions, 
         maxn_vmods_per_pep = maxn_vmods_per_pep, 
         maxn_sites_per_vmod = maxn_sites_per_vmod, 
         maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
         .path_fasta = .path_fasta, 
         digits = digits)
  
  invisible(NULL)
}


#' Helper of \link{batch_ms2ions} by time directory.
#' 
#' @param time A cached time (directory).
#' @param types The types of modification, e.g., "amods- tmod+ vnl- fnl-" etc.
#' @inheritParams batch_ms2ions
#' @inheritParams ms2match
hbatch_ms2ions <- function (time = NULL, types = NULL, mod_indexes = NULL, 
                            type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                            maxn_sites_per_vmod = 3L, 
                            maxn_vmods_sitescombi_per_pep = 32L, 
                            .path_fasta = NULL, digits = 5L) {
  
  path_ms1masses <- file.path(.path_fasta, "ms1masses")
  path_ms1time <- file.path(path_ms1masses, time)
  file_ms1time <- list.files(path = path_ms1time, 
                             pattern = paste0("^pepmasses_", "\\d+\\.rds$"))
  
  idxes <- gsub("^pepmasses_(\\d+)\\.rds$", "\\1", file_ms1time)
  idxes <- as.integer(idxes)
  ords <- order(idxes)
  file_ms1time <- file_ms1time[ords]
  
  rm(list = c("ords", "idxes"))
  
  path_ms2time <- create_dir(file.path(.path_fasta, "ms2masses", time))
  
  ## Targets 
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl-", 
                             "amods- tmod+ vnl- fnl-"))
  
  # inds <- NULL
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(.path_fasta, "ms1masses", time, 
                                     paste0("pepmasses_", i, ".rds")))
      data_ms1 <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
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
      
      n_cores <- detect_cores()
      
      # data_ms1 <- data_ms1[1:128]
      data_ms1 <- chunksplit(data_ms1, n_cores, "list")
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      # parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
      
      parallel::clusterExport(
        cl,
        c("mgen_ms2ions_base", 
          "gen_ms2ions_base", 
          "ms2ions_by_type", 
          "byions", "czions", "axions"), 
        envir = environment(proteoM:::gen_ms2ions_base)
      )

      ans <- parallel::clusterApply(
        cl, data_ms1, 
        mgen_ms2ions_base, 
        aa_masses = aa_masses, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      )
      
      parallel::stopCluster(cl)
      
      ans <- purrr::flatten(ans)
      saveRDS(ans, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      rm(list = "ans")
      gc()
    }
  }
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl+", 
                             "amods- tmod+ vnl- fnl+"))
  
  # inds <- NULL
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(.path_fasta, "ms1masses", time, 
                                     paste0("pepmasses_", i, ".rds")))
      data_ms1 <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
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
      
      n_cores <- detect_cores()
      
      # data_ms1 <- data_ms1[1:128]
      data_ms1 <- chunksplit(data_ms1, n_cores, "list")
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(
        cl,
        c("mgen_ms2ions_a0_vnl0_fnl1", 
          "gen_ms2ions_a0_vnl0_fnl1", 
          "gen_ms2ions_base", 
          "ms2ions_by_type", 
          "byions", "czions", "axions"), 
        envir = environment(proteoM:::gen_ms2ions_a0_vnl0_fnl1)
      )

      ans <- parallel::clusterApply(
        cl, data_ms1, 
        mgen_ms2ions_a0_vnl0_fnl1, 
        aa_masses = aa_masses, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      )
      
      parallel::stopCluster(cl)
      
      ans <- purrr::flatten(ans)
      saveRDS(ans, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      rm(list = "ans")
      gc()
    }
  }
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #        (ALL amods are vnl-)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl-", 
                             "amods+ tmod+ vnl- fnl-"))
  
  # inds <- NULL
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(.path_fasta, "ms1masses", time, 
                                     paste0("pepmasses_", i, ".rds")))
      data_ms1 <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
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
      
      n_cores <- detect_cores()
      
      # data_ms1 <- data_ms1[1:128]
      data_ms1 <- chunksplit(data_ms1, n_cores, "list")
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(
        cl,
        c("mgen_ms2ions_a1_vnl0_fnl0", 
          "gen_ms2ions_a1_vnl0_fnl0", 
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

      ans <- parallel::clusterApply(
        cl, data_ms1, 
        mgen_ms2ions_a1_vnl0_fnl0, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      )
      
      parallel::stopCluster(cl)
      
      ans <- purrr::flatten(ans)
      saveRDS(ans, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      rm(list = "ans")
      gc()
    }
  }
  
  # (9, 10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #         (ANY amod is vnl+)
  
  inds <- which(types %in% c("amods+ tmod- vnl+ fnl-", 
                             "amods+ tmod+ vnl+ fnl-"))
  
  if (length(inds)) {
    for (i in inds) {
      aa_masses <- readRDS(file.path(.path_fasta, "ms1masses", time, 
                                     paste0("pepmasses_", i, ".rds")))
      data_ms1 <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
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
      
      n_cores <- detect_cores()
      
      # data_ms1 <- data_ms1[1:128]
      data_ms1 <- chunksplit(data_ms1, n_cores, "list")

      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(
        cl,
        c("mgen_ms2ions_a1_vnl1_fnl0", 
          "gen_ms2ions_a1_vnl1_fnl0", 
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

      ans <- parallel::clusterApply(
        cl, data_ms1, 
        mgen_ms2ions_a1_vnl1_fnl0, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        vmods_nl = vmods_nl, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      )
      
      parallel::stopCluster(cl)
      
      ans <- purrr::flatten(ans)
      saveRDS(ans, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      rm(list = "ans")
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
      aa_masses <- readRDS(file.path(.path_fasta, "ms1masses", time, 
                                     paste0("pepmasses_", i, ".rds")))
      data_ms1 <- attr(aa_masses, "data")
      attr(aa_masses, "data") <- NULL
      
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
      
      n_cores <- detect_cores()
      
      # data_ms1 <- data_ms1[1:128]
      data_ms1 <- chunksplit(data_ms1, n_cores, "list")
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(
        cl,
        c("mgen_ms2ions_a1_vnl0_fnl1", 
          "gen_ms2ions_a1_vnl0_fnl1", 
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

      ans <- parallel::clusterApply(
        cl, data_ms1, 
        mgen_ms2ions_a1_vnl0_fnl1, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = 
          maxn_vmods_sitescombi_per_pep, 
        digits = digits
      )
      
      parallel::stopCluster(cl)
      
      ans <- purrr::flatten(ans)
      saveRDS(ans, file.path(path_ms2time, paste0("ms2masses_", i, ".rds")))
      rm(list = "ans")
      gc()
    }
  }
  
  ## Decoys...
  
  invisible(NULL)
}


#' Multiple \link{gen_ms2ions_base}.
#'
#' @param data Lists of named values. Peptide sequences in names and precursor
#'   masses in values.
#' @inheritParams gen_ms2ions_base
mgen_ms2ions_base <- function (data = NULL, aa_masses = NULL, 
                               ntmass = NULL, ctmass = NULL, mod_indexes = NULL, 
                               type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               digits = 4L) {
  
  seqs <- names(data)
  masses <- unname(data)
  
  ans <- mapply(
    gen_ms2ions_base, 
    seqs, 
    masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmass = NULL, ctmass = NULL, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  names(ans) <- seqs
  
  invisible(ans)
}


#' Multiple \link{gen_ms2ions_a0_vnl0_fnl1}.
#' 
#' @inheritParams gen_ms2ions_a0_vnl0_fnl1
#' @rdname mgen_ms2ions_base
mgen_ms2ions_a0_vnl0_fnl1 <- function (data = NULL, aa_masses = NULL, 
                                       ntmass = NULL, ctmass = NULL, fmods_nl = NULL, 
                                       mod_indexes = NULL, type_ms2ions = "by", 
                                       maxn_vmods_per_pep = 5L, 
                                       maxn_sites_per_vmod = 3L, 
                                       maxn_vmods_sitescombi_per_pep = 32L, 
                                       digits = 4L) {
  
  seqs <- names(data)
  masses <- unname(data)
  
  ans <- mapply(
    gen_ms2ions_a0_vnl0_fnl1, 
    seqs, 
    masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  names(ans) <- seqs
  
  invisible(ans)
}


#' Multiple \link{gen_ms2ions_a1_vnl0_fnl0}.
#' 
#' @inheritParams gen_ms2ions_a1_vnl0_fnl0
#' @rdname mgen_ms2ions_base
mgen_ms2ions_a1_vnl0_fnl0 <- function (data = NULL, aa_masses = NULL, 
                                       ntmod = NULL, ctmod = NULL, 
                                       ntmass = NULL, ctmass = NULL, amods = NULL, 
                                       mod_indexes = NULL, type_ms2ions = "by", 
                                       maxn_vmods_per_pep = 5L, 
                                       maxn_sites_per_vmod = 3L, 
                                       maxn_vmods_sitescombi_per_pep = 32L, 
                                       digits = 4L) {
  
  seqs <- names(data)
  masses <- unname(data)
  
  ans <- mapply(
    gen_ms2ions_a1_vnl0_fnl0, 
    seqs, 
    masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  names(ans) <- seqs
  
  invisible(ans)
}


#' Multiple \link{gen_ms2ions_a1_vnl1_fnl0}.
#' 
#' @inheritParams gen_ms2ions_a1_vnl1_fnl0
#' @rdname mgen_ms2ions_base
mgen_ms2ions_a1_vnl1_fnl0 <- function (data = NULL, aa_masses = NULL, 
                                       ntmod = NULL, ctmod = NULL, 
                                       ntmass = NULL, ctmass = NULL, 
                                       amods = NULL, vmods_nl = NULL, 
                                       mod_indexes = NULL, type_ms2ions = "by", 
                                       maxn_vmods_per_pep = 5L, 
                                       maxn_sites_per_vmod = 3L, 
                                       maxn_vmods_sitescombi_per_pep = 32L, 
                                       digits = 4L) {
  
  seqs <- names(data)
  masses <- unname(data)
  
  ans <- mapply(
    gen_ms2ions_a1_vnl1_fnl0, 
    seqs, 
    masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      vmods_nl = vmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  names(ans) <- seqs
  
  invisible(ans)
}


#' Multiple \link{gen_ms2ions_a1_vnl0_fnl1}.
#' 
#' @inheritParams gen_ms2ions_a1_vnl0_fnl1
#' @rdname mgen_ms2ions_base
mgen_ms2ions_a1_vnl0_fnl1 <- function (data = NULL, aa_masses = NULL, 
                                       ntmod = NULL, ctmod = NULL, 
                                       ntmass = NULL, ctmass = NULL, 
                                       amods = NULL, fmods_nl = NULL, 
                                       mod_indexes = NULL, type_ms2ions = "by", 
                                       maxn_vmods_per_pep = 5L, 
                                       maxn_sites_per_vmod = 3L, 
                                       maxn_vmods_sitescombi_per_pep = 32L, 
                                       digits = 4L) {
  
  seqs <- names(data)
  masses <- unname(data)
  
  ans <- mapply(
    gen_ms2ions_a1_vnl0_fnl1, 
    seqs, 
    masses, 
    MoreArgs = list(
      aa_masses = aa_masses, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      fmods_nl = fmods_nl, 
      mod_indexes = mod_indexes, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = 
        maxn_vmods_sitescombi_per_pep, 
      digits = digits
    ), 
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  names(ans) <- seqs
  
  invisible(ans)
}


