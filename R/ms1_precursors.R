#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param index_mods Logical; if TRUE, converts the names of modifications to
#'   indexes. Not currently used.
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from proteins. The enzyme
#'   is currently \code{trypsin}.
#' @param maxn_fasta_seqs Integer; the maximum number of protein sequences in
#'   fasta files.
#' @param min_len Integer; the minimum length of peptides. Shorter peptides will
#'   be excluded.
#' @param max_len Integer; the maximum length of peptides. Longer peptides will
#'   be excluded.
#' @param max_miss The maximum number of mis-cleavages per peptide sequence.
#' @param digits Integer; the number of decimal places to be used.
#' @inheritParams load_fasta
#' @inheritParams binTheoPeps
#' @inheritParams calc_aamasses
#' @inheritParams mcalc_monopep
#' @inheritParams calc_monopep
#' @inheritParams load_fasta2
#' @import parallel
#' @importFrom rlang current_env
#' @examples
#' \donttest{
#' res <- calc_pepmasses()
#'
#' library(purrr)
#' library(magrittr)
#'
#' res_attrs <- map(res, attributes)
#' map(res_attrs, names)
#' res_attrs %>% map(`[[`, "vmods")
#' res_mods <- map(res_attrs, `[`,
#'                 c("fmods", "fmods_ps", "fmods_neuloss",
#'                   "vmods", "vmods_ps", "vmods_neuloss"))
#'
#' res_data <- map(res_attrs, `[[`, "data")
#' peps_combi_1 <- res_data[[1]]
#'
#' # base: fixedmods without neulosses
#' unlist(res_data[[1]], use.names = FALSE) %>% length()
#'
#' # fixedmods, fixedmods + fixedmods_neulosses, varmods, varmods_neulosses
#' unlist(res_data, use.names = FALSE) %>% length()
#'
#' }
#'
#' @export
calc_pepmasses2 <- function (
  fasta = "~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
  acc_type = "uniprot_acc",
  acc_pattern = NULL,
  fixedmods = c("TMT6plex (K)",
                "Carbamidomethyl (C)"),
  varmods = c("TMT6plex (N-term)",
              "Acetyl (Protein N-term)",
              "Oxidation (M)",
              "Deamidated (N)",
              "Gln->pyro-Glu (N-term = Q)"),
  include_insource_nl = FALSE,
  index_mods = FALSE,
  enzyme = c("trypsin"),
  maxn_fasta_seqs = 50000L,
  maxn_vmods_setscombi = 64L,
  maxn_vmods_per_pep = 5L,
  maxn_sites_per_vmod = 3L,
  min_len = 7L, max_len = 100L, max_miss = 2L,
  out_path = NULL,
  digits = 4L,
  parallel = TRUE) {

  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)

  on.exit(
    if (exists(".savecall", envir = current_env())) {
      if (.savecall) {
        save_call2(path = file.path(.path_cache), fun = "calc_pepmasses",
                   time = .time_stamp)
      }
    } else {
      # warning("terminated abnormally.", call. = TRUE)
    },
    add = TRUE
  )

  # ---
  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, logical(1L)))
  stopifnot(min_len >= 0L, max_len >= min_len, max_miss <= 100L)

  # ---
  .path_cache <- create_dir("~/proteoM/.MSearch/Cache/Calls")

  .time_stamp <- match_calltime(path = .path_cache,
                                fun = "calc_pepmasses",
                                nms = c("fasta", "acc_type", "acc_pattern",
                                        "fixedmods", "varmods",
                                        "include_insource_nl", "enzyme",
                                        "maxn_fasta_seqs", "maxn_vmods_setscombi",
                                        "maxn_vmods_per_pep", "maxn_sites_per_vmod",
                                        "min_len", "max_len", "max_miss"))

  .path_fasta <- fasta %>%
    gsub("\\\\", "/", .) %>%
    gsub("(.*/)[^/]*", "\\1", .) %>%
    purrr::map_chr(find_dir) %>%
    `[[`(1)

  # ---
  if (length(.time_stamp) > 0L) {
    message("Loading peptide masses from cache.")

    aa_masses <- find_aa_masses(out_path = out_path,
                                fixedmods = fixedmods,
                                varmods = varmods,
                                maxn_vmods_setscombi = maxn_vmods_setscombi,
                                index_mods = index_mods)

    files <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp),
                        pattern = "pepmasses_\\d+\\.rds$")

    if (! length(files) == length(aa_masses)) {
      stop("Not all precursor masses were found: ", paste0("\n", files), ".\n",
           "Remove cache file: ", file.path(.path_cache, "...", .time_stamp),
           " and try again.",
           call. = FALSE)
    }

    rm(list = c("aa_masses", "files"))
    gc()

    fwd_peps <- NULL
    rev_peps <- NULL

    .savecall <- FALSE
  } else {
    delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$",
                                       "\\.xls$", "\\.csv$", "\\.txt$",
                                       "^mgf$", "^mgfs$"))

    .time_stamp <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")

    aa_masses <- find_aa_masses(out_path = out_path,
                                fixedmods = fixedmods,
                                varmods = varmods,
                                maxn_vmods_setscombi = maxn_vmods_setscombi,
                                index_mods = index_mods)

    aa_masses_1 <- aa_masses[[1]]

    gc()


    # --- Forward and reversed sequences  ---
    # (not yet concatenation by the number of missed cleavages)
    seqs_0 <- split_fastaseqs(fasta = fasta,
                              acc_type = acc_type,
                              acc_pattern = acc_pattern,
                              maxn_fasta_seqs = maxn_fasta_seqs,
                              max_miss = max_miss)


    # --- Special case of FIXED Protein [NC]-term modification(s) ---
    #
    # Dispatches `pep_seq`s by `fixedmod`s
    # (only incur at the rare case of Protein terminals being `fixedmods`)
    # (note that `fixed Protein N-term` -> no `variable Protein N-term`)

    is_fixed_protnt <- any(grepl("Protein N-term", fixedmods))
    is_fixed_protct <- any(grepl("Protein C-term", fixedmods))
    seqs_0 <- distri_fpeps(seqs_0, max_miss, is_fixed_protnt, is_fixed_protct)

    # Calculates terminal mass (after `distri_fpeps`)
    # (otherwise, fixed Protein N/C-terminal masses in aa_masses["N|C-term"]
    #   will be applied to both N|C-term and Protein N|C term)
    #
    # if ftmass is other than 18.010565 -> no variable Protein N-term etc.
    # (mostly ftmass will be 18.010565 unless fixed [NC] terminal modifications)

    ftmass <- unname(aa_masses_1["N-term"] + aa_masses_1["C-term"])


    # --- Masses of sequences: fixed mods + terminals ---
    message("Calculating bare peptide masses ...")

    # (most often ftmass == 18.010565, unless fixed [NC] terminal modifications)
    # typically equivalent to fixedmods = NULL, varmods = NULL
    
    fwd_peps <- ms1masses_bare(seqs = seqs_0,
                               aa_masses = aa_masses_1,
                               ftmass = ftmass,
                               max_miss = max_miss,
                               min_len = min_len,
                               max_len = max_len,
                               maxn_vmods_per_pep = maxn_vmods_per_pep,
                               maxn_sites_per_vmod = maxn_sites_per_vmod,
                               is_fixed_protnt = is_fixed_protnt,
                               is_fixed_protct = is_fixed_protct,
                               parallel = parallel,
                               digits = digits)

    message("\tCompleted bare peptide masses.")

        
    # --- Distribution ---
    message("Distributing peptides by variable modifications.")

    # Note one-to-multiple expansion: 
    #   `length(fw_peps) == length(aa_masses)` after this
    fwd_peps <-  distri_peps(aa_masses = aa_masses,
                             peps = fwd_peps,
                             max_miss = max_miss)

    message("\tCompleted bare peptides distributions.")

    rm(list = c("seqs_0"))
    gc()

    
    # --- Delta masses of `variable` terminals  ---
    len <- length(aa_masses)
    types <- purrr::map_chr(aa_masses, attr, "type", exact = TRUE)
    
    # (e.g., on top of the `fixed` 18.010565)
    message("Adding terminal masses (variable modifications) ...")

    inds <- grep("tmod+", types, fixed = TRUE)

    if (length(inds)) {
      nt_inds <- which(types %in% c("amods- tmod+ vnl- fnl-",
                                    "amods- tmod+ vnl- fnl+"))

      for (i in inds) {
        fwd_peps[[i]] <- add_term_mass(aa_masses[[i]], fwd_peps[[i]])

        if (i %in% nt_inds) {
          message("\tCompleted peptide terminal masses: ",
                  paste(attributes(aa_masses[[i]])$fmods,
                        attributes(aa_masses[[i]])$vmods,
                        collapse = ", "))
        }

        gc()
      }

      rm(list = "nt_inds")
    }

    
    # --- Mass of variable mods and/or NLs ---
    message("Calculating peptide masses (variable modifications + neutral losses) ...")

    fmods_ps <- lapply(aa_masses, attr, "fmods_ps", exact = TRUE)
    vmods_ps <- lapply(aa_masses, attr, "vmods_ps", exact = TRUE)
    fmods_nl <- lapply(aa_masses, attr, "fmods_nl", exact = TRUE)
    vmods_nl <- lapply(aa_masses, attr, "vmods_nl", exact = TRUE)
    amods <- lapply(aa_masses, attr, "amods", exact = TRUE)
    tmod <- lapply(aa_masses, attr, "tmod", exact = TRUE)

    # `amods-` and `fnl+` (must be vnl- since amods-)
    #
    # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"

    if (include_insource_nl) {
      n_cores <- detect_cores()

      inds <- which(types %in% c("amods- tmod- vnl- fnl+",
                                 "amods- tmod+ vnl- fnl+"))

      if (length(inds)) {
        for (i in inds) {
          fmods_ps_i <- fmods_ps[[i]]
          vmods_ps_i <- vmods_ps[[i]]
          fmods_nl_i <- fmods_nl[[i]]
          vmods_nl_i <- vmods_nl[[i]]
          amods_i <- amods[[i]]
          tmod_i <- tmod[[i]]

          aa_masses_i <- aa_masses[[i]]
          ntmod_i <- attr(aa_masses_i, "ntmod", exact = TRUE)
          ctmod_i <- attr(aa_masses_i, "ctmod", exact = TRUE)

          fwd_peps_i <- fwd_peps[[i]]
          fnl_combi_i <- expand.grid(fmods_nl_i)
          
          cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
          
          parallel::clusterExport(
            cl,
            c("ms1_a0_fnl1_byprot", 
              "ms1_a0_fnl1_bypep", 
              "delta_ms1_a0_fnl1"), 
            envir = environment(proteoM:::ms1_a0_fnl1_byprot))

          fwd_peps[[i]] <- parallel::clusterApply(
            cl = cl, 
            x = chunksplit(fwd_peps_i, n_cores, "list"), 
            fun = lapply, 
            FUN = "ms1_a0_fnl1_byprot", 
            fnl_combi = fnl_combi_i, 
            aa_masses = aa_masses_i,
            digits = digits)

          parallel::stopCluster(cl)
          gc()
        }
      }
    }
    

    # `amods+`; (9-14) nested under (7-8)
    #
    # (7-8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
    #   (9-10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
    #   (11-12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
    #   (13-14) "amods+ tmod- vnl+ fnl+", "amods+ tmod+ vnl+ fnl+"
    
    n_cores <- detect_cores()

    inds <- which(types %in% c("amods+ tmod- vnl- fnl-",
                               "amods+ tmod+ vnl- fnl-",
                               "amods+ tmod- vnl+ fnl-",
                               "amods+ tmod+ vnl+ fnl-",
                               "amods+ tmod- vnl- fnl+",
                               "amods+ tmod+ vnl- fnl+",
                               "amods+ tmod- vnl+ fnl+",
                               "amods+ tmod+ vnl+ fnl+"))

    if (length(inds)) {
      for (i in inds) {
        amods_i <- amods[[i]]
        aa_masses_i <- aa_masses[[i]]

        fwd_peps_i <- fwd_peps[[i]]

        vmods_nl_i = vmods_nl[[i]]
        fmods_nl_i = fmods_nl[[i]]
        
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("ms1_a1_vnl0_fnl0_byprot", 
            "ms1_a1_vnl0_fnl0_bypep", 
            "unique_mvmods", 
            "vmods_elements", 
            "find_intercombi", 
            "delta_ms1_a0_fnl1", 
            "find_unique_sets", 
            "recur_flatten"), 
          envir = environment(proteoM:::ms1_a1_vnl0_fnl0_byprot))

        # avoids parLapply: may cause `rho` errors
        fwd_peps[[i]] <- parallel::clusterApply(
          cl = cl, 
          x = chunksplit(fwd_peps_i, n_cores, "list"), 
          fun = lapply, 
          FUN = "ms1_a1_vnl0_fnl0_byprot", 
          amods = amods_i, 
          aa_masses = aa_masses_i,
          vmods_nl = vmods_nl_i, 
          fmods_nl = fmods_nl_i,
          include_insource_nl = include_insource_nl,
          maxn_vmods_per_pep = maxn_vmods_per_pep,
          maxn_sites_per_vmod = maxn_sites_per_vmod,
          digits = digits) %>% 
          purrr::flatten()
        
        parallel::stopCluster(cl)
        gc()

        message("\tCompleted peptide masses: ",
                paste(attributes(aa_masses_i)$fmods,
                      attributes(aa_masses_i)$vmods,
                      collapse = ", "))
      }
    }

    rm(list = c("amods_i", "fmods_nl", "fmods_ps", "fwd_peps_i", 
                "vmods_nl_i", "aa_masses_1", "aa_masses_i"))
    gc()


    # === Outputs ===
    path_masses <- create_dir(file.path(.path_fasta, "pepmasses", .time_stamp))

    fwd_peps <- purrr::map2(aa_masses, fwd_peps, ~ {
      attr(.x, "data") <- .y
      .x
    })

    names(fwd_peps) <- seq_along(aa_masses)
    gc()

    purrr::walk2(seq_along(fwd_peps), fwd_peps, ~ {
      saveRDS(.y, file.path(path_masses, paste0("pepmasses_", .x, ".rds")))
    })
    gc()

    # With large datasets and combinations of aa_masses
    mem_used <- memory.size(max = TRUE)/memory.limit()
    is_large <- (length(aa_masses) > 48L && mem_used > .9)

    if (is_large) {
      rm(list = "fwd_peps")
      gc()

      fwd_peps <- NULL
    } else {
      gc()
    }

    # ---
    .savecall <- TRUE
  }

  assign(".path_cache", .path_cache, envir = .GlobalEnv)
  assign(".time_stamp", .time_stamp, envir = .GlobalEnv)
  assign(".path_fasta", .path_fasta, envir = .GlobalEnv)

  # invisible(list(fwd = fwd_peps, rev = rev_peps))
  invisible(fwd_peps)
}


#' Finds the existence of \code{aa_masses_all.rds}.
#'
#' @inheritParams calc_pepmasses2
find_aa_masses  <- function(out_path = NULL, fixedmods = NULL, varmods = NULL,
                            maxn_vmods_setscombi = 64L, index_mods = FALSE) {

  if (!file.exists(file.path(out_path, "temp", "aa_masses_all.rds"))) {
    message("Computing the combinations of fixed and variable modifications.")

    dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)

    if (index_mods) {
      mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
        as.hexmode() %>%
        `names<-`(c(fixedmods, varmods))
    } else {
      mod_indexes <- NULL
    }

    aa_masses <- calc_aamasses(fixedmods = fixedmods,
                               varmods = varmods,
                               maxn_vmods_setscombi = maxn_vmods_setscombi,
                               mod_indexes = mod_indexes,
                               add_varmasses = FALSE,
                               add_nlmasses = FALSE) %T>%
      saveRDS(file.path(out_path, "temp", "aa_masses_all.rds"))
  } else {
    aa_masses <- readRDS(file.path(out_path, "temp", "aa_masses_all.rds"))
  }

  invisible(aa_masses)
}


#' Finds the existence of \code{max_mass.rds}.
#' 
#' Not used. 
#' 
#' @param .path_cache The file path to cache.
#' @param .time_stamp The time stamp to cache.
#' @inheritParams calc_pepmasses2
find_max_mass <- function (out_path = NULL, .path_fasta = NULL, .time_stamp = NULL) {

  if (!file.exists(file.path(out_path, "temp", "max_mass.rds"))) {

    dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)

    file <- file.path(.path_fasta, "pepmasses", .time_stamp, "pepmasses.rds")

    if (file.exists(file)) {
      fwd_peps <- readRDS(file.path(file))
    } else {
      stop("Cache results not found `", file, "`.\n",
           "Remove `.MSearch/Cache/Calls/calc_pepmasses/", .time_stamp, "` ",
           "and try again.",
           call. = FALSE)
    }

    max_mass <- fwd_peps %>%
      lapply(attr, "data") %>%
      unlist(use.names = FALSE) %>%
      max(na.rm = TRUE) %T>%
      saveRDS(file.path(out_path, "temp", "max_mass.rds"))

    rm(list = c("fwd_peps"))
    gc()
  } else {
    max_mass <- readRDS(file.path(out_path, "temp", "max_mass.rds"))
  }

  invisible(max_mass)
}


#' Finds the site of an AA residue.
#'
#' e.g. 'N-term' in `aa_masses`, not 'Q', for Gln-> pyro-Glu (N-term = Q).
#'
#' @param pos_site A named value. Position in name and site in value.
find_aa_site <- function (pos_site) {

  if (grepl("[NC]{1}-term", names(pos_site))) {
    site <- names(pos_site) %>%
      gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
  } else {
    site <- pos_site
  }
}


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param maxn_vmods_setscombi Integer; the maximum number of combinatorial variable
#'   modifications and neutral losses.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @inheritParams parse_aamasses
#' @inheritParams add_fixvar_masses
#' @examples
#' \donttest{
#' library(purrr)
#'
#' x <- calc_aamasses()
#' x_att <- map(x, attributes)
#' names(x_att[[1]])
#' x_vmods <- map(x_att, `[`, c("vmods"))
#'
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (C)"), c("Acetyl (N-term)",
#'                      "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))
#' 
#' stopifnot(length(x) == 6L)
#'
#' x <- calc_aamasses(fixedmods = NULL)
#' stopifnot(length(x) == 16L)
#' 
#' x <- calc_aamasses(fixedmods = NULL, varmods = NULL)
#' stopifnot(length(x) == 1L)
#'
#' # Fixed mod, no NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
#'                                  "Carbamidomethyl (. = C)"), varmods = NULL)
#' 
#' stopifnot(length(x) == 1L)
#' 
#' # Fixed mod + NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
#'                                  "Carbamidomethyl (. = M)"), varmods = NULL)
#' 
#' stopifnot(length(x) == 1L)
#' 
#' # Fixed mod, no NL; var mod, no NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
#'                                  "Carbamidomethyl (. = C)"),
#'                    varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)"))
#' 
#' stopifnot(length(x) == 3L)
#' 
#' # Fixed mod + NL; var mod + NL
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (. = M)",
#'                      "Deamidated (. = R)"),
#'                    c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)",
#'                      "Hex(5)HexNAc(2) (N)"))
#' 
#' stopifnot(length(x) == 6L)
#' 
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (. = M)", "Deamidated (. = R)"),
#'                    c("Acetyl (N-term)", "Carbamyl (. = M)",
#'                      "Gln->pyro-Glu (N-term = Q)", "Hex(5)HexNAc(2) (N)"))
#' 
#' stopifnot(length(x) == 18L)
#' 
#' ## Coercion       
#' # No fixed terminal or fixed anywhere coercion
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)"),
#'                    c("Carbamidomethyl (M)"))
#' 
#' stopifnot(length(x) == 2L)
#' 
#' x <- calc_aamasses(c("TMT6plex (K)", "Carbamidomethyl (C)"), 
#'                    c("Acetyl (Protein N-term)", "TMT6plex (N-term)", 
#'                      "Oxidation (M)", "Carbamidomethyl (M)"))
#' 
#' stopifnot(length(x) == 12L)
#' 
#' # Fixed terminal coercion
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)"),
#'                    c("Acetyl (Protein N-term)", "Oxidation (M)"))
#'                    
#' stopifnot(length(x) == 4L)
#' 
#' # Fixed anywhere coercion
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
#'                    c("Oxidation (M)"))
#'                    
#' stopifnot(length(x) == 3L)
#'                    
#' # Both fixed terminal and fixed anywhere coercion
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
#'                    c("Acetyl (Protein N-term)", "Oxidation (M)"))
#' 
#' stopifnot(length(x) == 6L)
#' 
#' }
#' \dontrun{
#' # conflicts
#' x <- calc_aamasses(c("Carbamidomethyl (N-term)", "TMT2plex (N-term)"), NULL)
#'
#' # need separate S and T
#' x <- calc_aamasses(NULL, "Phospho (ST)")
#' }
#' @export
calc_aamasses <- function (fixedmods = c("TMT6plex (K)",
                                         "Carbamidomethyl (. = C)"),
                           varmods = c("TMT6plex (N-term)",
                                       "Acetyl (Protein N-term)",
                                       "Oxidation (M)",
                                       "Deamidated (N)",
                                       "Gln->pyro-Glu (N-term = Q)"),
                           maxn_vmods_setscombi = 64,
                           mod_indexes = NULL,
                           add_varmasses = TRUE,
                           add_nlmasses = TRUE) {

  # title (position = site);
  # . stands for (a) anywhere in position or (b) any residue in site or both
  # Acetyl (Protein N-term) <-> Acetyl (Protein N-term = .)
  # Acetyl (N-term) <-> Acetyl (N-term = .)
  # Carbamidomethyl (C) <-> Carbamidomethyl (. = C)
  # Carbamidomethyl <-> Carbamidomethyl(. = .)
  # Gln->pyro-Glu (N-term Q) <-> Gln->pyro-Glu (N-term = Q)
  # anywhere, can be skipped: title (site)
  # N-term, C-term, Protein N-term, Protein C-term

  ## The same site but different mods
  #  (a) Among fixedmods; Failure (only 'one_of', prompt to devise of a joint Unimod)
  #      (a1) TMT6plex (N-term)
  #      (a2) Biotin (N-term)
  #  (b) Between fixedmods and varmods; Relaxation (from fixedmods to varmods)
  #      (b1) Fixed TMT (N-term) -> Variable TMT (N-term)
  #           Variable Acetyl (N-term)
  #      (b2) Fixed Oxidation (M) -> Variable Oxidation (M)
  #           variable Met->Ala (M)
  #  (c) Among varmods; OK if is 'one_of' in combination but not 'multiple'
  #      (c1) Oxidation (M)
  #      (c2) Met->Ala (M)
  #      excludes (c1) + (c2)
  #      OK if different mods to the same residue at different indexes of a peptide
  #      (c3) dHex(1)Hex(1) (S): VS(3)SALSPSK
  #      (c4) Phospho (S): VSS(4)ALSPSK

  options(digits = 9L)

  local({
    ## (0) Duplicated mods
    dup_mods <- intersect(fixedmods, varmods)

    if (length(dup_mods)) {
      stop("Modifications cannot be both 'fixed' and 'variable' at the same time: \n",
           paste(dup_mods, collapse = ", "), 
           call. = FALSE)
    }
  })

  new_mods <- local({
    # (a) Different fixedmods to the same site not allowed
    fmods_ps <- fixedmods %>%
      lapply(find_unimod) %>%
      lapply(`[[`, "position_site") %>%
      `names<-`(fixedmods) %>%
      purrr::flatten()
    
    dup_fixedmods <- fmods_ps %>%
      .[duplicated(.)]

    if (length(dup_fixedmods) > 0L) {
      stop("Multiple fixed modifications to the same site: \n",
           "'", purrr::reduce(dup_fixedmods, paste, sep = ", "), "'",
           call. = FALSE)
    }

    # (b) Coercion from fixedmods to varmods
    vmods_ps <- varmods %>%
      lapply(find_unimod) %>%
      lapply(`[[`, "position_site") %>%
      `names<-`(varmods) %>%
      purrr::flatten()
    
    dup_mods <- intersect(unlist(fmods_ps), unlist(vmods_ps))
    fV_coercion <- (length(dup_mods) > 0L)

    if (fV_coercion) {
      f_to_v <- local({
        idxes <- dup_mods %>% purrr::map(~ fmods_ps == .x)
        idxes %>% purrr::map_chr(~ fixedmods[.x])
      })

      varmods <- c(varmods, f_to_v)

      fixedmods <- local({
        idxes <- fixedmods %>% purrr::map_lgl(~ .x %in% f_to_v)
        fixedmods[!idxes]
      })

      warning("Coerce '",
              purrr::reduce(f_to_v, paste, sep = ", "), "'",
              " to conditional variable modifications.",
              call. = FALSE)
      
      fixednt_coerced <- any(grepl("N-term", f_to_v))
      fixedct_coerced <- any(grepl("C-term", f_to_v))
      anywhere_coreced_sites <- dup_mods %>% .[!grepl("[NC]-term", .)]
    } else {
      fixednt_coerced <- FALSE
      fixedct_coerced <- FALSE
      anywhere_coreced_sites <- NULL
    }

    default_mods <- c("initiator methionine from protein N-terminus")

    if (any(grepl(default_mods, c(fixedmods, varmods)))) {
      warning("Modifications defaulted and no need to specify: `\n",
              default_mods, "`.\n",
              call. = FALSE)

    }

    invisible(list(fixedmods = fixedmods, 
                   varmods = varmods, 
                   fixednt_coerced = fixednt_coerced, 
                   fixedct_coerced = fixedct_coerced, 
                   anywhere_coreced_sites = anywhere_coreced_sites, 
                   fV_coercion = fV_coercion))
  })

  fixedmods <- new_mods$fixedmods
  varmods <- new_mods$varmods
  fixednt_coerced <- new_mods$fixednt_coerced
  fixedct_coerced <- new_mods$fixedct_coerced
  anywhere_coreced_sites <- new_mods$anywhere_coreced_sites
  fV_coercion <- new_mods$fV_coercion
  rm(new_mods)

  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943,
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464,
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414,
    "N-term" = 1.007825, "C-term" = 17.002740,
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590,
    "-" = 0)

  ## (1) add fixed mods + NL
  aa_masses_fi2 <- add_fixvar_masses(mods = fixedmods,
                                     mod_type = "fmods",
                                     aa_masses = aa_masses,
                                     add_varmasses = add_varmasses)

  aa_masses_fi2 <- aa_masses_fi2 %>%
    purrr::map(~ {
      if (is.null(attr(.x, "fmods"))) {
        attr(.x, "fmods") <- ""
      }

      if (is.null(attr(.x, "fmods_neuloss"))) {
        attr(.x, "fmods_neuloss") <- ""
      }

      if (is.null(attr(.x, "fmods_mass"))) {
        attr(.x, "fmods_mass") <- 0
      }

      .x
    })

  ## (2) add variable mods + NL
  varmods_comb <- local({

    # (c) Remove entries with multiple terminal mods
    varmods_comb <- seq_along(varmods) %>%
      purrr::map(~ combn(varmods, .x, simplify = FALSE)) %>%
      purrr::flatten()

    # Check if multiple terminal varmods by names
    #   Gln->pyro Glu (N-term = Q) and Acetyl (Protein N-term = N-term)
    #   have different sites but 'N-term' in names; and also need to be excluded

    vmods_ps <- varmods %>%
      purrr::map(find_unimod) %>%
      purrr::map(`[[`, "position_site") %>%
      `names<-`(varmods) %>%
      purrr::flatten()

    vmods_ps_combi <- seq_along(vmods_ps) %>%
      purrr::map(~ combn(vmods_ps, .x, simplify = FALSE)) %>%
      purrr::flatten()

    dup_terms <-
      vmods_ps_combi %>% purrr::map_lgl(~ sum(grepl("N-term", names(.x))) >= 2L) |
      vmods_ps_combi %>% purrr::map_lgl(~ sum(grepl("C-term", names(.x))) >= 2L)

    # Not currently used
    dup_anywhere <- vmods_ps_combi %>%
      purrr::map_lgl(~ {
        .x <- .x %>%
          .[!grepl("[NC]{1}-term", names(.))]

        .x %>%
          duplicated() %>%
          any()
      })

    # Allow entries with different Anywhere mods to the same site
    #   (1) dHex(1)Hex(1) (S) and (2) Phospho (S)
    #   VS(1)S(2)ALSPSK

    varmods_comb <- varmods_comb  %>%
      .[!(dup_terms)] %>%
      purrr::map(unlist)
  })
  
  # respect users' choices: e.g., if `Fixed Anywhere N-term` by users 
  #   -> must have `N-term` in any realization of c(fixedmods, varmods)
  
  if (fixednt_coerced && length(varmods_comb)) {
    ok_nts <- sapply(varmods_comb, function (x) any(grepl("N-term", x)))
    varmods_comb <- varmods_comb[ok_nts]
    rm(list = c("ok_nts"))
  }
  
  if (fixedct_coerced && length(varmods_comb)) {
    ok_cts <- sapply(varmods_comb, function (x) any(grepl("C-term", x)))
    varmods_comb <- varmods_comb[ok_cts]
    rm(list = c("ok_cts"))
  }

  # return NULL if is.null(varmods_comb)
  aa_masses_var2 <- purrr::map(varmods_comb, ~ {
    varmods_i <- .x

    aa_masses_fi2 %>%
      purrr::map(~ add_fixvar_masses(mods = varmods_i,
                                     mod_type = "vmods",
                                     aa_masses = .x,
                                     add_varmasses = add_varmasses)) %>%
      purrr::flatten()
  }, aa_masses_fi2) %>%
    purrr::flatten()

  aa_masses_var2 <- aa_masses_var2 %>%
    purrr::map(~ {
      if (is.null(attr(.x, "vmods"))) {
        attr(.x, "vmods") <- ""
      }

      if (is.null(attr(.x, "vmods_neuloss"))) {
        attr(.x, "vmods_neuloss") <- ""
      }

      if (is.null(attr(.x, "vmods_mass"))) {
        attr(.x, "vmods_mass") <- 0
      }

      .x
    })

  ## (3) complete 'vmods', 'vmods_neuloss' and 'vmods_ps' to fixedmods
  aa_masses_fi2 <- aa_masses_fi2 %>%
    purrr::map(~ {
      if (is.null(attr(.x, "vmods"))) {
        attr(.x, "vmods") <- ""
      }

      if (is.null(attr(.x, "vmods_neuloss"))) {
        attr(.x, "vmods_neuloss") <- ""
      }

      if (is.null(attr(.x, "vmods_ps"))) {
        attr(.x, "vmods_ps") <- ""
      }

      if (is.null(attr(.x, "vmods_mass"))) {
        attr(.x, "vmods_mass") <- 0
      }

      .x
    })

  if (length(aa_masses_var2) >= maxn_vmods_setscombi) {
    warning("The ways of combinatorial variable modifications are ",
            length(aa_masses_var2), ".\n",
            "Dropping combinations at indexes greater than `maxn_vmods_setscombi = ", 
            maxn_vmods_setscombi, "`.",
            call. = FALSE)
    aa_masses_var2 <- aa_masses_var2[1:maxn_vmods_setscombi]
  }

  if (!is.null(mod_indexes)) {
    aa_masses_var2 <- aa_masses_var2 %>%
      purrr::map(~ {
        nms <- mod_indexes %>%
          .[names(.) %in% names(.x)]

        idxes <- which(names(.x) %in% names(nms))
        names(.x)[idxes] <- nms

        .x
      }, mod_indexes)

    aa_masses_fi2 <- aa_masses_fi2 %>%
      purrr::map(~ {
        names(attr(.x, "fmods_ps")) <- mod_indexes[fixedmods]
        names(attr(.x, "fmods_mass")) <- mod_indexes[fixedmods]

        .x
      })

    aa_masses_var2 <- aa_masses_var2 %>%
      purrr::map(~ {
        v_nms <- names(attr(.x, "vmods_ps"))
        names(attr(.x, "vmods_ps")) <- mod_indexes[v_nms]
        names(attr(.x, "vmods_mass")) <- mod_indexes[v_nms]

        f_nms <- names(attr(.x, "fmods_ps"))
        names(attr(.x, "fmods_ps")) <- mod_indexes[f_nms]
        names(attr(.x, "fmods_mass")) <- mod_indexes[f_nms]

        .x
      })
  }

  ## Coerced sites, including "[NC]-term", need to be present in a combination
  # 
  # with coercion, `aa_masses_fi2` becomes something not specified by users:
  #   TMT6plex (N-term) fixed -> variable; the `aa_masses_fi2` 
  #   c(""TMT6plex (K)", "Carbamidomethyl (C)") without `TMT6plex (N-term)` 
  #   is not a combination intended by users)

  if (fV_coercion) {
    aa_masses_all <- aa_masses_var2
  } else {
    aa_masses_all <- c(aa_masses_fi2, aa_masses_var2)
  }

  if (length(anywhere_coreced_sites)) {
    vmods_ps <- purrr::map(aa_masses_all, attr, "vmods_ps")
    
    for (i in seq_along(anywhere_coreced_sites)) {
      aa_masses_all <- check_anywhere_fmods_coercion(anywhere_coreced_sites[i], 
                                                     aa_masses_all, vmods_ps)
    }
    
    rm(list = c("vmods_ps"))
  }

  aa_masses_all <- purrr::map(aa_masses_all, parse_aamasses, add_nlmasses)
}


#' Checks fixed modifications after coercion to variable modifications.
#' 
#' The coerced site needs to be present in final \code{aa_masses_all}.
#' 
#' @param site An amino acid site; e.g. site = "M".
#' @param aa_masses_all All the amino acid lookup tables.
#' @param vmods_ps All of the "vmods_ps" attributes from "aa_masses_all".
check_anywhere_fmods_coercion <- function (site, aa_masses_all, vmods_ps) {
  
  oks <- purrr::map_lgl(vmods_ps, ~ any(grepl(site, .x)))
  aa_masses_all <- aa_masses_all[oks]
  
  if (!length(aa_masses_all)) {
    stop("Zero combination of AA tables ", 
         "Check `fixedmods` and `varmods`.")
  }
  
  aa_masses_all
}


#' Helper to add modification masses to amino-acid residues.
#'
#' It adds the masses of fixed, variable, and neutral-loss modifications to
#' amino-acid residues.
#'
#' @param mods A list of modifications.
#' @param mod_type The type of modification in one of \code{c("fmods", "vmods")}
#'   where \code{fmods}: fixed modifications and \code{vmods}: variable
#'   modifications.
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param add_varmasses Logical; if TRUE, add mass delta to the corresponding
#'   base mass for variable modifications. The argument at TRUE is for
#'   compatibility for MS1 precursor mass calculations before the approach of
#'   "rolling sum + fallthrough".
#' @return Lists of of amino-acid residues with modified mono-isotopic masses
#'   being incorporated.
add_fixvar_masses <- function (mods, mod_type, aa_masses, add_varmasses = TRUE) {

  stopifnot(mod_type %in% c("fmods", "vmods"),
            length(mod_type) == 1L)

  if (length(mods)) {
    all_mods <- paste(mods, collapse = ", ")
  } else {
    all_mods <- ""
  }

  res <- mods %>%
    purrr::map(find_unimod) %>%
    `names<-`(mods)

  mod_masses <- res %>% purrr::map(`[[`, 1)
  positions_sites <- res %>% purrr::map(`[[`, 2)
  neulosses <- res %>% purrr::map(`[[`, 3)
  rm(res)

  # the same `site` with different fixedmods
  local({
    if (mod_type == "fmods" && length(positions_sites) > 1L) {
      dups <- purrr::reduce(positions_sites, `c`) %>%
        .[duplicated(.)]

      if (length(dups) > 0L) {
        dups_in_each <- positions_sites %>%
          purrr::map(~ .x[.x == dups])

        warning("Conflicts in fixed modifications: \n",
                purrr::reduce(names(dups_in_each), paste, sep = "\n"), "\n",
                "May consider change from fixed to variable modifications(s); \n",
                "or create a new Unimod for joint modifications.",
                call. = FALSE)
      }
    }
  })

  if (mod_type == "fmods") {
    # Add mod_masses of fixed mods
    purrr::walk2(positions_sites, mod_masses, ~ {
      site <- find_aa_site(.x)
      m <- aa_masses[site]
      aa_masses[site] <<- m + .y
    })
  } else {
    #  Add mod_masses of variable mods (multiple lists)
    aas <- purrr::map2(positions_sites, mod_masses, ~ {
      site <- find_aa_site(.x)

      if (add_varmasses) {
        aa_masses[site] <- aa_masses[site] + .y
      } else {
        aa_masses[site] <- .y
      }

      aa_masses
    }, aa_masses)

    # Flatten the lists (with attributes being kept)
    aa_masses <- local({
      masses <- purrr::map2_dbl(positions_sites, aas, ~ .y[find_aa_site(.x)])
      attrs <- attributes(aa_masses)
      aa_masses <- c(aa_masses, masses)
      attrs$names <- names(aa_masses)
      attributes(aa_masses) <- attrs

      aa_masses
    })

    rm(list = c("aas"))
  }

  attr(aa_masses, mod_type) <- all_mods
  attr(aa_masses, paste0(mod_type, "_ps")) <- positions_sites
  attr(aa_masses, paste0(mod_type, "_mass")) <- mod_masses

  # add mod_masses of neutral losses

  no_nls <- neulosses %>%
    purrr::map_lgl(~ all(.x == 0)) %>%
    all()

  if (no_nls) {
    return(list(aa_masses))
  }

  attr(aa_masses, paste0(mod_type, "_neuloss")) <- neulosses

  list(aa_masses)
}


#' Parses \code{aa_masses}.
#'
#' @inheritParams add_fixvar_masses
#' @param add_nlmasses Logical; if TRUE, add an original mass to the
#'   corresponding (negative) neutral-loss mass.
parse_aamasses <- function (aa_masses, add_nlmasses = TRUE) {

  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)

  fmods_nl <- local({
    neulosses <- attr(aa_masses, "fmods_neuloss", exact = TRUE)

    if (all(neulosses == "")) return(character())

    # add `0` if absent
    no_zero <- purrr::map_lgl(neulosses, ~ !any(.x == 0)) %>%
      which()

    if (length(no_zero) > 0L) {
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])
    }

    ## entries with NL = 0 also kept
    ## (beneficial when calling `find_intercombi`)
    # idx <- map_lgl(neulosses, ~ length(.x) > 1)
    # neulosses <- neulosses[idx]

    ## In fixedmods: `M` instead of `Oxidation (M)`
    # fmods_ps <- fmods_ps[idx]

    names(neulosses) <- fmods_ps

    if (add_nlmasses) {
      resids <- as.list(aa_masses)

      fmods_nl <- purrr::imap(neulosses, ~ {
        resids[[.y]] - .x
      })
    } else {
      fmods_nl <- neulosses
    }

    invisible(fmods_nl)
  })

  vmods_nl <- local({
    neulosses <- attr(aa_masses, "vmods_neuloss", exact = TRUE)

    if (all(neulosses == "")) return(character())

    # add `0` if absent
    no_zero <- purrr::map_lgl(neulosses, ~ !any(.x == 0)) %>% which()

    if (length(no_zero)) {
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])
    }

    ## entries of "NL = 0" kept
    # idx <- map_lgl(neulosses, ~ length(.x) > 1)
    # neulosses <- neulosses[idx]

    if (add_nlmasses) {
      resids <- as.list(aa_masses)

      vmods_nl <- purrr::imap(neulosses, ~ {
        resids[[.y]] - .x
      })
    } else {
      vmods_nl <- neulosses
    }

    invisible(vmods_nl)
  })

  ## variable mods
  # multiple mods to [NC]-term already excluded from aa_masses

  amods <- local({
    sites <- vmods_ps %>%
      purrr::map(~ .x[grepl("Anywhere", names(.x))])

    empties <- sites %>% purrr::map_lgl(purrr::is_empty)

    sites <- sites[!empties]

    sites
  })

  # `TMT6plex (N-term)`
  # `Amidated (Protein C-term)`

  tmod <- vmods_ps %>% .[! . %in% amods]
  if (length(tmod) == 0L) {
    tmod <- NULL
  } else if (all(tmod == "")) {
    tmod <- NULL
  }

  # variable N-term, C-term
  #
  # $`Gln->pyro-Glu (N-term = Q)`
  # Any N-term
  # "Q"
  #
  # `tmod` is a length(1) list; R `==` and `grepl` works too
  # (OK to use either tmod or tmod[[1]])

  if (length(tmod) <= 1L) {
    ntmod <- tmod %>% .[. == "N-term" || grepl("N-term", names(.))]
    ctmod <- tmod %>% .[. == "C-term" || grepl("C-term", names(.))]
  } else if (length(tmod) == 2L) {

    # `TMT6plex (N-term)`
    # `Amidated (Protein C-term)`
    ntmod <- local({
      x <- purrr::map(tmod, ~ .x %>% .[. == "N-term" || grepl("N-term", names(.))])
      rows <- purrr::map_lgl(x, purrr::is_empty)
      x[!rows]
    })

    ctmod <- local({
      x <- purrr::map(tmod, ~ .x %>% .[. == "C-term" || grepl("C-term", names(.))])
      rows <- purrr::map_lgl(x, purrr::is_empty)
      x[!rows]
    })
  } else {
    stop("cannot have more than two terminal modifications.",
         call. = FALSE)
  }

  ## fixed mods
  famods <- local({
    sites <- fmods_ps %>%
      purrr::map(~ .x[grepl("Anywhere", names(.x))])

    empties <- sites %>% purrr::map_lgl(purrr::is_empty)

    sites <- sites[!empties]

    sites
  })

  ftmod <- fmods_ps %>% .[! . %in% famods]
  if (length(ftmod) == 0L) {
    ftmod <- NULL
  } else if (ftmod == "") {
    ftmod <- NULL
  }

  # fixed N-term, C-term
  fntmod <- ftmod %>% .[. == "N-term" || grepl("N-term", names(.))]
  fctmod <- ftmod %>% .[. == "C-term" || grepl("C-term", names(.))]

  # "amods- tmod- vnl- fnl-"
  if (length(fmods_nl) == 0L) {
    type <- "fnl-"
    if (length(vmods_nl) == 0L) {
      type <- paste("vnl-", type)
      if (length(tmod) == 0L) {
        type <- paste("tmod-", type)
        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 1
        } else {
          type <- paste("amods+", type) # 2
        }
      } else {
        type <- paste("tmod+", type)
        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 3
        } else {
          type <- paste("amods+", type) # 4
        }
      }
    } else {
      type <- paste("vnl+", type)
      if (length(tmod) == 0L) {
        type <- paste("tmod-", type)
        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 5
        } else {
          type <- paste("amods+", type) # 6
        }
      } else {
        type <- paste("tmod+", type)

        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 7
        } else {
          type <- paste("amods+", type) # 8
        }
      }
    }
  } else {
    type <- "fnl+"
    if (length(vmods_nl) == 0L) {
      type <- paste("vnl-", type)
      if (length(tmod) == 0L) {
        type <- paste("tmod-", type)
        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 1
        } else {
          type <- paste("amods+", type) # 2
        }
      } else {
        type <- paste("tmod+", type)

        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 3
        } else {
          type <- paste("amods+", type) # 4
        }
      }
    } else {
      type <- paste("vnl+", type)
      if (length(tmod) == 0L) {
        type <- paste("tmod-", type)
        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 5
        } else {
          type <- paste("amods+", type) # 6
        }
      } else {
        type <- paste("tmod+", type)

        if (length(amods) == 0L) {
          type <- paste("amods-", type) # 7
        } else {
          type <- paste("amods+", type) # 8
        }
      }
    }
  }

  attr(aa_masses, "type") <- type

  attr(aa_masses, "fmods_nl") <- fmods_nl
  attr(aa_masses, "famods") <- famods
  attr(aa_masses, "ftmod") <- ftmod
  attr(aa_masses, "fntmod") <- fntmod
  attr(aa_masses, "fctmod") <- fctmod

  attr(aa_masses, "vmods_nl") <- vmods_nl
  attr(aa_masses, "amods") <- amods
  attr(aa_masses, "tmod") <- tmod
  attr(aa_masses, "ntmod") <- ntmod
  attr(aa_masses, "ctmod") <- ctmod

  invisible(aa_masses)
}


#' Helper of \link{calc_pepmasses}.
#'
#' Prior to the calculation of peptide masses; for the base with fixed
#' modification only.
#'
#' @inheritParams calc_pepmasses
#' @inheritParams mcalc_monopep
#' @importFrom stringi stri_reverse
#' @return Two named list of "fwds" and "revs". List "fwds" contains peptide
#'   sequences split from forward fasta and "revs" from reversed fasta.
split_fastaseqs <- function (fasta, acc_type, acc_pattern, maxn_fasta_seqs,
                             max_miss = 2L) {

  # ---
  message("Loading fasta databases.")

  fasta_db <- load_fasta2(fasta, acc_type, acc_pattern)

  if (length(fasta_db) > maxn_fasta_seqs) {
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.\n",
         "  May consider a higher `maxn_fasta_seqs`.",
         call. = FALSE)
  }

  n_cores <- detect_cores()
  cl <- makeCluster(getOption("cl.cores", n_cores))

  clusterExport(cl, list("%>%"), 
                envir = environment(magrittr::`%>%`))
  clusterExport(cl, list("make_fastapeps0"), 
                envir = environment(proteoM:::make_fastapeps0))
  clusterExport(cl, list("keep_n_misses"), 
                envir = environment(proteoM:::keep_n_misses))

  fasta_db <- chunksplit(fasta_db, n_cores)

  # ---
  message("Splitting fasta sequences.")

  peps <- clusterApply(cl, fasta_db, make_fastapeps0, max_miss) %>%
    purrr::flatten()

  rm(list = c("fasta_db"))
  stopCluster(cl)

  invisible(peps)
}


#' Make peptide sequences from FASTA databases.
#'
#' Not yet concatenating peptides by the number of mis-cleavages.
#'
#' @param fasta_db Fasta database(s).
#' @inheritParams calc_pepmasses
make_fastapeps0 <- function (fasta_db, max_miss = 2L) {

  inds_m <- grep("^M", fasta_db)
  
  fasta_db <- lapply(fasta_db, function (x) {
    s <- gsub("([KR]{1})", paste0("\\1", "@"), x)
    paste0("-", s, "-")
  })
  
  fasta_dbm <- lapply(fasta_db[inds_m], function (x) {
    gsub("^-M", "-", x)
  })
  
  # --- with protein N-term (initiator) methionine ---
  peps <- lapply(fasta_db, function (x) {
    stringr::str_split(x, "@", simplify = TRUE)
  })
  
  # --- without protein N-term (initiator) methionine ---
  peps_m <- lapply(fasta_dbm, function (x) {
    s <- stringr::str_split(x, "@", simplify = TRUE)
    keep_n_misses(s, max_miss)
  })
  
  # Note: NA sequences -> NULL during mass calculations
  # (USE.NAMEs as they are prot_acc)
  peps[inds_m] <- mapply(list, peps_m, peps[inds_m], 
                         SIMPLIFY = FALSE)
  peps[-inds_m] <- lapply(peps[-inds_m], function (x) list(NA, x))

  rm(list = c("inds_m", "fasta_db", "fasta_dbm", "peps_m"))

  invisible(peps)
}


#' Helper in calculating peptide masses.
#'
#' (2) "amods- tmod+ vnl- fnl-".
#'
#' @inheritParams distri_peps
add_term_mass <- function (aa_masses, peps) {

  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)

  # No needs of is_empty(ntmod) && is_empty(ctmod)
  if (length(ntmod) && length(ctmod)) {
    delta <- aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  } else if (length(ntmod)) {
    delta <- aa_masses[names(ntmod)]
  } else if (length(ctmod)) {
    delta <- aa_masses[names(ctmod)]
  }

  out <- lapply(peps, `+`, delta)
}


#' Calculates mono-isotopic masses of peptide sequences.
#'
#' @param seqs Sequences of peptides from FASTAs by protein accessions. Each
#'   list contains two lists of sequences: (1) without and (2) with N-terminal
#'   methionine.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses
#' @inheritParams mcalc_monopep
#' @inheritParams distri_fpeps
ms1masses_bare <- function (seqs, aa_masses, ftmass = NULL,
                            max_miss = 2L, min_len = 7L, max_len = 100L,
                            maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                            is_fixed_protnt = FALSE, is_fixed_protct = FALSE,
                            parallel = TRUE, digits = 4L) {

  # (1) before rolling sum (not yet terminal H2O)
  # (1.1) without N-term methionine
  data_1 <- seqs %>%
    lapply(`[[`, 1) %>%
    ms1masses_noterm(aa_masses = aa_masses,
                     maxn_vmods_per_pep = maxn_vmods_per_pep,
                     maxn_sites_per_vmod = maxn_sites_per_vmod,
                     parallel = parallel,
                     digits = digits) %>%
    attr("data")

  # (1.2) with N-term methionine
  data_2 <- seqs %>%
    lapply(`[[`, 2) %>%
    ms1masses_noterm(aa_masses = aa_masses,
                     maxn_vmods_per_pep = maxn_vmods_per_pep,
                     maxn_sites_per_vmod = maxn_sites_per_vmod,
                     parallel = parallel,
                     digits = digits) %>%
    attr("data")

  # (2) rolling sum (not yet terminal H2O)
  n_cores <- detect_cores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

  parallel::clusterExport(cl, list("%>%"), 
                          envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("roll_sum"),
                          envir = environment(proteoM:::roll_sum))

  ms_1 <- parallel::clusterApply(
    cl = cl, 
    x = chunksplit(data_1, n_cores, "list"), 
    fun = lapply, 
    FUN = "roll_sum", 
    n = max_miss, 
    include_cts = FALSE
  ) %>% 
    purrr::flatten()

  if (is_fixed_protnt) {
    ms_2 <- parallel::clusterApply(
      cl = cl, 
      x = chunksplit(data_2, n_cores, "list"), 
      fun = lapply, 
      FUN = "roll_sum", 
      n = max_miss, 
      include_cts = FALSE
    ) %>% 
      purrr::flatten()
  } else {
    ms_2 <- parallel::clusterApply(
      cl = cl, 
      x = chunksplit(data_2, n_cores, "list"), 
      fun = lapply, 
      FUN = "roll_sum", 
      n = max_miss, 
      include_cts = TRUE
    ) %>% 
      purrr::flatten()
  }

  stopCluster(cl)

  # (3) putting together (+ terminal H2O)
  # (USE.NAMES of prot_acc)
  ms <- mapply(`c`, ms_1, ms_2, SIMPLIFY = FALSE)

  if (min_len > 1L && !is.infinite(max_len)) {
    ms <- lapply(ms, function (x) {
      x %>% .[str_exclude_count(names(.)) >= min_len &
                 str_exclude_count(names(.)) <= max_len]
    })
  }

  # "HD101_HUMAN" etc. has no tryptic peptides
  lens <- simplify2array(lapply(ms, length))

  # adding H2O or fixed N/C-term masses
  ms <- ms %>%
    .[lens > 0L] %>%
    lapply(function (x) x[!duplicated(names(x))]) %>%
    lapply(function (x) x + ftmass) # 18.010565

  invisible(ms)
}


#' Helper of \link{ms1masses_bare}.
#'
#' For either forward or reversed sequences.
#' 
#' @param aa_seqs Character string; a vector of peptide sequences with
#'   one-letter representation of amino acids.
#' @param parallel Logical; if TRUE, performs parallel computation. The default
#'   is TRUE.
#' @inheritParams matchMS
#' @inheritParams distri_peps
ms1masses_noterm <- function (aa_seqs, aa_masses,
                              maxn_vmods_per_pep = 5L,
                              maxn_sites_per_vmod = 3L,
                              parallel = TRUE,
                              digits = 4L) {

  options(digits = 9L)

  n_cores <- detect_cores()

  aa_seqs <- suppressWarnings(split(aa_seqs, seq_len(n_cores)))

  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

  parallel::clusterExport(cl, list("%>%"),
                          envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("map"),
                          envir = environment(purrr::map))
  parallel::clusterExport(cl, list("str_split"),
                          envir = environment(stringr::str_split))

  parallel::clusterExport(
    cl,
    c("calcms1mass_noterm", 
      "calcms1mass_noterm_byprot", 
      "calcms1mass_noterm_bypep"), 
    envir = environment(proteoM:::calcms1mass_noterm))

  out <- parallel::clusterApply(cl, aa_seqs, calcms1mass_noterm,
                                aa_masses = aa_masses,
                                maxn_vmods_per_pep = maxn_vmods_per_pep,
                                maxn_sites_per_vmod = maxn_sites_per_vmod,
                                digits = digits) %>%
    purrr::flatten()

  parallel::stopCluster(cl)

  attr(aa_masses, "data") <- out

  rm(list = c("aa_seqs", "out"))
  gc()

  invisible(aa_masses)
}


#' Helper function for parallel calculations peptide masses by proteins.
#'
#' For each split of multiple proteins; no terminal masses.
#'
#' @param aa_seqs Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calcms1mass_noterm <- function (aa_seqs, aa_masses,
                                maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                                digits = 4L) {

  lapply(aa_seqs, function (prot_peps) {
    calcms1mass_noterm_byprot(prot_peps = prot_peps,
                              aa_masses = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod,
                              digits = digits)
  })
}


#' Helper of \link{calcms1mass_noterm}.
#'
#' For single protein.
#'
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calcms1mass_noterm_byprot <- function (prot_peps, aa_masses,
                                       maxn_vmods_per_pep = 5L,
                                       maxn_sites_per_vmod = 3L,
                                       digits = 4L) {

  prot_peps %>%
    lapply(calcms1mass_noterm_bypep,
           aa_masses = aa_masses,
           maxn_vmods_per_pep = maxn_vmods_per_pep,
           maxn_sites_per_vmod = maxn_sites_per_vmod,
           digits = digits) %>% # by peptides
    unlist(use.names = TRUE)
}


#' Helper of \link{calcms1mass_noterm_byprot}.
#'
#' For single protein.
#'
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
calcms1mass_noterm_bypep <- function (aa_seq, aa_masses,
                                      maxn_vmods_per_pep = 5L,
                                      maxn_sites_per_vmod = 3L,
                                      digits = 4L) {

  if (is.na(aa_seq)) return(NULL)

  aa_seq %>%
    stringr::str_split("", simplify = TRUE) %>%
    aa_masses[.] %>%
    sum() %>%
    setNames(aa_seq) %>%
    round(digits = digits)
}



#' Distributes peptides by variable modifications.
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param peps Lists of peptides, either "fwds" or "revs", from
#'   \link{split_fastaseqs}.
#' @inheritParams calc_pepmasses
distri_peps <- function (aa_masses, peps, max_miss) {

  nms <- lapply(peps, names)

  out <- aa_masses %>%
    lapply(subpeps_by_vmods, nms)
  
  # USE.NAMEs of prot_acc
  out <- out %>%
    lapply(function (xs) {
      idxes <- mapply(fastmatch::fmatch, xs, nms, SIMPLIFY = FALSE)
      mapply(function (x, y) x[y], peps, idxes, SIMPLIFY = FALSE)
    })
  
  n2 <- ct_counts(max_miss)
  n1 <- (max_miss + 1L) * 2L
  
  out <- lapply(out, function(xs) {
    # ZN207_HUMAN: MGRKKKK (no N-term pep_seq at 2 misses and min_len >= 7L)
    len <- simplify2array(lapply(xs, length))
    xs <- xs[len > 0L]
    
    xs %>%
      lapply(rm_char_in_nfirst2, char = "^-", n = n1) %>%
      lapply(rm_char_in_nlast2, char = "-$", n = n2)
  })
}


#' Counts the number of trailing residues from C-term for the replacment of "-".
#'
#' n(i+1) = n(i) + (i+1)
#'
#' @param max_miss The maximum number of cleavages.
ct_counts <- function (max_miss = 2L) {

  ct <- integer(max_miss)

  if (max_miss > 0L) {
    ct[1] <- 2L

    for (i in 1:max_miss) {
      j <- i + 1L
      ct[j] <- ct[i] + j
    }
  } else {
    return(1L)
  }

  ct[max_miss]
}


#' Distributes peptides by fixed modifications.
#'
#' @inheritParams calc_pepmasses2
#' @param is_fixed_protnt Logical; is protein N-terminal modification fixed.
#' @param is_fixed_protct Logical; is protein C-terminal modification fixed.
#' @param data Lists of peptides by prot_accs.
distri_fpeps <- function (data, max_miss, is_fixed_protnt, is_fixed_protct) {

  if (is_fixed_protnt) {
    warning("At fixed `Protein N-term`, ",
            "non N-term peptides are removed.\n",
            "!!! Consider variable `Protein N-term` modifications. !!!",
            call. = FALSE)

    # peps: List of 2
    data <- lapply(data, function (peps) { 
      b <- peps[[2]]
      len <- min(max_miss + 1L, length(b))
      peps[[2]] <- b[1:len]

      peps
    })
  }

  if (is_fixed_protct) {
    warning("At fixed `Protein C-term`, ",
            "non C-term peptides are removed.\n",
            "!!! Consider variable `Protein C-term` modifications. !!!",
            call. = FALSE)

    data <- lapply(data, function (peps) {
      
      # Without N-term methionine
      a <- peps[[1]]
      end <- length(a)

      if (end > 0L && grepl("-$", a[end])) {
        start <- max(1L, end - max_miss)
        peps[[1]] <- a[start:end]
      }

      # With N-term methonine
      b <- peps[[2]]
      end <- length(b)

      if (end > 0L) {
        start <- max(1L, end - max_miss)
        peps[[2]] <- b[start:end]
      }

      peps
    })
  }

  invisible(data)
}




#' Helper of \link{calcms1mass_noterm}.
#'
#' By individual proteins.
#'
#' @param fnl_combi A data.frame of combinations of neutral losses for fixed
#'   modifications. Each row corresponds to a set of neutral loss. The first row
#'   corresponds to the combination without NLs (all zeros).
#' @param prot_peps Lists of named peptide masses under a protein.
#' @inheritParams calc_monopep
ms1_a0_fnl1_byprot <- function (prot_peps, fnl_combi, aa_masses, digits = 4L) {

  # Not yet tested USE.NAMES
  out <- mapply(ms1_a0_fnl1_bypep, prot_peps, names(prot_peps), 
                MoreArgs = list(
                  fnl_combi = fnl_combi,
                  aa_masses = aa_masses,
                  digits = digits, 
                ), SIMPLIFY = FALSE)
  
  # Not yet tested USE.NAMES
  nms <- mapply(function (x, y) rep(y, length(x)), out, names(out), 
                SIMPLIFY = FALSE) %>%
    unlist(recursive = FALSE, use.names = FALSE)
  
  out <- out %>%
    unlist(recursive = FALSE, use.names = FALSE) %>%
    `names<-`(nms)
}


#' Helper of \link{calcms1mass_noterm_byprot}.
#'
#' By individual peptides.
#'
#' @param mass The mass of a peptide.
#' @inheritParams ms1_a0_fnl1_byprot
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
ms1_a0_fnl1_bypep <- function (mass, aa_seq, fnl_combi, aa_masses, digits = 4L) {

  aas <- str_split::str_split(aa_seq, "", simplify = TRUE)

  delta <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)

  round(mass - delta, digits = digits)
}


#' Helper of peptide-mass calculation..
#'
#' (5) "amods- tmod+ vnl- fnl+"; (6) "amods- tmod- vnl- fnl+".
#'
#' The calculation goes through the rows in \code{fnl_combi}.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams ms1_a0_fnl1_byprot
#' @return A numeric vector
delta_ms1_a0_fnl1 <- function (fnl_combi, aas, aa_masses) {

  # not an option but always: include_insource_nl = TRUE
  # bring back the option later...

  nms <- colnames(fnl_combi)
  oks <- aas[aas %in% nms]

  if (!length(oks)) return (0L)

  len <- nrow(fnl_combi)
  out <- vector("numeric", len)
  out[[1]] <- 0L

  for (i in 2:len) {
    aa_masses[nms] <- unlist(fnl_combi[i, ])
    oks <- aas[aas %in% nms]
    out[[i]] <- sum(aa_masses[oks])
  }

  out <- unique(out)
}


#' Helper by individual proteins.
#'
#' (7) "amods+ tmod- vnl- fnl-"; (8) "amods+ tmod-+ vnl- fnl-".
#'
#' @param amods \code{Anywhere} variable modifications.
#' @inheritParams ms1_a0_fnl1_byprot
#' @inheritParams calc_monopep
#'
#' @examples
#' \donttest{
#' m0 <- calc_monopeptide("HQGVMCNVGMGQKMNSC", NULL, NULL)
#' stopifnot(unlist(m0$mass, use.names = FALSE) - 1822.7405 < 1e-4)
#'
#' m1 <- calc_monopeptide("HQGVMCNVGMGQKMNSC", "TMT6plex (N-term)", NULL)
#' stopifnot(unlist(m1$mass, use.names = FALSE) - 2051.9035 < 1e-4)
#'
#' # (7) "amods+ tmod- vnl- fnl-"
#' aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"),
#'                                varmods = c("Deamidated (N)",
#'                                            "Carbamidomethyl (C)"),
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' pep <- c("HQGVMCNVGMGQKMNSC" = 2051.90346)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"),
#'               `Carbamidomethyl (C)` = c("Anywhere" = "C"))
#'
#' x <- ms1_a1_vnl0_fnl0_byprot(pep, amods, aa_masses_all[[4]])
#'
#' stopifnot(x[[1]] - 2109.9089 < 1e-4,
#'           x[[2]] - 2166.9304 < 1e-4,
#'           x[[3]] - 2110.8930 < 1e-4,
#'           x[[4]] - 2167.9144 < 1e-4)
#'
#'
#' # (8) "amods+ tmod+ vnl- fnl-"
#' aa_masses_all <- calc_aamasses(fixedmods = "TMT6plex (K)",
#'                                varmods = c("Deamidated (N)",
#'                                            "Carbamidomethyl (S)",
#'                                            "Acetyl (Protein N-term)"),
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' pep <- c("HQGVMNVGMGQKSMNS" = 1932.9171) # + TMT6plex (K)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"),
#'               `Carbamidomethyl (S)` = c("Anywhere" = "S"))
#'
#' x <- ms1_a1_vnl0_fnl0_byprot(pep, amods, aa_masses_all[[8]])
#'
#' stopifnot(x[[1]] - 1990.9226 < 1e-4,
#'           x[[2]] - 1991.9066 < 1e-4,
#'           x[[3]] - 2047.9440 < 1e-4,
#'           x[[4]] - 2048.9281 < 1e-4)
#'
#' }
ms1_a1_vnl0_fnl0_byprot <- function (prot_peps, amods, aa_masses,
                                     vmods_nl = NULL, fmods_nl = NULL,
                                     include_insource_nl = FALSE,
                                     maxn_vmods_per_pep = 5L,
                                     maxn_sites_per_vmod = 3L,
                                     digits = 4L) {
  
  # USE.NAMES of prot_acc
  out <- mapply(ms1_a1_vnl0_fnl0_bypep, prot_peps, names(prot_peps), 
                MoreArgs = list (
                  amods = amods, 
                  aa_masses = aa_masses,
                  vmods_nl = vmods_nl, 
                  fmods_nl = fmods_nl,
                  include_insource_nl = include_insource_nl,
                  maxn_vmods_per_pep = maxn_vmods_per_pep,
                  maxn_sites_per_vmod = maxn_sites_per_vmod,
                  digits = digits
                ), SIMPLIFY = FALSE)
  
  nms <- mapply(function (x, y) rep(y, length(x)), out, names(out), 
                SIMPLIFY = FALSE, USE.NAMES = FALSE) %>%
    unlist(recursive = FALSE, use.names = FALSE)
  
  out <- out %>%
    unlist(recursive = FALSE, use.names = FALSE) %>%
    `names<-`(nms)
}


#' Helper by individual peptides.
#'
#' (7) "amods+ tmod- vnl- fnl-"; (8) "amods+ tmod-+ vnl- fnl-".
#'
#' @param mass The mass of a peptide.
#' @inheritParams ms1_a1_vnl0_fnl0_byprot
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
ms1_a1_vnl0_fnl0_bypep <- function (mass, aa_seq, amods, aa_masses,
                                    vmods_nl = NULL, fmods_nl = NULL,
                                    include_insource_nl = FALSE,
                                    maxn_vmods_per_pep = 5L,
                                    maxn_sites_per_vmod = 3L,
                                    digits = 4L) {

  aas <- stringr::str_split(aa_seq, "", simplify = TRUE)

  vmods_combi <- unique_mvmods(amods = amods, ntmod = NULL, ctmod = NULL,
                               aa_masses = aa_masses, aas = aas,
                               maxn_vmods_per_pep = maxn_vmods_per_pep,
                               maxn_sites_per_vmod = maxn_sites_per_vmod,
                               digits = digits) %>%
    find_intercombi()

  deltas <- lapply(vmods_combi, function (x) sum(aa_masses[x]))

  out <- simplify2array(lapply(deltas, function (x) round(mass + x, digits = digits)))

  if (include_insource_nl) {
    if (length(vmods_nl)) {
      vnl_combi <- lapply(vmods_combi, function (x) expand.grid(vmods_nl[x]))
      deltas_vnls <- lapply(vnl_combi, function (x) unique(rowSums(x)))
      out <- mapply(`-`, out, deltas_vnls, SIMPLIFY = FALSE)
    }

    if (length(fmods_nl)) {
      fnl_combi <- expand.grid(fmods_nl)
      deltas_fnls <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
      out <- mapply(`-`, out, deltas_fnls, SIMPLIFY = FALSE)
    }
  }

  invisible(out)
}

