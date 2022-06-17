#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
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
#' @param maxn_sites_per_vmod Integer; the maximum number of combinatorial
#'   variable modifications per site in a per peptide sequence.
#' @param maxn_vmods_per_pep The maximum number of variable modifications per
#'   peptide.
#' @param .path_ms1masses The file path to the theoretical masses of MS1
#'   precursors.
#' @param digits Integer; the number of decimal places to be used.
#' @inheritParams calc_aamasses
#' @inheritParams matchMS
#' @import parallel
#' @examples
#' \donttest{
#' res <- calc_pepmasses2()
#'
#' library(purrr)
#' library(magrittr)
#'
#' res_attrs <- lapply(res, attributes)
#' lapply(res_attrs, names)
#' lapply(res_attrs, `[[`, "vmods")
#' res_mods <- lapply(res_attrs, `[`,
#'                    c("fmods", "fmods_ps", "fmods_neuloss",
#'                      "vmods", "vmods_ps", "vmods_neuloss"))
#'
#' res_data <- lapply(res_attrs, `[[`, "data")
#' peps_combi_1 <- res_data[[1]]
#'
#' # base: fixedmods without neulosses
#' length(unlist(res_data[[1]], use.names = FALSE))
#'
#' # fixedmods, fixedmods + fixedmods_neulosses, varmods, varmods_neulosses
#' length(unlist(res_data, use.names = FALSE))
#'
#' }
calc_pepmasses2 <- function (
  fasta = "~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
  acc_type = "uniprot_acc",
  acc_pattern = NULL,
  fixedmods = c("TMT6plex (N-term)", 
                "TMT6plex (K)", 
                "Carbamidomethyl (C)"),
  varmods = c("Acetyl (Protein N-term)", 
              "Oxidation (M)", 
              "Deamidated (N)",
              "Gln->pyro-Glu (N-term = Q)"),
  mod_motifs = NULL, 
  enzyme = c("trypsin_p"),
  custom_enzyme = c(Cterm = NULL, Nterm = NULL), 
  maxn_fasta_seqs = 50000L,
  maxn_vmods_setscombi = 64L,
  maxn_vmods_per_pep = 5L,
  maxn_sites_per_vmod = 3L,
  min_len = 7L, max_len = 40L, max_miss = 2L,
  min_mass = 700L, max_mass = 4500L, 
  n_13c = 0L,
  out_path = NULL,
  digits = 4L,
  use_ms1_cache = TRUE, 
  .path_cache = NULL, 
  .path_fasta = NULL, 
  .path_ms1masses = NULL) 
{
  ## Enzymatic and Semi-enzymatic
  # (1) split_fastaseqs: 
  #     splits FASTA sequences by full-enzyme specificity
  # (2) ms1masses_bare: 
  # (2.1) ms1masses_noterm: 
  #     calculates the bare masses of (1) without terminal masses (e.g. H2O) 
  # (2.2) roll_sum
  #     concatenates sequences and masses according to the value of `max_miss`
  # (2.3) Adds FIXED terminal mass: 
  #     H2O 18.010565, FIXED N-term TMT modification 230.170757 + 17.002740 ...
  #     (which are reflected on aa_masses["N-term"] and aa_masses["C-term"])
  # 
  # Note:
  # For "historical" reasons, terminal "tmods+" refers to "VARIABLE tmods+". 
  # Thus, if N-term modification is FIXED, all the types of `aa_masses` must be 
  #   "tmods-". The same is true for anywhere "amods+/-". 
  # 
  # (3) distri_peps:
  #     distributes peptides by modifications;
  #     removes terminal tags of "-"
  # (4) hsemipeps_byprots:
  #     semi-enzyme only
  # (5) tbl_prots_peps:
  #     records/sets aside the peptide-and-protein associations
  # (6) flat_pepseqs:
  #     removes the information of protein
  # (7) Twelve types of VARIABLE modifications/masses:
  #     adds variable terminal masses (tmod+), variable anywhere masses (vmods+) 
  #     neutral losses (fnl+) etc.
  
  
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)

  on.exit(
    if (exists(".savecall", envir = fun_env)) {
      if (.savecall) {
        save_call2(path = .path_cache, fun = fun, time = .time_stamp)
      }
    },
    add = TRUE
  )

  # ---
  fun <- as.character(match.call()[[1]])
  fun_env <- environment()
  
  # argument_name-default_value pair
  new_args <- local({
    args <- c("mod_motifs")
    fmls <- formals(fun)[args]
    
    # for example `custom_enzyme` (should be NULL after eval)
    nargs <- lapply(fmls, function (x) if (is.call(x)) eval(x) else x)
    
    # unlisting with NUlls being kept
    are_nulls <- lapply(nargs, is.null)
    def_nulls <- unlist(are_nulls)
    c(unlist(nargs), nargs[def_nulls])
  })
  
  .time_stamp <- match_calltime(
    path = .path_cache,
    fun = fun,
    
    # `nms` must be matched in order to retrieve cached results
    nms = c("fasta", "acc_type", "acc_pattern",
            "fixedmods", "varmods", "mod_motifs", 
            "enzyme", "custom_enzyme",
            "maxn_fasta_seqs", "maxn_vmods_setscombi",
            "maxn_vmods_per_pep", "maxn_sites_per_vmod",
            "min_len", "max_len", "max_miss", 
            "min_mass", "max_mass", "n_13c"), 
    
    # exception: new arguments need matches but not defined in earlier versions
    new_args = new_args)

  # ---
  len_ts <- length(.time_stamp)
  
  if (len_ts && use_ms1_cache) {
    # can have multiple matches with use_ms1_cache on/off
    .time_stamp <- .time_stamp[len_ts]
    
    message("Loading peptide masses from cache.")

    aa_masses_all <- find_aa_masses(
      out_path = file.path(.path_fasta, "ms1masses", .time_stamp),
      fixedmods = fixedmods,
      varmods = varmods,
      mod_motifs = mod_motifs, 
      maxn_vmods_setscombi = maxn_vmods_setscombi)

    files <- list.files(path = file.path(.path_ms1masses, .time_stamp),
                        pattern = "pepmasses_\\d+\\.rds$")

    if (length(files) != length(aa_masses_all)) 
      stop("Not all precursor masses were found: ", 
           paste0("\n", files), ".\n",
           "Remove cache file: \n", 
           file.path(.path_cache, fun, paste0(.time_stamp, ".rda")),
           " and try again.")

    rm(list = c("aa_masses_all", "files"))
    gc()

    fwd_peps <- NULL
    rev_peps <- NULL

    .savecall <- FALSE
  } 
  else {
    # `mgf_quries.rds` kept 
    # (which only affected by min_mass, max_mass and ppm_ms1)
    delete_files(out_path, 
                 ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$",
                             "\\.xls$", "\\.csv$", "\\.txt$",
                             "^mgf$", "^mgfs$"))
    
    .time_stamp <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")

    aa_masses_all <- find_aa_masses(
      out_path = file.path(.path_fasta, "ms1masses", .time_stamp),
      fixedmods = fixedmods,
      varmods = varmods,
      mod_motifs = mod_motifs, 
      maxn_vmods_setscombi = maxn_vmods_setscombi)
    
    qs::qsave(aa_masses_all, file.path(out_path, "aa_masses_all.rds"), preset = "fast")

    ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
                           maxn_vmods_per_pep = maxn_vmods_per_pep,
                           maxn_sites_per_vmod = maxn_sites_per_vmod)

    len <- length(aa_masses_all)
    types <- purrr::map_chr(aa_masses_all, attr, "type", exact = TRUE)

    # By design, variable modifications, including variable [NC]-term, 
    # are appended in parallel to unmodified residues in `aa_masses`, e.g., 
    # "M", "Oxidation (M)" are two separate entries in `aa_masses`
    # 
    # (1) `aa_masses_1` is for base-mass calculations (without variable masses):
    #   aas <- strsplit("PEPTIDE")); mass <- aa_masses_1[aas]
    # (2) "N-term" and "C-term" masses are identical across `aa_masses_all`,
    #     e.g. 18 or 247 etc.
    # Thus, `aa_masses_1` can be any `aa_masses_all[[i]]` 
    # but uses `1` for tidiness
    
    aa_masses_1 <- aa_masses_all[[1]]
    gc()

    # --- Forward sequences  ---
    if (isTRUE(enzyme == "noenzyme")) {
      if (max_len > 25L) 
        warning("May be out of RAM at `max_len = ", max_len, "`.\n",
                "Consider a smaller value, e.g., `max_len = 25`")

      is_fixed_protnt <- any(grepl("Protein N-term", fixedmods))
      is_fixed_protct <- any(grepl("Protein C-term", fixedmods))

      if (any(is_fixed_protnt, is_fixed_protct))
        stop("Not yet support FIXED protein terminal modifications for ", 
             "noenzyme searches.\n", 
             "Change to VARIABLE protein terminal modifications.")

      seqs_0 <- NULL
      
      ftmass <- unname(aa_masses_1["N-term"] + aa_masses_1["C-term"])
      
      fwd_peps <- split_fastaseqs_noenz(fasta = fasta, 
                                        acc_type = acc_type,
                                        acc_pattern = acc_pattern,
                                        maxn_fasta_seqs = maxn_fasta_seqs, 
                                        min_len = min_len, 
                                        max_len = max_len, 
                                        aa_masses = aa_masses_1, 
                                        ftmass = ftmass, 
                                        digits = digits)
      
      gc()
    }
    else {
      # (not yet concatenation by the number of missed cleavages)
      seqs_0 <- split_fastaseqs(fasta = fasta,
                                enzyme = enzyme, 
                                custom_enzyme = custom_enzyme, 
                                acc_type = acc_type,
                                acc_pattern = acc_pattern,
                                maxn_fasta_seqs = maxn_fasta_seqs,
                                max_miss = max_miss)
      
      ### Special case of FIXED Protein [NC]-term modification(s) ---
      # 
      # (i) Dispatches `pep_seq`s by `fixedmod`s
      # (only incur at the rare case of Protein terminals being `fixedmods`)
      # (`fixed Protein [NC]-term` -> no `variable Protein [NC]-term`)
      # 
      # (ii) Calculates terminal mass (after `distri_fpeps`)
      # (otherwise, fixed "Protein N|C-terminal" masses in aa_masses["N|C-term"]
      #   will be applied to both "N|C-term" and "Protein N|C term")
      #
      # if `ftmass` is other than 18.010565 -> FIXED [NC]-term 
      #   -> NO variable "Protein N-term" etc.
      # (ftmass will be 18.010565 or plus fixed [NC] terminal modifications)
      #   aa_masses["N-term"] = 1.007825, aa_masses["C-term"] = 17.002740
      
      is_fixed_protnt <- any(grepl("Protein N-term", fixedmods))
      is_fixed_protct <- any(grepl("Protein C-term", fixedmods))
      seqs_0 <- distri_fpeps(data = seqs_0, max_miss = max_miss, 
                             is_fixed_protnt = is_fixed_protnt, 
                             is_fixed_protct = is_fixed_protct)

      # --- Masses of sequences: fixed mods + terminals ---
      message("Calculating bare peptide masses...")
      
      # e.g. if "TMT6plex (N-term)" is a fixedmod -> ftmass = (229 + 1) + 17;
      # if "TMT6plex (N-term)" coerced to a varmod -> ftmass = 1 + 17
      # see also calc_aamasses for details
      ftmass <- unname(aa_masses_1["N-term"] + aa_masses_1["C-term"])
      
      fwd_peps <- ms1masses_bare(seqs = seqs_0,
                                 aa_masses = aa_masses_1,
                                 ftmass = ftmass,
                                 max_miss = max_miss,
                                 min_len = min_len,
                                 max_len = max_len,
                                 min_mass = min_mass, 
                                 max_mass = max_mass, 
                                 maxn_vmods_per_pep = maxn_vmods_per_pep,
                                 maxn_sites_per_vmod = maxn_sites_per_vmod,
                                 is_fixed_protnt = is_fixed_protnt,
                                 is_fixed_protct = is_fixed_protct,
                                 digits = digits)

      message("\tCompleted bare peptide masses.")
    }

    
    # --- Semi-enzymatic and distribution ---
    # Note 1-to-n expansion: 
    #   `length(fwd_peps) == length(aa_masses_all)` after the step.
    n_cores <- detect_cores(16L)
    
    if (isTRUE(enzyme == "noenzyme"))
      n_cores <- floor(max(1L, n_cores/4L))

    fwd_peps <- chunksplit(fwd_peps, n_cores, "list")
    
    # (a) Optional semi-enzymatic peptides
    if (grepl("^semi", enzyme)) {
      message("Generating semi-enzymatic peptides.")
      
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      parallel::clusterExport(
        cl,
        c("hsemipeps_byprots", 
          "semipeps_byprots", 
          "calc_semipepmasses"), 
        envir = environment(proteoM:::calc_semipepmasses)
      )
      
      fwd_peps <- parallel::clusterApply(
        cl, 
        fwd_peps, 
        hsemipeps_byprots, 
        min_len = min_len , 
        aa_masses = aa_masses_1
      )
      
      parallel::stopCluster(cl)
    }
    
    # (b) Distribution
    message("Distributing peptides by variable modifications.")
    
    motifs_all <- lapply(aa_masses_all, find_motif_pat)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(
      cl,
      c("distri_peps", 
        "ct_counts", 
        "rm_char_in_nfirst", 
        "rm_char_in_nlast"), 
      envir = environment(proteoM:::distri_peps)
    )

    fwd_peps <- parallel::clusterApply(
      cl, 
      fwd_peps, 
      distri_peps, 
      aa_masses_all = aa_masses_all, 
      motifs_all = motifs_all, 
      max_miss = max_miss, 
      max_len = max_len, 
      enzyme = enzyme
    )
    
    parallel::stopCluster(cl)
    
    fwd_peps <- lapply(seq_len(len), function (i) {
      # by i-th aa_masses from each node
      fwd_peps_i <- lapply(fwd_peps, `[[`, i)
      
      # combines i-th results across nodes
      purrr::flatten(fwd_peps_i)
    })
    
    message("\tCompleted bare peptides distributions.")
    
    rm(list = c("seqs_0"))
    gc()

    # (c) Protein-peptide associations
    tbl_prots_peps(fwd_peps[[1]], file.path(.path_ms1masses, .time_stamp))
    gc()
    
    # (d) Flattened peptide lists (ripping off prot_acc's)
    fwd_peps <- lapply(fwd_peps, flat_pepseqs)
    gc()


    # --- Delta masses of `variable` terminals  ---
    # (e.g., on top of the `fixed` 18.010565)
    message("Adding terminal masses (variable modifications) ...")

    inds <- grep("tmod+", types, fixed = TRUE)

    if (length(inds)) {
      nt_inds <- which(types %in% c("amods- tmod+ vnl- fnl-",
                                    "amods- tmod+ vnl- fnl+"))

      for (i in inds) {
        fwd_peps[[i]] <- add_term_mass2(aa_masses_all[[i]], fwd_peps[[i]], 
                                        max_mass)

        if (i %in% nt_inds) 
          message("\tCompleted peptide terminal masses: ",
                  paste(attributes(aa_masses_all[[i]])$fmods,
                        attributes(aa_masses_all[[i]])$vmods,
                        collapse = ", "))

        gc()
      }

      rm(list = "nt_inds")
    }

    
    # --- Mass of variable mods and/or NLs ---
    message("Calculating peptide masses (variable modifications + neutral losses) ...")

    fmods_ps <- lapply(aa_masses_all, attr, "fmods_ps", exact = TRUE)
    vmods_ps <- lapply(aa_masses_all, attr, "vmods_ps", exact = TRUE)
    fmods_nl <- lapply(aa_masses_all, attr, "fmods_nl", exact = TRUE)
    vmods_nl <- lapply(aa_masses_all, attr, "vmods_nl", exact = TRUE)
    amods <- lapply(aa_masses_all, attr, "amods", exact = TRUE)
    tmod <- lapply(aa_masses_all, attr, "tmod", exact = TRUE)

    # `amods-` and `fnl+` (must be vnl- since amods-)
    #
    # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
    
    # in-source fragmentation
    if (FALSE) {
      n_cores <- detect_cores(32L)

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

          aa_masses_i <- aa_masses_all[[i]]
          ntmod_i <- attr(aa_masses_i, "ntmod", exact = TRUE)
          ctmod_i <- attr(aa_masses_i, "ctmod", exact = TRUE)

          fwd_peps_i <- fwd_peps[[i]]
          
          if (!length(fwd_peps_i))
            next
          
          fnl_combi_i <- expand_grid_rows(fmods_nl_i, use.names = TRUE)
          
          cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
          
          parallel::clusterExport(
            cl,
            c("hms1_a0_vnl0_fnl1", 
              "ms1_a0_vnl0_fnl1", 
              "expand_grid_rows", 
              "delta_ms1_a0_fnl1"), 
            envir = environment(proteoM:::ms1_a0_vnl0_fnl1))
          
          fwd_peps[[i]] <- parallel::clusterApply(
            cl, 
            chunksplit(fwd_peps_i, n_cores, "list"), 
            hms1_a0_vnl0_fnl1, 
            fnl_combi = fnl_combi_i, 
            aa_masses = aa_masses_i,
            max_mass = max_mass, 
            digits = digits
          ) %>% 
            purrr::flatten() %>% 
            unlist(recursive = FALSE, use.names = TRUE)

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
    
    # i = 4
    # 16L: 1.45542643 mins
    # 32L: 1.0508078 mins
    # 64L: 1.06206823 mins
    # 96L: 1.25min
    # 128L: 1.44mins
    
    # i = 7L
    # 16L: 16.572324 secs
    # 32L: 22.6
    # 64L: 39 

    n_cores <- detect_cores(32L)

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
        aa_masses_i <- aa_masses_all[[i]]
        ms1vmods_i <- ms1vmods_all[[i]]

        fwd_peps_i <- fwd_peps[[i]]
        
        if (!length(fwd_peps_i))
          next

        vmods_nl_i = vmods_nl[[i]]
        fmods_nl_i = fmods_nl[[i]]
        
        cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
        
        parallel::clusterExport(
          cl,
          c("hms1_a1_vnl0_fnl0", 
            "ms1_a1_vnl0_fnl0", 
            "match_mvmods", 
            "expand_grid_rows", 
            "recur_flatten", 
            "delta_ms1_a0_fnl1"), 
          envir = environment(proteoM:::ms1_a1_vnl0_fnl0))
        
        fwd_peps[[i]] <- parallel::clusterApply(
          cl, 
          chunksplit(fwd_peps_i, n_cores, "list"), 
          hms1_a1_vnl0_fnl0, 
          amods = amods_i, 
          aa_masses = aa_masses_i,
          vmods_nl = vmods_nl_i, 
          fmods_nl = fmods_nl_i,
          maxn_vmods_per_pep = maxn_vmods_per_pep,
          maxn_sites_per_vmod = maxn_sites_per_vmod,
          ms1vmods = ms1vmods_i, 
          max_mass = max_mass, 
          digits = digits
        ) %>% 
          purrr::flatten() %>% 
          unlist(recursive = FALSE, use.names = TRUE)

        parallel::stopCluster(cl)
        gc()

        message("\tCompleted peptide masses: ",
                paste(attributes(aa_masses_i)$fmods,
                      attributes(aa_masses_i)$vmods,
                      collapse = ", "))
      }
    }

    suppressWarnings(
      rm(list = c("amods_i", "fmods_nl", "fmods_ps", "fwd_peps_i", 
                  "vmods_nl_i", "aa_masses_1", "aa_masses_i"))
    )
    
    gc()

    fwd_peps <- lapply(fwd_peps, add_ms1_13c, n_13c)
    
    # lapply(fwd_peps, max)

    # === Outputs ===
    path_masses <- create_dir(file.path(.path_ms1masses, .time_stamp))

    fwd_peps <- purrr::map2(aa_masses_all, fwd_peps, ~ {
      attr(.x, "data") <- .y
      .x
    })

    names(fwd_peps) <- seq_along(aa_masses_all)
    gc()

    for (i in seq_along(fwd_peps)) {
      qs::qsave(fwd_peps[[i]], 
                file.path(path_masses, paste0("pepmasses_", i, ".rds")), 
                preset = "fast")
    }
    gc()

    # ---
    .savecall <- TRUE
    message("\n=== Completed MS1 precursor masses. ===\n")
  }

  assign(".path_cache", .path_cache, envir = .GlobalEnv)
  assign(".path_fasta", .path_fasta, envir = .GlobalEnv)
  assign(".path_ms1masses", .path_ms1masses, envir = .GlobalEnv)
  assign(".time_stamp", .time_stamp, envir = .GlobalEnv)
  
  local({
    .cache_info <- list(
      .path_cache = .path_cache, 
      .path_fasta = .path_fasta, 
      .path_ms1masses = .path_ms1masses, 
      .time_stamp = .time_stamp
    )
    
    path <- create_dir(file.path(out_path, "Calls"))
    qs::qsave(.cache_info, file.path(path, ".cache_info.rds"), preset = "fast")
  })

  invisible(fwd_peps)
}


#' Finds the existence of \code{aa_masses_all.rds}.
#'
#' @inheritParams calc_pepmasses2
find_aa_masses <- function(out_path = NULL, fixedmods = NULL, varmods = NULL,
                           mod_motifs = NULL, maxn_vmods_setscombi = 64L) 
{
  if (!file.exists(file.path(out_path, "aa_masses_all.rds"))) {
    message("Computing the combinations of fixed and variable modifications.")

    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                                   varmods = varmods,
                                   mod_motifs = mod_motifs, 
                                   maxn_vmods_setscombi = maxn_vmods_setscombi,
                                   out_path = out_path) 
    
    qs::qsave(aa_masses_all, file.path(out_path, "aa_masses_all.rds"), preset = "fast")
  } else {
    aa_masses_all <- qs::qread(file.path(out_path, "aa_masses_all.rds"))
  }

  invisible(aa_masses_all)
}


#' Finds the motifs of modifications.
#' 
#' @param aa_masses A list of amino acid lookups.
find_motif_pat <- function (aa_masses)
{
  mms <- attr(aa_masses, "mod_motifs", exact = TRUE)
  mms <- .Internal(unlist(mms, recursive = FALSE, use.names = FALSE))
  
  if (length(mms)) {
    motifs <- unique(mms)
    .Internal(paste0(list(motifs), collapse = "|", recycle0 = FALSE))
  }
  else {
    NULL
  }
}


#' A lookup table between prot_acc and pep_seqs.
#'
#' For both target and decoy peptides.
#'
#' Columns: \code{is_pnt}, is Protein N-term; \code{is_pct}, is Protein C-term.
#' For example, MEYEWKPDEQGLQQILQLLK: NP_002261, pep_start 9; NP_694858:
#' pep_start 1. Without the information of protein terminals, the peptide
#' sequence with the modification of Acetyl Protein N-term will match to both
#' protein accessions, with the match to NP_002261 being inconsistent.
#'
#' @param seqs Results from \link{distri_peps}.
#' @param path A file path.
#' @examples
#' \donttest{
#' tbl_prots_peps(fwd_peps[[1]])
#' }
tbl_prots_peps <- function (seqs, path) 
{
  message("Tabling the association of proteins and peptides.")
  
  path <- create_dir(path)
  
  pnt_idxes <- lapply(seqs, attr, "pnt_idxes", exact = TRUE)
  pct_idxes <- lapply(seqs, attr, "pct_idxes", exact = TRUE)
  
  seqs <- lapply(seqs, names)
  lens <- lapply(seqs, length)
  
  # Protein [NC]-term
  len <- length(seqs)
  cts <- nts <- vector("list", len)
  
  for (i in seq_len(len)) {
    seq_i <- seqs[[i]]
    pnt_i <- pnt_idxes[[i]]
    pct_i <- pct_idxes[[i]]
    ci <- ni <- rep(FALSE, lens[[i]])
    ni[pnt_i] <- TRUE
    ci[pct_i] <- TRUE
    
    nts[[i]] <- ni
    cts[[i]] <- ci
  }
  
  rm(list = c("seq_i", "pnt_i", "pct_i", "ni", "ci"))
  
  nts <- .Internal(unlist(nts, recursive = FALSE, use.names = FALSE))
  cts <- .Internal(unlist(cts, recursive = FALSE, use.names = FALSE))
  
  # Duplicated sequences within proteins
  dups <- lapply(seqs, duplicated.default)
  dups <- .Internal(unlist(dups, recursive = FALSE, use.names = FALSE))

  # prot_accs
  nms <- mapply(function (x, y) rep(x, y), 
                names(seqs), lens,
                SIMPLIFY = FALSE, USE.NAMES = FALSE) 
  
  nms <- .Internal(unlist(nms, recursive = FALSE, use.names = FALSE))
  
  # pep_seqs
  seqs <- unlist(seqs, recursive = FALSE, use.names = FALSE)
  
  ans <- data.frame(pep_seq = seqs, 
                    prot_acc = nms, 
                    is_pnt = nts, 
                    is_pct = cts, 
                    is_unique = !dups)
  
  # removals of duplicated sequences
  ans <- ans[with(ans, is_unique), ]
  ans$is_unique <- NULL
  
  ## Should be all clean after `distri_peps`
  # if (enzyme == "noenzyme") ans$pep_seq <- with(ans, gsub("-", "", pep_seq))
  
  qs::qsave(ans, file.path(path, "prot_pep_annots.rds"), 
            preset = "fast")
  qs::qsave(reverse_peps_in_frame(ans), file.path(path, "prot_pep_annots_rev.rds"), 
            preset = "fast")
  
  invisible(NULL)
}


#' Flattens pep_seqs with the removals of prot_accs.
#' 
#' @param x Lists of pep_seqs by prot_accs.
flat_pepseqs <- function (x) 
{
  x <- purrr::flatten(x)
  x <- .Internal(unlist(x, recursive = FALSE, use.names = TRUE))
  
  x[!duplicated.default(names(x))]
}


#' Finds the site of an AA residue.
#'
#' e.g. 'N-term' in `aa_masses`, not 'Q', for Gln-> pyro-Glu (N-term = Q).
#'
#' @param pos_site A named value. Position in name and site in value.
find_aa_site <- function (pos_site) 
{
  nms <- names(pos_site)
  
  site <- if (grepl("[NC]{1}-term", nms)) 
    gsub("(Protein|Any) ([NC]{1}-term)", "\\2", nms)
  else 
    pos_site
}


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param maxn_vmods_setscombi Integer; the maximum number of combinatorial variable
#'   modifications and neutral losses.
#' @param out_path An output path.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' x <- calc_aamasses()
#' x_att <- lapply(x, attributes)
#' names(x_att[[1]])
#' x_vmods <- lapply(x_att, `[`, c("vmods"))
#' x_fmods <- lapply(x_att, `[`, c("fmods"))
#' 
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (C)"), c("Acetyl (N-term)",
#'                      "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))
#' 
#' stopifnot(length(x) == 6L)
#'
#' # Fixed N-term mod (no coercion to variable mod)
#' x <- calc_aamasses("TMT6plex (N-term)", NULL)
#' x[[1]][["N-term"]]
#' 
#' # Fixed N-term mod (coerced to variable mod)
#' x <- calc_aamasses("TMT6plex (N-term)", "Acetyl (Protein N-term)")
#' lapply(x, `[[`, "N-term")
#' x[[1]][["TMT6plex (N-term)"]]
#' x[[2]][["Acetyl (Protein N-term)"]]
#' 
#' # No fixed mod
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
                           mod_motifs = NULL, 
                           maxn_vmods_setscombi = 64L,
                           out_path = NULL) 
{
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

  ## (0) Duplicated mods
  local({
    dup_mods <- intersect(fixedmods, varmods)

    if (length(dup_mods)) 
      stop("Modifications cannot be simultaneously 'fixed' and 'variable': \n",
           paste(dup_mods, collapse = ", "), "\n", 
           "Hint: the default \"fixedmods\" and \"varmods\" are not NULL.", 
           call. = FALSE)
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

    if (length(dup_fixedmods)) 
      stop("Multiple fixed modifications to the same site: \n",
           "'", paste(dup_fixedmods, collapse = ", "), "'",
           call. = FALSE)

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
        idxes <- lapply(dup_mods, function (x) fmods_ps == x)
        unlist(lapply(idxes, function (i) fixedmods[i]))
      })

      varmods <- c(f_to_v, varmods)

      fixedmods <- local({
        idxes <- unlist(lapply(fixedmods, function (x) x %in% f_to_v))
        fixedmods[!idxes]
      })

      warning("Coerce '",
              paste(f_to_v, collapse = ", "), "'",
              " to conditional variable modifications.",
              call. = FALSE)
      
      fixednt_coerced <- any(grepl("N-term", f_to_v))
      fixedct_coerced <- any(grepl("C-term", f_to_v))
      anywhere_coreced_sites <- dup_mods[!grepl("[NC]-term", dup_mods)]
    } 
    else {
      f_to_v <- NULL
      fixednt_coerced <- FALSE
      fixedct_coerced <- FALSE
      anywhere_coreced_sites <- NULL
    }

    default_mods <- c("initiator methionine from protein N-terminus")

    if (any(grepl(default_mods, c(fixedmods, varmods)))) 
      warning("Modifications defaulted and no need to specify: `\n",
              default_mods, "`.\n",
              call. = FALSE)

    invisible(list(fixedmods = fixedmods, 
                   varmods = varmods, 
                   f_to_v = f_to_v, 
                   fixednt_coerced = fixednt_coerced, 
                   fixedct_coerced = fixedct_coerced, 
                   anywhere_coreced_sites = anywhere_coreced_sites, 
                   fV_coercion = fV_coercion))
  })

  fixedmods <- new_mods$fixedmods
  varmods <- new_mods$varmods
  f_to_v <- new_mods$f_to_v
  fixednt_coerced <- new_mods$fixednt_coerced
  fixedct_coerced <- new_mods$fixedct_coerced
  anywhere_coreced_sites <- new_mods$anywhere_coreced_sites
  fV_coercion <- new_mods$fV_coercion
  rm(list = "new_mods")
  
  if (!is.null(mod_motifs)) {
    local({
      bads <- mod_motifs[!names(mod_motifs) %in% c(varmods, fixedmods)]
      
      if (length(bads))
        stop("\"mod_motifs\" not found in \"varmods\" or \"fixedmods\": ", 
             paste(bads, collapse = ", "))
    })
  }

  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943,
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464,
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414,
    "N-term" = 1.007825, "C-term" = 17.002740,
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590,
    "-" = 0)
  
  # At the end, the mass delta from fixed [NC]-term will be added to 
  # each aa_masses["N-term"] and aa_masses["C-term"]
  # e.g. with a TMT N-term tag, aa_masses["N-term"] + aa_masses["C-term"] = 
  # TMT + H2O
  # 
  # Variable [NC]-term will be a parallel entry in an aa_masses
  # e.g., aa_masses["Acetyl (N-term)"]
  # 
  # In other words, masses of fixed modification will be reflected in residues,
  # including "N-term" and "C-term"; 
  # Masses of variable modifications will be additive to the base vector of list.
  # E.g., if "Oxidated (M)" is a fixed modification -> M = 147; otherwise if
  # it is a variable modification, M = 131 and there is a second entry of 
  # "Oxidated (M)" = 147.

  ## (1) add fixed mods + NL
  
  fmod_motifs <- mod_motifs[names(mod_motifs) %in% fixedmods]
  vmod_motifs <- mod_motifs[names(mod_motifs) %in% varmods]

  aa_masses_fi2 <- add_fixvar_masses(mods = fixedmods,
                                     mod_type = "fmods",
                                     aa_masses = aa_masses, 
                                     mod_motifs = fmod_motifs)

  aa_masses_fi2 <- lapply(aa_masses_fi2, function (x) {
    if (is.null(attr(x, "fmods"))) 
      attr(x, "fmods") <- ""
    
    if (is.null(attr(x, "fmods_neuloss"))) 
      attr(x, "fmods_neuloss") <- ""
    
    if (is.null(attr(x, "fmods_mass"))) 
      attr(x, "fmods_mass") <- 0
    
    x
  })

  ## (2) add variable mods + NL
  varmods_comb <- local({
    # (c) Remove entries with multiple terminal mods
    varmods_comb <- lapply(seq_along(varmods), 
                           function (x) sim_combn(varmods, x)) 
    varmods_comb <- unlist(varmods_comb, recursive = FALSE)
    
    # Check if multiple terminal varmods by names
    #   Gln->pyro Glu (N-term = Q) and Acetyl (Protein N-term = N-term)
    #   have different sites but 'N-term' in names; and also need to be excluded

    vmods_ps <- varmods %>%
      lapply(find_unimod) %>%
      lapply(`[[`, "position_site") %>%
      `names<-`(varmods) %>% 
      purrr::flatten()

    vmods_ps_combi <- seq_along(vmods_ps) %>%
      lapply(function (x) sim_combn(vmods_ps, x)) %>%
      purrr::flatten()
    
    dup_terms <-
      purrr::map_lgl(vmods_ps_combi, ~ sum(grepl("N-term", names(.x))) >= 2L) |
      purrr::map_lgl(vmods_ps_combi, ~ sum(grepl("C-term", names(.x))) >= 2L)

    # Not currently used
    dup_anywhere <- 
      purrr::map_lgl(vmods_ps_combi, function (x) {
        x <- x[!grepl("[NC]{1}-term", names(x))]
        any(duplicated(x))
      })

    # Allow entries with different Anywhere mods to the same site
    #   (1) dHex(1)Hex(1) (S) and (2) Phospho (S)
    #   VS(1)S(2)ALSPSK

    varmods_comb  %>%
      .[!(dup_terms)] %>%
      lapply(unlist)
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
  aa_masses_var2 <- lapply(varmods_comb, function (vi) {
    lapply(aa_masses_fi2, function (x) 
      add_fixvar_masses(mods = vi,
                        mod_type = "vmods",
                        aa_masses = x, 
                        mod_motifs = vmod_motifs)) %>%
      purrr::flatten()
  }) %>%
    purrr::flatten()

  aa_masses_var2 <- lapply(aa_masses_var2, function (x) {
    if (is.null(attr(x, "vmods"))) 
      attr(x, "vmods") <- ""
    if (is.null(attr(x, "vmods_neuloss"))) 
      attr(x, "vmods_neuloss") <- ""
    if (is.null(attr(x, "vmods_mass"))) 
      attr(x, "vmods_mass") <- 0
    
    x
  })

  ## (3) complete 'vmods', 'vmods_neuloss' and 'vmods_ps' to fixedmods
  aa_masses_fi2 <- aa_masses_fi2 %>%
    lapply(function (x) {
      if (is.null(attr(x, "vmods"))) 
        attr(x, "vmods") <- ""
      if (is.null(attr(x, "vmods_neuloss"))) 
        attr(x, "vmods_neuloss") <- ""
      if (is.null(attr(x, "vmods_ps"))) 
        attr(x, "vmods_ps") <- ""
      if (is.null(attr(x, "vmods_mass"))) 
        attr(x, "vmods_mass") <- 0
      
      x
    })
  
  if (length(aa_masses_var2) >= maxn_vmods_setscombi) {
    warning("The ways of combinatorial variable modifications are ",
            length(aa_masses_var2), ".\n",
            "Dropping combinations at indexes greater than `maxn_vmods_setscombi = ", 
            maxn_vmods_setscombi, "`.",
            call. = FALSE)
    aa_masses_var2 <- aa_masses_var2[1:maxn_vmods_setscombi]
  }

  ## Coerced sites, including "[NC]-term", need to be present in a combination
  # 
  # with coercion, `aa_masses_fi2` becomes something not specified by users:
  #   TMT6plex (N-term) fixed -> variable; the `aa_masses_fi2` 
  #   c(""TMT6plex (K)", "Carbamidomethyl (C)") without `TMT6plex (N-term)` 
  #   is not a combination intended by users)

  aa_masses_all <- if (fV_coercion) 
    aa_masses_var2
  else 
    c(aa_masses_fi2, aa_masses_var2)

  if (length(anywhere_coreced_sites)) {
    vmods_ps <- lapply(aa_masses_all, attr, "vmods_ps")
    
    for (i in seq_along(anywhere_coreced_sites)) 
      aa_masses_all <- check_anywhere_fmods_coercion(anywhere_coreced_sites[i], 
                                                     aa_masses_all, vmods_ps)
    
    rm(list = c("vmods_ps"))
  }

  # Indexes of modifications
  if (!is.null(out_path)) {
    mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
      as.hexmode() %>%
      `names<-`(c(fixedmods, varmods))
    
    is_coerced <- if (length(f_to_v)) 
      names(mod_indexes) %in% f_to_v
    else 
      rep(FALSE, length(mod_indexes))

    Desc <- if (length(mod_indexes)) 
      names(mod_indexes)
    else 
      character()

    ## At NULL fixedmods and varmods: 
    # [1] Abbr    Desc    Type    Coerced
    # <0 rows> (or 0-length row.names)
    
    df_mods <- data.frame(Abbr = as.character(mod_indexes),
                          Desc = Desc, 
                          Type = c(rep("fixed", length(fixedmods)), 
                                   rep("variable", length(varmods))), 
                          Coerced = is_coerced) 
    
    readr::write_tsv(df_mods, file.path(out_path, "mod_indexes.txt"))
    
    rm(list = c("mod_indexes", "is_coerced", "Desc", "df_mods"))
  }

  aa_masses_all <- lapply(aa_masses_all, parse_aamasses)
  
  invisible(aa_masses_all)
}


#' Checks fixed modifications after coercion to variable modifications.
#'
#' The coerced site needs to be present in the final \code{aa_masses_all},
#' as either a fixed or a variable modification.
#'
#' @param site An amino acid site; e.g. site = "M".
#' @param aa_masses_all All the amino acid lookup tables.
#' @param vmods_ps All of the "vmods_ps" attributes from "aa_masses_all".
check_anywhere_fmods_coercion <- function (site, aa_masses_all, vmods_ps) 
{
  oks <- purrr::map_lgl(vmods_ps, ~ any(grepl(site, .x)))
  aa_masses_all <- aa_masses_all[oks]
  
  if (!length(aa_masses_all)) 
    stop("Zero combination of AA tables ", 
         "Check `fixedmods` and `varmods`.", 
         call. = FALSE)
  
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
#' @inheritParams matchMS
#' @return Lists of of amino-acid residues with modified mono-isotopic masses
#'   being incorporated.
add_fixvar_masses <- function (mods, mod_type, aa_masses, mod_motifs = NULL) 
{
  all_mods <- if (length(mods)) paste(mods, collapse = ", ") else ""

  res <- mods %>%
    lapply(find_unimod) %>%
    `names<-`(mods)
  
  local({
    if (length(res)) {
      x <- res[[1]]
      
      nm_seqs <- c("title", "monomass", "position_site", "nl")
      ok <- identical(names(x), nm_seqs)
      
      if (!ok)
        stop("The structures from `find_unimod` is not in the order of: ", 
             paste(nm_seqs, collapse = ", "))
    }
  })

  mod_masses <- lapply(res, `[[`, "monomass")
  positions_sites <- lapply(res, `[[`, "position_site")
  neulosses <- lapply(res, `[[`, "nl")
  rm(list = c("res"))

  # the same `site` with different fixedmods
  local({
    if (isTRUE(mod_type == "fmods") && (length(positions_sites) > 1L)) {
      dups <- purrr::reduce(positions_sites, `c`) %>%
        .[duplicated(.)]

      if (length(dups)) {
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
  } 
  else {
    #  Add mod_masses of variable mods (multiple lists)
    aas <- purrr::map2(positions_sites, mod_masses, ~ {
      site <- find_aa_site(.x)
      aa_masses[site] <- .y
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
  
  # Nee update subset_by_prps and subset_anysite if changed from NULL to ""
  nms <- names(positions_sites)
  pmod_motifs <- mod_motifs[nms]
  
  # pmod_motifs <- lapply(pmod_motifs, function (x) if (is.null(x)) "" else x)
  # names(pmod_motifs) <- nms
  if (!is.null(pmod_motifs)) names(pmod_motifs) <- nms
  
  attr(aa_masses, "mod_motifs") <- pmod_motifs
  rm(list = c("nms", "pmod_motifs"))
  
  ### 
  # Variable mods: need both the "original" and the "delta" forms
  # whereas fixed mods only has one form.
  # 
  # (1) 
  # fixedmods = c("Carbamidomethyl (C)"),
  # varmods = c("K8 (C-term)", "TMT6plex (N-term)"),
  # 
  # N-term      C-term      K8 (C-term)   TMT6plex (N-term) 
  # 1.007825    17.002740   8.014200      229.162932
  # 
  # (2)
  # fixedmods = c("Carbamidomethyl (C)"),
  # varmods = c("Oxidation (M)"), 
  # M           Oxidation (M) 
  # 131.040485  15.994915 
  ### 

  # add mod_masses of neutral losses
  no_nls <- neulosses %>%
    purrr::map_lgl(~ all(.x == 0)) %>%
    all()

  if (no_nls) 
    return(list(aa_masses))

  attr(aa_masses, paste0(mod_type, "_neuloss")) <- neulosses

  list(aa_masses)
}


#' Parses \code{aa_masses}.
#'
#' @inheritParams add_fixvar_masses
parse_aamasses <- function (aa_masses) 
{
  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)

  fmods_nl <- local({
    neulosses <- attr(aa_masses, "fmods_neuloss", exact = TRUE)

    if (all(neulosses == "")) 
      return(character())

    # add `0` if absent
    no_zero <- which(purrr::map_lgl(neulosses, ~ !any(.x == 0)))

    if (length(no_zero)) 
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])

    # Entries with NL = 0 also kept (beneficial when expand.grid).
    # In fixedmods: `M` instead of `Oxidation (M)`.

    names(neulosses) <- fmods_ps
    neulosses
  })

  vmods_nl <- local({
    neulosses <- attr(aa_masses, "vmods_neuloss", exact = TRUE)

    if (all(neulosses == "")) 
      return(character())

    no_zero <- which(purrr::map_lgl(neulosses, ~ !any(.x == 0)))

    if (length(no_zero)) 
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])

    neulosses
  })

  ## variable mods
  # multiple mods to [NC]-term already excluded from aa_masses

  amods <- local({
    sites <- purrr::map(vmods_ps, ~ .x[grepl("Anywhere", names(.x))])
    empties <- purrr::map_lgl(sites, purrr::is_empty)
    sites[!empties]
  })
  
  min_n_res <- amods %>% 
    unlist(recursive = FALSE, use.names = FALSE) %>% 
    count_elements()
  
  # Is "same Anywhere mod existed"
  is_same <- any(length(min_n_res) > 1L)

  # `TMT6plex (N-term)` and `Amidated (Protein C-term)`
  tmod <- vmods_ps %>% .[! . %in% amods]
  
  if (!length(tmod)) 
    tmod <- NULL
  else if (all(tmod == "")) 
    tmod <- NULL

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
  } 
  else if (length(tmod) == 2L) {
    # e.g., `TMT6plex (N-term)` + `Amidated (Protein C-term)`
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
  } 
  else {
    stop("cannot have more than two terminal modifications.",
         call. = FALSE)
  }

  ## fixed mods
  famods <- local({
    sites <- lapply(fmods_ps, function (x) x[grepl("Anywhere", names(x))])
    empties <- purrr::map_lgl(sites, purrr::is_empty)
    sites[!empties]
  })

  ftmod <- fmods_ps %>% .[! . %in% famods]
  
  if (!length(ftmod)) {
    ftmod <- NULL
  }
  # length(ftmod) > 1L: `TMT6plex (N-term)` and `K8 (C-term)`
  else if (all(ftmod == "")) {
    ftmod <- NULL
  }

  # fixed N-term, C-term
  fntmod <- ftmod %>% .[. == "N-term" || grepl("N-term", names(.))]
  fctmod <- ftmod %>% .[. == "C-term" || grepl("C-term", names(.))]

  # "amods- tmod- vnl- fnl-"
  if (!length(fmods_nl)) {
    type <- "fnl-"
    
    if (!length(vmods_nl)) {
      type <- paste("vnl-", type)
      
      if (!length(tmod)) {
        type <- paste("tmod-", type)
        
        type <- if (!length(amods)) 
          paste("amods-", type) # 1
        else 
          paste("amods+", type) # 2
      } else {
        type <- paste("tmod+", type)
        
        type <- if (!length(amods)) 
          paste("amods-", type) # 3
        else 
          paste("amods+", type) # 4
      }
    } else {
      type <- paste("vnl+", type)
      
      if (!length(tmod)) {
        type <- paste("tmod-", type)
        
        type <- if (!length(amods)) 
          paste("amods-", type) # 5
        else 
          paste("amods+", type) # 6
      } else {
        type <- paste("tmod+", type)

        type <- if (!length(amods)) 
          paste("amods-", type) # 7
        else 
          paste("amods+", type) # 8
      }
    }
  } else {
    type <- "fnl+"
    
    if (!length(vmods_nl)) {
      type <- paste("vnl-", type)
      
      if (!length(tmod)) {
        type <- paste("tmod-", type)
        
        type <- if (!length(amods)) 
          paste("amods-", type) # 1
        else 
          paste("amods+", type) # 2
      } else {
        type <- paste("tmod+", type)

        type <- if (!length(amods)) 
          paste("amods-", type) # 3
        else 
          paste("amods+", type) # 4
      }
    } else {
      type <- paste("vnl+", type)
      
      if (!length(tmod)) {
        type <- paste("tmod-", type)
        
        type <- if (!length(amods)) 
          paste("amods-", type) # 5
        else 
          paste("amods+", type) # 6
      } else {
        type <- paste("tmod+", type)

        type <- if (!length(amods)) 
          paste("amods-", type) # 7
        else 
          paste("amods+", type) # 8
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
  
  attr(aa_masses, "min_n_res") <- min_n_res
  attr(aa_masses, "is_same") <- is_same

  invisible(aa_masses)
}


#' Helper of \link{calc_pepmasses2}.
#'
#' Prior to the calculation of peptide masses; for the base with fixed
#' modification only.
#'
#' @inheritParams calc_pepmasses2
#' @importFrom stringi stri_reverse
#' @return Two named list of "fwds" and "revs". List "fwds" contains peptide
#'   sequences split from forward fasta and "revs" from reversed fasta.
#' @examples 
#' x <- split_fastaseqs("~/proteoM/dbs/fasta/uniprot/uniprot_mm_2020_11.fasta")
split_fastaseqs <- function (fasta = NULL, enzyme = "trypsin_p", 
                             custom_enzyme = c(Cterm = NULL, Nterm = NULL), 
                             acc_type = "uniprot_acc", acc_pattern = NULL, 
                             maxn_fasta_seqs = 200000L, max_miss = 2L) 
{
  message("Loading fasta databases.")

  fasta_db <- load_fasta2(fasta, acc_type, acc_pattern)

  if (length(fasta_db) > maxn_fasta_seqs) 
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.\n",
         "  May consider a higher `maxn_fasta_seqs`.")
  
  n_cores <- detect_cores(16L)

  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

  parallel::clusterExport(
    cl,
    c("make_fastapeps0", 
      "keep_n_misses"), 
    envir = environment(proteoM:::make_fastapeps0))

  # ---
  message("Splitting fasta sequences.")

  peps <- parallel::clusterApply(cl, chunksplit(fasta_db, n_cores), 
                                 make_fastapeps0, 
                                 enzyme, custom_enzyme, max_miss) %>%
    purrr::flatten()

  parallel::stopCluster(cl)
  
  rm(list = c("fasta_db"))
  gc()

  invisible(peps)
}


#' Makes peptide sequences from FASTA databases.
#'
#' A step before concatenating peptides by the number of mis-cleavages.
#'
#' The uses of "@" tags is faster than lookahead and lookbehind
#' (https://www.r-bloggers.com/2018/04/strsplit-but-keeping-the-delimiter/).
#'
#' @param fasta_db Fasta database(s).
#' @inheritParams calc_pepmasses2
make_fastapeps0 <- function (fasta_db, enzyme = "trypsin_p", custom_enzyme = NULL, 
                             max_miss = 2L) 
{
  inds_m <- grep("^M", fasta_db)

  if (is.null(enzyme)) {
    patc <- custom_enzyme[["Cterm"]]
    patn <- custom_enzyme[["Nterm"]]
    
    if (!is.null(patc)) {
      if (grepl("\\^", patc)) {
        fasta_db <- lapply(fasta_db, function (x) 
          .Internal(gsub(patc, paste0("\\1", "@", "\\2"), x, 
                         ignore.case = FALSE, perl = FALSE, 
                         fixed = FALSE, useBytes = FALSE)))
      }
      else {
        fasta_db <- lapply(fasta_db, function (x) 
          .Internal(gsub(patc, paste0("\\1", "@"), x, 
                         ignore.case = FALSE, perl = FALSE, 
                         fixed = FALSE, useBytes = FALSE)))
      }
    }
    
    if (!is.null(patn)) {
      if (grepl("\\^", patn)) {
        fasta_db <- lapply(fasta_db, function (x) 
          .Internal(gsub(patn, paste0("\\1", "@", "\\2"), x, 
                         ignore.case = FALSE, perl = FALSE, 
                         fixed = FALSE, useBytes = FALSE)))
      }
      else {
        fasta_db <- lapply(fasta_db, function (x) 
          .Internal(gsub(patn, paste0("@", "\\1"), x, 
                         ignore.case = FALSE, perl = FALSE, 
                         fixed = FALSE, useBytes = FALSE)))
      }
    }
    
    fasta_db <- lapply(fasta_db, function (x) gsub("@+", "@", x))
    
    fasta_db <- lapply(fasta_db, function (x) paste0("-", x, "-"))
  }
  else if (enzyme == "trypsin_p" || enzyme == "semitrypsin_p") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([KR]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "trypsin" || enzyme == "semitrypsin") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([KR]{1})([^P]{1})", paste0("\\1", "@", "\\2"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "lysc" || enzyme == "semilysc") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([K]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "lysn" || enzyme == "semilysn") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([K]{1})", paste0("@", "\\1"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "argc"|| enzyme == "semiargc") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([R]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "lysc_p"|| enzyme == "semilysc_p") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([K]{1})([^P]{1})", paste0("\\1", "@", "\\2"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "chymotrypsin"|| enzyme == "semichymotrypsin") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([FWY]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "gluc"|| enzyme == "semigluc") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([E]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "glun"|| enzyme == "semiglun") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([E]{1})", paste0("@", "\\1"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "aspc"|| enzyme == "semiaspc") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([D]{1})", paste0("\\1", "@"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "aspn"|| enzyme == "semiaspn") {
    fasta_db <- lapply(fasta_db, function (x) 
      paste0("-", 
             .Internal(gsub("([D]{1})", paste0("@", "\\1"), x, 
                            ignore.case = FALSE, perl = FALSE, 
                            fixed = FALSE, useBytes = FALSE)), 
             "-"))
  }
  else if (enzyme == "nodigest") {
    fasta_db <- lapply(fasta_db, function (x) paste0("-", x, "-"))
  }
  else {
    stop("Unknown enzyme.")
  }

  fasta_dbm <- lapply(fasta_db[inds_m], function (x) gsub("^-M", "-", x))
  
  # --- with protein N-term (initiator) methionine ---
  peps <- lapply(fasta_db, function (x) {
    s <- .Internal(strsplit(x, "@", fixed = FALSE, perl = FALSE, useBytes = FALSE))
    s <- .Internal(unlist(s, recursive = FALSE, use.names = FALSE))
  })
  
  # --- without protein N-term (initiator) methionine ---
  peps_m <- lapply(fasta_dbm, function (x) {
    s <- .Internal(strsplit(x, "@", fixed = FALSE, perl = FALSE, useBytes = FALSE))
    s <- .Internal(unlist(s, recursive = FALSE, use.names = FALSE))
    
    keep_n_misses(s, max_miss)
  })
  
  # Note: NA sequences -> NULL during mass calculations
  # (USE.NAMEs as they are prot_acc)
  if (length(inds_m)) {
    peps[inds_m] <- mapply(list, peps_m, peps[inds_m], SIMPLIFY = FALSE)
    peps[-inds_m] <- lapply(peps[-inds_m], function (x) list(NA, x))
  }
  else {
    peps <- lapply(peps, function (x) list(NA, x))
  }

  invisible(peps)
}


#' Splits fasta sequences by proteins.
#' 
#' For \code{no enzyme} searches.
#' 
#' @param aa_masses A lookup table of the masses of amino-acid residues.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
split_fastaseqs_noenz <- function (fasta = NULL, acc_type = "uniprot_acc", 
                                   acc_pattern = NULL, maxn_fasta_seqs = 200000L, 
                                   min_len = 7L, max_len = 40L, aa_masses = NULL, 
                                   ftmass = 18.010565, digits = 5L) 
{
  message("Loading fasta databases.")
  
  fasta_db <- load_fasta2(fasta, acc_type, acc_pattern)
  
  if (length(fasta_db) > maxn_fasta_seqs) 
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.\n",
         "  May consider a higher `maxn_fasta_seqs`.")

  n_cores <- detect_cores(16L)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(
    cl,
    c("make_noenzpeps", 
      "mmake_noenzpeps", 
      "hmake_noenzpeps", 
      "ms1masses_bare_noenz"), 
    envir = environment(proteoM:::make_noenzpeps))
  
  message("Splitting fasta sequences.")
  
  peps <- parallel::clusterApply(cl, chunksplit(fasta_db, n_cores), 
                                 mmake_noenzpeps, 
                                 min_len = min_len, 
                                 max_len = max_len, 
                                 aa_masses = aa_masses, 
                                 ftmass = ftmass, 
                                 digits = digits) %>%
    purrr::flatten()
  
  parallel::stopCluster(cl)
  
  invisible(peps)
}


#' Helper: multiple applications of \link{make_noenzpeps}.
#' 
#' For \code{no enzyme} searches.
#' 
#' @param fasta_db Fasta database(s).
#' @param aa_masses A lookup table of the masses of amino-acid residues.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
mmake_noenzpeps <- function (fasta_db = NULL, min_len = 7L, max_len = 40L, 
                             aa_masses = NULL, ftmass = 18.010565, digits = 5L) 
{
  lapply(fasta_db, make_noenzpeps, min_len = min_len, max_len = max_len, 
         aa_masses = aa_masses, ftmass = ftmass, digits = digits)
}


#' Makes peptides and masses for a protein.
#'
#' For \code{no enzyme} searches.
#'
#' @param prot A FASTA entry of protein.
#' @param aa_masses A lookup table of the masses of amino-acid residues.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
#' @examples 
#' \donttest{
#' prot <- paste0(
#'   "MSSKQHCVKLNDGHLIPALGFGTYKPKEVPKSKSLEAACLA", 
#'   "LDVGYRHVDTAYAYQVEEEIGQAIQSKIKAGVVKREDLFIT", 
#'   "TKLWCTCFRPELVKPALEKSLKKLQLDYVDLYIMHYPVPMK", 
#'   "SGDNDFPVNEQGKSLLDTVDFCDTWERLEECKDAGLVKSIG", 
#'   "VSNFNHRQLERILNKPGLKYKPVCNQVECHLYLNQRKLLDY",
#'   "CESKDIVLVAYGALGTQRYKEWVDQNSPVLLNDPVLCDVAK", 
#'   "KNKRSPALIALRYLIQRGIVPLAQSFKENEMRENLQVFGFQ", 
#'   "LSPEDMKTLDGLNKNFRYLPAEFLVDHPEYPFVEEY")
#' 
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'               "Carbamidomethyl (C)")
#' 
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", 
#'             "Deamidated (N)","Gln->pyro-Glu (N-term = Q)")
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' aa_masses <- aa_masses_all[[1]]
#' 
#' ans <- make_noenzpeps(prot, 7, 40, aa_masses)
#' 
#' # short FASTA (both N-term and C-term on a sequence)
#' prot <- substring(prot, 1, 10)
#' ans <- make_noenzpeps(prot, 7, 40, aa_masses)
#' }
make_noenzpeps <- function (prot = NULL, min_len = 7L, max_len = 40L, 
                            aa_masses = NULL, ftmass = 18.010565, digits = 5L) 
{
  len <- nchar(prot)
  
  if (len < min_len)
    return(NULL)
  
  if (len == min_len)
    return(prot[[1]])
  
  max_len <- min(max_len, len)
  starts <- 1:(len - min_len + 1L)
  
  # (1) Finds sub-sequence with a `start` value
  # (2) Calculates masses
  # (3) Adds C-term tag "-" after (2); no effect on masses
  ans <- lapply(starts, hmake_noenzpeps, prot, min_len, max_len, len, 
                aa_masses, ftmass, digits)

  # N-term peptides
  nms_1 <- names(ans[[1]])
  names(ans[[1]]) <- paste0("-", nms_1)

  if (grepl("^M", nms_1[[1]]))
    names(ans[[2]]) <- paste0("-", names(ans[[2]]))
  
  .Internal(unlist(ans, recursive = FALSE, use.names = TRUE))
}


#' Helper of \link{make_noenzpeps} (one start positions).
#'
#' For \code{no enzyme} searches.
#' 
#' C-term tagged with "-" after the mass calculations (but not yet N-term).
#'
#' @param prot A FASTA entry of protein.
#' @param start The staring position of an amino-acid in a \code{prot}.
#' @param len The number of amino-acid residues in a \code{prot}.
#' @param aa_masses A lookup table of the masses of amino-acid residues.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
#' 
#' @examples 
#' \donttest{
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'               "Carbamidomethyl (C)")
#' 
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", 
#'             "Deamidated (N)","Gln->pyro-Glu (N-term = Q)")
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' aa_masses <- aa_masses_all[[1]]
#' 
#' aas <- LETTERS[LETTERS %in% names(aa_masses)]
#' prot <- paste0(aas, collapse = "")
#' len <- nchar(prot)
#' masses <- hmake_noenzpeps(1, prot, 7, 40, len, aa_masses)
#' }
hmake_noenzpeps <- function (start = 1L, prot = NULL, min_len = 7L, max_len = 40L, 
                             len = NULL, aa_masses = NULL, ftmass = 18.010565, 
                             digits = 5L) 
{
  end_fi <- min_len + start  - 1L
  end_la <- min(max_len + start - 1L, len)
  ends <- end_fi:end_la
  
  peps <- substring(prot, start, ends)
  masses <- ms1masses_bare_noenz(peps, aa_masses, ftmass, digits)

  if (end_la == len) {
    len_a <- length(peps)
    names(masses)[len_a] <- paste0(peps[len_a], "-")
  }
  
  masses
}


#' Calculates masses of peptides.
#' 
#' For no enzyme workflow.
#' 
#' @param x A vector of peptide sequences.
#' @param aa_masses A lookup table of the masses of amino-acid residues.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @param digits The number of decimal places.
#' @examples
#' \donttest{
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'               "Carbamidomethyl (C)")
#' 
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", 
#'             "Deamidated (N)","Gln->pyro-Glu (N-term = Q)")
#' 
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' aa_masses <- aa_masses_all[[1]]
#' 
#' x <- c("MSSKQHC", "MSSKQHCV", "MSSKQHCVK", "MSSKQHCVKL")
#' masses1 <- ms1masses_bare_noenz(x, aa_masses)
#' 
#' x <- c("MSSKQHC", "MSSKQHCV", "MSSKQHCVK", "MSSKQHCVKL-")
#' masses2 <- ms1masses_bare_noenz(x, aa_masses)
#' 
#' stopifnot(identical(unname(masses1), unname(masses2)))
#' }
ms1masses_bare_noenz <- function (x, aa_masses, ftmass = 18.010565, digits = 5L) 
{
  len <- length(x)
  aas <- .Internal(strsplit(x[len], "", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  
  len_a <- length(aas)
  
  # No need: paste of N-term "-" after `hmake_noenzpeps`
  # if (aas[1] == "-") aas <- aas[2:len_a]
  
  ## No need: paste of C-term "-" after `ms1masses_bare_noenz`
  # if (aas[len_a] == "-") aas <- aas[1:(len_a - 1L)]
  
  masses <- cumsum(aa_masses[aas])
  masses <- tail(masses, len) + ftmass
  
  names(masses) <- x
  
  masses
}


#' Finds mis-cleavages in a vector.
#'
#' A convenience utility may be used to extract the first \eqn{n+1} peptides
#' from 0 to n mis-cleavages. It also assumes that the data were already sorted
#' in a desirable way.
#'
#' @param x A vector of data.
#' @param n Integer. The number of mis-cleavages.
keep_n_misses <- function (x, n) 
{
  len <- length(x)
  
  if (n < 0L) 
    stop("`n` cannot be nagative integers: ", n)
  
  if (!len) 
    stop("Length of `x` cannot be zero.")
  
  x[1:min(n + 1, len)]
}


#' Excludes mis-cleavages in a vector.
#'
#' @inheritParams keep_n_misses
#' @seealso keep_n_misses
exclude_n_misses <- function (x, n) 
{
  len <- length(x)
  
  if (n < 0L) 
    stop("`n` cannot be nagative integers: ", n)
  
  if (!len) 
    stop("Length of `x` cannot be zero.")
  
  x[-(1:min(n + 1, len))]
}


#' Excludes a character in string counting
#'
#' @param x A character string
#' @param char A character to be excluded for counting.
#' @importFrom stringi stri_length stri_count_fixed
str_exclude_count <- function (x, char = "-") 
{
  stringi::stri_length(x) - stringi::stri_count_fixed(x, char)
}


#' Removes a starting character from the first \code{n} entries.
#' 
#' @param x A list of character strings. Peptide sequences in names and masses
#'   in values.
#' @param char A starting character to be removed.
#' @param n The number of beginning entries to be considered.
#' @inheritParams matchMS
rm_char_in_nfirst <- function (x, char = "-", n = (max_miss + 1L) * 2L, 
                               max_len = 40L) 
{
  nms <- names(x)
  len <- length(nms)
  n <- min(len, n)
  seqs <- seq_len(n)
  
  # the possible space
  nms2 <- nms[seqs]
  idxes <- which(stringi::stri_startswith_fixed(nms2, char))
  
  # the exact space
  nms3 <- nms2[idxes]
  nms3 <- substr(nms3, 2L, max_len)
  
  # update the possible space
  nms2[idxes] <- nms3
  
  # update the full space
  nms[seqs] <- nms2
  names(x) <- nms
  
  attr(x, "pnt_idxes") <- idxes
  
  x
}


#' Removes a trailing character from the last \code{n} entries.
#' 
#' The default value of \code{n} is for full-enzymatic peptides.
#' 
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst
rm_char_in_nlast <- function (x, char = "-", n = (max_miss + 1L) * 2L) 
{
  nms <- names(x)
  len <- length(nms)
  n <- min(len, n)
  seqs <- (len - n + 1L):len
  
  # possible space
  nms2 <- nms[seqs]
  idxes <- which(stringi::stri_endswith_fixed(nms2, char))
  
  # exact space
  nms3 <- nms2[idxes]
  stops <- nchar(nms3) - 1L
  nms3 <- substr(nms3, 1L, stops)
  
  # update possible space
  nms2[idxes] <- nms3
  
  # update full space
  nms[seqs] <- nms2
  names(x) <- nms
  
  attr(x, "pct_idxes") <- len - length(seqs) + idxes
  
  x
}


#' Helper in calculating peptide masses.
#'
#' (2) "amods- tmod+ vnl- fnl-".
#' 
#' @param peps A list of peptide sequences.
#' @inheritParams add_fixvar_masses
#' @inheritParams distri_peps
#' @inheritParams matchMS
add_term_mass2 <- function (aa_masses, peps, max_mass = 4500L) 
{
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  # No needs of is_empty(ntmod) && is_empty(ctmod)
  delta <- if (length(ntmod) && length(ctmod)) 
    aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  else if (length(ntmod)) 
    aa_masses[names(ntmod)]
  else if (length(ctmod)) 
    aa_masses[names(ctmod)]

  ans <- peps + delta
  
  ans[ans <= max_mass]
}


#' Calculates mono-isotopic masses of peptide sequences.
#'
#' @param seqs Sequences of peptides from FASTAs by protein accessions. Each
#'   list contains two lists of sequences: (1) without and (2) with N-terminal
#'   methionine.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
#' @inheritParams add_fixvar_masses
#' @inheritParams distri_fpeps
ms1masses_bare <- function (seqs = NULL, aa_masses = NULL, ftmass = NULL,
                            max_miss = 2L, min_len = 7L, max_len = 40L,
                            min_mass = 700L, max_mass = 4500L, 
                            maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                            is_fixed_protnt = FALSE, is_fixed_protct = FALSE, 
                            digits = 4L) 
{
  # (1) before rolling sum (not yet terminal H2O)
  # (1.1) without N-term methionine
  data_1 <- lapply(seqs, `[[`, 1)
  data_1 <- ms1masses_noterm(data_1, aa_masses = aa_masses,
                             maxn_vmods_per_pep = maxn_vmods_per_pep,
                             maxn_sites_per_vmod = maxn_sites_per_vmod,
                             digits = digits)
  data_1 <- attr(data_1, "data")

  # (1.2) with N-term methionine
  data_2 <- lapply(seqs, `[[`, 2)
  data_2 <- ms1masses_noterm(data_2, aa_masses = aa_masses,
                             maxn_vmods_per_pep = maxn_vmods_per_pep,
                             maxn_sites_per_vmod = maxn_sites_per_vmod,
                             digits = digits)
  data_2 <- attr(data_2, "data")

  # (2) rolling sum (not yet terminal H2O)
  n_cores <- detect_cores(16L)
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

  parallel::clusterExport(cl, c("roll_sum", "accumulate_char"), 
                          envir = environment(proteoM:::roll_sum))

  ms_1 <- parallel::clusterApply(
    cl = cl, 
    x = chunksplit(data_1, n_cores, "list"), 
    fun = lapply, 
    FUN = "roll_sum", 
    n = max_miss, 
    include_cts = FALSE
  ) 
  
  ms_1 <- purrr::flatten(ms_1)
  
  # ms_1 <- lapply(ms_1, function (x) x[x >= min_mass & x <= max_mass])

  if (is_fixed_protnt) {
    ms_2 <- parallel::clusterApply(
      cl = cl, 
      x = chunksplit(data_2, n_cores, "list"), 
      fun = lapply, 
      FUN = "roll_sum", 
      n = max_miss, 
      include_cts = FALSE
    ) 
    
    ms_2 <- purrr::flatten(ms_2)
  } 
  else {
    ms_2 <- parallel::clusterApply(
      cl = cl, 
      x = chunksplit(data_2, n_cores, "list"), 
      fun = lapply, 
      FUN = "roll_sum", 
      n = max_miss, 
      include_cts = TRUE
    ) 
    
    ms_2 <- purrr::flatten(ms_2)
  }
  
  parallel::stopCluster(cl)
  
  # ms_2 <- lapply(ms_2, function (x) x[x >= min_mass & x <= max_mass])

  # (3) putting together (+ terminal H2O)
  # (USE.NAMES of prot_acc)
  ms <- mapply(`c`, ms_1, ms_2, SIMPLIFY = FALSE, USE.NAMES = TRUE)

  if (min_len > 1L && !is.infinite(max_len)) {
    ms <- lapply(ms, function (x) {
      cts <- str_exclude_count(names(x), "-")
      x[cts >= min_len & cts <= max_len]
    })
  }

  # "HD101_HUMAN" etc. has no tryptic peptides
  lens <- unlist(lapply(ms, length), recursive = FALSE, use.names = FALSE)

  # adding H2O or fixed N/C-term masses
  ms <- ms[lens > 0L]
  ms <- lapply(ms, function (x) x[!duplicated.default(names(x))] + ftmass)
  ms <- lapply(ms, function (x) x[x >= min_mass & x <= max_mass])

  invisible(ms)
}


#' Adds Carbon-13 masses.
#'
#' @param peps A named vector of peptide sequences. Sequences in names and
#'   masses in values.
#' @inheritParams matchMS
add_ms1_13c <- function (peps, n_13c = 1L, max_mass = 4500L) 
{
  # consider -1L too...
  
  if (n_13c <= 0L) 
    return(peps)
  
  len <- n_13c + 1L
  out <- vector("list", len)
  out[[1]] <- peps
  
  for (i in 2:len) 
    out[[i]] <- out[[i-1]] + 1.00335483
  
  out <- .Internal(unlist(out, recursive = FALSE, use.names = TRUE))
  
  out[out <= max_mass]
}


#' Helper of \link{ms1masses_bare}.
#'
#' For either forward or reversed sequences.
#' 
#' @param aa_seqs Character string; a vector of peptide sequences with
#'   one-letter representation of amino acids.
#' @inheritParams matchMS
#' @inheritParams add_fixvar_masses
#' @inheritParams distri_peps
ms1masses_noterm <- function (aa_seqs, aa_masses, maxn_vmods_per_pep = 5L,
                              maxn_sites_per_vmod = 3L, digits = 4L) 
{
  options(digits = 9L)

  n_cores <- detect_cores(16L)

  aa_seqs <- suppressWarnings(split(aa_seqs, seq_len(n_cores)))

  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

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
                                digits = digits) 
  
  parallel::stopCluster(cl)
  
  out <- purrr::flatten(out)

  attr(aa_masses, "data") <- out

  rm(list = c("aa_seqs", "out"))
  gc()

  invisible(aa_masses)
}


#' Helper function for parallel calculations peptide masses by proteins.
#'
#' For each split of multiple proteins; no terminal masses.
#'
#' @inheritParams ms1masses_noterm
calcms1mass_noterm <- function (aa_seqs, aa_masses, maxn_vmods_per_pep = 5L, 
                                maxn_sites_per_vmod = 3L, digits = 4L) 
{
  lapply(aa_seqs, function (x) 
    calcms1mass_noterm_byprot(prot_peps = x,
                              aa_masses = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod,
                              digits = digits))
}


#' Helper of \link{calcms1mass_noterm}.
#'
#' For single protein.
#'
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams ms1masses_noterm
calcms1mass_noterm_byprot <- function (prot_peps, aa_masses, 
                                       maxn_vmods_per_pep = 5L,
                                       maxn_sites_per_vmod = 3L, digits = 4L) 
{
  # by peptides
  ans <- lapply(prot_peps, calcms1mass_noterm_bypep,
                aa_masses = aa_masses,
                maxn_vmods_per_pep = maxn_vmods_per_pep,
                maxn_sites_per_vmod = maxn_sites_per_vmod,
                digits = digits)
  
  .Internal(unlist(ans, recursive = FALSE, use.names = TRUE))
}


#' Helper of \link{calcms1mass_noterm_byprot}.
#'
#' For single protein.
#'
#' @param aa_seq Character string; a peptide sequence with one-letter
#'   representation of amino acids.
#' @inheritParams ms1masses_noterm
#' @importFrom stringr str_split
calcms1mass_noterm_bypep <- function (aa_seq, aa_masses, maxn_vmods_per_pep = 5L,
                                      maxn_sites_per_vmod = 3L, digits = 4L) 
{
  if (is.na(aa_seq)) return(NULL)
  
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))

  aas <- aa_masses[aas]
  aas <- sum(aas)
  names(aas) <- aa_seq
  
  round(aas, digits = digits)
}


#' Distributes peptides by variable modifications.
#'
#' @param prps Lists of peptide sequences with a one-letter representation of
#'   amino acid residues. Each list is named by protein accession.
#' @param aa_masses_all All the amino acid look-up tables.
#' @param motifs_all Lists of motifs of modifications.
#' @inheritParams calc_pepmasses2
#' @param enzyme A character string of enzyme. The information is used for
#'   faster replacement of "-" in peptide sequences from protein C-terminals.
distri_peps <- function (prps, aa_masses_all, motifs_all, max_miss = 2L, 
                         max_len = 40L, enzyme = "trypsin_p") 
{
  # proteins without applicable peptide sequences
  # bads <- lapply(prps, is.null)
  # bads <- .Internal(unlist(bads, use.names = FALSE, recursive = FALSE))
  # prps <- prps[!bads]
  # rm(list = "bads")
  
  nms <- lapply(prps, names)
  
  out <- mapply(function (aa_masses, motifs) {
    nms_sub <- subpeps_by_vmods(aa_masses, nms, motifs = motifs)
    
    # USE.NAMEs of prot_acc
    idxes <- mapply(fastmatch::fmatch, nms_sub, nms, 
                    SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    mapply(function (x, y) x[y], prps, idxes, 
           SIMPLIFY = FALSE, USE.NAMES = TRUE)
  }, aa_masses_all, motifs_all, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  # semi enzymes: don't know the maximum number of 
  # C-term peptides that can end with "-"
  # (N-term remains the same)
  n1 <- if (enzyme == "noenzyme")
    Inf
  else 
    (max_miss + 1L) * 2L

  n2 <- if (grepl("^semi", enzyme) || enzyme == "noenzyme")
    Inf
  else
    ct_counts(max_miss)

  # ZN207_HUMAN: MGRKKKK (no N-term pep_seq at 2 misses and min_len >= 7L)
  
  out <- lapply(out, function(xs) {
    len <- .Internal(unlist(lapply(xs, length), recursive = FALSE, use.names = FALSE))
    xs <- xs[len > 0L]
    xs <- lapply(xs, rm_char_in_nfirst, char = "-", n = n1, max_len = max_len)
    xs <- lapply(xs, rm_char_in_nlast, char = "-", n = n2)
  })
}


#' Counts the number of trailing residues from C-term for the replacment of "-".
#'
#' For full-enzymes: n(i+1) = n(i) + (i+1). Not applicable for semi-enzymes.
#'
#' @param max_miss The maximum number of cleavages.
ct_counts <- function (max_miss = 2L) 
{
  ct <- integer(max_miss)

  if (max_miss > 0L) {
    ct[1] <- 2L

    for (i in 1:max_miss) {
      j <- i + 1L
      ct[j] <- ct[i] + j
    }
  } 
  else
    return(1L)

  ct[max_miss]
}


#' Distributes peptides by fixed modifications.
#'
#' @inheritParams calc_pepmasses2
#' @inheritParams matchMS
#' @param is_fixed_protnt Logical; is protein N-terminal modification fixed.
#' @param is_fixed_protct Logical; is protein C-terminal modification fixed.
#' @param data Lists of peptides by prot_accs.
distri_fpeps <- function (data = NULL, max_miss = 2L, is_fixed_protnt = FALSE, 
                          is_fixed_protct = FALSE) 
{
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


#' Concatenates adjacent peptides in a list (with mass).
#' 
#' @param peps A list of peptide sequences with a one-letter representation of
#'   amino acid residues.
#' @param n The number of mis-cleavages for consideration.
#' @param include_cts Logical; the list, \code{peps}, includes the protein
#'   C-terminal sequence or not. At the default of TRUE, mis-cleaved peptides at
#'   the end of the protein C-terms will be added as they should. The arguments
#'   would be typically at FALSE, for example, when used for generating
#'   mis-cleaved peptides from the N-terminal of peptides with the removal of a
#'   starting residue \code{M}.
#' @examples
#' \donttest{
#' peps <- 1:26
#' names(peps) <- LETTERS
#' res1 <- roll_sum(peps, 2)
#' res2 <- roll_sum(peps, 2, FALSE)
#' 
#' # length shorter than n
#' peps <- c(a = 1)
#' res1 <- roll_sum(peps, 2)
#' res2 <- roll_sum(peps, 2, FALSE)
#' 
#' peps <- c(a = 1, b = 2, c = 3)
#' res1 <- roll_sum(peps, 4)
#' res2 <- roll_sum(peps, 4, FALSE)
#' }
roll_sum <- function (peps = NULL, n = 2L, include_cts = TRUE) 
{
  len <- length(peps)
  
  if (!len) 
    return(NULL)
  
  if (n >= len) 
    n <- len - 1L
  
  res <- lapply(seq_len(len - n), function (x) {
    psub <- peps[x:(x + n)]
    
    vals <- cumsum(psub)
    nms <- accumulate_char(names(psub), paste0)
    names(vals) <- nms
    
    vals
  }) 
  
  res <- .Internal(unlist(res, recursive = FALSE, use.names = TRUE))
  
  if (include_cts && n >= 1L) {
    ends <- peps[(len - n + 1L):len]
    
    res2 <- lapply(n:1L, function (x) {
      psub <- tail(ends, x)
      
      vals <- cumsum(psub)
      nms <- accumulate_char(names(psub), paste0)
      names(vals)  <- nms
      
      vals
    })
    
    res2 <- .Internal(unlist(res2, recursive = FALSE, use.names = TRUE))
    ans <- c(res, res2)
  } 
  else
    ans <- res
  
  invisible(ans)
}


#' Helper of \link{semipeps_byprots}.
#'
#' @param prots Lists of proteins with full-enzymatic sequences. For each entry
#'   under a protein, the value is a mass and the name is a peptide sequence.
#' @param min_len The minimum length of peptide sequences for consideration.
#' @inheritParams add_fixvar_masses
hsemipeps_byprots <- function (prots, min_len = 7L, aa_masses = NULL) 
{
  lapply(prots, semipeps_byprots, min_len, aa_masses)
}


#' Finds the semi-enzymatic peptides for a proteins.
#'
#' Redundancy such as peptides from N-term methionine cleavage or semi-tryptic
#' ladders are handled.
#'
#' @param vals A list of full-enzymatic peptides under a protein.
#' @inheritParams matchMS
#' @inheritParams add_fixvar_masses
#' @return A vector of full- and semi-enzymatic peptides. Sequences in names and
#'   masses in values.
semipeps_byprots <- function (vals, min_len = 7L, aa_masses = NULL) 
{
  peps <- names(vals)
  peps <- gsub("^-", "", peps)
  len <- length(peps)
  
  cts <- grepl("-$", peps)
  ots <- !cts
  peps_ct <- peps[cts]
  peps_ot <- peps[ots]
  vals_ct <- vals[cts]
  vals_ot <- vals[ots]
  
  semis_ot <- mapply(calc_semipepmasses, vals_ot, peps_ot, 
                     MoreArgs = list(min_len = min_len, 
                                     aa_masses = aa_masses, 
                                     ct_offset = 0L), 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
  semis_ot <- .Internal(unlist(semis_ot, recursive = FALSE, use.names = TRUE))
  
  semis_ct <- mapply(calc_semipepmasses, vals_ct, peps_ct, 
                     MoreArgs = list(min_len = min_len, 
                                     aa_masses = aa_masses, 
                                     ct_offset = 1L), 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
  semis_ct <- .Internal(unlist(semis_ct, recursive = FALSE, use.names = TRUE))
  
  ans <- c(vals, semis_ot, semis_ct)
  ans <- ans[!duplicated.default(ans)]
  
  invisible(ans)
}


#' Finds and calculates the masses of semi-enzymatic sequences.
#'
#' Semi-enzymatic sequences are built on full-enzymatic where the N-term
#' residues are removed sequentially.
#'
#' N-term "-" has no effect on the semi-enzymatic generation since
#' semi-enzymatic sequences are always NONE N-term (original enzymatic sequences
#' keep separately and concatenated later). Thus, the N-term tag of "-" was
#' removed from \code{pep} in \link{semipeps_byprots}, before calling this
#' function.
#'
#' @param val A mass of a peptide.
#' @param pep A character string of peptide.
#' @param ct_offset Zero or one to account for the "-" in the C-term of pep.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
#'               "Carbamidomethyl (C)")
#'
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'             "Deamidated (N)","Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' aa_masses <- aa_masses_all[[1]]
#'
#' val <- 4237.89756
#' pep <- "ALELNQSAEYYYEENEMNYTHDYSQYEVICIK"
#' ans <- calc_semipepmasses(val, pep, 7, aa_masses)
#'
#' val <- 2423.1017
#' pep <- "QNVEEIPFDSEGPTEPTSSFTI-"
#' ans <- calc_semipepmasses(val, pep, 7, aa_masses, 1L)
#' }
calc_semipepmasses <- function (val, pep, min_len = 7L, aa_masses = NULL, 
                                ct_offset = 0L) 
{
  options(digits = 9L)
  
  len <- nchar(pep)
  len2 <- len - ct_offset
  
  if (len2 <= min_len) 
    return(NULL)
  
  span <- len2 - min_len
  
  semipeps <- substring(pep, 2:(span + 1L), len)
  aas <- .Internal(strsplit(pep, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  aas <- aas[1:span]
  aas2 <- aa_masses[aas]
  delta <- cumsum(aas2)
  
  ans <- val - delta
  names(ans) <- semipeps
  
  invisible(ans)
}


#' Helper of peptide-mass calculation..
#'
#' (5) "amods- tmod+ vnl- fnl+"; (6) "amods- tmod- vnl- fnl+".
#'
#' The calculation goes through the rows in \code{fnl_combi}.
#' 
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @inheritParams hms1_a0_vnl0_fnl1
#' @return A numeric vector
delta_ms1_a0_fnl1 <- function (fnl_combi, aas, aa_masses) 
{
  nms <- lapply(fnl_combi, names)
  nms <- .Internal(unlist(nms, recursive = FALSE, use.names = FALSE))
  oks <- aas[aas %in% nms]

  if (!length(oks)) 
    return (0L)

  len <- length(fnl_combi)
  out <- vector("numeric", len)
  out[[1]] <- 0L

  for (i in 2:len) {
    row <- fnl_combi[[i]]
    aa_masses[nms] <- .Internal(unlist(row, recursive = FALSE, use.names = FALSE))
    oks <- aas[aas %in% nms]
    out[[i]] <- sum(aa_masses[oks])
  }

  out[!duplicated.default(out)]
}


#' Helper of \link{ms1_a0_vnl0_fnl1}.
#' 
#' @param masses A named list of peptide masses.
#' @param fnl_combi A data.frame of combinations of neutral losses for fixed
#'   modifications. Each row corresponds to a set of neutral loss. The first row
#'   corresponds to the combination without NLs (all zeros).
#' @inheritParams matchMS
#' @inheritParams add_fixvar_masses
hms1_a0_vnl0_fnl1 <- function (masses, fnl_combi, aa_masses, max_mass = 4500L, 
                               digits = 4L) 
{
  mapply(ms1_a0_vnl0_fnl1, 
         masses, names(masses), 
         MoreArgs = list (
           fnl_combi = fnl_combi, 
           aa_masses = aa_masses,
           max_mass = max_mass, 
           digits = digits
         ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


#' Helper by individual peptides. 
#'
#' (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
#'
#' @param mass The mass of a peptide.
#' @param aa_seq Character string; a peptide sequence with one-letter
#'   representation of amino acids.
#' @inheritParams hms1_a0_vnl0_fnl1
#' @importFrom stringr str_split
ms1_a0_vnl0_fnl1 <- function (mass, aa_seq, fnl_combi, aa_masses, 
                              max_mass = 4500L,digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))

  delta <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
  out <- round(mass - delta, digits = digits)
  
  out <- out[out <= max_mass]
  names(out) <- rep(aa_seq, length(out))
  
  invisible(out)
}


#' Helper of \link{ms1_a1_vnl0_fnl0}.
#' 
#' @param masses A named list of peptide masses.
#' @param amods \code{Anywhere} variable modifications.
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @param vmods_nl The attribute of \code{vmods_nl} from an \code{aa_masses}.
#' @param ms1vmods The set of all possible MS1 vmod labels at a given aa_masses.
#' @inheritParams matchMS
#' @inheritParams add_fixvar_masses
hms1_a1_vnl0_fnl0 <- function (masses, amods, aa_masses,
                               vmods_nl = NULL, fmods_nl = NULL,
                               maxn_vmods_per_pep = 5L,
                               maxn_sites_per_vmod = 3L,
                               ms1vmods = NULL, 
                               max_mass = 4500L, 
                               digits = 4L) 
{
  mapply(ms1_a1_vnl0_fnl0, 
         masses, names(masses), 
         MoreArgs = list (
           amods = amods, 
           aa_masses = aa_masses,
           vmods_nl = vmods_nl, 
           fmods_nl = fmods_nl,
           maxn_vmods_per_pep = maxn_vmods_per_pep,
           maxn_sites_per_vmod = maxn_sites_per_vmod,
           ms1vmods = ms1vmods, 
           max_mass = max_mass, 
           digits = digits
         ), SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


#' Helper by individual peptides.
#'
#' (7) "amods+ tmod- vnl- fnl-"; (8) "amods+ tmod+ vnl- fnl-".
#'
#' @param mass The mass of a peptide.
#' @param aa_seq Character string; a peptide sequence with one-letter
#'   representation of amino acids.
#' @inheritParams hms1_a1_vnl0_fnl0
#' @importFrom stringr str_split
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
#'                                            "Carbamidomethyl (C)"))
#'
#' pep <- c("HQGVMCNVGMGQKMNSC" = 2051.90346)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"),
#'               `Carbamidomethyl (C)` = c("Anywhere" = "C"))
#'
#' x <- ms1_a1_vnl0_fnl0(pep, names(pep), amods, aa_masses_all[[4]])
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
#'                                            "Acetyl (Protein N-term)"))
#'
#' pep <- c("HQGVMNVGMGQKSMNS" = 1932.9171) # + TMT6plex (K)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"),
#'               `Carbamidomethyl (S)` = c("Anywhere" = "S"))
#'
#' x <- ms1_a1_vnl0_fnl0(pep, names(pep), amods, aa_masses_all[[8]])
#'
#' stopifnot(x[[1]] - 1990.9226 < 1e-4,
#'           x[[2]] - 1991.9066 < 1e-4,
#'           x[[3]] - 2047.9440 < 1e-4,
#'           x[[4]] - 2048.9281 < 1e-4)
#' 
#' 
#' # (8-b)
#' .ms1_vmodsets <- make_ms1_vmodsets(aa_masses_all = aa_masses_all, 
#'                                    maxn_vmods_per_pep = 5L, 
#'                                    maxn_sites_per_vmod = 3L)
#' .base_ent <- lapply(.ms1_vmodsets, `[[`, 1)
#' 
#'  x <- ms1_a1_vnl0_fnl0(pep, names(pep), amods, aa_masses_all[[8]], 
#'                       .ms1_vmodsets = .ms1_vmodsets, 
#'                       .base_ent = .base_ent)
#' 
#' # (8-c)
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' 
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'             "Gln->pyro-Glu (N-term = Q)")
#' 
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                out_path = NULL)
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i)
#' 
#' i <- 10
#' aa_masses <- aa_masses_all[[i]]
#' ms1vmods <- ms1vmods_all[[i]]
#' 
#' pep <- c("HQGVMCNVGMGQKMNSC" = 2051.90346)
#' amods <- attr(aa_masses, "amods")
#' 
#' x <- ms1_a1_vnl0_fnl0(pep, names(pep), amods, aa_masses_all[[i]], 
#'                       ms1vmods = ms1vmods)
#' 
#' # y <- ms1_a1_vnl0_fnl0(pep, names(pep), amods, aa_masses_all[[i]], 
#' #                       ms1vmods = NULL)
#' 
#' # identical(x, y)
#' }
ms1_a1_vnl0_fnl0 <- function (mass, aa_seq, amods, aa_masses,
                              vmods_nl = NULL, fmods_nl = NULL,
                              maxn_vmods_per_pep = 5L,
                              maxn_sites_per_vmod = 3L,
                              ms1vmods = NULL, 
                              max_mass = 4500L, 
                              digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  
  vmods_combi <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)$ms1
  
  # "Carbamidomethyl (M)", "Carbamyl (M)" requires two "M"s
  # aas may only have one "M"
  if (!length(vmods_combi)) return(NULL)

  deltas <- lapply(vmods_combi, function (x) sum(aa_masses[x]))
  
  masses <- lapply(deltas, function (x) round(mass + x, digits = digits))
  out <- .Internal(unlist(masses, recursive = FALSE, use.names = FALSE))

  if (FALSE) {
    if (length(vmods_nl)) {
      vnl_combi <- lapply(vmods_combi, function (x) expand_grid_rows(vmods_nl[x]))

      deltas_vnls <- lapply(vnl_combi, function (x) {
        ss <- lapply(x, sum)
        ss <- .Internal(unlist(ss, recursive = FALSE, use.names = FALSE))
        ss <- unique(ss)
      })

      out <- mapply(`-`, out, deltas_vnls, SIMPLIFY = FALSE)
    }
    
    if (length(fmods_nl)) {
      fnl_combi <- expand_grid_rows(fmods_nl)
                               
      deltas_fnls <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
      out <- mapply(`-`, out, deltas_fnls, SIMPLIFY = FALSE)
    }
  }
  
  # should be NULL if called from matchMS(max_mass = ...)
  out <- out[out <= max_mass]
  names(out) <- rep(aa_seq, length(out))
  
  invisible(out)
}


