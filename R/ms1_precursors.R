#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param aa_masses An amino acid mass look-up.
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
calc_pepmasses2 <- function (aa_masses = NULL, 
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
  fixedlabs = NULL, 
  varlabs = NULL, 
  mod_motifs = NULL, 
  enzyme = c("trypsin_p"),
  custom_enzyme = c(Cterm = NULL, Nterm = NULL), 
  noenzyme_maxn = 0L, 
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
  # (2.4) Adds coerced fixed masses
  #     E.g., with the original fixed TMT6plex (K) and variable Acetyl (K),
  #     TMT6plex (K) is back-coerced to fixedmod under some combinations.
  # (2.5) Adds variable masses with the look-up of aa_masses_all
  # 
  # Note:
  # For "historical" reasons, terminal "tmods+" refers to "VARIABLE tmods+". 
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
    args <- c("noenzyme_maxn")
    fmls <- formals(fun)[args]
    
    # for example `custom_enzyme` (should be NULL after eval)
    nargs <- lapply(fmls, function (x) if (is.call(x)) eval(x) else x)
    
    # unlisting with NUlls being kept
    are_nulls <- lapply(nargs, is.null)
    def_nulls <- unlist(are_nulls)
    c(unlist(nargs), nargs[def_nulls])
  })
  
  ### currently not new argument to bypass
  # new_args <- NULL
  ###
  
  .time_stamp <- match_calltime(
    path = .path_cache,
    fun = fun,
    
    # `nms` must be matched in order to retrieve cached results
    nms = c("fasta", "acc_type", "acc_pattern",
            "fixedmods", "varmods", "mod_motifs", 
            "fixedlabs", "varlabs", 
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
      aa_masses = aa_masses, 
      out_path = file.path(.path_fasta, "ms1masses", .time_stamp),
      fixedmods = fixedmods,
      varmods = varmods,
      varlabs = varlabs, 
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
    # `mgf_quries.rds` kept (only affected by min_mass, max_mass and ppm_ms1)
    delete_files(out_path, 
                 ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$",
                             "\\.xls$", "\\.csv$", "\\.txt$",
                             "^mgf$", "^mgfs$"))
    
    .time_stamp <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")
    path_tstamp <- file.path(.path_fasta, "ms1masses", .time_stamp)
    file_ms1 <- file.path(path_tstamp, "aa_masses_ms1.rds")

    aa_masses_all <- find_aa_masses(
      aa_masses = aa_masses, 
      out_path = path_tstamp,
      fixedmods = fixedmods,
      varmods = varmods,
      varlabs = varlabs,
      mod_motifs = mod_motifs, 
      maxn_vmods_setscombi = maxn_vmods_setscombi)
    
    if (file.exists(file_ms1))
      aa_masses_ms1 <- qs::qread(file_ms1)
    else
      stop("File not found: ", file_ms1)
    
    aa_masses_0 <- aa_masses_ms1[[1]]
    aa_masses_1 <- aa_masses_all[[1]]
    
    ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
                           maxn_vmods_per_pep = maxn_vmods_per_pep,
                           maxn_sites_per_vmod = maxn_sites_per_vmod)

    # Design:
    #   variable modifications, including [NC]-term, were appended in parallel 
    #   to unmodified residues in `aa_masses`, e.g., "M", "Oxidation (M)" are 
    #   two separate entries in `aa_masses`
    # 
    # aa_masses_0 - coerced fixedmods: the base for varmod mass calculations
    # aa_masses_1 - original fixedmods: for all-fixed mass calculations
    # 
    # The first entry in aa_masses_all is aa_masses_1, not aa_masses_0: 
    #   Fixed "TMT (K)" coerced to variable under the presence of "Acetyl (K)"; 
    #   If use aa_masses_0, sequences without K will drop when dispatching 
    #   (must contain K as varmods were considered realized).
    # 
    # "tmods-" refers to VARIABLE terminal modifications
    #   aa_masses_0 and aa_masses_1 are always "tmods-"
    
    is_fixed_protnt <- any(grepl("Protein N-term", fixedmods))
    is_fixed_protct <- any(grepl("Protein C-term", fixedmods))

    # --- Forward sequences  ---
    if (isTRUE(enzyme == "noenzyme")) {
      if (max_len > 25L) 
        warning("May be out of RAM at `max_len = ", max_len, "`.\n",
                "Consider a sectional search, e.g., `noenzyme_maxn = 10`")

      if (any(is_fixed_protnt, is_fixed_protct))
        stop("Not yet support FIXED protein terminal modifications for ", 
             "noenzyme searches.\n", 
             "Change to VARIABLE protein terminal modifications.")

      seqs_0 <- NULL
      ftmass <- unname(aa_masses_0["N-term"] + aa_masses_0["C-term"])
      
      fwd_peps <- split_fastaseqs_noenz(fasta = fasta, 
                                        acc_type = acc_type,
                                        acc_pattern = acc_pattern,
                                        maxn_fasta_seqs = maxn_fasta_seqs, 
                                        min_len = min_len, 
                                        max_len = max_len, 
                                        aa_masses = aa_masses_0, 
                                        ftmass = ftmass)
      
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
      
      ### Special case of FIXED Protein [NC]-term modification(s)
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
      
      seqs_0 <- distri_fpeps(data = seqs_0, max_miss = max_miss, 
                             is_fixed_protnt = is_fixed_protnt, 
                             is_fixed_protct = is_fixed_protct)
      ###

      
      # --- Masses of sequences: fixed mods + terminals ---
      message("Calculating bare peptide masses...")
      
      # e.g. if "TMT6plex (N-term)" is a fixedmod -> ftmass = (229 + 1) + 17;
      # if "TMT6plex (N-term)" coerced to a varmod -> ftmass = 1 + 17
      ftmass <- unname(aa_masses_0["N-term"] + aa_masses_0["C-term"])
      
      fwd_peps <- ms1masses_bare(seqs = seqs_0,
                                 aa_masses = aa_masses_0,
                                 ftmass = ftmass,
                                 max_miss = max_miss,
                                 min_len = min_len,
                                 max_len = max_len,
                                 maxn_vmods_per_pep = maxn_vmods_per_pep,
                                 maxn_sites_per_vmod = maxn_sites_per_vmod,
                                 is_fixed_protnt = is_fixed_protnt,
                                 is_fixed_protct = is_fixed_protct)

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
        max_len = max_len, 
        aa_masses = aa_masses_0
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

    # aa_masses_all[[1]] is for the original all-fixed mode not for the coerced,
    # otherwise, e.g. fixed to variable coercion of "TMT (K)" with a conflicting 
    # "Acetyl (K)", sequences without "K" will be dropped.

    fwd_peps <- parallel::clusterApply(
      cl, 
      fwd_peps, 
      distri_peps, 
      aa_masses_all = aa_masses_all, 
      motifs_all = motifs_all, 
      max_miss = max_miss, 
      max_len = max_len, # different purpose
      enzyme = enzyme
    )
    
    parallel::stopCluster(cl)
    
    fwd_peps <- lapply(seq_along(aa_masses_all), function (i) {
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
    
    # (d) Flattened peptide lists (prot_acc's removed)
    fwd_peps <- lapply(fwd_peps, flat_pepseqs)
    gc()

    # (e) Adjusted base masses if with fixed-to-variable coercion
    # (filtered by min_mass and max_mass since it is a final)
    fwd_peps[[1]] <- adj_base_masses(fwd_peps[[1]], aa_masses_0, aa_masses_1, 
                                     min_mass = min_mass, max_mass = max_mass, 
                                     digits = digits)
    
    gc()
    
    # --- Delta masses of `variable` terminals  ---
    # (e.g., on top of the `fixed` 18.010565)
    message("Adding terminal masses...")
    fwd_peps <- mapply(add_term_mass, fwd_peps, aa_masses_ms1, 
                       MoreArgs = list(min_mass = min_mass, max_mass = max_mass))
    
    message("Adding coerced fixed masses...")
    fwd_peps <- mapply(adj_anywhere_masses, fwd_peps, aa_masses_ms1)

    # --- Mass of variable mods and/or NLs ---
    message("Adding variable masses...")

    # (switch to aa_masses_all)
    types <- purrr::map_chr(aa_masses_all, attr, "type", exact = TRUE)

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
            min_mass = min_mass, 
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

    n_cores <- detect_cores(16L)

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
          min_mass = min_mass, 
          max_mass = max_mass, 
          digits = digits
        ) %>% 
          purrr::flatten() %>% 
          unlist(recursive = FALSE, use.names = TRUE)

        parallel::stopCluster(cl)
        gc()

        message("\tCompleted peptide masses: ",
                paste(attributes(aa_masses_i)$fmods, "|", 
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
find_aa_masses <- function(aa_masses = NULL, out_path = NULL, fixedmods = NULL, 
                           varmods = NULL, varlabs = NULL, 
                           mod_motifs = NULL, maxn_vmods_setscombi = 64L) 
{
  file <- file.path(out_path, "aa_masses_all.rds")
  
  if (file.exists(file))
    return(qs::qread(file))
  
  message("Computing the combinations of fixed and variable modifications.")
  
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                                 varmods = varmods,
                                 aa_masses = aa_masses, 
                                 varlabs = varlabs, 
                                 mod_motifs = mod_motifs, 
                                 maxn_vmods_setscombi = maxn_vmods_setscombi,
                                 out_path = out_path)

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
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"), 
#'                    varmods   = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))
#' 
#' stopifnot(length(x) == 6L)
#'
#' # Fixed N-term mod (no coercion to variable mod)
#' x <- calc_aamasses(fixedmods = "TMT6plex (N-term)", varmods = NULL)
#' x[[1]][["N-term"]]
#' 
#' # Fixed N-term mod (coerced to variable mod)
#' x <- calc_aamasses(fixedmods = "TMT6plex (N-term)", varmods = "Acetyl (Protein N-term)")
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
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (. = M)",
#'                      "Deamidated (. = R)"),
#'                    varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)",
#'                      "Hex(5)HexNAc(2) (N)"))
#' 
#' stopifnot(length(x) == 6L)
#' 
#' x <- calc_aamasses(c(fixedmods = "TMT6plex (N-term)", "TMT6plex (K)",
#'                      "Carbamidomethyl (. = M)", "Deamidated (. = R)"),
#'                    varmods = c("Acetyl (N-term)", "Carbamyl (. = M)",
#'                      "Gln->pyro-Glu (N-term = Q)", "Hex(5)HexNAc(2) (N)"))
#' 
#' stopifnot(length(x) == 18L)
#' 
#' ## Coercion       
#' # No fixed terminal or fixed anywhere coercion
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)"),
#'                    varmods = c("Carbamidomethyl (M)"))
#' 
#' stopifnot(length(x) == 2L)
#' 
#' x <- calc_aamasses(fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"), 
#'                    varmods = c("Acetyl (Protein N-term)", "TMT6plex (N-term)", 
#'                      "Oxidation (M)", "Carbamidomethyl (M)"))
#' 
#' stopifnot(length(x) == 12L)
#' 
#' # Fixed terminal coercion
#' x <- calc_aamasses(fixedmos = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)"),
#'                    varmods = c("Acetyl (Protein N-term)", "Oxidation (M)"))
#'                    
#' stopifnot(length(x) == 4L)
#' 
#' # Fixed anywhere coercion
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
#'                    varmods = c("Oxidation (M)"))
#'                    
#' stopifnot(length(x) == 3L)
#'                    
#' # Both fixed terminal and fixed anywhere coercion
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)", "Carbamidomethyl (M)"),
#'                    varmods = c("Acetyl (Protein N-term)", "Oxidation (M)"))
#' 
#' stopifnot(length(x) == 6L)
#' 
#' }
#' \dontrun{
#' # conflicts
#' x <- calc_aamasses(fixedmods = c("Carbamidomethyl (N-term)", "TMT2plex (N-term)"), 
#'                    varmods = NULL)
#'
#' # need separate S and T
#' x <- calc_aamasses(fixedmods = NULL, varmods = "Phospho (ST)")
#' }
#' @export
calc_aamasses <- function (fixedmods = c("TMT6plex (K)",
                                         "Carbamidomethyl (. = C)"),
                           varmods = c("TMT6plex (N-term)",
                                       "Acetyl (Protein N-term)",
                                       "Oxidation (M)",
                                       "Deamidated (N)",
                                       "Gln->pyro-Glu (N-term = Q)"),
                           aa_masses = NULL, 
                           varlabs = NULL, mod_motifs = NULL, 
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

  ## (0) Prep
  check_dupfvmods(fixedmods, varmods)
  
  fixedmods_orig <- fixedmods
  new_mods <- coerce_fvmods(fixedmods, varmods)
  fixedmods <- new_mods$fixedmods
  varmods <- new_mods$varmods
  f_to_v <- new_mods$f_to_v
  anywhere_coerce_sites <- new_mods$anywhere_coerce_sites
  nt_coerce_site <- new_mods$nt_coerce_site
  ct_coerce_site <- new_mods$ct_coerce_site
  rm(list = "new_mods")
  
  if (!is.null(mod_motifs)) check_mod_motifs(mod_motifs, c(fixedmods, varmods))
  fmod_motifs <- mod_motifs[names(mod_motifs) %in% fixedmods]
  vmod_motifs <- mod_motifs[names(mod_motifs) %in% varmods]
  
  if (is.null(aa_masses)) {
    aa_masses <- c(
      A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943,
      C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464,
      H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
      M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
      T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414,
      "N-term" = 1.007825, "C-term" = 17.002740,
      U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590,
      "-" = 0)
  }

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

  ## (1, 2) add fixed mods + NL, coerced fixed mods + NL
  aa_masses_fi <- add_fixed_masses(fixedmods_orig, aa_masses, fmod_motifs)
  
  aa_masses_fc <- if (is.null(f_to_v))
    aa_masses_fi
  else
    add_fixed_masses(fixedmods, aa_masses, fmod_motifs)
  
  attr(aa_masses_fc, "nt_coerce_site") <- nt_coerce_site
  attr(aa_masses_fc, "ct_coerce_site") <- ct_coerce_site
  attr(aa_masses_fc, "anywhere_coerce_sites") <- anywhere_coerce_sites

  ## (3) add variable mods + NL
  varmods_comb <- find_aamasses_vmodscombi(varmods, f_to_v, anywhere_coerce_sites)

  aa_masses_var <- lapply(varmods_comb, add_var_masses, aa_masses = aa_masses_fc, 
                          varlabs = varlabs, mod_motifs = vmod_motifs, 
                          anywhere_coerce_sites = anywhere_coerce_sites, 
                          nt_coerce_site = nt_coerce_site, 
                          ct_coerce_site = ct_coerce_site)

  if (length(aa_masses_var) >= maxn_vmods_setscombi) {
    warning("The ways of combinatorial variable modifications are ",
            length(aa_masses_var), ".\n",
            "Dropping combinations at indexes greater than `maxn_vmods_setscombi = ", 
            maxn_vmods_setscombi, "`.",
            call. = FALSE)
    
    aa_masses_var <- aa_masses_var[1:maxn_vmods_setscombi]
  }
  
  aa_masses_ms1 <- lapply(c(list(aa_masses_fc), aa_masses_var), parse_aamasses)
  
  aa_masses_all <- lapply(c(list(aa_masses_fi), aa_masses_var), finalize_aamasses, 
                          aa_masses = aa_masses, varlabs = varlabs, 
                          mod_motifs = mod_motifs)

  if (!is.null(out_path)) {
    save_mod_indexes(out_path, fixedmods, varmods, f_to_v)
    qs::qsave(aa_masses_ms1, file.path(out_path, "aa_masses_ms1.rds"), preset = "fast")
    qs::qsave(aa_masses_all, file.path(out_path, "aa_masses_all.rds"), preset = "fast")
  }
  
  invisible(aa_masses_all)
}


#' Finalizes \code{aa_masses_all}
#'
#' Replaces interim fixed and variable modifications with the finals. Results in
#' correct attributes in aa_masses_i. 
#'
#' @param aa_masses_i The i-th aa_masses_all.
#' @param aa_masses The original look-ups of AA masses.
#' @inheritParams matchMS
finalize_aamasses <- function (aa_masses_i, aa_masses, varlabs = NULL, mod_motifs = NULL)
{
  fixedmods <- names(attr(aa_masses_i, "fmods_ps", exact = TRUE))
  varmods <- names(attr(aa_masses_i, "vmods_ps", exact = TRUE))
  
  anywhere_excepts <- attr(aa_masses_i, "anywhere_excepts", exact = TRUE)
  nt_except <- attr(aa_masses_i, "nt_except", exact = TRUE)
  ct_except <- attr(aa_masses_i, "ct_except", exact = TRUE)
  excepts <- c(anywhere_excepts, nt_except, ct_except)
  cmods <- names(excepts)
  
  # coerced mods
  if (length(cmods)) {
    fixedmods <- c(fixedmods, cmods)
    varmods <- varmods[!varmods %in% cmods]
  }
  
  fmod_motifs <- mod_motifs[names(mod_motifs) %in% fixedmods]
  vmod_motifs <- mod_motifs[names(mod_motifs) %in% varmods]
  
  # may skip the mass updates if no coerced mods
  aa_masses_i <- add_fixed_masses(mods = fixedmods, aa_masses = aa_masses, 
                                  mod_motifs = fmod_motifs)
  
  aa_masses_i <- add_var_masses(mods = varmods, aa_masses = aa_masses_i, 
                                varlabs = varlabs, mod_motifs = vmod_motifs)
  
  aa_masses_i <- parse_aamasses(aa_masses_i)
}



#' Saves mod_indexes.txt
#' 
#' @param out_path An output path
#' @param fixedmods Fixed modifications
#' @param varmods Variable modifications
#' @param f_to_v Coerced fixed to variable modifications
save_mod_indexes <- function (out_path = NULL, fixedmods, varmods, 
                              f_to_v)
{
  if (is.null(out_path))
    return(NULL)
  
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
}


#' Checks duplicated modifications between fixedmods and varmods
#' 
#' @inheritParams calc_aamasses
check_dupfvmods <- function (fixedmods, varmods)
{
  dup_mods <- intersect(fixedmods, varmods)
  
  if (length(dup_mods)) 
    stop("Modifications cannot be simultaneously 'fixed' and 'variable': \n",
         paste(dup_mods, collapse = ", "), "\n", 
         "Hint: the default \"fixedmods\" and \"varmods\" are not NULL.", 
         call. = FALSE)
}


#' Coerces conflicting fixedmods to varmods
#' 
#' @inheritParams calc_aamasses
coerce_fvmods <- function (fixedmods, varmods)
{
  fmods_ps <- find_modps(fixedmods)
  vmods_ps <- find_modps(varmods)
  
  local({
    dup_fsites <- fmods_ps[duplicated(fmods_ps)]
    
    if (length(dup_fsites)) {
      stop("Multiple fixed modifications to the same site: \n",
           "'", paste(dup_fsites, collapse = ", "), "'")
    }
  })
  
  coerce_sites <- find_f_to_v(fixedmods, fmods_ps, vmods_ps)
  fV_coercion <- (length(coerce_sites) > 0L)

  if (fV_coercion) {
    f_to_v <- names(coerce_sites)
    varmods <- c(f_to_v, varmods)
    fixedmods <- fixedmods[!fixedmods %in% f_to_v]
    
    warning("Coerce '",
            paste(f_to_v, collapse = ", "), "'",
            " to conditional variable modifications.",
            call. = FALSE)

    oks_a <- !grepl("[NC]-term", f_to_v)
    oks_n <- grepl("N-term", f_to_v)
    oks_c <- grepl("C-term", f_to_v)
    
    # ok since the parity between coerce_sites and f_to_v
    anywhere_coerce_sites <- coerce_sites[oks_a]
    nt_coerce_site <- coerce_sites[oks_n]
    ct_coerce_site <- coerce_sites[oks_c]
    if (!length(anywhere_coerce_sites)) anywhere_coerce_sites <- NULL
    if (!length(nt_coerce_site)) nt_coerce_site <- NULL
    if (!length(ct_coerce_site)) ct_coerce_site <- NULL
  } 
  else {
    fixedmods <- fixedmods
    varmods <- varmods
    f_to_v <- NULL
    anywhere_coerce_sites <- NULL
    nt_coerce_site <- NULL
    ct_coerce_site <- NULL
  }

  default_mods <- c("initiator methionine from protein N-terminus")
  
  if (any(grepl(default_mods, c(fixedmods, varmods)))) 
    warning("Modifications defaulted and no need to specify: `\n",
            default_mods, "`.\n",
            call. = FALSE)
  
  list(fixedmods = fixedmods, 
       varmods = varmods, 
       f_to_v = f_to_v, 
       anywhere_coerce_sites = anywhere_coerce_sites, 
       nt_coerce_site = nt_coerce_site, 
       ct_coerce_site = ct_coerce_site)
}


#' Helper of finding the coercion sites.
#'
#' @param fmods_ps Positions and sites of fixed modifications.
#' @param vmods_ps Positions and sites of variable modifications.
#' @inheritParams calc_aamasses
#' @return A named vector, e.g., \code{c("TMT6plex (K) = K", "TMT6plex
#'   (N-term)" = "N-term")}
find_f_to_v <- function (fixedmods, fmods_ps, vmods_ps)
{
  f_nms <- names(fmods_ps)
  v_nms <- names(vmods_ps)
  
  ok_ft <- grepl("[NC]-term", f_nms)
  ok_vt <- grepl("[NC]-term", v_nms)
  fmods_ps_any <- fmods_ps[!ok_ft]
  vmods_ps_any <- vmods_ps[!ok_vt]
  fmods_ps_term <- fmods_ps[ok_ft]
  vmods_ps_term <- vmods_ps[ok_vt]

  coerce_asites <- fmods_ps_any[fmods_ps_any %in% vmods_ps_any]
  coerce_tsites <- fmods_ps_term[fmods_ps_term %in% vmods_ps_term]
  
  # e.g. "N-term" can be matched by both site and position
  # (no guarantee in the order of coerce_sites; so match names one at a time)
  coerce_sites <- unique(c(coerce_asites, coerce_tsites)) %>% 
    lapply(function (x) {
      names(x) <- fixedmods[fmods_ps == x]
      x
    })
  
  unlist(coerce_sites)
}

#' Checks mod_motifs
#' 
#' @param mods A concatenated list of fixed and variable modifications.
#' @inheritParams calc_aamasses
check_mod_motifs <- function (mod_motifs, mods)
{
  bads <- mod_motifs[!names(mod_motifs) %in% mods]
  
  if (length(bads))
    stop("\"mod_motifs\" not found in \"varmods\" or \"fixedmods\": ", 
         paste(bads, collapse = ", "))
}


#' Finds the combination of varmods
#' 
#' @param f_to_v Coerced fixed to variable modifications
#' @param anywhere_coerce_sites \code{Anywhere} sites coerced from fixed to
#'   variable modifications.
#' @inheritParams calc_aamasses
find_aamasses_vmodscombi <- function (varmods = NULL, f_to_v = NULL, 
                                      anywhere_coerce_sites = NULL) 
{
  if (is.null(varmods))
    return(NULL)
  
  varmods_comb <- lapply(seq_along(varmods), function (x) sim_combn(varmods, x)) 
  varmods_comb <- unlist(varmods_comb, recursive = FALSE)
  vmods_ps <- find_modps(varmods)

  vmods_ps_combi <- seq_along(vmods_ps) %>%
    lapply(function (x) sim_combn(vmods_ps, x)) %>%
    purrr::flatten()

  ## Remove the combinations without anywhere_coerce_sites
  # [x] e.g. "TMT6plex (K)" coerced from fixedmod to varmod, 
  #     the combination after coercion must contain "K"
  # [x] in the (less frequent) case of multiple varmod coercions, 
  #     apply a strong condition of "all" sites
  
  # if no K -> need everything else: TMT6plex (N-term), Acetyl (K), Deamidated (N)
  
  len_a <- length(anywhere_coerce_sites)
  
  if (len_a > 1L) {
    warning("Multiple coercions from fixed to variable modifications: ", 
            paste(anywhere_coerce_sites, collapse = ", "), "\n", 
            "Suggest change (some of) them to variable modifications.")
  }

  if (len_a) {
    ok_coerced_sites <- 
      lapply(vmods_ps_combi, function (v) all(anywhere_coerce_sites %in% unlist(v)))
    ok_coerced_sites <- unlist(ok_coerced_sites)
  }
  else {
    ok_coerced_sites <- NULL
  }

  ## Remove the combinations with multiple terminal mods
  #  by names: Gln->pyro Glu (N-term = Q) and Acetyl (Protein N-term = N-term)
  #  have different sites but 'N-term' in names; and also need to be excluded
  dup_terms <-
    purrr::map_lgl(vmods_ps_combi, ~ sum(grepl("N-term", names(.x))) >= 2L) |
    purrr::map_lgl(vmods_ps_combi, ~ sum(grepl("C-term", names(.x))) >= 2L)
  
  # Allow entries with different Anywhere mods to the same site
  #   (1) dHex(1)Hex(1) (S) and (2) Phospho (S)
  #   VS(1)S(2)ALSPSK
  
  varmods_comb <- if (is.null(ok_coerced_sites))
    varmods_comb[!dup_terms]
  else 
    varmods_comb[!dup_terms & ok_coerced_sites]

  ## Remove the combinations without terminal mods
  # e.g., if `Fixed Anywhere N-term` by users 
  #   -> must have `N-term` in any realization of varmods
  fixednt_coercion <- any(grepl("N-term", f_to_v))
  fixedct_coercion <- any(grepl("C-term", f_to_v))

  if (fixednt_coercion && length(varmods_comb)) {
    ok_nts <- lapply(varmods_comb, function (x) any(grepl("N-term", x)))
    varmods_comb <- varmods_comb[unlist(ok_nts)]
    rm(list = c("ok_nts"))
  }
  
  if (fixedct_coercion && length(varmods_comb)) {
    ok_cts <- lapply(varmods_comb, function (x) any(grepl("C-term", x)))
    varmods_comb <- varmods_comb[unlist(ok_cts)]
    rm(list = c("ok_cts"))
  }

  if (length(f_to_v)) {
    is_base <- lapply(varmods_comb, function (x) all(x %in% f_to_v))
    is_base <- unlist(is_base)
    varmods_comb <- varmods_comb[!is_base]
  }
    
  varmods_comb
}


#' Helper to add modification masses to amino-acid residues.
#'
#' It adds the masses of fixed, variable, and neutral-loss modifications to
#' amino-acid residues.
#'
#' @param mods A list of modifications.
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param anywhere_f_to_v Anywhere modifications coerced from fixed to variable;
#'   for example, \code{TMT6plex (K)}.
#' @param anywhere_coerce_sites The sites of coerced Anywhere modifications; for
#'   example, \code{K}.
#' @inheritParams matchMS
#' @return Lists of of amino-acid residues with modified mono-isotopic masses
#'   being incorporated. Returns NULL if \code{is.null(varmods_comb)}.
add_var_masses <- function (mods, aa_masses, varlabs = NULL, mod_motifs = NULL, 
                            anywhere_coerce_sites = NULL, nt_coerce_site = NULL, 
                            ct_coerce_site = NULL) 
{
  mod_type <- "vmods"
  all_mods <- if (length(mods)) paste(mods, collapse = ", ") else ""
  
  res <- extract_umods(mods)
  mod_masses <- lapply(res, `[[`, "monomass")
  positions_sites <- lapply(res, `[[`, "position_site")
  neulosses <- lapply(res, `[[`, "nl")
  positions <- unlist(lapply(positions_sites, names), use.names = FALSE)
  sites <- unlist(lapply(positions_sites, find_aa_site), use.names = FALSE)
  rm(list = c("res"))

  if (!is.null(varlabs)) {
    # subsets positions_sites2 by positions_sites
    res2 <- extract_umods(varlabs)
    mod_masses2 <- lapply(res2, `[[`, "monomass")
    positions_sites2 <- lapply(res2, `[[`, "position_site")
    neulosses2 <- lapply(res2, `[[`, "nl")
    positions2 <- unlist(lapply(positions_sites2, names), use.names = FALSE)
    sites2 <- unlist(lapply(positions_sites2, find_aa_site), use.names = FALSE)
    rm(list = "res2")

    oks2 <- (positions2 %in% positions) & (sites2 %in% sites)
    positions_sites2 <- positions_sites2[oks2]
    mod_masses2 <- mod_masses2[oks2]
    neulosses2 <- neulosses2[oks2]
    positions2 <- positions2[oks2]
    sites2 <- sites2[oks2]
    rm(list = c("oks2"))

    # to ensure the same order
    ps <- paste(positions, sites, sep = ".")
    ps2 <- paste(positions2, sites2, sep = ".")
    idxes <- match(ps2, ps)
    rm(list = c("ps", "ps2"))
    
    # updates mass deltas and NLs
    mod_masses[idxes] <- mod_masses[idxes] %+% mod_masses2

    if (length(neulosses2))
      neulosses[idxes] <- mapply(function (x, y) x + y, neulosses[idxes], neulosses2, 
                                 USE.NAMES = TRUE, SIMPLIFY = FALSE)
    
    rm(list = c("mod_masses2", "positions_sites2", "neulosses2", "idxes", 
                "positions2", "sites2"))
  }

  anywhere_excepts <- find_except_sites(anywhere_coerce_sites, positions_sites)
  nt_except <- find_except_sites(nt_coerce_site, positions_sites)
  ct_except <- find_except_sites(ct_coerce_site, positions_sites)

  for (i in seq_along(positions_sites))
    aa_masses[names(positions_sites[i])] <- mod_masses[[i]]

  attr(aa_masses, mod_type) <- all_mods
  attr(aa_masses, paste0(mod_type, "_ps")) <- positions_sites
  attr(aa_masses, paste0(mod_type, "_mass")) <- mod_masses
  aa_masses <- add_aamasses_motifs(aa_masses, mod_motifs, positions_sites)
  aa_masses <- add_aamasses_neulosses(aa_masses, neulosses, mod_type)
  
  if (is.null(attr(aa_masses, "vmods"))) attr(aa_masses, "vmods") <- ""
  if (is.null(attr(aa_masses, "vmods_neuloss"))) attr(aa_masses, "vmods_neuloss") <- ""
  if (is.null(attr(aa_masses, "vmods_mass"))) attr(aa_masses, "vmods_mass") <- 0
  attr(aa_masses, "anywhere_excepts") <- anywhere_excepts
  attr(aa_masses, "nt_except") <- nt_except
  attr(aa_masses, "ct_except") <- ct_except

  ### 
  # Variable mods: need both the "original" and the "delta" forms
  # (in combinatorial some are modified and some are not) 
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

  aa_masses
}


#' Adds fixed masses
#' 
#' @inheritParams add_var_masses
add_fixed_masses <- function (mods, aa_masses, mod_motifs = NULL) 
{
  mod_type <- "fmods"
  all_mods <- if (length(mods)) paste(mods, collapse = ", ") else ""
  
  res <- extract_umods(mods)
  mod_masses <- lapply(res, `[[`, "monomass")
  positions_sites <- lapply(res, `[[`, "position_site")
  neulosses <- lapply(res, `[[`, "nl")
  rm(list = c("res"))
  
  check_fmods_pos_site(positions_sites)
  
  for (i in seq_along(positions_sites)) {
    p <- positions_sites[[i]]
    s <- find_aa_site(p)
    aa_masses[s] <- aa_masses[s] + mod_masses[[i]]
  }

  attr(aa_masses, mod_type) <- all_mods
  attr(aa_masses, paste0(mod_type, "_ps")) <- positions_sites
  attr(aa_masses, paste0(mod_type, "_mass")) <- mod_masses
  aa_masses <- add_aamasses_motifs(aa_masses, mod_motifs, positions_sites)
  aa_masses <- add_aamasses_neulosses(aa_masses, neulosses, mod_type)
  
  if (is.null(attr(aa_masses, "fmods"))) attr(aa_masses, "fmods") <- ""
  if (is.null(attr(aa_masses, "fmods_neuloss"))) attr(aa_masses, "fmods_neuloss") <- ""
  if (is.null(attr(aa_masses, "fmods_mass"))) attr(aa_masses, "fmods_mass") <- 0
  if (is.null(attr(aa_masses, "vmods"))) attr(aa_masses, "vmods") <- ""
  if (is.null(attr(aa_masses, "vmods_neuloss"))) attr(aa_masses, "vmods_neuloss") <- ""
  if (is.null(attr(aa_masses, "vmods_ps"))) attr(aa_masses, "vmods_ps") <- ""
  if (is.null(attr(aa_masses, "vmods_mass"))) attr(aa_masses, "vmods_mass") <- 0

  aa_masses
}


#' Finds the coerced fixed sites that can be coerced back to fixed sites.
#'
#' @param positions_sites Lists of positions and sites.
#' @param fixed_sites Sites coerced from fixed to variable modifications.
#' @examples
#' \donttest{
#' mods <- c("TMT6plex (K)", "Carbamyl (M)", "Oxidation (M)", "Gln->pyro-Glu (N-term = Q)")
#' res <- extract_umods(mods)
#' positions_sites <- lapply(res, `[[`, "position_site")
#' fixed_sites <- c("TMT6plex (K)" = "K", "Carbamyl (M)" = "M")
#' }
find_except_sites <- function (fixed_sites, positions_sites)
{
  fixed_mods <- names(fixed_sites)
  csites <- fixed_sites[fixed_sites %in% positions_sites]
  
  except_sites <- lapply(csites, function (csite) {
    vmods_at_csite <- names(positions_sites[positions_sites == csite])
    other_csites <- setdiff(vmods_at_csite, fixed_mods)
    if (length(other_csites)) NULL else csite
  })
  
  except_sites <- unlist(except_sites, recursive = FALSE)
  except_vmods <- fixed_mods[fixed_sites %in% except_sites]
  
  if (!is.null(except_sites)) 
    names(except_sites) <- except_vmods
  
  invisible(except_sites)
}


#' Finds the positions and sites of modifications.
#'
#' @param mods A vector of modifications, e.g., \code{c("TMT6plex (K)", "Acetyl
#'   (K)")}
find_modps <- function (mods)
{
  ans <- lapply(mods, find_unimod) 
  ps <- lapply(ans, `[[`, "position_site")
  purrr::flatten(ps)
}


#' Extracts Unimod results.
#' 
#' @param mods A vector of modifications. 
extract_umods <- function (mods)
{
  res <- lapply(mods, find_unimod)
  names(res) <- mods
  check_resunimod(res)
}


#' Checks the results from find_unimod
#' 
#' Shows warning if the same \code{site} with different fixedmods.
#' 
#' @param res Results from \link{find_unimod}
check_resunimod <- function (res)
{
  if (length(res)) {
    x <- res[[1]]
    
    nm_seqs <- c("title", "monomass", "position_site", "nl")
    ok <- identical(names(x), nm_seqs)
    
    if (!ok)
      stop("The structures from `find_unimod` is not in the order of: ", 
           paste(nm_seqs, collapse = ", "))
  }
  
  invisible(res)
}


#' Checks the positions_sites in fixedmods
#' 
#' @param positions_sites Lists of positions and sites.
check_fmods_pos_site <- function (positions_sites)
{
  if (length(positions_sites) > 1L) {
    dups <- purrr::reduce(positions_sites, `c`) %>%
      .[duplicated(.)]
    
    if (length(dups)) {
      dups_in_each <- lapply(positions_sites, function (x) x == dups)
      dup_mods <- names(positions_sites[unlist(dups_in_each)])
      
      warning("Conflicts in fixed modifications: \n",
              paste(dup_mods, collapse = ", "), "\n",
              "May consider change from fixed to variable modifications(s); \n",
              "or create a new Unimod for joint modifications.",
              call. = FALSE)
    }
  }
}


#' Adds the attribute of neuloss.
#' 
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param neulosses Neutral losses
#' @param mod_type The type of modification
add_aamasses_neulosses <- function (aa_masses, neulosses, mod_type) 
{
  no_nls <- neulosses %>%
    purrr::map_lgl(~ all(.x == 0)) %>%
    all()
  
  # if (no_nls) return(list(aa_masses))
  if (no_nls) return(aa_masses)
  
  attr(aa_masses, paste0(mod_type, "_neuloss")) <- neulosses
  
  aa_masses
}


#' Adds the motifs attribute to aa_masses
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param positions_sites Named list of positions (in names) and sites (in
#'   values).
#' @inheritParams matchMS
add_aamasses_motifs <- function (aa_masses, mod_motifs, positions_sites)
{
  # Need to update subset_by_prps and subset_anysite if changed from NULL to ""
  nms <- names(positions_sites)
  pmod_motifs <- mod_motifs[nms]
  
  # pmod_motifs <- lapply(pmod_motifs, function (x) if (is.null(x)) "" else x)
  # names(pmod_motifs) <- nms
  if (!is.null(pmod_motifs)) 
    names(pmod_motifs) <- nms
  
  attr(aa_masses, "mod_motifs") <- pmod_motifs
  
  aa_masses
}


#' Parses \code{aa_masses}.
#'
#' @inheritParams add_var_masses
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
  
  min_n_res <- local({
    anywhere_backs <- attr(aa_masses, "anywhere_backs", exact = TRUE)
    amods2 <- if (is.null(anywhere_backs)) amods else amods[!names(amods) %in% anywhere_backs]
    count_elements(unlist(amods2, recursive = FALSE, use.names = FALSE))
    # count_elements(unlist(amods, recursive = FALSE, use.names = FALSE))
  })

  # Is "same Anywhere mod existed"
  is_same <- any(length(min_n_res) > 1L)

  # `TMT6plex (N-term)` and `Amidated (Protein C-term)`
  tmod <- vmods_ps[!vmods_ps %in% amods]
  tmod <- if (!length(tmod)) NULL else if (all(tmod == "")) NULL else tmod

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

  # length(ftmod) > 1L: `TMT6plex (N-term)` and `K8 (C-term)`
  ftmod <- fmods_ps %>% .[! . %in% famods]
  ftmod <- if (!length(ftmod)) NULL else if (all(ftmod == "")) NULL
  
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
                                   ftmass = 18.010565) 
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
                                 ftmass = ftmass) 
  peps <- purrr::flatten(peps)
  
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
                             aa_masses = NULL, ftmass = 18.010565) 
{
  lapply(fasta_db, make_noenzpeps, min_len = min_len, max_len = max_len, 
         aa_masses = aa_masses, ftmass = ftmass)
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
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods, varmods = varmods)
#' aa_masses <- aa_masses_all[[1]]
#' 
#' ans <- make_noenzpeps(prot, 7, 40, aa_masses)
#' 
#' # short FASTA (both N-term and C-term on a sequence)
#' prot <- substring(prot, 1, 10)
#' ans <- make_noenzpeps(prot, 7, 40, aa_masses)
#' }
make_noenzpeps <- function (prot = NULL, min_len = 7L, max_len = 40L, 
                            aa_masses = NULL, ftmass = 18.010565) 
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
                aa_masses, ftmass)

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
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods, varmods = varmods)
#' aa_masses <- aa_masses_all[[1]]
#' 
#' aas <- LETTERS[LETTERS %in% names(aa_masses)]
#' prot <- paste0(aas, collapse = "")
#' len <- nchar(prot)
#' masses <- hmake_noenzpeps(1, prot, 7, 40, len, aa_masses)
#' }
hmake_noenzpeps <- function (start = 1L, prot = NULL, min_len = 7L, max_len = 40L, 
                             len = NULL, aa_masses = NULL, ftmass = 18.010565) 
{
  end_fi <- min_len + start  - 1L
  end_la <- min(max_len + start - 1L, len)
  ends <- end_fi:end_la
  
  peps <- substring(prot, start, ends)
  masses <- ms1masses_bare_noenz(peps, aa_masses, ftmass)

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
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods, varmods = varmods)
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
ms1masses_bare_noenz <- function (x, aa_masses, ftmass = 18.010565) 
{
  len <- length(x)
  aas <- .Internal(strsplit(x[len], "", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
  len_a <- length(aas)
  
  masses <- cumsum(aa_masses[aas])
  len_m <- length(masses)
  masses <- masses[(len_m - len + 1L):len_m] + ftmass
  # masses <- tail(masses, len) + ftmass
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


#' Adjusts the masses for the all-fixed mode
#' 
#' The \code{fwd_peps_1} was at first based on the coerced aa_masses_0;
#'
#' @param fwd_peps_1 The first list of forward peptides with masses.
#' @param aa_masses_0 The amino-acid masses look-up table with the coercion of
#'   fixed to variable modification.
#' @param aa_masses_1 The amino-acid masses look-up table without the coercion
#'   of fixed to variable modification.
#' @inheritParams matchMS
adj_base_masses <- function (fwd_peps_1, aa_masses_0, aa_masses_1, 
                             min_mass = 700L, max_mass = 4500L, digits = 4L)
{
  nt_coerce_site <- attr(aa_masses_0, "nt_coerce_site", exact = TRUE)
  ct_coerce_site <- attr(aa_masses_0, "ct_coerce_site", exact = TRUE)
  anywhere_coerce_sites <- attr(aa_masses_0, "anywhere_coerce_sites", exact = TRUE)
  
  if (length(nt_coerce_site) > 1L)
    stop("Cannot have more than one N-term site coercion.")
  
  if (length(ct_coerce_site) > 1L)
    stop("Cannot have more than one C-term site coercion.")
  
  ok <- is.null(nt_coerce_site) && 
    is.null(ct_coerce_site) && 
    is.null(anywhere_coerce_sites)
  
  if (ok)
    return(fwd_peps_1[fwd_peps_1 <= max_mass])

  delta_nt <- if (length(nt_coerce_site))
    aa_masses_1[[nt_coerce_site]] - aa_masses_0[[nt_coerce_site]]
  else
    0
  
  delta_ct <- if (length(ct_coerce_site))
    aa_masses_1[[ct_coerce_site]] - aa_masses_0[[ct_coerce_site]]
  else 
    0
  
  len_a <- length(anywhere_coerce_sites)

  if (len_a) {
    ds <- lapply(anywhere_coerce_sites, function (s) {
      aa_masses_1[s] - aa_masses_0[s]
    })
    
    ns <- lapply(anywhere_coerce_sites, function (s) {
      .Call(stringi:::C_stri_count_fixed, str = names(fwd_peps_1), 
            pattern = s, opts_fixed = NULL)
    })
    
    dns <- mapply(function (d, n) {
      d * n
    }, ds, ns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    delta_anywhere <- dns[[1]]
    
    if (len_a > 1L)
      for (i in 2:len_a) delta_anywhere <- delta_anywhere + dns[[i]]
  }
  else {
    delta_anywhere <- 0
  }

  delta <- delta_nt + delta_ct + delta_anywhere
  fwd_peps_1 <- fwd_peps_1 + delta
  fwd_peps_1 <- fwd_peps_1[fwd_peps_1 >= min_mass & fwd_peps_1 <= max_mass]
  
  round(fwd_peps_1, digits = digits)
}


#' Adjusts the masses for the coerced sites
#' 
#' @param peps The list of forward peptides with masses.
#' @param aa_masses An amino-acid masses look-up table.
adj_anywhere_masses <- function (peps, aa_masses)
{
  anywhere_excepts <- attr(aa_masses, "anywhere_excepts", exact = TRUE)
  len_a <- length(anywhere_excepts)
  
  if (!len_a)
    return(peps)
  
  ds <- aa_masses[names(anywhere_excepts)]
  
  ns <- lapply(anywhere_excepts, function (s) {
    .Call(stringi:::C_stri_count_fixed, str = names(peps), 
          pattern = s, opts_fixed = NULL)
  })
  
  dns <- mapply(function (d, n) {
    d * n
  }, ds, ns, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  for (i in seq_along(dns)) {
    peps <- peps + dns[[i]]
  }
  
  invisible(peps)
}


#' Helper in calculating peptide masses.
#'
#' (2) "amods- tmod+ vnl- fnl-".
#' 
#' No needs of checking \code{is_empty(ntmod) && is_empty(ctmod)}
#' 
#' @param peps A list of peptide sequences.
#' @inheritParams add_var_masses
#' @inheritParams distri_peps
#' @inheritParams matchMS
add_term_mass <- function (peps, aa_masses, min_mass = 700L, max_mass = 4500L) 
{
  type <- attr(aa_masses, "type", exact = TRUE)
  
  if (!grepl("tmod+", type, fixed = TRUE))
    return(peps)
  
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  delta <- if (length(ntmod) && length(ctmod)) 
    aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  else if (length(ntmod)) 
    aa_masses[names(ntmod)]
  else if (length(ctmod)) 
    aa_masses[names(ctmod)]

  ans <- peps + delta
  
  ans[ans >= min_mass & ans <= max_mass]
}


#' Calculates mono-isotopic masses of peptide sequences.
#'
#' @param seqs Sequences of peptides from FASTAs by protein accessions. Each
#'   list contains two lists of sequences: (1) without and (2) with N-terminal
#'   methionine.
#' @param ftmass The sum of masses of \code{fixed} N-term and C-term
#'   modifications.
#' @inheritParams calc_pepmasses2
#' @inheritParams add_var_masses
#' @inheritParams distri_fpeps
ms1masses_bare <- function (seqs = NULL, aa_masses = NULL, ftmass = NULL,
                            max_miss = 2L, min_len = 7L, max_len = 40L,
                            maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                            is_fixed_protnt = FALSE, is_fixed_protct = FALSE) 
{
  # (1) before rolling sum (not yet terminal H2O)
  # (1.1) without N-term methionine
  data_1 <- lapply(seqs, `[[`, 1)
  data_1 <- ms1masses_noterm(data_1, aa_masses = aa_masses,
                             maxn_vmods_per_pep = maxn_vmods_per_pep,
                             maxn_sites_per_vmod = maxn_sites_per_vmod)
  data_1 <- attr(data_1, "data")

  # (1.2) with N-term methionine
  data_2 <- lapply(seqs, `[[`, 2)
  data_2 <- ms1masses_noterm(data_2, aa_masses = aa_masses,
                             maxn_vmods_per_pep = maxn_vmods_per_pep,
                             maxn_sites_per_vmod = maxn_sites_per_vmod)
  data_2 <- attr(data_2, "data")

  # (2) rolling sum (not yet terminal masses, e.g, H2O)
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

  # (3) putting together (+ terminal masses)
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
#' @inheritParams add_var_masses
#' @inheritParams distri_peps
ms1masses_noterm <- function (aa_seqs, aa_masses, maxn_vmods_per_pep = 5L,
                              maxn_sites_per_vmod = 3L) 
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
                                maxn_sites_per_vmod = maxn_sites_per_vmod) 
  
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
                                maxn_sites_per_vmod = 3L) 
{
  lapply(aa_seqs, function (x) 
    calcms1mass_noterm_byprot(prot_peps = x,
                              aa_masses = aa_masses,
                              maxn_vmods_per_pep = maxn_vmods_per_pep,
                              maxn_sites_per_vmod = maxn_sites_per_vmod))
}


#' Helper of \link{calcms1mass_noterm}.
#'
#' For single protein.
#'
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams ms1masses_noterm
calcms1mass_noterm_byprot <- function (prot_peps, aa_masses, 
                                       maxn_vmods_per_pep = 5L,
                                       maxn_sites_per_vmod = 3L) 
{
  # by peptides
  ans <- lapply(prot_peps, calcms1mass_noterm_bypep,
                aa_masses = aa_masses,
                maxn_vmods_per_pep = maxn_vmods_per_pep,
                maxn_sites_per_vmod = maxn_sites_per_vmod)
  
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
                                      maxn_sites_per_vmod = 3L) 
{
  if (is.na(aa_seq)) return(NULL)
  
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))

  aas <- aa_masses[aas]
  aas <- sum(aas)
  names(aas) <- aa_seq

  invisible(aas)
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
#' @inheritParams add_var_masses
hsemipeps_byprots <- function (prots, min_len = 7L, max_len = 40L, aa_masses) 
{
  lapply(prots, semipeps_byprots, min_len, max_len, aa_masses)
}


#' Finds the semi-enzymatic peptides for a proteins.
#'
#' Redundancy such as peptides from N-term methionine cleavage or semi-tryptic
#' ladders are handled.
#'
#' @param vals A list of full-enzymatic peptides under a protein.
#' @inheritParams matchMS
#' @inheritParams add_var_masses
#' @return A vector of full- and semi-enzymatic peptides. Sequences in names and
#'   masses in values.
semipeps_byprots <- function (vals, min_len = 7L, max_len = 40L, aa_masses) 
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
                                     max_len = max_len, 
                                     aa_masses = aa_masses, 
                                     ct_offset = 0L), 
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)
  semis_ot <- .Internal(unlist(semis_ot, recursive = FALSE, use.names = TRUE))
  
  semis_ct <- mapply(calc_semipepmasses, vals_ct, peps_ct, 
                     MoreArgs = list(min_len = min_len, 
                                     max_len = max_len, 
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
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods, varmods = varmods)
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
calc_semipepmasses <- function (val, pep, min_len = 7L, max_len = 40L, aa_masses, 
                                ct_offset = 0L) 
{
  options(digits = 9L)
  
  len <- nchar(pep)
  len2 <- len - ct_offset
  
  if (len2 < min_len || len2 > max_len) 
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
#' @inheritParams add_var_masses
hms1_a0_vnl0_fnl1 <- function (masses, fnl_combi, aa_masses, 
                               min_mass = 700, max_mass = 4500L, digits = 4L) 
{
  mapply(ms1_a0_vnl0_fnl1, 
         masses, names(masses), 
         MoreArgs = list (
           fnl_combi = fnl_combi, 
           aa_masses = aa_masses,
           min_mass = min_mass, 
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
                              min_mass = 700, max_mass = 4500L,digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, 
                            useBytes = FALSE))
  aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))

  delta <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
  out <- round(mass - delta, digits = digits)
  
  out <- out[out >= min_mass & out <= max_mass]
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
#' @inheritParams add_var_masses
hms1_a1_vnl0_fnl0 <- function (masses, amods, aa_masses,
                               vmods_nl = NULL, fmods_nl = NULL,
                               maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L,
                               ms1vmods = NULL, min_mass = 700L, max_mass = 4500L, 
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
           min_mass = min_mass, 
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
                              min_mass = 700L, max_mass = 4500L, 
                              digits = 4L) 
{
  aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
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
  out <- out[out >= min_mass & out <= max_mass]
  names(out) <- rep(aa_seq, length(out))
  
  invisible(out)
}


