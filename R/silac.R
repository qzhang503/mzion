#' SILAC
#'
#' Searches against mixed SILAC groups (heave and light mixed into one sample).
#'
#' @param this_call An expression from match.call.
#' @param aa_masses An amino acid look-ups.
#' @inheritParams matchMS
matchMS_silac_mix <- function (silac_mix = list(base = NULL, heavy = c("K8 (K)", "R10 (R)")), 
                               this_call, out_path, mgf_path, aa_masses) 
{
  message("Searches against SILAC groups ", 
          "(heavy, light etc. mixed into one sample)")
  message("[x] For reprocessing (with new score function) ", 
          "delete cached 'pepscores_[...]' and 'calc_pepscores.rda'.\n")
  
  lapply(c("silac_mix", "mgf_path", "this_call", "out_path"), function (x) {
    if (is.null(x)) stop("`", x, "` cannot be NULL.")
  })
  
  if (!is.list(silac_mix)) {
    stop("Supply silac_mix groups as lists, e.g., \n\n", 
         "  # unlabelled base (see ?add_unimod for more)\n", 
         "  silac_mix = list(base = c(fixedlabs = NULL, varlabs = NULL), \n", 
         "                   grpC = c(fixedlabs = c(\"Label:13C(6) (R)\", ...), \n", 
         "                            varlabs   = c(\"Label:13C(2) (Protein N-term)\", ...), \n", 
         "                   grpN = c(fixedlabs = c(\"Label:15N(2) (K)\", ...), \n", 
         "                            varlabs   = c(\"Label:15N(-1) (N)\", ...), \n\n", 
         "  # labelled base\n", 
         "  silac_mix = list(base = c(fixedlabs = ..., varlabs = ...), \n", 
         "                   grpC = c(fixedlabs = ..., varlabs = ...), \n", 
         "                   grpN = c(fixedlabs = ..., varlabs = ...), \n\n")
  }
  
  nms <- names(silac_mix)
  
  if (!"base" %in% nms) {
    warning("Required `base` group not found; assume: base = NULL")
    nms <- c("base", nms)
    silac_mix <- c(list(base = NULL), silac_mix)
  }
  
  if (all(unlist(lapply(silac_mix, is.null))))
    stop("`silac_mix` groups cannot be all NULL.")
  
  if ((len <- length(nms)) < 2L)
    stop("Need at least two `silac_mix` groups.")
  
  out_paths <- vector("list", len)
  
  for (i in seq_len(len)) {
    sub_call <- this_call
    sub_nm <- nms[i]
    sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))

    if (file.exists(file.path(sub_path, "psmQ.txt"))) 
      next
    
    if (i > 1L) {
      if (file.exists(mgf_call <- file.path(out_paths[[1]], "Calls", "load_mgfs.rda"))) {
        sub_call_path <- create_dir(file.path(sub_path, "Calls"))
        file.copy(mgf_call, file.path(sub_call_path, "load_mgfs.rda"), overwrite = TRUE)
      }
    }
    
    sub_mods <- silac_mix[[sub_nm]]
    
    if (is.null(sub_mods)) {
      this_aa <- NULL
    }
    else {
      this_aa <- aa_masses

      # fixedlabs and varlabs
      lab_nms <- names(sub_mods)
      oks <- grepl("fixedlabs\\d*", lab_nms)
      fixedlabs <- unname(sub_mods[oks])
      varlabs <- unname(sub_mods[!oks])

      if (!length(fixedlabs)) fixedlabs <- NULL
      if (!length(varlabs)) varlabs <- NULL
      
      sub_call$fixedlabs <- fixedlabs
      sub_call$varlabs <- varlabs
      
      rm(list = c("lab_nms", "oks", "fixedlabs", "varlabs"))
    }
    
    rm(list = "sub_mods")
    
    sub_call$fixedmods <- eval(sub_call$fixedmods)
    sub_call$out_path <- sub_path
    sub_call$mgf_path <- mgf_path
    sub_call$bypass_silac_mix <- TRUE
    sub_call$silac_mix <- NULL
    sub_call$bypass_from_pepscores <- FALSE # flowthrough: silac + noenzyme
    sub_call$aa_masses <- this_aa

    df <- tryCatch(eval(sub_call), error = function (e) NULL)
    
    if (is.null(df)) 
      warning("No results at `silac_mix = ", sub_nm, "`.")
    
    rm(list = c("df"))
    gc()
  }
  
  message("Combine mixed SILAC results.")
  comine_PSMsubs(sub_paths = out_paths, groups = nms, out_path = out_path)
  gc()
  
  enzyme <- as.character(this_call[["enzyme"]])
  
  if (isTRUE(enzyme == "noenzyme"))
    return(NULL)
  
  message("Done (mixed SILAC search).")
  options(show.error.messages = FALSE)
  stop()
}


#' Searches by sets of parameters.
#'
#' @param grp_args The names of arguments in \code{par_groups}.
#' @param mgf_paths The paths to MGF (with group searches).
#' @param this_call An expression from match.call.
#' @inheritParams matchMS
matchMS_par_groups <- function (par_groups = NULL, grp_args = NULL, 
                                mgf_paths = NULL, this_call = NULL, 
                                out_path = NULL) 
{
  message("Multiple searches by parameter groups...")
  message("[x] For reprocessing (with new score function) ", 
          "delete cached 'pepscores_[...]' and 'calc_pepscores.rda'.\n", 
          "[x] To bypass cached results, use \"use_ms1_cache = FALSE\".\n")
  
  lapply(c("par_groups", "grp_args", "this_call", "out_path"), function (x) {
    if (is.null(x))
      stop("`", x, "` cannot be NULL.")
  })
  
  if (!is.list(par_groups)) {
    stop("Supply `par_groups` as list, e.g., \n\n", 
         "par_groups = list(\n", 
         "  list(mgf_path  = \"~/mzion/my_proj/mgfs/grp_1\"", ",\n", 
         "       fixedmods = c(\"Carbamidomethyl (C)\")", "),\n", 
         "  list(mgf_path  = \"~/mzion/my_proj/mgfs/grp_2\"", ",\n", 
         "       fixedmods = c(\"Carbamidomethyl (C)\", \"K8 (K)\", \"R10 (R)\")", ")\n", 
         "  )")
  }
  
  if (all(unlist(lapply(par_groups, is.null))))
    stop("`par_groups` cannot be all NULL.")
  
  nms <- names(par_groups)
  len <- length(nms)
  
  if (len < 2L)
    stop("Need at least two `par_groups`.")
  
  ans <- out_paths <- vector("list", len)
  
  for (i in seq_len(len)) {
    sub_call <- this_call
    sub_nm <- nms[i]
    sub_pars <- par_groups[[i]]
    sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))
    
    file_peploc <- file.path(sub_path, "temp", "peploc.rds")
    
    if (file.exists(file_peploc)) {
      ans[[i]] <- qs::qread(file_peploc)
      next
    }
    
    for (arg in grp_args) 
      sub_call[[arg]] <- sub_pars[[arg]]
    
    sub_call$out_path <- sub_path
    sub_call$bypass_par_groups <- TRUE
    sub_call$bypass_from_protacc <- TRUE
    sub_call$par_groups <- NULL
    
    local({
      mgf_path <- eval(sub_call$mgf_path)
      checkMGF(mgf_path, grp_args = NULL, error = "stop")
    })
    
    df <- tryCatch(eval(sub_call), error = function (e) NULL)
    
    if (!file.exists(file_peploc)) 
      stop("File not found: ", file_peploc)
    
    ans[[i]] <- if (is.null(df)) qs::qread(file_peploc) else df
    rm(list = "df")
    
    if (is.null(ans[[i]])) 
      warning("No results at `par_groups = ", sub_nm, "`.")
    
    gc()
  }
  
  # map `raw_file` and `scan_title` before data combination
  if (!is.null(mgf_paths)) {
    for (i in seq_along(ans)) {
      ans[[i]] <- map_raw_n_scan(ans[[i]], mgf_paths[[i]])
      ans[[i]]$pep_group <- nms[i]
    }
  }
  
  out <- dplyr::bind_rows(ans)
  rm(list = "ans")
  dir.create(file.path(out_path, "temp"), showWarnings = FALSE, recursive = FALSE)
  qs::qsave(out, file.path(out_path, "temp", "peploc.rds"), preset = "fast")
  qs::qsave(out_paths, file.path(out_path, "temp", "out_paths.rds"), preset = "fast")
  rm(list = "out")
  
  gc()
  
  this_call$bypass_pepmasses <- TRUE
  this_call$bypass_bin_ms1 <- TRUE
  this_call$bypass_mgf <- TRUE
  this_call$bypass_ms2match <- TRUE
  this_call$bypass_pepscores <- TRUE
  this_call$bypass_pepfdr <- TRUE
  this_call$bypass_peploc <- TRUE
  this_call$bypass_par_groups <- TRUE
  this_call$bypss_mgf_checks <- TRUE
  this_call$from_group_search <- TRUE
  this_call$par_groups <- NULL
  
  out <- tryCatch(eval(this_call), error = function (e) NULL)
  
  message("Done (group search).")
  
  invisible(out)
}


#' Adds fixedlab masses to aa_masses
#' 
#' @param aa_masses A look-up of amino-acid residue masses.
#' @inheritParams matchMS
add_fixedlab_masses <- function (fixedlabs, aa_masses)
{
  # not yet check the validity of fixedlabs
  # ...
  # not yet test multiple sites, e.g., Oxidation (M) and Carbamyl (M)
  # ...
  
  umods <- lapply(fixedlabs, find_unimod)
  monos <- unlist(lapply(umods, `[[`, "monomass"), use.names = FALSE)
  sites <- unlist(lapply(umods, `[[`, "position_site"), use.names = FALSE)
  names(monos) <- sites
  
  # e.g. the same site to both N and C labels
  us <- unique(sites)
  ds <- unlist(lapply(us, function (x) sum(monos[names(monos) == x])))
  names(ds) <- us
  
  idxes <- match(us, names(aa_masses))
  aa_masses[idxes] <- aa_masses[idxes] + ds
  
  aa_masses
}


#' Noenzyme search.
#'
#' @param this_call An expression from match.call.
#' @param silac_noenzyme Logical; is the search a combination of SILAC and
#'   noenzyme.
#' @param groups_noenzyme Logical; is the search a combination of group search
#'   and noenzyme.
#' @inheritParams matchMS
matchMS_noenzyme <- function (this_call = NULL, min_len = 7L, max_len = 40L, 
                              fasta = NULL, out_path = NULL, mgf_path = NULL, 
                              noenzyme_maxn = 0L, quant = "none", 
                              silac_noenzyme = FALSE, 
                              groups_noenzyme = FALSE) 
{
  if (groups_noenzyme)
    stop("Not yet support group searches with no enzyme specificity")
  
  message("Searches with no enzyme specificity...")
  message("[x] For reprocessing (with new score function) ", 
          "delete cached 'pepscores_[...]' and 'calc_pepscores.rda'.\n")
  
  size <- local({
    if (noenzyme_maxn) 
      return (noenzyme_maxn)
    
    mouse_fasta_size <- 11 
    fasta_size <- sum(unlist(lapply(fasta, file.size)))/1024^2
    
    # large RAM -> large `fct_mem`
    fct_mem <- local({
      mgf_files <- list.files(mgf_path, pattern = "\\.mgf$", full.names = TRUE)
      fct_mgf <- max(1, sum(unlist(lapply(mgf_files, file.size)))/1024^3)
      
      ans <- tryCatch(
        find_free_mem()/1024/fct_mgf,
        error = function (e) NA)
      
      if (is.na(ans))
        ans <- fct_mgf
      
      ans
    })
    
    # large fasta -> large `fct_fasta` (>= 1)
    fct_fasta <- max(1, fasta_size/mouse_fasta_size)
    
    # ^1.5, 0.6: 90% RAM aggressiveness with uniprot fasta human + mouse
    max(1L, floor(fct_mem/(fct_fasta^1.5) * .5))
  })
  
  len <- length(min_len:max_len)
  
  if (len > size) {
    message("Split data by peptide lengths. To overrule, may consider ", 
            "\"noenzyme_maxn = ", len + 1, "\".")

    if (size == 1L) {
      n_chunks <- len
      spans <- split(min_len:max_len, 1:len)
    }
    else {
      n_chunks <- ceiling(len/size)
      spans <- chunksplit(min_len:max_len, n_chunks, rightmost.closed = TRUE)
    }
    
    sub_nms <- out_paths <- vector("list", n_chunks)
    
    for (i in seq_len(n_chunks)) {
      sub_call <- this_call
      span <- spans[[i]]
      start <- span[1]
      end <- span[length(span)]
      sub_nm <- sub_nms[[i]] <- paste0("sub", i, "_", start, "_", end)
      sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))
      
      ok <- file.exists(file.path(sub_path, "psmQ.txt"))
      
      if (ok) next
      
      if (i > 1L) {
        mgf_call <- file.path(out_paths[[1]], "Calls", "load_mgfs.rda")
        
        if (file.exists(mgf_call)) {
          sub_call_path <- create_dir(file.path(sub_path, "Calls"))
          file.copy(mgf_call, file.path(sub_call_path, "load_mgfs.rda"), overwrite = TRUE)
        }
      }
      
      sub_call$min_len <- start
      sub_call$max_len <- end
      sub_call$out_path <- sub_path
      sub_call$mgf_path <- mgf_path
      sub_call$bypass_noenzyme <- TRUE
      sub_call$bypass_from_pepscores <- TRUE
      sub_call$silac_noenzyme <- silac_noenzyme # silac + noenzyme
      
      ans <- tryCatch(eval(sub_call), error = function (e) NULL)
      
      message("Completed `min_len = ", start, "` to `max_len = ", end, "`.")
      gc()
    }
    
    # not necessary for silac_noenzyme
    file.copy(file.path(out_paths[[1]], "Calls"), out_path, recursive = TRUE)
    combine_ion_matches(out_path, out_paths, type = "ion_matches_")
    combine_ion_matches(out_path, out_paths, type = "reporters_")

    this_call$bypass_noenzyme <- TRUE
    this_call$bypass_pepmasses <- TRUE
    this_call$bypass_bin_ms1 <- TRUE
    this_call$bypass_mgf <- TRUE
    this_call$bypass_ms2match <- TRUE
    
    # (a) nested silac + noenzyme: flow through to psmQ -> combine -> done
    # 
    # (b) noenzyme only: early return bypass_from_pepscores -> next length range 
    #     ... -> all lengths -> combine ion_matches ... -> final psmQ
    
    silac_mix <- eval(this_call$silac_mix)
    
    # if (a) nested or (b) noenzyme only
    if (length(silac_mix) > 1L) {
      lapply(c("psmC.txt", "psmQ.txt", "psmT2.txt", "psmT3.txt"), function (file) {
        dfs <- lapply(out_paths, function (sub_path) {
          fi <- file.path(sub_path, file)
          
          if (file.exists(fi)) 
            df <- readr::read_tsv(fi, show_col_types = FALSE)
          else 
            df <- NULL
        })
        
        dfs <- dplyr::bind_rows(dfs)
        readr::write_tsv(dfs, file.path(out_path, file))
        
        invisible(NULL)
      })
      
      this_call$silac_mix <- NULL # just in case 
    }
    else {
      ans <- tryCatch(eval(this_call), error = function (e) NULL)
    }
    
    message("Done (noenzyme search).")
    options(show.error.messages = FALSE)
    stop()
  }
  else {
    this_call$bypass_noenzyme <- TRUE
    this_call$silac_noenzyme <- TRUE # for regular add_protacc
    ans <- tryCatch(eval(this_call), error = function (e) NULL)
  }
  
  invisible(NULL)
}


#' Combines the results of ion matches.
#' 
#' @param out_path A parent output path.
#' @param out_paths Sub output pathes.
#' @param type The type of data for combining.
combine_ion_matches <- function (out_path, out_paths, type = "ion_matches_") 
{
  out_path_temp  <- create_dir(file.path(out_path, "temp"))
  out_paths_temp <- lapply(out_paths, function(x) file.path(x, "temp"))
  
  pat  <- paste0(type, "[0-9]+\\.rds$")
  pat2 <- paste0(type, "([0-9]+)\\.rds$")
  xs   <- list.files(out_paths_temp[[1]], pattern = pat)
  
  files <- if (length(xs)) 
    paste0(type, sort(as.integer(gsub(pat2, "\\1", xs))), ".rds")
  else 
    NULL
  
  if (!length(files)) {
    warning("Files not found: ", type)
    return(NULL)
  }
  
  for (i in seq_along(files)) {
    ans <- lapply(out_paths_temp, function (x) qs::qread(file.path(x, files[i])))
    ans <- dplyr::bind_rows(ans)
    qs::qsave(ans, file.path(out_path_temp, files[i]), preset = "fast")
  }
  
  invisible(NULL)
}


#' Combines PSMs from sub folders.
#'
#' @param sub_paths A list of sub paths.
#' @param groups A character vector of sample groups. The group names will also
#'   be applied to the \code{pep_group} in PSM tables.
#' @param out_path An output path.
comine_PSMsubs <- function (sub_paths, groups, out_path) 
{
  lapply(c("psmC.txt", "psmQ.txt", "psmT2.txt", "psmT3.txt"), function (file) {
    df <- lapply(groups, function (group) {
      fi <- file.path(out_path, group, file)
      
      if (file.exists(fi)) {
        df <- readr::read_tsv(fi, show_col_types = FALSE)
        df$pep_group <- group
        nms <- names(df)
        
        df <- dplyr::bind_cols(
          df[, grepl("^prot_", nms)], 
          df[, grepl("^pep_", nms)], 
          df[, !grepl("^(prot|pep)_", nms)], 
        )
      }
      else {
        df <- NULL
      }
    }) 
    
    df <- dplyr::bind_rows(df)
    readr::write_tsv(df, file.path(out_path, file))
    
    invisible(NULL)
  })
  
  create_dir(file.path(out_path, "Calls"))
  file.copy(file.path(sub_paths[[1]], "Calls"), out_path, recursive = TRUE)
  combine_ion_matches(out_path, sub_paths, type = "ion_matches_")
  suppressWarnings(combine_ion_matches(out_path, sub_paths, type = "reporters_"))

  invisible(NULL)
}


#' Wrapper of \link{matchMS}.
#' 
#' With the calibration of precursor masses.
#' 
#' @param ... Arguments for \link{matchMS}.
matchMS_ms1calib <- function (...)
{
  matchMS(...)
  
  this_call <- match.call()
  fun <- as.character(this_call[1])
  this_fml <- formals()
  
  out_path <- this_call[["out_path"]]
  mgf_path <- this_call[["mgf_path"]]
  calib_masses <- this_call[["calib_masses"]]
  
  if (passed_ms1calib <- check_ms1calib(out_path, calib_masses)) {
    file <- file.path(mgf_path, "ppm_ms1calib.rds")
    ppm_ms1_calib <- if (file.exists(file)) qs::qread(file) else stop()
    ppm_ms1_af <- ppm_ms1_calib[["ppm_ms1_af"]]
  }
  else 
    stop()
  
  matchMS(..., 
          ppm_ms1 = ppm_ms1_af, 
          bypass_mgf = TRUE, 
          calib_masses = FALSE, 
  )
}


