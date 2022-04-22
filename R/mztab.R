#' Makes an mzTab file.
#'
#' With proteoM searches and proteoQ preprocessing.
#'
#' @param out_path A parent path where the outputs of \code{PSM}, \code{Peptide}
#'   and \code{Protein} files and folders are.
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @import dplyr
make_mztab <- function (out_path = stop("Provide the path.", call. = FALSE)) 
{
  ## MTD
  load(file.path(out_path, "Calls", "matchMS.rda"))

  # Header
  hdrs <- local({
    nms <- c("mzTab-version", "mzTab-mode", "mzTab-type", "title", "description") 
    vals <- c("1", "Complete", "Identification", "null", "null")
    data.frame(nm = nms, val = vals)
  })

  # Instrument and MGF format
  ans_mgfs <- local({
    mgf_path <- call_pars$mgf_path
    info_mgfs <- qs::qread(file.path(mgf_path, "info_format.rds"))
    
    data_format <- info_mgfs$data_format
    val_data_format <- paste0("[MS, , ", data_format, ", ]")
    mgf_format <- info_mgfs$mgf_format
    val_mgf_format <- paste0("[MS, , ", mgf_format, ", ]")
    
    mgf_queries <- qs::qread(file.path(mgf_path, "mgf_queries.rds"))
    raw_files <- names(qs::qread(file.path(mgf_path, "raw_indexes.rds")))
    
    ans_mgfs <- vector("list", length(raw_files))
    
    for (i in seq_along(ans_mgfs)) {
      nm_format <- paste0("ms_run[", i, "]-format")
      nm_location <- paste0("ms_run[", i, "]-location")
      nm_id_format <- paste0("ms_run[", i, "]-id_format")
      
      val_location <- raw_files[i]
      
      nms <- c(nm_format, nm_location, nm_id_format)
      vals <- c(val_data_format, val_location, val_mgf_format)
      
      ans_mgfs[[i]] <- data.frame(nm = nms, val = vals)
    }
    
    data.frame(do.call(rbind, ans_mgfs))
  })
  
  # Software settings
  load(file.path(out_path, "Calls", "proteoM.rda"))
  
  proteom_info <- session_info$otherPkgs[[1]]
  proteom_ver <- proteom_info$Version

  ans_software_1 <- local({
    ln_software_1 <- data.frame(nm = "software[1]", 
                                val = paste0("[MS, MS:0000000, proteoM,", 
                                             proteom_ver, "]"))
    
    idxes <- which(unlist(lapply(call_pars, is.null)))
    call_pars[[idxes]] <- "NULL"
    rm(list = "idxes")
    
    fixedmods <- call_pars$fixedmods
    varmods <- call_pars$varmods
    fixedmods <- lapply(fixedmods, proteoM::find_unimod)
    varmods <- lapply(varmods, proteoM::find_unimod)
    
    call_pars <- lapply(call_pars, paste, collapse = ", ")
    call_pars <- data.frame(do.call(rbind, call_pars))
    
    vals <- mapply(paste, rownames(call_pars), call_pars[, 1], sep = "=", 
                   USE.NAMES = FALSE)
    nms <- rep("software[1]-setting", length(call_pars))
    ans_call_pars <- data.frame(nm = nms, val = vals)

    ans_fixedmods <- local({
      ans_fixedmods <- vector("list", length(fixedmods))
      
      for (i in seq_along(fixedmods)) {
        fixedmod <- fixedmods[[i]]
        title <- fixedmod$title
        monomass <- fixedmod$monomass
        site <- fixedmod$position_site
        position <- names(site)
        
        nm_t <- paste0("fixedmods[", i, "]")
        nm_p <- paste0(nm_t, "-position")
        nm_s <- paste0(nm_t, "-site")
        nms <- c(nm_t, nm_p, nm_s)
        
        vals <- c(paste0("[CHEMMOD, CHEMMOD:", monomass, ", ", title, ",]"), 
                  position, unname(site))
        
        ans_fixedmods[[i]] <- data.frame(nm = nms, val = vals)
      }
      
      data.frame(do.call(rbind, ans_fixedmods))
    })
    
    ans_varmods <- local({
      ans_varmods <- vector("list", length(varmods))
      
      for (i in seq_along(varmods)) {
        varmod <- varmods[[i]]
        title <- varmod$title
        monomass <- varmod$monomass
        site <- varmod$position_site
        position <- names(site)

        nm_t <- paste0("varmods[", i, "]")
        nm_p <- paste0(nm_t, "-position")
        nm_s <- paste0(nm_t, "-site")
        nms <- c(nm_t, nm_p, nm_s)

        vals <- c(paste0("[CHEMMOD, CHEMMOD:", monomass, ", ", title, ",]"), 
                  position, unname(site))
        
        ans_varmods[[i]] <- data.frame(nm = nms, val = vals)
      }
      
      data.frame(do.call(rbind, ans_varmods))
    })
    
    do.call(rbind, list(ln_software_1, ans_call_pars, 
                        ans_fixedmods, ans_varmods))
  })
  
  mtd <- do.call(rbind, list(hdrs, ans_mgfs, ans_software_1))
  mtd <- cbind(field = "MTD", mtd)

  ## Proteins
  df_prots <- readr::read_tsv(file.path(out_path, "Protein", "Protein.txt"), 
                              show_col_types = FALSE) 

  df_peps <- readr::read_tsv(file.path(out_path, "Peptide", "Peptide.txt"), 
                             show_col_types = FALSE) 

  df_shared_prot_accs <- local({
    df <- unique(df_peps[, c("prot_acc", "shared_prot_accs")]) %>% 
      dplyr::mutate(len = stringr::str_count(shared_prot_accs, ",")) %>% 
      dplyr::arrange(prot_acc, -len) %>% 
      dplyr::group_by(prot_acc) %>% 
      dplyr::mutate(group_count = dplyr::n())
    
    df_1 <- df %>% 
      dplyr::filter(group_count == 1L) %>% 
      dplyr::mutate(shared_prot_accs = "null")
    
    df_n <- df %>% 
      dplyr::filter(group_count > 1L) %>% 
      dplyr::filter(row_number() == 1L)
    
    dplyr::bind_rows(df_1, df_n) %>% 
      dplyr::select(-len, -group_count)
  })

  cols_prt <- c(
    "PRH", "accession", "description", "taxid", "species", "database", 
    "database_version", "search_engine", "best_search_engine_score[1]", 
    "ambiguity_members", "modifications", "protein_coverage", 
    "search_engine_score[1]_ms_run[1]", 
    "num_psms_ms_run[1]", "num_psms_unique_ms_run[1]", 
    "num_peptides_distinct_ms_run[1]", "num_peptides_unique_ms_run[1]",
    "opt_global_mass",
    "reliability"
  )

  prt <- matrix(ncol = length(cols_prt), 
                nrow = nrow(df_prots)) %>% 
    data.frame(check.names = FALSE) %>% 
    setNames(cols_prt)
  
  for (i in seq_along(prt)) prt[[i]] <- "null"

  prt <- prt %>% 
    dplyr::mutate(PRH = "PRT", 
                  accession = df_prots$prot_acc, 
                  description = df_prots$prot_desc, 
                  taxid = "null", 
                  species = df_prots$species, 
                  database = "null", 
                  database_version = "null", 
                  search_engine = "proteoM", 
                  "best_search_engine_score[1]" = df_prots$prot_es, 
                  ) %>% 
    dplyr::left_join(df_shared_prot_accs, by = c("accession" = "prot_acc")) %>% 
    dplyr::mutate(ambiguity_members = shared_prot_accs) %>% 
    dplyr::select(-shared_prot_accs)
  
  prt <- prt %>% 
    dplyr::mutate(modifications = "null") %>% 
    dplyr::left_join(df_prots[, c("prot_acc", "prot_cover")], 
                     by = c("accession" = "prot_acc")) %>% 
    dplyr::mutate(protein_coverage = prot_cover) %>% 
    dplyr::select(-prot_cover)
  
  prt <- prt %>% 
    dplyr::left_join(df_prots[, c("prot_acc", "prot_mass")], 
                     by = c("accession" = "prot_acc")) %>% 
    dplyr::mutate(opt_global_mass = prot_mass) %>% 
    dplyr::select(-prot_mass)
  
  prt <- prt %>% 
    dplyr::mutate(reliability = 1L)
  
  prt <- prt %>% 
    dplyr::left_join(df_prots[, c("prot_acc", "prot_n_psm", "prot_n_uniqpsm", 
                              "prot_n_pep", "prot_n_uniqpep")], 
                     by = c("accession" = "prot_acc")) %>% 
    dplyr::mutate("num_psms_distinct_ms_run[1]" = prot_n_psm, 
                  "num_psms_unique_ms_run[1]" = prot_n_uniqpsm, 
                  "num_peptides_distinct_ms_run[1]" = prot_n_pep, 
                  "num_peptides_unique_ms_run[1]" = prot_n_uniqpep) %>% 
    dplyr::select(-c("prot_n_psm", "prot_n_uniqpsm", 
                     "prot_n_pep", "prot_n_uniqpep"))
  
  prt <- local({
    df <- df_prots[, grepl("^I[0-9]+", names(df_prots))]
    colnames(df) <- paste0("protein_abundance_study_variable[", 1:ncol(df), "]")
    
    df <- cbind(prot_acc = df_prots$prot_acc, 
                df) %>% 
      data.frame(check.names = FALSE)
    
    prt <- prt %>% 
      dplyr::left_join(df, by = c("accession" = "prot_acc"))
  })
  
  ## Peptides
  cols_pep <- c(
    "PEH", "sequence", "accession", "unique", "database", "database_version", 
    "search_engine", "best_search_engine_score[1]", 
    "modifications", "retention_time", "charge", "mass_to_charge",
    # "num_psms_ms_run[1]", "num_psms_unique_ms_run[1]", 
    "opt_global_mass", "opt_global_missed_cleavages",
    "reliability"
  )
  
  pep <- matrix(ncol = length(cols_pep), 
                nrow = nrow(df_peps)) %>% 
    data.frame(check.names = FALSE) %>% 
    setNames(cols_pep)
  
  for (i in seq_along(pep)) pep[[i]] <- "null"
  
  if ("pep_seq" %in% names(df_peps)) 
    df_peps$sequence <- df_peps$pep_seq
  else if ("pep_seq_mod" %in% names(df_peps))
    df_peps$sequence <- df_peps$pep_seq_mod
  
  pep <- pep %>% 
    dplyr::mutate(PEH = "PEP", 
                  sequence = df_peps$sequence,
                  accession = df_peps$prot_acc, 
                  unique = df_peps$pep_isunique, 
                  database = "null", 
                  database_version = "null", 
                  search_engine = "proteoM", 
                  "best_search_engine_score[1]" = df_peps$pep_score, 
                  modifications = df_peps$pep_vmod, 

                  opt_global_missed_cleavages = df_peps$pep_miss, 
                  reliability = 1L, )
  
  pep <- local({
    df <- df_peps[, grepl("^I[0-9]+", names(df_peps))]
    colnames(df) <- paste0("peptide_abundance_study_variable[", 1:ncol(df), "]")
    
    cbind(pep, df) %>% 
      data.frame(check.names = FALSE)
  })
  
  ## PSMs
  psm_files <- list.files(path = file.path(out_path, "PSM"),
                          pattern = "TMTset[0-9]+_LCMSinj[0-9]+_PSM_N\\.txt$",
                          all.files = TRUE)

  df_psms <- lapply(psm_files, 
                    function (x) readr::read_tsv(file.path(out_path, "PSM", x), 
                                                 show_col_types = FALSE)) %>% 
    dplyr::bind_rows()

  cols_psm <- c(
    "PSH", "sequence", "PSM_ID", "accession", "unique", "database", "database_version", 
    "search_engine", "best_search_engine_score[1]", "modifications", "retention_time", 
    "charge", "exp_mass_to_charge", "calc_molecular_weight", "opt_global_mass", 
    "opt_global_missed_cleavages", "opt_global_spectrum_file", 
    "opt_global_scan_number", "pre", "post", "start", "end",
    "reliability")

  local({
    # redundancy kept at "prot_acc" and "pep_ivmod"
    uniq_id <- c("raw_file", "pep_scan_num", "prot_acc", "pep_ivmod")
  })
  
  
  psm <- matrix(ncol = length(cols_psm), 
                nrow = nrow(df_psms)) %>% 
    data.frame(check.names = FALSE) %>% 
    setNames(cols_psm)
  
  for (i in seq_along(psm)) psm[[i]] <- "null"

  psm <- psm %>% 
    dplyr::mutate(PSH = "PSM", 
                  sequence = df_psms$pep_seq,
                  # PSM_ID = seq_len(nrow(df_psms)), 
                  accession = df_psms$prot_acc, 
                  unique = df_psms$pep_isunique, 
                  database = "null", 
                  database_version = "null", 
                  search_engine = "proteoM", 
                  "best_search_engine_score[1]" = df_psms$pep_score, 
                  modifications = df_psms$pep_vmod, 
                  retention_time = df_psms$pep_ret_range, 
                  charge = df_psms$pep_exp_z, 
                  exp_mass_to_charge = df_psms$pep_exp_mz, 
                  calc_molecular_weight = df_psms$pep_calc_mr, 
                  opt_global_mass = df_psms$pep_exp_mr,
                  opt_global_missed_cleavages = df_psms$pep_miss, 
                  opt_global_spectrum_file = df_psms$raw_file, 
                  opt_global_scan_number = df_psms$pep_scan_num, 
                  pre = df_psms$pep_res_before,
                  post = df_psms$pep_res_after,
                  start = df_psms$pep_start,
                  end = df_psms$pep_end,
                  reliability = 1L, ) 

  psm <- local({
    df <- df_psms[, grepl("^I[0-9]+", names(df_psms))]
    colnames(df) <- paste0("peptide_abundance_study_variable[", 1:ncol(df), "]")
    
    cbind(psm, df) %>% 
      data.frame(check.names = FALSE)
  })
  
  ## Outputs
  lines_mtd <- capture.output(write.table(mtd, stdout(), sep = "\t", quote = FALSE, 
                                          row.names = FALSE, col.names = FALSE))
  lines_mtd <- paste(lines_mtd, collapse = "\n")
  lines_mtd <- paste(lines_mtd, "\n")
  
  lines_prt <- capture.output(write.table(prt, stdout(), sep = "\t", quote = FALSE, 
                                          row.names = FALSE))
  lines_prt <- paste(lines_prt, collapse = "\n")
  lines_prt <- paste(lines_prt, "\n")
  
  lines_pep <- capture.output(write.table(pep, stdout(), sep = "\t", quote = FALSE, 
                                          row.names = FALSE))
  lines_pep <- paste(lines_pep, collapse = "\n")
  lines_pep <- paste(lines_pep, "\n")
  
  lines_psm <- capture.output(write.table(psm, stdout(), sep = "\t",quote = FALSE,  
                                          row.names = FALSE))
  lines_psm <- paste(lines_psm, collapse = "\n")
  lines_psm <- paste(lines_psm, "\n")

  dir.create(file.path(out_path, "mzTab"), showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_path, "mzTab", "mztab.mzTab")
  
  out <- Reduce(append, list(lines_mtd, lines_prt, lines_pep, lines_psm)) %T>% 
    writeLines(out_file)
}


