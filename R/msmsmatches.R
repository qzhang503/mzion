#' Finds the indexes of top-n entries without re-ordering.
#'
#' At length(x) >= n, the length of output may be shorter than n with ties.
#'
#' @param x A numeric vector.
#' @param n The number of top entries to keep.
#' @return The indexes of the top-n entries.
#' @examples
#' \donttest{
#' which_topx(c(1:5), 50)
#'
#' length(which_topx(sample(5000, 500), 100))
#'
#' length(which_topx(sample(100, 100, replace = TRUE), 100))
#' }
which_topx <- function(x, n = 50L, ...) {

  len <- length(x)
  p <- len - n

  if (p  <= 0L) return(seq_along(x))

  xp <- sort(x, partial = p, ...)[p]

  which(x > xp)
}


#' Finds the indexes of top-n entries without re-ordering.
#'
#' @param x A numeric vector.
#' @param n The number of top entries to keep.
#' @return The indexes of the top-n entries.
which_topx2 <- function(x, n = 50L, ...) {

  len <- length(x)
  p <- len - n

  if (p  <= 0L) return(seq_along(x))

  xp <- sort(x, partial = p, ...)[p]

  ans <- which(x > xp)

  # in case of ties -> length(ans) < n
  # detrimental e.g. ms2_n = 500 and n = 100
  #   -> expect 100 `ms2_moverzs` guaranteed but may be only 99
  #
  # MGF `ms2_moverzs` is increasing
  # `ans2` goes first to ensure non-decreasing index for `ms2_moverzs`

  d <- n - length(ans)

  if (d > 0L) {
    ans2 <- which(x == xp)
    ans <- c(ans2[1:d], ans)
    ans <- sort(ans)
  }

  invisible(ans)
}


#' Finds the top-n entries without re-ordering.
#'
#' @inheritParams which_topx
#' @return The top-n entries.
topx <- function(x, n = 50L, ...) {

  len <- length(x)
  p <- len - n

  if (p  <= 0L) return(x)

  xp <- sort(x, partial = p, ...)[p]

  x[x > xp]
}


#' Finds the numeric difference in ppm.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The difference between \eqn{x} and \eqn{y} in ppm.
find_ppm_error <- function (x = 1000, y = 1000.01) {
  (y - x)/y * 1E6
}


#' Finds the error range of a number.
#'
#' Assumes \eqn{x} is positive without checking.
#'
#' @param x A numeric value.
#' @param ppm Numeric; the ppm allowed from \code{x}.
#' @return The lower and the upper bound to \eqn{x} by \eqn{ppm}.
find_mass_error_range <- function (x = 500L, ppm = 20L) {
  d <- x * ppm/1E6
  c(x-d, x+d)
}


#' Splits data by groups then into chunks.
#' 
#' Not currently used: groupProts <- map_pepprot <- chunk_groupsplit
#' 
#' @inheritParams chunksplit
#' @param f A factor; see also base \code{split}.
chunk_groupsplit <- function (data, f, n_chunks) {

  if (n_chunks <= 1L) return(data)

  data <- split(data, f)
  len <- length(data)

  labs <- levels(cut(1:len, n_chunks))

  x <- cbind(
    lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
    upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))

  grps <- findInterval(1:len, x[, 1])

  data <- split(data, grps)

  lapply(data, function (x) do.call(rbind, x))
}


#' Splits data into chunks by length.
#'
#' @param data Input data.
#' @param n_chunks The number of chunks.
#' @param type The type of data for splitting.
chunksplit <- function (data, n_chunks = 5L, type = "list") {

  stopifnot(type %in% c("list", "row"))

  if (n_chunks <= 1L) return(data)

  if (type == "list") {
    len <- length(data)
  } else if (type == "row") {
    len <- nrow(data)
  } else {
    stop("Unknown type.", call. = TRUE)
  }

  if (len == 0L) return(data)

  labs <- levels(cut(1:len, n_chunks))

  x <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
             upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))

  grps <- findInterval(1:len, x[, 1])
  split(data, grps)
}


#' Splits data into chunks with approximately equal sizes.
#'
#' @param nx Positive integer; an arbitrarily large number for data to be split
#'   into for estimating the cumulative sizes.
#' @inheritParams chunksplit
chunksplitLB <- function (data, n_chunks = 5L, nx = 100L, type = "list") {

  stopifnot(type %in% c("list", "row"))

  if (n_chunks <= 1L) return(data)

  if (type == "list") {
    len <- length(data)
  } else if (type == "row") {
    len <- nrow(data)
  } else {
    stop("Unknown type.", call. = TRUE)
  }

  if (len == 0L) return(data)

  # The finer groups by 'nx'
  grps_nx <- local({
    labsx <- levels(cut(1:len, nx))

    xx <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1",
                                              labsx))),
                upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1",
                                                labsx))))

    findInterval(1:len, xx[, 1])
  })

  # The equated size for a chunk
  size_chunk <- local({
    size_nx <- data %>%
      split(., grps_nx) %>%
      lapply(object.size) %>%
      cumsum()

    size_nx[length(size_nx)]/n_chunks
  })

  #  Intervals
  grps <- local({
    size_data <- data %>%
      lapply(object.size) %>%
      cumsum()

    # the position indexes
    ps <- purrr::map_dbl(1:(n_chunks-1), function (x) {
      which(size_data < size_chunk * x) %>% `[`(length(.))
    })

    grps <- findInterval(1:len, ps)
  })

  split(data, grps)
}


#' Searches for MS ions.
#'
#' Database searches of MSMS data.
#' @param out_path A file path of outputs.
#' @param mgf_path The file path to a list of MGF files. The experimenters need
#'   to supply the files. Note that the supported MGFs are from MSConvert or
#'   Proteome Discoverer.
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. The experimenters need to supply the files.
#' @param acc_type Character string(s); the types of protein accessions in one
#'   of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
#'   the corresponding regular expressions need to be supplied via argument
#'   \code{acc_pattern}.
#' @param acc_pattern Regular expression(s) describing the patterns to separate
#'   the header lines of fasta entries. At the \code{NULL} default, the pattern
#'   will be automated when \code{acc_type} are among c("uniprot_acc",
#'   "uniprot_id", "refseq_acc", "other"). See also \link{load_fasta2} for
#'   custom examples.
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param exclude_phospho_nl If TRUE, excludes neutral losses in the MS2 ion
#'   searches against variable modifications of phospho sites.
#' @param include_insource_nl Logical Logical; if TRUE, includes MS1 precursor
#'   masses with the losses of neutral species prior to MS2 fragmentation. The
#'   default is FALSE. The setting at TRUE remains experimenting by allowing
#'   additional masses in the universe of MS1 precursors. The offsets in NLs (by
#'   adding back precursor masses) have not yet taken into account in MS2 ion
#'   searches. A more systemically approach such as \code{open} MS1 masses might
#'   be developed in the future.
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from proteins. The enzyme
#'   is currently \code{trypsin}.
#' @param maxn_fasta_seqs Integer; the maximum number of protein sequences in
#'   fasta files.
#' @param maxn_vmods_setscombi Integer; the maximum number of combinatorial
#'   variable modifications and neutral losses.
#' @param maxn_vmods_per_pep The maximum number of variable modifications per
#'   peptide.
#' @param maxn_sites_per_vmod Integer; the maximum number of combinatorial
#'   variable modifications per site in a per peptide sequence.
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @param min_len Integer; the minimum length of peptides. Shorter peptides will
#'   be excluded.
#' @param max_len Integer; the maximum length of peptides. Longer peptides will
#'   be excluded.
#' @param max_miss The maximum number of mis-cleavages per peptide sequence.
#' @param min_mass The minimum precursor mass for interrogation.
#' @param max_mass The maximum precursor mass for interrogation.
#' @param min_ms2mass The minimum MS2 mass for interrogation.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @param topn_ms2ions A non-negative integer; the top-n species for uses in MS2
#'   ion searches. The default is to use the top-100 ions in an MS2 event.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
#' @param quant A quantitation method. The default is "none". Additional choices
#'   include \code{tmt6} etc. For other multiplicities of \code{tmt}, use the
#'   compatible higher plexes, for example, \code{tmt16} for \code{tmt12} etc.
#'   and \code{tmt10} for \code{tmt8} etc.
#' @param target_fdr Numeric; a targeted false-discovery rate (FDR) at the
#'   levels of PSM, peptide or protein. See also argument \code{fdr_type}.
#' @param fdr_type Character string; the type of FDR control. The value is in
#'   one of c("psm", "peptide", "protein"). Note that \code{fdr_type = protein}
#'   is equivalent to \code{fdr_type = peptide} with the additional filtration
#'   of data at \code{prot_tier == 1}. A variant is to set \code{fdr_type =
#'   psm}, followed by a data filtration at \code{prot_tier == 1}.
#' @param combine_tier_three Logical; if TRUE, combines all protein results to
#'   the output of \code{psmQ.txt}. Outputs under the option TRUE are often
#'   comparable to Mascot outputs with FDR controls at the levels of PSMs or
#'   peptides.
#' @param .path_cache The file path of cached search parameters. At the NULL
#'   default, the path is \code{"~/proteoM/.MSearches/Cache/Calls/"}. The
#'   parameter is for users' awareness of the structure of file folders and the
#'   default is suggested.
#' @param .path_fasta The parent file path to the theoretical masses of MS1
#'   precursors. At the NULL default, the path is \code{gsub("(.*)\\.[^\\.]*$",
#'   "\\1", get("fasta", envir = environment())[1])}. The parameter is for
#'   users' awareness of the structure of file folders and the default is
#'   suggested.
#' @param digits Integer; the number of decimal places to be used.
#' @seealso \link{load_fasta2} for setting the values of \code{acc_type} and
#'   \code{acc_pattern}. \link{parse_unimod} for the grammar of Unimod.
#' @return A list of complete PSMs in \code{psmC.txt}; a list of quality PSMs in
#'   \code{psmQ.txt}.
#' @examples
#' \donttest{
#' matchMS(
#'   fasta = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   max_miss = 2,
#'   quant = "tmt10",
#'   fdr_type = "protein",
#'   out_path = "~/proteoM/examples",
#' )
#'
#' \dontrun{
#' # Supposed phosphopeptides and TMTpro
#' matchMS(
#'   fasta = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("TMTpro (N-term)", "TMTpro (K)", "Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'               "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'              "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   max_miss = 2,
#'   quant = "tmt16",
#'   fdr_type = "protein",
#'   out_path = "~/proteoM/examples",
#' )
#' }
#'
#' }
#' @export
matchMS <- function (out_path = "~/proteoM/outs",
                     mgf_path = file.path(out_path, "mgf"),
                     fasta = c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
                               "~/proteoM/dbs/fasta/crap/crap.fasta"),
                     acc_type = c("uniprot_acc", "other"),
                     acc_pattern = NULL,
                     fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                   "Carbamidomethyl (C)"),
                     varmods = c("Acetyl (Protein N-term)",
                                 "Oxidation (M)", "Deamidated (N)",
                                 "Gln->pyro-Glu (N-term = Q)"),
                     include_insource_nl = FALSE,
                     exclude_phospho_nl = TRUE, 
                     enzyme = c("trypsin"),
                     maxn_fasta_seqs = 200000L,
                     maxn_vmods_setscombi = 64L,
                     maxn_vmods_per_pep = 5L,
                     maxn_sites_per_vmod = 3L,
                     maxn_vmods_sitescombi_per_pep = 64L,
                     min_len = 7L, max_len = 50L, max_miss = 2L,
                     min_mass = 500L, max_mass = 6000L, min_ms2mass = 110L, 
                     type_ms2ions = "by",
                     topn_ms2ions = 100L,
                     minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L,
                     ppm_reporters = 10L,
                     quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"),
                     target_fdr = 0.01,
                     fdr_type = c("psm", "peptide", "protein"),
                     combine_tier_three = FALSE,
                     .path_cache = NULL, 
                     .path_fasta = NULL,
                     digits = 4L) {

  options(digits = 9L)

  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"),
                   fun = as.character(match.call()[[1]]))
      }
    },
    add = TRUE
  )

  ## Preparation
  # accession pattern
  if ((!is.null(acc_pattern)) && acc_pattern == "") {
    acc_pattern <- NULL
  }
  
  # numeric types 
  stopifnot(vapply(c(maxn_fasta_seqs, maxn_vmods_setscombi, maxn_vmods_per_pep, 
                     maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                     min_len, max_len, max_miss, topn_ms2ions, minn_ms2, 
                     min_mass, max_mass, min_ms2mass, 
                     ppm_ms1, ppm_ms2, ppm_reporters, digits, 
                     target_fdr), 
                   is.numeric, logical(1L)))

  # (a) integers casting for parameter matching when calling cached)
  maxn_fasta_seqs <- as.integer(maxn_fasta_seqs)
  maxn_vmods_setscombi <- as.integer(maxn_vmods_setscombi)
  maxn_vmods_per_pep <- as.integer(maxn_vmods_per_pep)
  maxn_sites_per_vmod <- as.integer(maxn_sites_per_vmod)
  maxn_vmods_sitescombi_per_pep <- as.integer(maxn_vmods_sitescombi_per_pep)
  min_len <- as.integer(min_len)
  max_len <- as.integer(max_len)
  max_miss <- as.integer(max_miss)
  topn_ms2ions <- as.integer(topn_ms2ions)
  minn_ms2 <- as.integer(minn_ms2)
  min_mass <- as.integer(min_mass)
  max_mass <- as.integer(max_mass)
  min_ms2mass <- as.integer(min_ms2mass)
  ppm_ms1 <- as.integer(ppm_ms1)
  ppm_ms2 <- as.integer(ppm_ms2)
  ppm_reporters <- as.integer(ppm_reporters)
  digits <- as.integer(digits)
  
  stopifnot(min_len >= 0L, max_len >= min_len, max_miss <= 10L, 
            min_mass >= 0L, max_mass >= min_mass, min_ms2mass >= 0L)

  # (b) doubles
  target_fdr <- as.double(target_fdr)
  
  stopifnot(target_fdr < .5)
  
  # fdr_type
  fdr_type <- rlang::enexpr(fdr_type)
  oks <- eval(formals()[["fdr_type"]])

  if (length(fdr_type) > 1L) {
    fdr_type <- oks[[1]]
  } else {
    fdr_type <- rlang::as_string(fdr_type)
  }

  stopifnot(fdr_type %in% oks)
  rm(list = c("oks"))

  # Quantitation method
  quant <- rlang::enexpr(quant)
  oks <- eval(formals()[["quant"]])

  if (length(quant) > 1L) {
    quant <- oks[[1]]
  } else {
    quant <- rlang::as_string(quant)
  }

  stopifnot(quant %in% oks, length(quant) == 1L)
  rm(list = c("oks"))

  # Output path
  out_path <- create_dir(out_path)

  filelist <- list.files(path = file.path(mgf_path), pattern = "\\.mgf$")

  if (!length(filelist)) {
    stop("No `.mgf` files under ", mgf_path, call. = FALSE)
  }

  rm(list = c("filelist"))
  
  # file paths
  if (is.null(.path_cache)) {
    .path_cache <- "~/proteoM/.MSearches/Cache/Calls/"
  }
  
  if (is.null(.path_fasta)) {
    .path_fasta <- file.path(gsub("(.*)\\.[^\\.]*$", "\\1", fasta[1]))
  }
  
  .path_cache <- create_dir(.path_cache)
  .path_ms1masses <- create_dir(file.path(.path_fasta, "ms1masses"))

  ## Theoretical MS1 masses
  res <- calc_pepmasses2(
    fasta = fasta,
    acc_type = acc_type,
    acc_pattern = acc_pattern,
    fixedmods = fixedmods,
    varmods = varmods,
    include_insource_nl = include_insource_nl,
    exclude_phospho_nl = exclude_phospho_nl, 
    enzyme = enzyme,
    maxn_fasta_seqs = maxn_fasta_seqs,
    maxn_vmods_setscombi = maxn_vmods_setscombi,
    maxn_vmods_per_pep = maxn_vmods_per_pep,
    maxn_sites_per_vmod = maxn_sites_per_vmod,
    min_len = min_len,
    max_len = max_len,
    max_miss = max_miss,
    out_path = out_path,
    digits = digits,
    .path_cache = .path_cache, 
    .path_ms1masses = .path_ms1masses
  )

  ## Bin theoretical peptides
  bin_ms1masses(res = res, 
                min_mass = min_mass, 
                max_mass = max_mass, 
                ppm_ms1 = ppm_ms1, 
                .path_cache = .path_cache, 
                .path_ms1masses = .path_ms1masses)
  
  rm(list = c("res"))
  gc()

  ## MGFs
  load_mgfs(mgf_path = mgf_path,
            min_mass = min_mass,
            max_mass = max_mass, 
            min_ms2mass = min_ms2mass,
            topn_ms2ions = topn_ms2ions,
            ppm_ms1 = ppm_ms1,
            ppm_ms2 = ppm_ms2,
            index_ms2 = FALSE)

  ## MSMS matches
  ms2match(mgf_path = mgf_path,
           aa_masses_all = readRDS(file.path(out_path, "temp", "aa_masses_all.rds")),
           out_path = out_path,
           mod_indexes = find_mod_indexes(out_path),
           type_ms2ions = type_ms2ions,
           maxn_vmods_per_pep = maxn_vmods_per_pep,
           maxn_sites_per_vmod = maxn_sites_per_vmod,
           maxn_vmods_sitescombi_per_pep =
             maxn_vmods_sitescombi_per_pep,
           minn_ms2 = minn_ms2,
           ppm_ms1 = ppm_ms1,
           ppm_ms2 = ppm_ms2,
           min_ms2mass = min_ms2mass,
           quant = quant,
           ppm_reporters = ppm_reporters,
           
           # dummy for argument matching
           fasta = fasta,
           acc_type = acc_type,
           acc_pattern = acc_pattern,
           topn_ms2ions = topn_ms2ions,
           fixedmods = fixedmods,
           varmods = varmods,
           include_insource_nl = include_insource_nl,
           enzyme = enzyme,
           maxn_fasta_seqs = maxn_fasta_seqs,
           maxn_vmods_setscombi = maxn_vmods_setscombi,
           min_len = min_len,
           max_len = max_len,
           max_miss = max_miss,
           
           digits = digits)

  ## Peptide scores
  out <- calc_pepscores(topn_ms2ions = topn_ms2ions,
                        type_ms2ions = type_ms2ions,
                        target_fdr = target_fdr,
                        fdr_type = fdr_type,
                        min_len = min_len,
                        max_len = max_len,
                        penalize_sions = TRUE,
                        ppm_ms2 = ppm_ms2,
                        out_path = out_path,
                        digits = digits)
  
  ## Peptide ranks and score deltas between `pep_ivmod`
  out <- calc_peploc(out)

  gc()

  ## Protein accessions, score cut-offs and optional reporter ions
  out <- add_prot_acc(out, out_path)
  out <- calc_protfdr(out, target_fdr)
  out <- add_rptrs(out, quant, out_path)

  gc()

  ## Clean-ups
  out <- out %>%
    dplyr::mutate(pep_ms1_delta = ms1_mass - theo_ms1) %>%
    dplyr::rename(pep_scan_title = scan_title,
                  pep_exp_mz = ms1_moverz,
                  pep_exp_mr = ms1_mass,
                  pep_exp_z = ms1_charge,
                  pep_calc_mr = theo_ms1,
                  pep_delta = pep_ms1_delta,
                  pep_tot_int = ms1_int,
                  pep_ret_time = ret_time,
                  pep_scan_num = scan_num,
                  pep_ms2_n = ms2_n,
                  pep_frame = frame)

  out <- dplyr::bind_cols(
    out %>% .[grepl("^prot_", names(.))],
    out %>% .[grepl("^pep_", names(.))],
    out %>% .[grepl("^psm_", names(.))],
    out %>% .[!grepl("^prot_|^pep_|^psm_", names(.))],
  ) %>%
    reloc_col_after("pep_exp_z", "pep_exp_mr") %>%
    reloc_col_after("pep_calc_mr", "pep_exp_z") %>%
    reloc_col_after("pep_delta", "pep_calc_mr") %T>%
    readr::write_tsv(file.path(out_path, "psmC.txt"))

  ## psmC to psmQ
  out <- try_psmC2Q(out, out_path = out_path,
                    fdr_type = fdr_type,
                    combine_tier_three = combine_tier_three)

  .savecall <- TRUE

  invisible(out)
}


#' Helper of \link{ms2match}.
#' 
#' @param aa_masses_all All the amino acid lookup tables.
#' 
#' @inheritParams matchMS
#' @inheritParams calc_aamasses
#' @inheritParams load_mgfs
#' @examples
#' try_ms2match(mgf_path = mgf_path,
#'   aa_masses_all = aa_masses_all,
#'   out_path = out_path,
#'   mod_indexes = find_mod_indexes(out_path),
#'   type_ms2ions = type_ms2ions,
#'   maxn_vmods_per_pep = maxn_vmods_per_pep,
#'   maxn_sites_per_vmod = maxn_sites_per_vmod,
#'   maxn_vmods_sitescombi_per_pep =
#'     maxn_vmods_sitescombi_per_pep,
#'   minn_ms2 = minn_ms2,
#'   ppm_ms1 = ppm_ms1,
#'   ppm_ms2 = ppm_ms2,
#'   min_ms2mass = min_ms2mass,
#'   quant = quant,
#'   ppm_reporters = ppm_reporters,
#'
#'   # dummy for argument matching
#'   fasta = fasta,
#'   acc_type = acc_type,
#'   acc_pattern = acc_pattern,
#'   topn_ms2ions = topn_ms2ions,
#'   fixedmods = fixedmods,
#'   varmods = varmods,
#'   include_insource_nl = include_insource_nl,
#'   enzyme = enzyme,
#'   maxn_fasta_seqs = maxn_fasta_seqs,
#'   maxn_vmods_setscombi = maxn_vmods_setscombi,
#'   min_len = min_len,
#'   max_len = max_len,
#'   max_miss = max_miss,
#'
#'   target_fdr = target_fdr,
#'   fdr_type = fdr_type,
#'   combine_tier_three = combine_tier_three,
#'
#'   digits = digits)
try_ms2match <- function (mgf_path, aa_masses_all, out_path, mod_indexes, 
                          type_ms2ions, maxn_vmods_per_pep, maxn_sites_per_vmod,
                          maxn_vmods_sitescombi_per_pep, minn_ms2, ppm_ms1, 
                          ppm_ms2, min_ms2mass, quant, ppm_reporters,
                          fasta, acc_type, acc_pattern, topn_ms2ions,
                          fixedmods, varmods, include_insource_nl, enzyme, 
                          maxn_fasta_seqs, maxn_vmods_setscombi, 
                          min_len, max_len, max_miss, 
                          target_fdr, fdr_type, combine_tier_three, 
                          digits) {
  
  ans <- tryCatch(
    ms2match(mgf_path = mgf_path,
             aa_masses_all = aa_masses_all,
             out_path = out_path,
             mod_indexes = find_mod_indexes(out_path),
             type_ms2ions = type_ms2ions,
             maxn_vmods_per_pep = maxn_vmods_per_pep,
             maxn_sites_per_vmod = maxn_sites_per_vmod,
             maxn_vmods_sitescombi_per_pep =
               maxn_vmods_sitescombi_per_pep,
             minn_ms2 = minn_ms2,
             ppm_ms1 = ppm_ms1,
             ppm_ms2 = ppm_ms2,
             min_ms2mass = min_ms2mass,
             quant = quant,
             ppm_reporters = ppm_reporters,
             
             # dummy for argument matching
             fasta = fasta,
             acc_type = acc_type,
             acc_pattern = acc_pattern,
             topn_ms2ions = topn_ms2ions,
             fixedmods = fixedmods,
             varmods = varmods,
             include_insource_nl = include_insource_nl,
             enzyme = enzyme,
             maxn_fasta_seqs = maxn_fasta_seqs,
             maxn_vmods_setscombi = maxn_vmods_setscombi,
             min_len = min_len,
             max_len = max_len,
             max_miss = max_miss,
             
             digits = digits),
    error = function(e) NA
  )
  
  if (!is.null(ans)) {
    message(
      "Retry `matchMS()` with a new R session...\n", 
      "\"Error in save_call2\" at the end will not affect the research results."
    )
    
    fileConn <- file(file.path("~/matchMS.R"))
    
    lines <- c(
      "library(proteoM)\n",
      "proteoM::matchMS(",
      paste0("  out_path = \"", out_path, "\","),
      paste0("  mgf_path = \"", mgf_path, "\","),
      paste0("  fasta = c(\"", paste(fasta, collapse = "\", \""), "\"),"),
      paste0("  acc_type = c(\"", paste(acc_type, collapse = "\", \""), "\"),"),
      paste0("  acc_pattern = \"", acc_pattern, "\","),
      paste0("  fixedmods = c(\"", paste(fixedmods, collapse = "\", \""), "\"),"),
      paste0("  varmods = c(\"", paste(varmods, collapse = "\", \""), "\"),"),
      paste0("  include_insource_nl = ", include_insource_nl, ","),
      paste0("  enzyme = c(\"", paste(enzyme, collapse = "\", \""), "\"),"),
      paste0("  maxn_fasta_seqs = ", maxn_fasta_seqs, "L,"),
      paste0("  maxn_vmods_setscombi = ", maxn_vmods_setscombi, "L,"),
      paste0("  maxn_vmods_per_pep = ", maxn_vmods_per_pep, "L,"),
      paste0("  maxn_sites_per_vmod = ", maxn_sites_per_vmod, "L,"),
      paste0("  maxn_vmods_sitescombi_per_pep = ", maxn_vmods_sitescombi_per_pep, "L,"),
      paste0("  min_len = ", min_len, "L,"),
      paste0("  max_len = ", max_len, "L,"),
      paste0("  max_miss = ", max_miss, "L,"),
      paste0("  type_ms2ions = \"", type_ms2ions, "\","),
      paste0("  topn_ms2ions = ", topn_ms2ions, "L,"),
      paste0("  minn_ms2 = ", minn_ms2, "L,"),
      paste0("  ppm_ms1 = ", ppm_ms1, "L,"),
      paste0("  ppm_ms2 = ", ppm_ms2, "L,"),
      paste0("  ppm_reporters = ", ppm_reporters, "L,"),
      paste0("  quant = \"", quant, "\","),
      
      paste0("  target_fdr = ", target_fdr, ","),
      paste0("  fdr_type = \"", fdr_type, "\","),
      paste0("  combine_tier_three = ", combine_tier_three, ","),
      
      paste0("  digits = ", digits, "L"),
      ")\n",
      "unlink(\"~/matchMS.R\")"
    )
    
    writeLines(lines, fileConn)
    close(fileConn)
    
    rstudioapi::restartSession(command='source("~/matchMS.R")')
  } else {
    rm(list = "ans")
    gc()
  }
  
  invisible(NULL)
}


#' Helper of \link{psmC2Q}.
#'
#' @inheritParams psmC2Q
#' @importFrom magrittr %>% %T>%
try_psmC2Q <- function (out = NULL, out_path = NULL, fdr_type = "protein",
                        combine_tier_three = FALSE) {

  n_peps <- length(unique(out$pep_seq))
  n_prots <- length(unique(out$prot_acc))

  # `n_peps` and `n_prots` including both targets and decoys:
  # `n_prots` about 1:1
  # `n_peps` about 1.8:1

  if (n_peps > 1000000L && n_prots > 100000L) {
    out <- NA
  } else {
    out <- tryCatch(
      psmC2Q(out,
             out_path = out_path,
             fdr_type = fdr_type,
             combine_tier_three = combine_tier_three),
      error = function(e) NA
    )
  }

  if (length(out) == 1L && is.na(out)) {
    message("Retry with a new R session: \n\n",
            "proteoM:::reproc_psmC(\n",
            "  out_path = \"", out_path, "\",\n",
            "  fdr_type = \"", fdr_type, "\",\n",
            "  combine_tier_three  = ", combine_tier_three, "\n",
            ")")

    fileConn <- file(file.path("~/post_psmC.R"))

    lines <- c(
      "library(proteoM)\n",
      "proteoM:::reproc_psmC(",
      paste0("  out_path = \"", out_path, "\","),
      paste0("  fdr_type = \"", fdr_type, "\","),
      paste0("  combine_tier_three = ", combine_tier_three),
      ")\n",
      "unlink(\"~/post_psmC.R\")"
    )

    writeLines(lines, fileConn)
    close(fileConn)

    rstudioapi::restartSession(command='source("~/post_psmC.R")')
  } else {
    try(rm(list = c(".path_cache", ".path_ms1masses", ".time_stamp"),
           envir = .GlobalEnv))

    message("Done.")
  }

  invisible(out)
}



#' Reprocessing of \code{psmC.txt}.
#'
#' Protein grouping from \code{psmC.txt} to \code{psmQ.txt}.
#'
#' May solve some memory shortage issues for large data sets (e.g., over a
#' million peptide sequences * 35000 proteins from \code{psmC.txt}).
#'
#' @inheritParams matchMS
reproc_psmC <- function (out_path = NULL, fdr_type = "protein",
                         combine_tier_three = FALSE) {

  if (is.null(out_path)) {
    stop("`out_path` cannot be NULL.", call. = FALSE)
  }

  message("Please wait for the `Search completed` message...")

  readr::read_tsv(file.path(out_path, "psmC.txt"),
                  show_col_types = FALSE) %>%
    psmC2Q(out_path = out_path,
           fdr_type = fdr_type,
           combine_tier_three = combine_tier_three)

  message("Done.")
}


#' From \code{psmC.txt} to \code{psmQ.txt}.
#'
#' @param out A result of \code{psmC.txt}.
#' @inheritParams matchMS
psmC2Q <- function (out = NULL, out_path = NULL, fdr_type = "protein",
                    combine_tier_three = FALSE) {

  message("\n=================================\n",
          "prot_tier  prot_issig  prot_n_pep \n",
          "    1          [y]          \n",
          "    2          [n]          > 1\n",
          "    3          [n]          = 1\n",
          "=================================\n")

  out <- out %>%
    dplyr::filter(pep_issig, !pep_isdecoy, !grepl("^-", prot_acc))

  # Set aside one-hit wonders
  out3 <- out %>%
    dplyr::filter(!prot_issig, prot_n_pep == 1L) %>%
    dplyr::mutate(prot_tier = 3L)

  out <- dplyr::bind_rows(
    out %>% dplyr::filter(prot_issig),
    out %>% dplyr::filter(!prot_issig, prot_n_pep >= 2L)
  ) %>%
    dplyr::mutate(prot_tier = ifelse(prot_issig, 1L, 2L))

  gc()

  # Protein groups
  message("Building protein-peptide maps.")

  if (length(unique(out$prot_acc)) > 35000L) {
    if (fdr_type != "protein") {
      warning("Coerce to `fdr_type = protein` ",
              "and saved peptides of tier-2 proteins to `psmT2.txt`.",
              call. = FALSE)

      # dummy
      fdr_type <- "protein"
    }

    out2 <- dplyr::filter(out, prot_tier == 2L)
    out <- dplyr::filter(out, prot_tier == 1L)
  } else {
    if (fdr_type == "protein") {
      out2 <- dplyr::filter(out, prot_tier == 2L)
      out <- dplyr::filter(out, prot_tier == 1L)
    } else {
      out2 <- out[0, ]
      out <- out
    }
  }

  out <- grp_prots(out, file.path(out_path, "temp1"))

  if (nrow(out2)) {
    out2 <- out2 %>% grp_prots(file.path(out_path, "temp2"))
  } else {
    out2 <- out[0, ]
  }

  out3 <- grp_prots(out3, file.path(out_path, "temp3"))

  # Cleanup
  out <- dplyr::bind_cols(
    out %>% .[grepl("^prot_", names(.))],
    out %>% .[grepl("^pep_", names(.))],
    out %>% .[grepl("^psm_", names(.))],
    out %>% .[!grepl("^prot_|^pep_|^psm_", names(.))],
  ) %>%
    reloc_col_after("prot_es", "prot_family_member") %>%
    reloc_col_after("prot_es_co", "prot_es") %>%
    reloc_col_after("prot_tier", "prot_isess")

  # Three-tier combines
  max <- max(out$prot_hit_num, na.rm = TRUE)

  if (fdr_type == "protein") {
    if (combine_tier_three) {
      warning("Coerce to `combine_tier_three = FALSE` at `fdr_type = protein`.",
              call. = FALSE)
      combine_tier_three <- FALSE
    }
  }

  if (combine_tier_three) {
    out <- list(out, out2, out3) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(prot_acc, pep_seq) %T>%
      readr::write_tsv(file.path(out_path, "psmQ.txt"))
  } else {
    out <- out %>%
      dplyr::arrange(prot_acc, pep_seq) %T>%
      readr::write_tsv(file.path(out_path, "psmQ.txt"))

    if (nrow(out2)) {
      out2 <- out2[names(out)] %>%
        dplyr::mutate(prot_hit_num = prot_hit_num + max)  %T>%
        readr::write_tsv(file.path(out_path, "psmT2.txt"))

      max <- max(out2$prot_hit_num, na.rm = TRUE)
    }

    out3 <- out3[names(out)] %>%
      dplyr::mutate(prot_hit_num = prot_hit_num + max)  %T>%
      readr::write_tsv(file.path(out_path, "psmT3.txt"))
  }

  invisible(out)
}


