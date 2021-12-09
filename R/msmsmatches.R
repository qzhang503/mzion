#' Searches for MS ions.
#'
#' Database searches of MSMS data.
#'
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
#' @param varmods A character vector of variable modifications. Multiple
#'   modifications to the same residue are allowed, for example, both a less
#'   common \code{Carbamyl (M)} and a common \code{Oxidation (M)}.
#' @param exclude_phospho_nl If TRUE, excludes neutral losses in the MS2 ion
#'   searches against variable modifications of phospho sites. May toggle it to
#'   \code{FALSE} if search speed against phospho data is not a primary concern.
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
#' @param maxn_vmods_per_pep The maximum number of \code{Anywhere}
#'   (non-terminal) variable modifications per peptide.
#' @param maxn_sites_per_vmod Integer; the maximum number of combinatorial
#'   \code{Anywhere} (non-terminal) variable modifications per site in a per
#'   peptide sequence.
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence. The default is
#'   64. May consider a smaller value, i.e. 32, when searching against
#'   phosphopeptide data.
#' @param min_len A positive integer; the minimum length of peptides. Shorter
#'   peptides will be excluded.
#' @param max_len A positive integer; the maximum length of peptides. Longer
#'   peptides will be excluded.
#' @param max_miss A non-negative integer; the maximum number of mis-cleavages
#'   per peptide sequence.
#' @param min_mass A non-negative integer; the minimum precursor mass for
#'   interrogation.
#' @param max_mass A non-negative integer; the maximum precursor mass for
#'   interrogation.
#' @param min_ms2mass A non-negative integer; the minimum MS2 mass for
#'   interrogation.
#' @param n_13c A non-negative integer; the maximum number of 13C off-sets for
#'   consideration in MS1 masses. The default is 0 with no off-sets.
#'   Peak-pickings by various MGF conversion tools may have attempted to adjust
#'   precursor masses to the corresponding mono-isotopic masses in isotope
#'   envelopes. Nevertheless, by setting \code{n_13c = 1}, another 1% increase
#'   in the number of PSM may be readily achieved at a relatively small cost of
#'   search time.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @param topn_ms2ions A non-negative integer; the top-n species for uses in MS2
#'   ion searches. The default is to use the top-100 ions in an MS2 event.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 A positive integer; the mass tolerance of MS1 species.
#' @param ppm_ms2 A positive integer; the mass tolerance of MS2 species.
#' @param ppm_reporters A positive integer; the mass tolerance of MS2 reporter
#'   ions.
#' @param quant A quantitation method. The default is "none". Additional choices
#'   include \code{tmt6, tmt10, tmt11, tmt16 and tmt18}. For other
#'   multiplicities of \code{tmt}, use the compatible higher plexes. For
#'   example, apply \code{tmt16} for \code{tmt12} provided a set of 12-plexes
#'   being constructed from a 16-plex TMTpro (7 * 13C + 2 * 15N). It is also
#'   possible that an experimenter may construct a \code{tmt12} from a 18-plex
#'   TMTpro (8 *13C + 1 * 15N) and thus \code{quant = tmt18} is suitable.
#' @param target_fdr Numeric; a targeted false-discovery rate (FDR) at the
#'   levels of PSM, peptide or protein. See also argument \code{fdr_type}.
#' @param fdr_type Character string; the type of FDR control. The value is in
#'   one of c("psm", "peptide", "protein"). Note that \code{fdr_type = protein}
#'   is equivalent to \code{fdr_type = peptide} with the additional filtration
#'   of data at \code{prot_tier == 1}. A variant is to set \code{fdr_type =
#'   psm}, followed by a data filtration at \code{prot_tier == 1}.
#' @param max_pepscores_co Numeric; the upper limit in the cut-offs of peptide
#'   scores for discriminating significant and insignificant identities. The
#'   default is \code{Inf} without any restriction. Experimenters might consider
#'   to relax the restriction by a define threshold, i.e.,
#'   \code{max_pepscores_co = 50}. A graphic summary of experimentally
#'   determined cut-offs can be found at
#'   \code{`out_path`/temp/pepscore_len.pdf}.
#' @param max_protscores_co Numeric; the upper limit in the cut-offs of protein
#'   scores for discriminating significant and insignificant identities. The
#'   default is \code{Inf} without any restriction. Experimenters might consider
#'   to relax the restriction by a define threshold, i.e.,
#'   \code{max_protscores_co = 50}. A graphic summary of experimentally
#'   determined cut-offs can be found at
#'   \code{`out_path`/temp/protein_score_co.pdf}.
#' @param combine_tier_three Logical; if TRUE, combines search results at tiers
#'   1, 2 and 3 to the single output of \code{psmQ.txt}. The default is FALSE in
#'   that data will be segregated into the three quality tiers according to the
#'   choice of \code{fdr_type}. The (convenience) parameter matters since
#'   \href{http://github.com/qzhang503/proteoQ}{proteoQ} will only look for the
#'   inputs of \code{psmQ[...].txt}.
#'
#'   For instance, if the aim is to bypass the constraint by protein FDR and
#'   focus on PSMs that have met the cut-offs specified by \code{target_fdr},
#'   experimenters may set \code{combine_tier_three = TRUE} and hence pool all
#'   significant peptides in \code{psmQ.txt} for downstream proteoQ.
#'
#'   Tier-1: both proteins and peptides with scores above significance
#'   thresholds.
#'
#'   Tier-2: \eqn{\ge} 2 significant peptides but protein scores below
#'   significance thresholds.
#'
#'   Tier-3: one significant peptide and protein scores below significance
#'   thresholds.
#'
#' @param use_ms1_cache Logical; if TRUE, use cached precursor masses.
#'
#'   Set \code{use_ms1_cache = TRUE} for reprocessing of data, e.g., from
#'   \code{fdr_type = psm} to \code{fdr_type = protein}.
#' @param .path_cache The file path of cached search parameters. At the NULL
#'   default, the path is \code{"~/proteoM/.MSearches/Cache/Calls/"}. The
#'   parameter is for users' awareness of the structure of file folders and the
#'   default is suggested. Occasionally experimenters may remove the file folder
#'   for disk space or (hopefully not) in events of (improper) structural change
#'   incurred by the developer.
#' @param .path_fasta The parent file path to the theoretical masses of MS1
#'   precursors. At the NULL default, the path is \code{gsub("(.*)\\.[^\\.]*$",
#'   "\\1", get("fasta", envir = environment())[1])}. The parameter is for
#'   users' awareness of the structure of file folders and the default is
#'   suggested. Occasionally experimenters may remove the file folder for disk
#'   space or (hopefully) in rare events of structural change incurred by the
#'   developer.
#' @param digits Integer; the number of decimal places to be used.
#' @seealso \link{load_fasta2} for setting the values of \code{acc_type} and
#'   \code{acc_pattern}. \cr \link{table_unimods} summarizes
#'   \href{https://www.unimod.org/}{Unimod} into a table format. \cr
#'   \link{parse_unimod} for the grammar of Unimod.
#'   \href{https://proteoq.netlify.app/post/mixing-data-at-different-tmt-plexes/}{For
#'    example}, the name tag of "TMT6plex" is common among TMT-6, -10 and -11
#'   while "TMTpro" is specific to TMT-16.
#' @return A list of complete PSMs in \code{psmC.txt}; a list of quality PSMs in
#'   \code{psmQ.txt}.
#' @examples
#' \donttest{
#' # A hypothetical example
#' # (see also https://github.com/qzhang503/proteoM)
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
#' # Hypothetical phosphopeptides and 16-plex TMTpro
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
#'   maxn_vmods_sitescombi_per_pep = 32,
#'   quant = "tmt16",
#'   fdr_type = "protein",
#'   out_path = "~/proteoM/examples",
#' )
#'
#' # Hypothetical phosphopeptides and 18-plex TMTpro
#' matchMS(
#'   fasta = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'             "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("TMTpro18 (N-term)", "TMTpro18 (K)", "Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'               "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'              "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   max_miss = 2,
#'   maxn_vmods_sitescombi_per_pep = 32,
#'   quant = "tmt18",
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
                     min_mass = 500L, max_mass = 6000L, ppm_ms1 = 20L, 
                     n_13c = 0L, 

                     type_ms2ions = "by", topn_ms2ions = 100L, 
                     min_ms2mass = 110L, minn_ms2 = 6L, 
                     ppm_ms2 = 25L, ppm_reporters = 10L,
                     quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16", "tmt18"),
                     
                     target_fdr = 0.01,
                     fdr_type = c("psm", "peptide", "protein"),
                     max_pepscores_co = Inf, max_protscores_co = Inf, 
                     
                     combine_tier_three = FALSE,
                     use_ms1_cache = TRUE, 
                     .path_cache = NULL, 
                     .path_fasta = NULL,
                     digits = 4L) 
{
  options(digits = 9L)

  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"),
                   fun = as.character(match.call()[[1]]))
      }
    },
    add = TRUE
  )

  ## Tentative; not for users
  if (include_insource_nl) 
    stop("Currently only supports `include_insource_nl = TRUE`.")
  
  ## Preparation
  # modifications
  fixedmods <- sort(fixedmods)
  varmods <- sort(varmods)
  
  # accession pattern
  if ((!is.null(acc_pattern)) && acc_pattern == "") 
    acc_pattern <- NULL
  
  # numeric types 
  stopifnot(vapply(c(maxn_fasta_seqs, maxn_vmods_setscombi, maxn_vmods_per_pep, 
                     maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                     min_len, max_len, max_miss, topn_ms2ions, minn_ms2, 
                     min_mass, max_mass, min_ms2mass, n_13c, 
                     ppm_ms1, ppm_ms2, ppm_reporters, digits, 
                     target_fdr, max_pepscores_co, max_protscores_co), 
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
  n_13c <- as.integer(n_13c)
  ppm_ms1 <- as.integer(ppm_ms1)
  ppm_ms2 <- as.integer(ppm_ms2)
  ppm_reporters <- as.integer(ppm_reporters)
  digits <- as.integer(digits)
  
  stopifnot(min_len >= 0L, max_len >= min_len, max_miss <= 10L, 
            min_mass >= 0L, max_mass >= min_mass, min_ms2mass >= 0L, 
            n_13c >= 0L, 
            maxn_vmods_per_pep >= maxn_sites_per_vmod)

  # (b) doubles
  target_fdr <- as.double(target_fdr)
  target_fdr <- round(target_fdr, digits = 2L)
  
  if (target_fdr > .25) 
    stop("Choose a smaller `target_fdr`.", call. = FALSE)
  
  max_pepscores_co <- round(max_pepscores_co, digits = 2L)
  max_protscores_co <- round(max_protscores_co, digits = 2L)
  
  stopifnot(max_pepscores_co >= 0, max_protscores_co >= 0)
  
  # fdr_type
  fdr_type <- rlang::enexpr(fdr_type)
  oks <- eval(formals()[["fdr_type"]])

  if (length(fdr_type) > 1L) 
    fdr_type <- oks[[1]]
  else 
    fdr_type <- rlang::as_string(fdr_type)

  stopifnot(fdr_type %in% oks)
  rm(list = c("oks"))

  # Quantitation method
  quant <- rlang::enexpr(quant)
  oks <- eval(formals()[["quant"]])

  quant <- if (length(quant) > 1L) 
    oks[[1]]
  else 
    rlang::as_string(quant)

  stopifnot(quant %in% oks, length(quant) == 1L)
  rm(list = c("oks"))
  
  check_tmt_pars(fixedmods, varmods, quant) 

  # Output path
  out_path <- create_dir(out_path)
  
  # MGF path
  mgf_path <- find_dir(mgf_path)
  
  if (is.null(mgf_path)) 
    stop("`mgf_path` not found.", call. = FALSE)

  filelist <- list.files(path = mgf_path, pattern = "\\.mgf$")

  if (!length(filelist)) 
    stop("No `.mgf` files under ", mgf_path, call. = FALSE)

  rm(list = c("filelist"))
  
  # system paths
  if (is.null(.path_cache)) {
    .path_cache <- "~/proteoM/.MSearches/Cache/Calls/"
  }
  .path_cache <- create_dir(.path_cache)
  
  if (is.null(.path_fasta)) {
    .path_fasta <- file.path(gsub("(.*)\\.[^\\.]*$", "\\1", fasta[1]))
  } 
  .path_fasta <- create_dir(.path_fasta)
  
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
    n_13c = n_13c,
    out_path = out_path,
    digits = digits,
    use_ms1_cache = use_ms1_cache, 
    .path_cache = .path_cache, 
    .path_fasta = .path_fasta, 
    .path_ms1masses = .path_ms1masses
  )

  ## Bin theoretical peptides
  bin_ms1masses(res = res, 
                min_mass = min_mass, 
                max_mass = max_mass, 
                ppm_ms1 = ppm_ms1, 
                use_ms1_cache = use_ms1_cache, 
                .path_cache = .path_cache, 
                .path_ms1masses = .path_ms1masses)
  
  rm(list = c("res"))
  gc()

  ## MGFs
  load_mgfs(out_path = out_path, 
            mgf_path = mgf_path,
            min_mass = min_mass,
            max_mass = max_mass, 
            min_ms2mass = min_ms2mass,
            topn_ms2ions = topn_ms2ions,
            ppm_ms1 = ppm_ms1,
            ppm_ms2 = ppm_ms2,
            index_ms2 = FALSE)

  ## MSMS matches
  ms2match(mgf_path = mgf_path,
           aa_masses_all = 
             readRDS(file.path(.path_ms1masses, .time_stamp, "aa_masses_all.rds")),
           out_path = out_path,
           mod_indexes = 
             find_mod_indexes(file.path(.path_ms1masses, .time_stamp, "mod_indexes.txt")),
           type_ms2ions = type_ms2ions,
           maxn_vmods_per_pep = maxn_vmods_per_pep,
           maxn_sites_per_vmod = maxn_sites_per_vmod,
           maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep,
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
  calc_pepscores(topn_ms2ions = topn_ms2ions,
                 type_ms2ions = type_ms2ions,
                 target_fdr = target_fdr,
                 fdr_type = fdr_type,
                 min_len = min_len,
                 max_len = max_len,
                 penalize_sions = TRUE,
                 ppm_ms2 = ppm_ms2,
                 out_path = out_path,
                 
                 # dummies
                 mgf_path = mgf_path,
                 maxn_vmods_per_pep = maxn_vmods_per_pep,
                 maxn_sites_per_vmod = maxn_sites_per_vmod,
                 maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep,
                 minn_ms2 = minn_ms2,
                 ppm_ms1 = ppm_ms1,
                 min_ms2mass = min_ms2mass,
                 quant = quant,
                 ppm_reporters = ppm_reporters,
                 fasta = fasta,
                 acc_type = acc_type,
                 acc_pattern = acc_pattern,
                 fixedmods = fixedmods,
                 varmods = varmods,
                 include_insource_nl = include_insource_nl,
                 enzyme = enzyme,
                 maxn_fasta_seqs = maxn_fasta_seqs,
                 maxn_vmods_setscombi = maxn_vmods_setscombi,
                 
                 digits = digits)

  ## Peptide FDR 
  out <- calc_pepfdr(target_fdr = target_fdr, 
                     fdr_type = fdr_type, 
                     min_len = min_len, 
                     max_len = max_len, 
                     max_pepscores_co = max_pepscores_co, 
                     out_path = out_path) %>% 
    post_pepfdr(out_path)

  ## Peptide ranks and score deltas between `pep_ivmod`
  out <- calc_peploc(out)
  gc()

  ## Protein accessions, score cut-offs and optional reporter ions
  out <- add_prot_acc(out, out_path)
  out <- calc_protfdr(df = out, 
                      target_fdr = target_fdr, 
                      max_protscores_co = max_protscores_co, 
                      out_path = out_path)
  out <- add_rptrs(out, quant, out_path)
  gc()

  ## Clean-ups
  out <- local({
    raws <- readRDS(file.path(mgf_path, "raw_indexes.rds"))
    raws2 <- names(raws)
    names(raws2) <- raws
    out$raw_file <- unname(raws2[out$raw_file])
    
    scans <- readRDS(file.path(mgf_path, "scan_indexes.rds"))
    scans2 <- names(scans)
    names(scans2) <- scans
    out$scan_title <- unname(scans2[out$scan_title])
    
    out
  })
  
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
#' Not yet used.
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
                          digits) 
{
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
    
    rstudioapi::restartSession(command = 'source("~/matchMS.R")')
  } 
  else {
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
                        combine_tier_three = FALSE) 
{
  n_peps <- length(unique(out$pep_seq))
  n_prots <- length(unique(out$prot_acc))

  # `n_peps` and `n_prots` including both targets and decoys:
  # `n_prots` about 1:1
  # `n_peps` about 1.8:1

  if (n_peps > 1000000L && n_prots > 100000L) {
    out <- NA
  } 
  else {
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
  } 
  else {
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
                         combine_tier_three = FALSE) 
{
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.", call. = FALSE)

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
                    combine_tier_three = FALSE) 
{
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

      fdr_type <- "protein"
    }

    out2 <- dplyr::filter(out, prot_tier == 2L)
    out <- dplyr::filter(out, prot_tier == 1L)
  } 
  else {
    if (fdr_type == "protein") {
      out2 <- dplyr::filter(out, prot_tier == 2L)
      out <- dplyr::filter(out, prot_tier == 1L)
    } 
    else {
      message("No tier-2 outputs at `fdr_type = ", fdr_type, "`.")
      out2 <- out[0, ]
      out <- out
    }
  }

  out <- grp_prots(out, file.path(out_path, "temp1"))

  if (nrow(out2)) 
    out2 <- grp_prots(out2, file.path(out_path, "temp2"))
  else 
    out2 <- out[0, ]
  
  if (nrow(out3)) 
    out3 <- grp_prots(out3, file.path(out_path, "temp3"))
  else 
    out3 <- out[0, ]

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
    
    local({
      file_t2 <- file.path(out_path, "psmT2.txt")
      file_t3 <- file.path(out_path, "psmT3.txt")
      
      if (file.exists(file_t2)) {
        message("Delete `psmT2.txt` at `combine_tier_three = TRUE`.")
        unlink(file_t2)
      }

      if (file.exists(file_t3)) {
        message("Delete `psmT3.txt` at `combine_tier_three = TRUE`.")
        unlink(file_t3)
      }
    })
  } 
  else {
    out <- out %>%
      dplyr::arrange(prot_acc, pep_seq) %T>%
      readr::write_tsv(file.path(out_path, "psmQ.txt"))

    if (nrow(out2)) {
      out2 <- out2[names(out)] %>%
        dplyr::mutate(prot_hit_num = prot_hit_num + max)  %T>%
        readr::write_tsv(file.path(out_path, "psmT2.txt"))

      max <- max(out2$prot_hit_num, na.rm = TRUE)
    }

    if (nrow(out3)) {
      out3 <- out3[names(out)] %>%
        dplyr::mutate(prot_hit_num = prot_hit_num + max)  %T>%
        readr::write_tsv(file.path(out_path, "psmT3.txt"))
    }
  }

  invisible(out)
}


#' Checks the compatibility of TMT names and plexes.
#' 
#' @inheritParams matchMS
check_tmt_pars <- function (fixedmods, varmods, quant) 
{
  if (TRUE) {
    # mono-isotopic
    H <- 1.007825035
    O <- 15.99491463
    C <- 12
    N <- 14.003074
    C13 <- 13.00335483
    N15 <- 15.00010897
    
    # average
    H_a <- 1.00794
    O_a <- 15.9994
    C_a <- 12.0107
    C13_a <- C13
    N_a <- 14.0067
    N15_a <- N15
    
    # TMTpro-16
    H * 25 + C * 8 + C13 * 7 + N * 1 + N15 * 2 + O * 3 # 304.207146
    H_a * 25 + C_a * 8 + C13_a * 7 + N_a * 1 + N15_a * 2 + O_a * 3 # 304.312702
    
    # TMTpro-18
    H*(25) + C*(7) + C13*(8) + N*(2) + N15*(1) + O*(3) # 304.213465
    H_a*(25) + C_a*(7) + C13_a*(8) + N_a*(2) + N15_a*(1) + O_a*(3) # 304.311948
    
    # TMTpro-zero
    H*(25) + C*(8) + C*(7) + N + N*(2) + O*(3) # 295.189592
    H_a*(25) + C_a*(8) + C_a*(7) + N_a*(1) + N_a*(2) + O_a*(3) # 295.3773
  }
  
  tmt_msg_1 <- "*** TMT6plex for tmt6, tmt10, tmt11 ***"
  tmt_msg_2 <- "*** TMTpro for tmt16 ***"
  tmt_msg_3 <- "*** TMTpro18 for tmt18 ***"
  
  fvmods <- c(fixedmods, varmods)
  
  if (grepl("^tmt[0-9]+", quant)) {
    possibles <- fvmods[grepl("^TMT", fvmods)]
    
    if (quant == "tmt18") {
      ok <- all(grepl("TMTpro18 ", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMTpro18` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    } else if (quant == "tmt16") {
      ok <- all(grepl("TMTpro ", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMTpro` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    } else {
      ok <- all(grepl("TMT6plex", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMT6plex` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    }
  }
  
  invisible(NULL)
}
