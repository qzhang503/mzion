#' Searches for MS ions.
#'
#' Database searches of MSMS data.
#'
#' The annotation of protein attributes, including percent coverages, will be
#' performed with \link[proteoQ]{normPSM} given that values will be affected
#' with the combination of multiple PSM tables.
#'
#' The search engine does not assume that variable peptides are descendants of
#' fixed peptides. In other words, each combination of variable and fixed
#' modifications is a set of \emph{realization} and applied freshly for searches
#' against all possible candidate sequences. Cares in search space were taken by
#' restricting only possible candidates at a given realization.
#'
#' The search is a two-way match: (a) a forward matching of theoretical values
#' to experiment ones and (b) a backward matching of the experimental values to
#' the theoretical ones. This allows the establishment of one-to-one
#' correspondences between experiments and theoreticals. The correspondences are
#' made available to users in files of \code{ion_matches_1.rds} etc., and, e.g.,
#' are used in \link{mapMS2ions} for rapid retrievals and visualizations of MS2
#' ion ladders.
#'
#' @param out_path A file path of outputs.
#' @param mgf_path A file path to a list of MGF files. The experimenter needs to
#'   supply the files.
#'
#'   The supported MGFs are in the formats of (1) MSConvert against \code{.raw}
#'   from Thermo's Orbitrap or \code{.d} from Bruker's timsTOF Pro, (2) Thermo's
#'   Proteome Discoverer or (3) Bruker's DataAnalysis.
#'
#'   With MSConvert, the default \code{titleMaker} is required for correct
#'   parsing (don't think it can be altered by users, but just in case).
#'
#'   Individuality in MGF files are slightly preferred to take advantage of
#'   parallel reading of the files.
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. The experimenter needs to supply the files.
#' @param acc_type Character string(s); the types of protein accessions in one
#'   of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
#'   the corresponding regular expressions need to be supplied via argument
#'   \code{acc_pattern}.
#' @param acc_pattern Regular expression(s) describing the patterns to separate
#'   the header lines of fasta entries. At the \code{NULL} default, the pattern
#'   will be automated when \code{acc_type} are among c("uniprot_acc",
#'   "uniprot_id", "refseq_acc", "other"). See also \link{load_fasta2} for
#'   custom examples.
#' @param fixedmods Character string(s) of fixed modifications.
#'
#' @param varmods Character string(s) of variable modifications. Multiple
#'   modifications to the same residue are allowed, for example, both
#'   \code{Carbamyl (M)} and \code{Oxidation (M)}.
#'
#'   For both \code{fixedmods} and \code{varmods}, the modification title,
#'   \code{TMT6plex}, applies to all of TMT-6, TMT-10, TMT-11. It is also
#'   possible to use aliased: (1) \code{TMT10plex} for TMT-10 and
#'   \code{TMT11plex} for TMT-11, (2) \code{TMT16plex} for TMTpro and (3)
#'   \code{TMT18plex} for TMTpro18. See also \link{parse_unimod} for grammars of
#'   modification \code{title}, \code{position} and \code{site}.
#' @param exclude_phospho_nl Logical; if TRUE, excludes neutral losses in the
#'   MS2 ion searches against variable modifications of phospho sites. The
#'   default if TRUE. May toggle it to \code{FALSE} if search speed against
#'   phospho data is not a primary concern.
#' @param include_insource_nl Not yet used. Logical; if TRUE, includes MS1
#'   precursor masses with the losses of neutral species prior to MS2
#'   fragmentation. The default is FALSE. The setting at TRUE remains
#'   experimenting by allowing additional masses in the universe of MS1
#'   precursors.
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from proteins. The
#'   default is \code{Trypsin_P}.
#'
#' @param custom_enzyme Regular expression(s) for custom enzyme specificity. The
#'   default is NULL. Uses of custom enzyme specificity is probably rather
#'   infrequent. Should there be such need, the argument \code{enzyme} will be
#'   ignored and the following may be applied:
#'
#'   \cr ## Examples \cr \cr # Equivalent to Trypsin \cr # at the Cterm of K or
#'   R but not followed by P \cr # (the quantifiers "\{1\}" can be skipped at a
#'   small cost of speed) \cr custom_enzyme = c(Cterm =
#'   "([KR]\{1\})([^P]\{1\})")
#'
#'   \cr # GluN again \cr custom_enzyme = c(Nterm = "([E]\{1\})")
#'
#'   \cr # Trypsin_P + GluN \cr custom_enzyme = c(Cterm = "([KR]\{1\})", Nterm =
#'   "([E]\{1\})")
#'
#'   \cr # Faked: Trypsin, proline not allowed on neither Nterm or Cterm \cr
#'   custom_enzyme = c(Cterm = "([KR]\{1\})([^P]\{1\})", Nterm =
#'   "([^P]\{1\})([KR]\{1\})")
#'
#' @param noenzyme_maxn Non-negative integer; the maximum number of peptide
#'   lengths for sectional searches at \code{noenzyme} specificity. The argument
#'   may be used to guard against RAM exhaustion. At the zero default, The
#'   peptide lengths from \code{min_len} to \code{max_len} will be broken
#'   automatically into continuous sections. At value 1, searches will be
#'   performed against individual peptide lengths; at value 2, two lengths will
#'   be taken at a time, etc.
#' @param maxn_fasta_seqs A positive integer; the maximum number of protein
#'   sequences in fasta files. The default is 200000.
#' @param maxn_vmods_setscombi A non-negative integer; the maximum number of
#'   sets of combinatorial variable modifications. The default is 64.
#' @param maxn_vmods_per_pep A non-negative integer; the maximum number of
#'   \code{Anywhere} (non-terminal) variable modifications per peptide. The
#'   default is 5.
#' @param maxn_sites_per_vmod A non-negative integer; the maximum number of
#'   combinatorial \code{Anywhere} (non-terminal) variable modifications per
#'   site in a peptide sequence. The default is 3.
#'
#'   For instance, variable modifications of \code{Carbamyl (M)} and
#'   \code{Oxidation (M)} both have site \code{M}. In order to have a
#'   combination of two \code{Carbamyl (M)} and two \code{Oxidation (M)} being
#'   considered, the value of \code{maxn_sites_per_vmod} needs to be four or
#'   greater.
#' @param maxn_vmods_sitescombi_per_pep A non-negative integer; the maximum
#'   number of combinatorial variable modifications per peptide sequence. The
#'   default is 64.
#' @param min_len A positive integer; the minimum length of peptide sequences
#'   for considerations. Shorter peptides will be excluded. The default is 7.
#' @param max_len A positive integer; the maximum length of peptide sequences
#'   for considerations. Longer peptides will be excluded. The default is 40.
#' @param max_miss A non-negative integer; the maximum number of mis-cleavages
#'   per peptide sequence for considerations. The default is 2.
#' @param min_mass A positive integer; the minimum precursor mass for
#'   interrogation. The default is 700.
#' @param max_mass A positive integer; the maximum precursor mass for
#'   interrogation. The default is 4500.
#' @param min_ms2mass A positive integer; the minimum MS2 mass for
#'   interrogation. The default is 110.
#' @param n_13c A non-negative integer; the maximum number of 13C off-sets for
#'   consideration in MS1 masses. The default is 0 with no off-sets.
#'   Peak-pickings by various MGF conversion tools may have attempted to adjust
#'   precursor masses to the corresponding mono-isotopic masses in isotope
#'   envelopes. Nevertheless, by setting \code{n_13c = 1}, some increases in the
#'   number of PSMs may be readily achieved at a relatively small cost of search
#'   time.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @param topn_ms2ions A positive integer; the top-n species for uses in MS2 ion
#'   searches. The default is to use the top-100 ions in an MS2 event.
#' @param minn_ms2 A positive integer; the minimum number of matched MS2 ions
#'   for consideration as a hit. The default is 6.
#' @param min_ms1_charge A positive integer; the minimum MS1 charge state for
#'   considerations. The default is 2.
#' @param max_ms1_charge A positive integer; the maximum MS1 charge state for
#'   considerations. The default is 6.
#' @param min_scan_num A positive integer; the minimum scan number for
#'   considerations. The default is 1. The setting only applies to MGFs with
#'   numeric scan numbers.
#' @param max_scan_num A positive integer; the maximum scan number for
#'   considerations. The default is the maximum machine integer. The setting
#'   only applies to MGFs with numeric scan numbers.
#' @param min_ret_time A non-negative numeric; the minimum retention time in
#'   seconds for considerations. The default is 0.
#' @param max_ret_time A non-negative numeric; the maximum retention time in
#'   seconds for considerations. The default is \code{Inf}.
#' @param ppm_ms1 A positive integer; the mass tolerance of MS1 species. The
#'   default is 20.
#' @param ppm_ms2 A positive integer; the mass tolerance of MS2 species. The
#'   default is 25.
#' @param ppm_reporters A positive integer; the mass tolerance of MS2 reporter
#'   ions. The default is 10.
#' @param quant A character string; the quantitation method. The default is
#'   "none". Additional choices include \code{tmt6, tmt10, tmt11, tmt16 and
#'   tmt18}. For other multiplicities of \code{tmt}, use the compatible higher
#'   plexes. For example, apply \code{tmt16} for \code{tmt12} provided a set of
#'   12-plexes being constructed from a 16-plex TMTpro (7 * 13C + 2 * 15N). It
#'   is also possible that an experimenter may construct a \code{tmt12} from a
#'   18-plex TMTpro (8 *13C + 1 * 15N) and thus \code{quant = tmt18} is
#'   suitable.
#' @param target_fdr A numeric; the targeted false-discovery rate (FDR) at the
#'   levels of PSM, peptide or protein. The default is 0.01. See also argument
#'   \code{fdr_type}.
#' @param fdr_type A character string; the type of FDR control. The value is in
#'   one of c("psm", "peptide", "protein"). The current default is \code{psm}.
#'
#'   Note that \code{fdr_type = protein} is equivalent to \code{fdr_type =
#'   peptide} with the additional filtration of data at \code{prot_tier == 1}. A
#'   variant is to set \code{fdr_type = psm}, followed by a data filtration at
#'   \code{prot_tier == 1}.
#' @param max_pepscores_co A positive numeric; the upper limit in the cut-offs
#'   of peptide scores for discriminating significant and insignificant
#'   identities. For higher quality and data-driven thresholds, choose the
#'   default \code{max_pepscores_co = Inf}.
#' @param max_protscores_co A positive numeric; the upper limit in the cut-offs
#'   of protein scores for discriminating significant and insignificant
#'   identities.  For higher quality and data-driven thresholds, choose the
#'   default \code{max_protscores_co = Inf}.
#' @param match_pepfdr Logical; if TRUE, matches empirically the highest peptide
#'   probability (corresponding to the lowest score) cut-offs to the pre-defined
#'   level of \code{target_fdr}. The default if TRUE. Choose \code{match_pepfdr
#'   = FALSE} for higher data quality (at a price of fewer hits).
#' @param combine_tier_three Logical; if TRUE, combines search results at tiers
#'   1, 2 and 3 to the single output of \code{psmQ.txt}. The default is FALSE in
#'   that data will be segregated into the three quality tiers according to the
#'   choice of \code{fdr_type}. The (convenience) parameter matters since
#'   \href{http://github.com/qzhang503/proteoQ}{proteoQ} will only look for the
#'   inputs of \code{psmQ[...].txt}.
#'
#'   For instance, if the aim is to bypass the constraint by protein FDR and
#'   focus on PSMs that have met the cut-offs specified by \code{target_fdr}, an
#'   experimenter may set \code{combine_tier_three = TRUE} and hence pool all
#'   significant peptides in \code{psmQ.txt} for downstream proteoQ.
#'
#'   Tier-1: both proteins and peptides with scores above significance
#'   thresholds.
#'
#'   Tier-2: \eqn{\ge} 2 significant peptides but protein scores below
#'   significance thresholds.
#'
#'   Tier-3: one significant peptide per protein and protein scores below
#'   significance thresholds.
#'
#' @param max_n_prots A positive integer to threshold the maximum number of
#'   protein entries before coercing \code{fdr_type} from \code{psm} or
#'   \code{peptide} to \code{protein}. The argument has no effect if
#'   \code{fdr_type} is already \code{protein}. In general, there is no need to
#'   change the default.
#'
#'   Note that for memory efficiency proteins at tiers 1, 2 and 3 are grouped
#'   separately. Further note that there is no tier-2 proteins at
#'   \code{fdr_type} of \code{psm} or \code{peptide}. For very large data sets,
#'   a lower value of \code{max_n_prots} can be used to reduce the chance of
#'   memory exhaustion by setting aside some protein entries from tier 1 to 2.
#' @param use_ms1_cache Logical; at the TRUE default, use cached precursor
#'   masses.
#'
#'   Set \code{use_ms1_cache = TRUE} for reprocessing of data, e.g., from
#'   \code{fdr_type = psm} to \code{fdr_type = protein}.
#' @param .path_cache The file path of cached search parameters. At the NULL
#'   default, the path is \code{"~/proteoM/.MSearches/Cache/Calls/"}. The
#'   parameter is for the users' awareness of the underlying structure of file
#'   folders and the use of default is suggested. Occasionally experimenters may
#'   remove the file folder for disk space or under infrequent events of
#'   modified framework incurred by the developer.
#' @param .path_fasta The parent file path to the theoretical masses of MS1
#'   precursors. At the NULL default, the path is \code{gsub("(.*)\\.[^\\.]*$",
#'   "\\1", get("fasta", envir = environment())[1])}. The parameter is for the
#'   users' awareness of the structure of file folders and the use of default is
#'   suggested. Occasionally experimenters may remove the file folder for disk
#'   space or under infrequent events of modified framework incurred by the
#'   developer.
#' @param digits A non-negative integer; the number of decimal places to be
#'   used. The default is 4.
#' @param ... Not currently used.
#' @section \code{FASTA}: \link{load_fasta2} sets the values of \code{acc_type}
#'   and \code{acc_pattern}. \cr
#' @section \code{Unimod}: \link{table_unimods} summarizes
#'   \href{https://www.unimod.org/}{Unimod} into a table format. \cr\cr
#'   \link{find_unimod} finds the mono-isotopic mass, position, site and neutral
#'   losses of a modification \cr\cr \link{parse_unimod} parses a Unimod.
#'   \href{https://proteoq.netlify.app/post/mixing-data-at-different-tmt-plexes/}{For
#'    example}, the name tag of "TMT6plex" is common among TMT-6, -10 and -11
#'   while "TMTpro" is for TMT-16 and "TMTpro18" for TMT-18. Experimenters may
#'   use aliases of "TMT10plex", "TMT11plex", "TMT16plex" and "TMT18plex".\cr\cr
#'   \link{calc_unimod_compmass} calculates the composition masses of a Unimod
#'   \cr\cr \link{add_unimod} adds a Unimod entry. \cr\cr \link{remove_unimod}
#'   removes a Unimod entry \cr\cr \link{remove_unimod_title} removes a Unimod
#'   entry by title.
#' @section \code{Visualization}: \link{mapMS2ions} visualizes the MS2 ion
#'   ladders.
#' @section \code{mzTab}: \link{make_mztab} converts outputs from the proteoM ->
#'   proteoQ pipeline to mzTab files.
#' @return A list of complete PSMs in \code{psmC.txt}; a list of quality PSMs in
#'   \code{psmQ.txt}.
#' @examples
#' \donttest{
#' # A hypothetical example
#' # (see also https://github.com/qzhang503/proteoM)
#' matchMS(
#'   fasta    = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   max_miss = 2,
#'   quant    = "tmt10",
#'   fdr_type = "protein",
#'   out_path = "~/proteoM/examples",
#' )
#'
#' \dontrun{
#' # Hypothetical phosphopeptides and 16-plex TMTpro
#' matchMS(
#'   fasta     = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type  = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("TMTpro (N-term)", "TMTpro (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'                 "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   max_miss  = 2,
#'   maxn_vmods_sitescombi_per_pep
#'             = 32,
#'   quant     = "tmt16",
#'   fdr_type  = "protein",
#'   out_path  = "~/proteoM/examples",
#' )
#'
#' # Hypothetical phosphopeptides and 18-plex TMTpro
#' matchMS(
#'   fasta     = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type  = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("TMTpro18 (N-term)", "TMTpro18 (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'                 "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   max_miss  = 2,
#'   maxn_vmods_sitescombi_per_pep
#'             = 32,
#'   quant     = "tmt18",
#'   fdr_type  = "protein",
#'   out_path  = "~/proteoM/examples",
#' )
#'
#' # Hypothetical Bruker's PASEF
#' matchMS(
#'   fasta     = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                 "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type  = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)"),
#'   ppm_ms1   = 25,
#'   ppm_ms2   = 40,
#'   quant     = "none",
#'   fdr_type  = "protein",
#'   out_path  = "~/proteoM/examples_pasef",
#' )
#'
#' # An exemplary custom Unimod (Oxi+Carbamidomethyl)
#' add_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'            full_name   = "Oxidation and iodoacetamide derivative"),
#'            specificity = c(site        = "M",
#'                            position    = "Anywhere"),
#'            delta       = c(mono_mass   = "73.016379",
#'                            avge_mass   = "73.0507",
#'                            composition = "H(3) C(2) N O(2)"),
#'            neuloss     = c(mono_mass   = "63.998285",
#'                            avge_mass   = "64.1069",
#'                            composition = "H(4) C O S"))
#'
#' matchMS(
#'   fasta    = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'               "Oxi+Carbamidomethyl (M)"),
#'   max_miss = 2,
#'   quant    = "tmt10",
#'   fdr_type = "protein",
#'   out_path = "~/proteoM/examples",
#' )
#'
#' }
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
                     enzyme = c("Trypsin_P", "Trypsin", "LysC", "LysN", "ArgC", 
                                "LysC_P", "Chymotrypsin", "GluC", "GluN", 
                                "AspC", "AspN", "SemiTrypsin_P", "SemiTrypsin", 
                                "SemiLysC", "SemiLysN", "SemiArgC", 
                                "SemiLysC_P", "SemiChymotrypsin", "SemiGluC", 
                                "SemiGluN", "SemiAspC", "SemiAspN", "Noenzyme", 
                                "Nodigest"),
                     custom_enzyme = c(Cterm = NULL, Nterm = NULL), 
                     noenzyme_maxn = 0L, 
                     maxn_fasta_seqs = 200000L,
                     maxn_vmods_setscombi = 64L,
                     maxn_vmods_per_pep = 5L,
                     maxn_sites_per_vmod = 3L,
                     maxn_vmods_sitescombi_per_pep = 64L,
                     min_len = 7L, max_len = 40L, max_miss = 2L, 
                     min_mass = 700L, max_mass = 4500L, 
                     ppm_ms1 = 20L, 
                     n_13c = 0L, 

                     type_ms2ions = "by", 
                     min_ms2mass = 115L, 
                     minn_ms2 = 6L, 
                     ppm_ms2 = 25L, 
                     
                     ppm_reporters = 10L,
                     quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16", "tmt18"),
                     
                     target_fdr = 0.01,
                     fdr_type = c("psm", "peptide", "protein"),
                     max_pepscores_co = Inf, max_protscores_co = Inf, 
                     match_pepfdr = TRUE, 
                     
                     combine_tier_three = FALSE,
                     max_n_prots = 40000L, 
                     use_ms1_cache = TRUE, 
                     .path_cache = NULL, 
                     .path_fasta = NULL,
                     
                     topn_ms2ions = 100L, 
                     min_ms1_charge = 2L, max_ms1_charge = 6L, 
                     min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                     min_ret_time = 0, max_ret_time = Inf, 

                     digits = 4L, ...) 
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
  
  ## Developer's dots
  dots <- as.list(substitute(...()))

  ## Tentative; not for users
  if (include_insource_nl) 
    stop("Currently only supports `include_insource_nl = FALSE`.")
  
  ## Preparation
  # modifications
  fixedmods <- sort(fixedmods)
  varmods <- sort(varmods)
  
  # accession pattern
  if ((!is.null(acc_pattern)) && (acc_pattern == "")) 
    acc_pattern <- NULL
  
  # logical types 
  stopifnot(vapply(c(include_insource_nl, exclude_phospho_nl, match_pepfdr, 
                     combine_tier_three, use_ms1_cache), 
                   is.logical, logical(1L)))

  # numeric types 
  stopifnot(vapply(c(maxn_fasta_seqs, maxn_vmods_setscombi, maxn_vmods_per_pep, 
                     maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                     min_len, max_len, max_miss, topn_ms2ions, minn_ms2, 
                     min_mass, max_mass, min_ms2mass, n_13c, 
                     ppm_ms1, ppm_ms2, ppm_reporters, max_n_prots, digits, 
                     target_fdr, max_pepscores_co, max_protscores_co), 
                   is.numeric, logical(1L)))

  # (a) integers casting for parameter matching when calling cached)
  max_integer <- .Machine$integer.max

  if (is.infinite(max_len)) max_len <- max_integer
  if (is.infinite(maxn_fasta_seqs)) maxn_fasta_seqs <- max_integer
  if (is.infinite(maxn_vmods_setscombi)) maxn_vmods_setscombi <- max_integer
  if (is.infinite(maxn_vmods_per_pep)) maxn_vmods_per_pep <- max_integer
  if (is.infinite(maxn_sites_per_vmod)) maxn_sites_per_vmod <- max_integer
  if (is.infinite(maxn_vmods_sitescombi_per_pep)) maxn_vmods_sitescombi_per_pep <- max_integer
  if (is.infinite(max_miss)) max_miss <- max_integer
  if (is.infinite(max_mass)) max_mass <- max_integer
  if (is.infinite(max_n_prots)) max_n_prots <- max_integer
  if (is.infinite(topn_ms2ions)) topn_ms2ions <- max_integer
  if (is.infinite(max_scan_num)) max_scan_num <- max_integer
  if (is.infinite(max_ret_time)) max_ret_time <- max_integer

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
  noenzyme_maxn <- as.integer(noenzyme_maxn)
  ppm_ms1 <- as.integer(ppm_ms1)
  ppm_ms2 <- as.integer(ppm_ms2)
  ppm_reporters <- as.integer(ppm_reporters)
  max_n_prots <- as.integer(max_n_prots)
  min_ms1_charge <- as.integer(min_ms1_charge)
  max_ms1_charge <- as.integer(max_ms1_charge)
  min_scan_num <- as.integer(min_scan_num)
  max_scan_num <- as.integer(max_scan_num)
  digits <- as.integer(digits)
  
  stopifnot(min_len >= 1L, max_len >= min_len, max_miss <= 10L, minn_ms2 >= 1L, 
            min_mass >= 1L, max_mass >= min_mass, min_ms2mass >= 1L, 
            n_13c >= 0L, noenzyme_maxn >= 0L, 
            maxn_vmods_per_pep >= maxn_sites_per_vmod, max_n_prots > 1000L, 
            min_ms1_charge >= 1L, max_ms1_charge >= min_ms1_charge, 
            min_scan_num >= 1L, max_scan_num >= min_scan_num)
            

  # (b) doubles
  target_fdr <- as.double(target_fdr)
  target_fdr <- round(target_fdr, digits = 2L)
  
  if (target_fdr > .25) 
    stop("Choose a smaller `target_fdr`.", call. = FALSE)
  
  min_ret_time <- round(min_ret_time, digits = 2L)
  max_ret_time <- round(max_ret_time, digits = 2L)
  max_pepscores_co <- round(max_pepscores_co, digits = 2L)
  max_protscores_co <- round(max_protscores_co, digits = 2L)

  stopifnot(max_pepscores_co >= 0, max_protscores_co >= 0, 
            min_ret_time >= 0, max_ret_time >= min_ret_time)
  
  # enzyme
  oks <- eval(formals()[["enzyme"]])
  oks_lwr <- tolower(oks)
  enzyme <- substitute(enzyme)
  
  if (length(enzyme) > 1L) {
    enzyme <- oks[1]
    enzyme_lwr <- tolower(enzyme)
  }
  else {
    enzyme <- as.character(enzyme)
    enzyme_lwr <- tolower(enzyme)
    
    if ((!enzyme_lwr %in% oks_lwr) && is.null(custom_enzyme)) 
      stop("Incorrect `enzyme = ", enzyme, "`.")
  }
  
  if (!is.null(custom_enzyme)) {
    warning("Overrule `enzyme` with `custom_enzyme`.", call. = FALSE)
    enzyme <- NULL
  }
  else
    enzyme <- enzyme_lwr
  
  if ((!is.null(enzyme)) && (enzyme == "noenzyme"))
    max_miss <- 0L
  
  if ((!is.null(enzyme)) && (enzyme == "nodigest")) {
    max_miss <- 0L
    max_len <- max_integer
  }

  rm(list = c("oks", "oks_lwr", "enzyme_lwr"))
  
  # fdr_type
  oks <- eval(formals()[["fdr_type"]])
  fdr_type <- substitute(fdr_type)
  
  if (length(fdr_type) > 1L)
    fdr_type <- oks[1]
  else {
    fdr_type <- as.character(fdr_type) 
    
    if (!fdr_type %in% oks)
      stop("Incorrect `fdr_type`.")
  }
  
  rm(list = c("oks"))

  # Quantitation method
  oks <- eval(formals()[["quant"]])
  quant <- substitute(quant)
  
  if (length(quant) > 1L)
    quant <- oks[1]
  else {
    quant <- as.character(quant) 
    
    if (!quant %in% oks)
      stop("Incorrect `quant`.")
  }
  
  rm(list = c("oks"))

  # TMT
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
  
  
  ## noenzyme specificity
  # set dots$recalled <- TRUE in `matchMS_noenzyme` and thus 
  # bypassing when calling matchMS again from `matchMS_noenzyme`
  recalled <- if (isTRUE(dots$recalled)) FALSE else TRUE
  
  if (enzyme == "noenzyme" && recalled) 
    matchMS_noenzyme(match.call(), min_len, max_len, fasta, out_path, mgf_path, 
                     noenzyme_maxn)
  
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
    custom_enzyme = custom_enzyme, 
    maxn_fasta_seqs = maxn_fasta_seqs,
    maxn_vmods_setscombi = maxn_vmods_setscombi,
    maxn_vmods_per_pep = maxn_vmods_per_pep,
    maxn_sites_per_vmod = maxn_sites_per_vmod,
    min_len = min_len,
    max_len = max_len,
    max_miss = max_miss,
    min_mass = min_mass, 
    max_mass = max_mass, 
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
            min_ms1_charge = min_ms1_charge, 
            max_ms1_charge = max_ms1_charge, 
            min_scan_num = min_scan_num, 
            max_scan_num = max_scan_num, 
            min_ret_time = min_ret_time,
            max_ret_time = max_ret_time, 
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
                     match_pepfdr = match_pepfdr, 
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
                    combine_tier_three = combine_tier_three, 
                    max_n_prots = max_n_prots)

  local({
    session_info <- sessionInfo()
    save(session_info, file = file.path(out_path, "Calls", "proteoM.rda"))
  })

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
#' @inheritParams ms2match
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
                        combine_tier_three = FALSE, max_n_prots = 40000L) 
{
  n_peps <- length(unique(out$pep_seq))
  n_prots <- length(unique(out$prot_acc))

  # `n_peps` and `n_prots` including both targets and decoys:
  # `n_prots` about 1:1
  # `n_peps` about 1.8:1

  if (n_prots == 1L) {
    message("No protein groups with the number of of proteins = ", n_prots, ".\n",
            "Search completed successfully.")
    options(show.error.messages = FALSE)
    stop()
  }
    
  if (n_peps > 1000000L && n_prots > 100000L) {
    out <- NA
  } 
  else {
    out <- tryCatch(
      psmC2Q(out,
             out_path = out_path,
             fdr_type = fdr_type,
             combine_tier_three = combine_tier_three, 
             max_n_prots = max_n_prots),
      error = function(e) NA
    )
  }

  if (length(out) == 1L && is.na(out)) {
    message("Retry with a new R session: \n\n",
            "proteoM:::reproc_psmC(\n",
            "  out_path = \"", out_path, "\",\n",
            "  fdr_type = \"", fdr_type, "\",\n",
            "  combine_tier_three  = ", combine_tier_three, ",\n",
            "  max_n_prots  = ", max_n_prots, "\n",
            ")")

    fileConn <- file(file.path("~/post_psmC.R"))

    lines <- c(
      "library(proteoM)\n",
      "proteoM:::reproc_psmC(",
      paste0("  out_path = \"", out_path, "\","),
      paste0("  fdr_type = \"", fdr_type, "\","),
      paste0("  combine_tier_three = ", combine_tier_three, ","),
      paste0("  max_n_prots = ", max_n_prots),
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
                         combine_tier_three = FALSE, max_n_prots = 40000L) 
{
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.", call. = FALSE)

  message("Leave the session open and wait for the `Search completed` message.")
  
  df <- readr::read_tsv(file.path(out_path, "psmC.txt"), show_col_types = FALSE)
  
  psmC2Q(df, out_path = out_path,
         fdr_type = fdr_type,
         combine_tier_three = combine_tier_three, 
         max_n_prots = max_n_prots)

  message("Done.")
}


#' From \code{psmC.txt} to \code{psmQ.txt}.
#'
#' @param out A result of \code{psmC.txt}.
#' @inheritParams matchMS
psmC2Q <- function (out = NULL, out_path = NULL, fdr_type = "protein",
                    combine_tier_three = FALSE, max_n_prots = 40000L) 
{
  options(warn = 1L)
  
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
  
  len_prots <- length(unique(out$prot_acc))

  if (len_prots > max_n_prots && fdr_type != "protein") {
    warning("Large number of proteins at ", len_prots, ".\n", 
            "Coerce to `fdr_type = protein` ",
            "and save peptide results of tier-2 proteins in `psmT2.txt`.",
            call. = FALSE)
    
    fdr_type <- "protein"

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
      
      if (len_prots > max_n_prots) {
        warning("The number of proteins is ", len_prots, ".\n", 
                "Consider `fdr_type = protein`.",
                call. = FALSE)
      }

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

  if (fdr_type == "protein" && combine_tier_three) {
    warning("Coerce to `combine_tier_three = FALSE` at `fdr_type = protein`.",
            call. = FALSE)
    combine_tier_three <- FALSE
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
      ok <- all(grepl("TMTpro18.* |TMT18plex.* ", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMTpro18` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    } 
    else if (quant == "tmt16") {
      ok <- all(grepl("TMTpro.* |TMT16plex.* ", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMTpro` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    } 
    else {
      ok <- all(grepl("TMT6plex.* |TMT10plex.* |TMT11plex.* ", possibles))
      
      if (!ok) 
        stop("All TMT modifications need to be `TMT6plex` at `", quant, "`.\n", 
             tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
             call. = FALSE)
    }
  }
  
  invisible(NULL)
}


#' Noenzyme search.
#' 
#' @param this_call An expression from match.call.
#' @inheritParams matchMS
matchMS_noenzyme <- function (this_call = NULL, min_len = 7L, max_len = 40L, 
                              fasta = NULL, out_path = NULL, mgf_path = NULL, 
                              noenzyme_maxn = 0L) 
{
  message("Searches with no enzyme specificity...")
  
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
    
    # ^1.5, 0.6: 90% RAM aggressiveness with fasta human + mouse
    max(1L, floor(fct_mem/(fct_fasta^1.5) * .5))
    # max(1L, floor(fct_mem/(fct_fasta^1.5) * .2))
  })
  
  len <- length(min_len:max_len)
  
  if (len > size) {
    if (size == 1L) {
      n_chunks <- len
      spans <- split(min_len:max_len, 1:len)
    }
    else {
      n_chunks <- ceiling(len/size)
      spans <- chunksplit(min_len:max_len, n_chunks, rightmost.closed = TRUE)
    }

    out_paths <- vector("list", n_chunks)
    
    for (i in seq_len(n_chunks)) {
      sub_call <- this_call
      span <- spans[[i]]
      start <- span[1]
      end <- span[length(span)]
      sub_nm <- paste0("sub", i, "_", start, "_", end)
      sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))
      
      ok <- file.exists(file.path(sub_path, "psmQ.txt"))
      if (ok) next
      
      sub_call$min_len <- start
      sub_call$max_len <- end
      sub_call$out_path <- sub_path
      sub_call$mgf_path <- mgf_path
      sub_call$recalled <- TRUE

      eval(sub_call, envir = environment())
      gc()
    }
    
    # psmC
    df <- lapply(out_paths, function (x) {
      file <- file.path(x, "psmC.txt")
      
      if (file.exists(file))
        readr::read_tsv(file, show_col_types = FALSE)
      else
        NULL
    })
    
    df <- dplyr::bind_rows(df)
    readr::write_tsv(df, file.path(out_path, "psmC.txt"))
    
    # psmQ
    df <- lapply(out_paths, function (x) {
      file <- file.path(x, "psmQ.txt")
      
      if (file.exists(file))
        readr::read_tsv(file, show_col_types = FALSE)
      else
        NULL
    })
    
    df <- dplyr::bind_rows(df)
    readr::write_tsv(df, file.path(out_path, "psmQ.txt"))
    
    # psmT2
    df <- lapply(out_paths, function (x) {
      file <- file.path(x, "psmT2.txt")
      
      if (file.exists(file))
        readr::read_tsv(file, show_col_types = FALSE)
      else
        NULL
    })
    
    df <- dplyr::bind_rows(df)
    readr::write_tsv(df, file.path(out_path, "psmT2.txt"))
    
    # psmT3
    df <- lapply(out_paths, function (x) {
      file <- file.path(x, "psmT3.txt")
      
      if (file.exists(file))
        readr::read_tsv(file, show_col_types = FALSE)
      else
        NULL
    })
    
    df <- dplyr::bind_rows(df)
    readr::write_tsv(df, file.path(out_path, "psmT3.txt"))
    
    message("Done (noenzyme search).")
    options(show.error.messages = FALSE)
    stop()
  }
  
  invisible(NULL)
}

