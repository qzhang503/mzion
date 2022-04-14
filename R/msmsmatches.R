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
#' made available to users in files of \code{ion_matches_1.rds...} (nested form)
#' and \code{list_table_1.rds...} etc. (flat form). A more self-contained output
#' can be made available at the TRUE of \code{add_ms2theos},
#' \code{add_ms2theos2}, \code{add_ms2moverzs} and \code{add_ms2ints}.
#'
#' When there is no evidence to distinguish, e.g. distinct primary sequences of
#' \code{P[EMPTY]EPTIDE} and \code{P[MTYPE]EPTIDE}, both will be reported by
#' proteoM and further kept by proteoQ. For peptides under the same primary
#' sequence, the redundancy in the positions and/or neutral losses of
#' \code{Anywhere} variable modifications are also kept in the outputs of
#' proteoM but removed with proteoQ.
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
#' @param acc_pattern Regular expression(s) describing the patterns in
#'   separating the header lines of fasta entries. At the \code{NULL} default,
#'   the pattern will be automated when \code{acc_type} are among
#'   c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). See also
#'   \link{load_fasta2} for custom examples.
#' @param fixedmods Character string(s) of fixed modifications.
#'
#' @param varmods Character string(s) of variable modifications. Multiple
#'   modifications to the same residue are allowed, for example, both
#'   \code{Carbamyl (M)} and \code{Oxidation (M)}.
#'
#'   For both \code{fixedmods} and \code{varmods}, the modification title,
#'   \code{TMT6plex}, applies to all of TMT-6, TMT-10, TMT-11. It is also
#'   possible to use aliased: (1) \code{TMT10plex} for TMT-10 and
#'   \code{TMT11plex} for TMT-11 and (2) \code{TMT16plex} for TMTpro. See also
#'   \link{parse_unimod} for grammars of modification \code{title},
#'   \code{position} and \code{site}.
#' @param include_insource_nl Not yet used. Logical; if TRUE, includes MS1
#'   precursor masses with the losses of neutral species prior to MS2
#'   fragmentation. The default is FALSE. The setting at TRUE remains
#'   experimenting by allowing additional masses in the universe of MS1
#'   precursors.
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from protein entries. The
#'   default is \code{Trypsin_P}. See also parameter \code{custom_enzyme}.
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
#'   performed against individual peptide lengths; at value 2, two adjacent
#'   lengths will be taken at a time, etc.
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
#' @param par_groups Parameter(s) of \code{matchMS} multiplied by sets of values
#'   in groups. Multiple searches will be performed separately against the
#'   parameter groups. For instance with one set of samples in SILAC light and
#'   the other in SILAC heavy, the experimenters may specify two arguments for
#'   parameter \code{mgf_path} and two arguments for parameter \code{fixedmods}
#'   that link to the respective samples. In this way, there is no need to
#'   search against, e.g. heavy-isotope-labeled K8R10 with the light samples and
#'   vice versa. Note that results will be combined at the end, with the group
#'   names indicated under column \code{pep_group}. The default is NULL without
#'   grouped searches. See the examples under SILAC and Group searches.
#' @param silac_mix A list of labels indicating SILAC groups in samples. The
#'   parameter is most relevant for SILAC experiments where peptides of heavy,
#'   light etc. were \emph{mixed} into one sample. The default is NULL
#'   indicating a none mixed-SILAC experiment. See also the examples under
#'   SILAC.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @param topn_ms2ions A positive integer; the top-n species for uses in MS2 ion
#'   searches. The default is to use the top-100 ions in an MS2 event.
#' @param minn_ms2 A positive integer; the minimum number of matched MS2 ions
#'   for consideration as a hit. The default is 6. Counts of secondary ions,
#'   e.g. b0, b* etc., are not part of the threshold.
#' @param min_ms1_charge A positive integer; the minimum MS1 charge state for
#'   considerations. The default is 2.
#' @param max_ms1_charge A positive integer; the maximum MS1 charge state for
#'   considerations. The default is 6.
#' @param min_scan_num A positive integer; the minimum scan number for
#'   considerations. The default is 1. The setting only applies to MGFs with
#'   numeric scan numbers. For example, it has no effects on Bruker's timsTOF
#'   data.
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
#'   12-plexes being constructed from a 16-plex TMTpro \eqn{(7 * 13C + 2 *
#'   15N)}. It is also possible that an experimenter may construct a
#'   \code{tmt12} from a 18-plex TMTpro \eqn{(8 *13C + 1 * 15N)} where
#'   \code{quant = tmt18} is suitable.
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
#' @param topn_seqs_per_query Positive integer; a threshold to discard peptide
#'   matches under the same MS query with scores beyond the top-n. The default
#'   is 3.
#'
#'   The same \code{MS query} refers to the identity in \code{MS scan number}
#'   and \code{MS raw file name}. Target and decoys matches are treated
#'   separately.
#' @param topn_mods_per_seq Positive integer; a threshold to discard variable
#'   modifications under the same peptide match with scores beyond the top-n.
#'   The default is 3.
#'
#'   The same \code{peptide match} refers to matches with identities in \code{MS
#'   scan number}, \code{MS raw file name} and \code{peptide sequence}. Target
#'   and decoys matches are treated separately.
#'
#'   For a variable modification with multiple neutral losses (NL), the
#'   best-scored NL will be used in the ranking.
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
#' @param add_ms2theos Logical. If true, adds the sequence of primary
#'   theoretical MS2 m/z values (\code{pep_ms2_theos}). The sequence order at a
#'   given \code{type_ms2ions} is:
#'
#'   \tabular{ll}{ \strong{Type}   \tab \strong{Sequence}\cr by \tab \eqn{b1,
#'   b2..., y1, y2...} \cr ax \tab \eqn{a1, a2..., x1, x2...} \cr cz \tab
#'   \eqn{c1, c2..., z1, z2...} \cr }
#'
#' @param add_ms2theos2 Logical. If true, adds the sequence of secondary
#'   theoretical MS2 m/z values (\code{pep_ms2_theos2}). The sequence order at a
#'   given \code{type_ms2ions} is:
#'
#'   \tabular{ll}{ \strong{Type}   \tab \strong{Order of sequences}\cr by \tab
#'   \eqn{b2, b*, b*2, b0, b02, y2, y*, y*2, y0, y02} \cr ax \tab \eqn{a2, a*,
#'   a*2, a0, a02, x2} \cr cz \tab \eqn{c2, z2} \cr }
#'
#' @param add_ms2moverzs Logical; if TRUE, adds the sequence of experimental
#'   \eqn{m/z} values (\code{pep_ms2_moverzs}).
#' @param add_ms2ints Logical; if TRUE, adds the sequence of experimental MS2
#'   intensity values (\code{pep_ms2_ints}).
#' @param .path_cache The file path of cached search parameters. The parameter
#'   is for the users' awareness of the underlying structure of file folders and
#'   the use of default is suggested. Occasionally experimenters may remove the
#'   file folder for disk space or under infrequent events of modified framework
#'   incurred by the developer.
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
#'   while "TMTpro" is for TMT-16. Experimenters may use aliases of "TMT10plex",
#'   "TMT11plex" and "TMT16plex.\cr\cr \link{calc_unimod_compmass} calculates
#'   the composition masses of a Unimod \cr\cr \link{add_unimod} adds a Unimod
#'   entry. \cr\cr \link{remove_unimod} removes a Unimod entry \cr\cr
#'   \link{remove_unimod_title} removes a Unimod entry by title.
#' @section \code{Visualization}: \link{mapMS2ions} visualizes the MS2 ion
#'   ladders.
#' @section \code{mzTab}: \link{make_mztab} converts outputs from the proteoM ->
#'   proteoQ pipeline to mzTab files.
#' @section \code{Output columns}: system.file("extdata", "column_keys.txt",
#'   package = "proteoM") \cr
#' @return A list of complete PSMs in \code{psmC.txt}; a list of quality PSMs in
#'   \code{psmQ.txt}.
#' @examples
#' \donttest{
#' ## All examples are hypothetical
#' ## (some real ones at https://github.com/qzhang503/proteoM)
#' 
#' # TMT10
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
#' # TMT16, phospho
#' matchMS(
#'   fixedmods = c("TMTpro (N-term)", "TMTpro (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'                 "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   quant     = "tmt16",
#'   fdr_type  = "psm",
#'   combine_tier_three = TRUE,
#'   out_path  = "~/proteoM/examples",
#' )
#'
#' # Bruker's PASEF
#' matchMS(
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
#' # Custom Unimod (Oxi+Carbamidomethyl)
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
#'   fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'               "Oxi+Carbamidomethyl (M)"),
#'   quant    = "tmt10",
#' )
#'
#' 
#' ## Stable isotope-labeled K and R
#' # K8, R10
#' matchMS(
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'               "K8 (C-term = K)", "R10 (C-term = R)"),
#'   quant    = "none",
#' )
#' 
#' # TMT+K8, TMT+R10
#' matchMS(
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'               "TMT+K8 (C-term = K)", "TMT+R10 (C-term = R)"),
#'   max_miss = 2,
#'   quant    = "tmt10",
#' )
#' 
#' 
#' #######################################
#' # SILAC
#' #######################################
#'
#' ## 1. heavy and light mixed into one sample (typical SILAC)
#'
#' # unlabeled base
#' matchMS(
#'   silac_mix = list(base = NULL, heavy = c("K8 (K)", "R10 (R)")),
#'   ...
#' )
#'
#' # labeled base
#' # (first add K4 Unimod if not yet present)
#' K4 <- calc_unimod_compmass("2H(4) H(-4)")
#' mono_mass <- K4$mono_mass
#' avge_mass <- K4$avge_mass
#'
#' add_unimod(header      = c(title       = "K4",
#'                            full_name   = "Heavy lysine 2H(4) H(-4)"),
#'            specificity = c(site        = "K",
#'                            position    = "Anywhere"),
#'            delta       = c(mono_mass   = "4.025108",
#'                            avge_mass   = "4.02464",
#'                            composition = "2H(4) H(-4)"),
#'            neuloss     = c(mono_mass   = "0",
#'                            avge_mass   = "0",
#'                            composition = "0"))
#'
#' matchMS(
#'   silac_mix = list(base   = c("K4 (K)"),
#'                    median = c("K6 (K)", "R6 (R)"),
#'                    heavy  = c("K8 (K)", "R10 (R)")),
#'   ...
#' )
#'
#' # custom groups
#' matchMS(
#'   silac_mix = list(base = NULL,
#'                    grp1 = c("K6 (K)", "R6 (R)"),
#'                    grp2 = c("K8 (K)", "R10 (R)")),
#'   ...
#' )
#'
#' matchMS(
#'   silac_mix = list(base = c("K6 (K)", "R6 (R)"),
#'                    grp1 = c("K8 (K)", "R10 (R)"),
#'   ...
#' )
#'
#'
#' ## 2. Heavy and light in separate samples
#' #  (toy examples assessing the technical quality of SILAC)
#'
#' # MGFs of light and heavy samples under separate folders;
#' # Heavy modifications being fixedmods
#'
#' # (i) SILAC but low throughput since no sample mixing
#' matchMS(
#'   par_groups = list(
#'     light = list(mgf_path  = "~/proteoM/my_project/mgf/light",
#'                  fixedmods = "Carbamidomethyl (C)"),
#'     heavy = list(mgf_path  = "~/proteoM/my_project/mgf/heavy",
#'                  fixedmods = c("Carbamidomethyl (C)", "K8 (K)", "R10 (R)"))
#'   ),
#'   quant = "none",
#'   ...
#' )
#'
#' # The results next processed by proteoQ just like LFQ.
#'
#' # (ii) SILAC at low-throughput + TMT
#' matchMS(
#'   par_groups = list(
#'     light = list(mgf_path  = "~/proteoM/my_project/mgf/light",
#'                  fixedmods = "Carbamidomethyl (C)"),
#'     heavy = list(mgf_path  = "~/proteoM/my_project/mgf/heavy",
#'                  fixedmods = c("Carbamidomethyl (C)", "K8 (K)", "R10 (R)"))
#'   ),
#'   quant = "TMT10",
#'   ...
#' )
#'
#' # Next processed by proteoQ just like TMT
#' # ...
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
                     
                     par_groups = NULL, 
                     silac_mix = NULL, 

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
                     
                     topn_mods_per_seq = 3L, 
                     topn_seqs_per_query = 3L, 
                     
                     combine_tier_three = FALSE,
                     max_n_prots = 40000L, 
                     use_ms1_cache = TRUE, 
                     .path_cache = "~/proteoM/.MSearches (1.1.0.0)/Cache/Calls", 
                     .path_fasta = NULL,
                     
                     topn_ms2ions = 100L, 
                     min_ms1_charge = 2L, max_ms1_charge = 6L, 
                     min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                     min_ret_time = 0, max_ret_time = Inf, 
                     
                     add_ms2theos = FALSE, add_ms2theos2 = FALSE, 
                     add_ms2moverzs = FALSE, add_ms2ints = FALSE,

                     digits = 4L, ...) 
{
  options(digits = 9L)

  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) {
        tryCatch(save_call2(path = file.path(out_path, "Calls"),
                            fun = fun), 
                 error = function(e) NA)
      }
    },
    add = TRUE
  )
  
  this_call <- match.call()
  fun <- as.character(this_call[1])
  
  suppressWarnings(
    rm(list = c(".path_cache", ".path_fasta", ".path_ms1masses", 
                ".time_stamp", ".time_bin", ".path_bin"), 
       envir = .GlobalEnv)
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
  stopifnot(vapply(c(include_insource_nl, match_pepfdr, combine_tier_three, 
                     use_ms1_cache, add_ms2theos, add_ms2theos2, add_ms2moverzs, 
                     add_ms2ints), 
                   is.logical, logical(1L)))

  # numeric types 
  stopifnot(vapply(c(maxn_fasta_seqs, maxn_vmods_setscombi, maxn_vmods_per_pep, 
                     maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                     min_len, max_len, max_miss, topn_ms2ions, minn_ms2, 
                     min_mass, max_mass, min_ms2mass, n_13c, 
                     ppm_ms1, ppm_ms2, ppm_reporters, max_n_prots, digits, 
                     target_fdr, max_pepscores_co, max_protscores_co, 
                     topn_mods_per_seq, topn_seqs_per_query), 
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
  if (is.infinite(topn_mods_per_seq)) topn_mods_per_seq <- max_integer
  if (is.infinite(topn_seqs_per_query)) topn_seqs_per_query <- max_integer

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
  topn_mods_per_seq <- as.integer(topn_mods_per_seq)
  topn_seqs_per_query <- as.integer(topn_seqs_per_query)
  digits <- as.integer(digits)
  
  stopifnot(min_len >= 1L, max_len >= min_len, max_miss <= 10L, minn_ms2 >= 2L, 
            min_mass >= 1L, max_mass >= min_mass, min_ms2mass >= 1L, 
            n_13c >= 0L, noenzyme_maxn >= 0L, 
            maxn_vmods_per_pep >= maxn_sites_per_vmod, max_n_prots > 1000L, 
            min_ms1_charge >= 1L, max_ms1_charge >= min_ms1_charge, 
            min_scan_num >= 1L, max_scan_num >= min_scan_num, 
            topn_mods_per_seq >= 1L, topn_seqs_per_query >= 1L)

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
  else {
    enzyme <- enzyme_lwr
  }

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
  
  # Output path
  out_path <- create_dir(out_path)
  dir.create(file.path(out_path, "Calls"), showWarnings = FALSE, recursive = FALSE)
  dir.create(file.path(out_path, "temp"), showWarnings = FALSE, recursive = FALSE)

  # grouped searches 
  # (this step before checking mgf_path)
  if (length(par_groups)) {
    if ("out_path" %in% names(par_groups))
      stop("Do not include `out_path` in `par_groups`.\n", 
           "The same parent `out_path` is assumed.", 
           call. = FALSE)
    
    if ("fasta" %in% names(par_groups))
      stop("Do not include `fasta` in `par_groups`.\n", 
           "The same set of `fasta` files is assumed.", 
           call. = FALSE)
    
    grp_args <- local({
      nms <- lapply(par_groups, names)
      all_nms <- sort(unique(unlist(nms, use.names = FALSE, recursive = FALSE)))
      nms_1 <- sort(nms[[1]])
      
      if (!identical(nms_1, all_nms))
        stop("Not all names are identical to those in the first group: ", 
             paste(nms_1, collapse = ", "), 
             call. = FALSE)
      
      fargs <- formalArgs(fun)
      bads <- nms_1[! nms_1 %in% fargs]
      
      if (length(bads)) 
        stop("Arguments in `par_groups` not defined in `", fun, "`:\n  ", 
             paste(bads, collapse = ", "), call. = FALSE)
      
      cargs <- names(this_call)
      cargs <- cargs[cargs != ""]
      dups <- nms_1[nms_1 %in% cargs]
      
      if (length(dups))
        stop("Arguments in `par_groups` already in the call", ":\n  ", 
             paste(dups, collapse = ", "), call. = FALSE)
      
      nms_1
    })
  }
  else {
    grp_args <- NULL
  }

  # MGF path
  if ("mgf_path" %in% grp_args) {
    mgf_path <- NULL
    
    mgf_paths <- lapply(par_groups, `[[`, "mgf_path")
    mgf_paths <- lapply(mgf_paths, checkMGF, grp_args, error = "warn")
    
    for (i in seq_along(mgf_paths)) 
      par_groups[[i]][["mgf_path"]] <- mgf_paths[[i]]
    
    rm(list = c("i"))
  }
  else {
    # (MGFs in sub-folders if group searches)
    mgf_path <- checkMGF(mgf_path, error = "warn")
    mgf_paths <- NULL
  }
  
  ## No-enzyme searches
  exec_noenzyme <- if (isTRUE(dots$bypass_noenzyme)) FALSE else TRUE
  
  if (enzyme == "noenzyme" && exec_noenzyme) {
    matchMS_noenzyme(this_call = this_call, min_len = min_len, max_len = max_len, 
                     fasta = fasta, out_path = out_path, mgf_path = mgf_path, 
                     noenzyme_maxn = noenzyme_maxn, quant = quant)
    
    return(NULL)
  }

  ## Mixed SILAC
  exec_silac_mix <- if (isTRUE(dots$bypass_silac_mix)) FALSE else TRUE
  
  if (length(silac_mix) && exec_silac_mix) {
    matchMS_silac_mix(silac = silac_mix, 
                      this_call = this_call, 
                      out_path = out_path, 
                      mgf_path = mgf_path)
    
    return(NULL)
  }
  
  # Searches by group (separate SILACs)
  exec_par_groups <- if (isTRUE(dots$bypass_par_groups)) FALSE else TRUE
  
  if (length(par_groups) && exec_par_groups) {
    df <- matchMS_par_groups(par_groups = par_groups, 
                             grp_args = grp_args,
                             mgf_paths = mgf_paths, 
                             this_call = this_call, 
                             out_path = out_path)
    
    return(df)
  }
  
  ## Theoretical MS1 masses
  bypass_pepmasses <- dots$bypass_pepmasses
  if (is.null(bypass_pepmasses)) bypass_pepmasses <- FALSE

  if (!bypass_pepmasses) {
    res <- calc_pepmasses2(
      fasta = fasta,
      acc_type = acc_type,
      acc_pattern = acc_pattern,
      fixedmods = fixedmods,
      varmods = varmods,
      include_insource_nl = include_insource_nl,
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
  }

  ## Bin theoretical peptides
  bypass_bin_ms1 <- dots$bypass_bin_ms1
  if (is.null(bypass_bin_ms1)) bypass_bin_ms1 <- FALSE
  
  if (!bypass_bin_ms1) {
    bin_ms1masses(res = res, 
                  min_mass = min_mass, 
                  max_mass = max_mass, 
                  ppm_ms1 = ppm_ms1, 
                  use_ms1_cache = use_ms1_cache, 
                  .path_cache = .path_cache, 
                  .path_ms1masses = .path_ms1masses, 
                  out_path = out_path)
    
    try(rm(list = "res"), silent = TRUE)
    gc()
  }

  ## MGFs
  bypass_mgf <- dots$bypass_mgf
  if (is.null(bypass_mgf)) bypass_mgf <- FALSE

  if (!bypass_mgf) {
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
              index_ms2 = FALSE, 
              enzyme = enzyme)
  }

  ## MSMS matches
  bypass_ms2match <- dots$bypass_ms2match
  if (is.null(bypass_ms2match)) bypass_ms2match <- FALSE
  
  use_first_rev <- dots$use_first_rev
  if (is.null(use_first_rev)) use_first_rev <- FALSE

  if (!bypass_ms2match) {
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
             use_first_rev = use_first_rev, 
             
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
  }

  ## Peptide scores
  bypass_from_pepscores <- dots$bypass_from_pepscores
  if (is.null(bypass_from_pepscores)) bypass_from_pepscores <- FALSE

  if (bypass_from_pepscores) 
    return(NULL)
  
  bypass_pepscores <- dots$bypass_pepscores
  if (is.null(bypass_pepscores)) bypass_pepscores <- FALSE
  
  if (!bypass_pepscores) {
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
                   add_ms2theos = add_ms2theos, 
                   add_ms2theos2 = add_ms2theos2, 
                   add_ms2moverzs = add_ms2moverzs, 
                   add_ms2ints = add_ms2ints,
                   digits = digits)
  }
  
  ## Peptide FDR 
  bypass_pepfdr <- dots$bypass_pepfdr
  if (is.null(bypass_pepfdr)) bypass_pepfdr <- FALSE
  
  if (!bypass_pepfdr) {
    prob_cos <- calc_pepfdr(target_fdr = target_fdr, 
                       fdr_type = fdr_type, 
                       min_len = min_len, 
                       max_len = max_len, 
                       max_pepscores_co = max_pepscores_co, 
                       match_pepfdr = match_pepfdr, 
                       out_path = out_path)
    
    post_pepfdr(prob_cos, out_path)
    rm(list = "prob_cos")
  }

  ## Peptide ranks and score deltas between `pep_ivmod`
  bypass_peploc <- dots$bypass_peploc
  if (is.null(bypass_peploc)) bypass_peploc <- FALSE
  
  if (!bypass_peploc) {
    calc_peploc(out_path = out_path, 
                topn_mods_per_seq = topn_mods_per_seq, 
                topn_seqs_per_query = topn_seqs_per_query)
    gc()
  }

  ## Protein accessions, score cut-offs and optional reporter ions
  bypass_from_protacc <- dots$bypass_from_protacc
  if (is.null(bypass_from_protacc)) bypass_from_protacc <- FALSE
  
  if (bypass_from_protacc) 
    return(NULL)

  if (enzyme != "noenzyme") {
    df <- add_prot_acc(out_path = out_path, 
                       .path_cache = .path_cache, 
                       .path_fasta = .path_fasta)
  }
  else {
    df <- add_prot_acc2(out_path = out_path, 
                        .path_cache = .path_cache, 
                        .path_fasta = .path_fasta)
  }

  df <- calc_protfdr(df = df, 
                     target_fdr = target_fdr, 
                     max_protscores_co = max_protscores_co, 
                     out_path = out_path)
  
  df <- add_rptrs(df, quant, out_path)
  gc()

  ## Clean-ups
  # (raw_file etc. already mapped if `from_group_search`)
  from_group_search <- dots$from_group_search
  if (!isTRUE(from_group_search)) df <- map_raw_n_scan(df, mgf_path)

  df$pep_ms1_delta <- df$ms1_mass - df$theo_ms1

  df <- dplyr::rename(df, 
    pep_scan_title = scan_title,
    pep_exp_mz = ms1_moverz,
    pep_exp_mr = ms1_mass,
    pep_exp_z = ms1_charge,
    pep_calc_mr = theo_ms1,
    pep_delta = pep_ms1_delta,
    pep_tot_int = ms1_int,
    pep_ret_range = ret_time,
    pep_scan_num = scan_num,
    pep_n_ms2 = ms2_n,
    pep_frame = frame)

  nms <- names(df)
  df <- dplyr::bind_cols(
    df[grepl("^prot_", nms)],
    df[grepl("^pep_", nms)],
    df[grepl("^psm_", nms)],
    df[!grepl("^prot_|^pep_|^psm_", nms)],
  )
  rm(list = "nms")
  
  df <- reloc_col_after(df, "pep_exp_z", "pep_exp_mr")
  df <- reloc_col_after(df, "pep_calc_mr", "pep_exp_z")
  df <- reloc_col_after(df, "pep_delta", "pep_calc_mr")
  readr::write_tsv(df, file.path(out_path, "psmC.txt"))
  
  ## psmC to psmQ
  df <- try_psmC2Q(df, 
                   out_path = out_path,
                   fdr_type = fdr_type, 
                   combine_tier_three = combine_tier_three, 
                   max_n_prots = max_n_prots)

  local({
    session_info <- sessionInfo()
    save(session_info, file = file.path(out_path, "Calls", "proteoM.rda"))
  })

  .savecall <- TRUE

  invisible(df)
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
    suppressWarnings(
      rm(list = c(".path_cache", ".path_ms1masses", ".time_stamp"),
         envir = .GlobalEnv)
    )

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

  out <- dplyr::filter(out, pep_issig, !pep_isdecoy, !grepl("^-", prot_acc))

  # Set aside one-hit wonders
  out3 <- out %>%
    dplyr::filter(!prot_issig, prot_n_pep == 1L) %>%
    dplyr::mutate(prot_tier = 3L)

  out <- dplyr::bind_rows(
    out %>% dplyr::filter(prot_issig),
    out %>% dplyr::filter(!prot_issig, prot_n_pep >= 2L)
  ) %>%
    dplyr::mutate(prot_tier = ifelse(prot_issig, 1L, 2L))

  # the same peptide can be present in all three protein tiers; 
  # steps up if pep_seq(s) in tier 3 also in tiers 1, 2
  if (FALSE) {
    rows <- out3$pep_seq %in% out$pep_seq
    
    out <- dplyr::bind_rows(out, out3[rows, ])
    out3 <- out3[!rows, ]
    
    rm(list = "rows")
    gc()
  }

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
  if (FALSE) {
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
        warning("All TMT modifications need to be `TMTpro18` or `TMT18plex` at `", 
                quant, "`.\n", 
                tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
                call. = FALSE)
    } 
    else if (quant == "tmt16") {
      ok <- all(grepl("TMTpro.* |TMT16plex.* ", possibles))
      
      if (!ok) 
        warning("All TMT modifications need to be `TMTpro` or `TMT16plex` at `", 
                quant, "`.\n", 
                tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3, 
                call. = FALSE)
    } 
    else {
      ok <- all(grepl("TMT6plex.* |TMT10plex.* |TMT11plex.* ", possibles))
      
      if (!ok) 
        warning("All TMT modifications need to be `TMT6plex`, `TMT10plex` or `TMT11plex` at `", 
                quant, "`.\n", 
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
                              noenzyme_maxn = 0L, quant = "none") 
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
    
    # ^1.5, 0.6: 90% RAM aggressiveness with uniprot fasta human + mouse
    max(1L, floor(fct_mem/(fct_fasta^1.5) * .5))
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
      
      if (i > 1L) {
        mgf_call <- file.path(out_paths[[1]], "Calls", "load_mgfs.rda")
        
        if (file.exists(mgf_call)) {
          sub_call_path <- create_dir(file.path(sub_path, "Calls"))
          
          file.copy(mgf_call, file.path(sub_call_path, "load_mgfs.rda"), 
                    overwrite = TRUE)
        }
      }

      sub_call$min_len <- start
      sub_call$max_len <- end
      sub_call$out_path <- sub_path
      sub_call$mgf_path <- mgf_path
      sub_call$bypass_noenzyme <- TRUE
      sub_call$bypass_from_pepscores <- TRUE
      sub_call$use_first_rev <- TRUE
      
      ans <- tryCatch(eval(sub_call), error = function (e) NULL)

      if (is.null(ans)) {
        message("No results from `min_len = ", start, "` to `max_len = ", end, "`.")
        # unlink(sub_path, recursive = TRUE)
      }
      
      file.copy(file.path(.path_ms1masses, .time_stamp, "prot_pep_annots.rds"), 
                file.path(sub_path))
      file.copy(file.path(.path_ms1masses, .time_stamp, "prot_pep_annots_rev.rds"), 
                file.path(sub_path))

      gc()
    }
    
    file.copy(file.path(out_paths[[1]], "Calls"), out_path, recursive = TRUE)
    combine_ion_matches(out_path, out_paths, type = "ion_matches_")
    combine_ion_matches(out_path, out_paths, type = "reporters_")
    combine_ion_matches(out_path, out_paths, type = "ion_matches_rev_")

    this_call$bypass_noenzyme <- TRUE
    this_call$bypass_pepmasses <- TRUE
    this_call$bypass_bin_ms1 <- TRUE
    this_call$bypass_mgf <- TRUE
    this_call$bypass_ms2match <- TRUE

    ans <- tryCatch(eval(this_call), error = function (e) NULL)

    message("Done (noenzyme search).")
    options(show.error.messages = FALSE)
    stop()
  }
  
  invisible(NULL)
}


#' SILAC
#'
#' Searches against mixed SILAC groups (heave and light mixed into one sample).
#' 
#' @param this_call An expression from match.call.
#' @inheritParams matchMS
matchMS_silac_mix <- function (silac_mix = list(base = NULL, 
                                            heavy = c("K8 (K)", "R10 (R)")), 
                               this_call, out_path, mgf_path) 
{
  message("Searches against SILAC groups ", 
          "(heavy, light etc. mixed into one sample)")
  
  lapply(c("silac_mix", "mgf_path", "this_call", "out_path"), function (x) {
    if (is.null(x)) stop("`", x, "` cannot be NULL.")
  })
  
  if (!is.list(silac_mix)) {
    stop("Supply silac_mix groups as list, e.g., \n\n", 
         "  # unlabelled base\n", 
         "  # (see ?add_unimod for custom silac_mix groups)\n", 
         "  silac_mix = list(base = NULL, heavy = c(\"K8 (K)\", \"R10 (R)\"))\n\n", 
         "  # labelled base \n", 
         "  silac_mix = list(base   = c(\"K4 (K)\", \"R4 (R)\"), \n", 
         "                   median = c(\"K6 (K)\", \"R6 (R)\"), \n", 
         "                   heavy  = c(\"K8 (K)\", \"R10 (R)\"))\n\n",  
         "# custom groups\n",
         "# (always require base)\n", 
         "  silac_mix = list(base    = NULL, \n", 
         "                   my_grp1 = ...)\n\n",  
         "  silac_mix = list(base    = c(\"K6 (K)\", \"R6 (R)\"), \n", 
         "                   my_grp1 = ...)\n\n",  
         call. = FALSE)
  }
  
  nms <- names(silac_mix)
  
  if (!"base" %in% nms) {
    warning("Required `base` group not found; assume: base = NULL")
    nms <- c("base", nms)
    silac_mix <- c(list(base = NULL), silac_mix)
  }

  if (all(unlist(lapply(silac_mix, is.null))))
    stop("`silac_mix` groups cannot be all NULL.", call. = FALSE)
  
  len <- length(nms)
  
  if (len < 2L)
    stop("Need at least two `silac_mix` groups.", call. = FALSE)
  
  out_paths <- vector("list", len)
  
  for (i in seq_len(len)) {
    sub_call <- this_call
    sub_nm <- nms[i]
    sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))
    
    ok <- file.exists(file.path(sub_path, "psmQ.txt"))
    
    if (ok) next
    
    if (i > 1L) {
      mgf_call <- file.path(out_paths[[1]], "Calls", "load_mgfs.rda")
      
      if (file.exists(mgf_call)) {
        sub_call_path <- create_dir(file.path(sub_path, "Calls"))
        
        file.copy(mgf_call, file.path(sub_call_path, "load_mgfs.rda"), 
                  overwrite = TRUE)
      }
    }
    
    sub_call$fixedmods <- c(eval(sub_call$fixedmods), silac_mix[[sub_nm]])
    sub_call$out_path <- sub_path
    sub_call$mgf_path <- mgf_path
    sub_call$bypass_silac_mix <- TRUE
    sub_call$silac_mix <- NULL

    df <- tryCatch(eval(sub_call), error = function (e) NULL)
    
    if (is.null(df)) {
      warning("No results at `silac_mix = ", sub_nm, "`.")
      # unlink(sub_path, recursive = TRUE)
    }

    rm(list = c("df"))
    gc()
  }
  
  message("Combine mixed SILAC results.")
  comine_PSMsubs(sub_paths = out_paths, groups = nms, out_path = out_path)

  gc()
  
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
  
  lapply(c("par_groups", "grp_args", "this_call", "out_path"), function (x) {
    if (is.null(x))
      stop("`", x, "` cannot be NULL.")
  })
  
  if (!is.list(par_groups)) {
    stop("Supply `par_groups` as list, e.g., \n\n", 
         "par_groups = list(\n", 
         "  list(mgf_path  = \"~/proteoM/my_proj/mgfs/grp_1\"", ",\n", 
         "       fixedmods = c(\"Carbamidomethyl (C)\")", "),\n", 
         "  list(mgf_path  = \"~/proteoM/my_proj/mgfs/grp_2\"", ",\n", 
         "       fixedmods = c(\"Carbamidomethyl (C)\", \"K8 (K)\", \"R10 (R)\")", ")\n", 
         "  )", 
         call. = FALSE)
  }
  
  nms <- names(par_groups)
  
  if (all(unlist(lapply(par_groups, is.null))))
    stop("`par_groups` cannot be all NULL.")
  
  len <- length(nms)
  
  if (len < 2L)
    stop("Need at least two `par_groups`.", call. = FALSE)
  
  ans <- out_paths <- vector("list", len)
  
  for (i in seq_len(len)) {
    sub_call <- this_call
    sub_nm <- nms[i]
    sub_pars <- par_groups[[i]]
    sub_path <- out_paths[[i]] <- create_dir(file.path(out_path, sub_nm))
    
    file_peploc <- file.path(sub_path, "temp", "peploc.rds")
    
    if (file.exists(file_peploc)) {
      ans[[i]] <- readRDS(file_peploc)
      next
    }

    for (arg in grp_args) 
      sub_call[[arg]] <- sub_pars[[arg]]

    sub_call$out_path <- sub_path
    sub_call$bypass_par_groups <- TRUE
    sub_call$bypass_from_protacc <- TRUE
    sub_call$par_groups <- NULL

    df <- tryCatch(eval(sub_call), error = function (e) NULL)
    
    if (!file.exists(file_peploc)) 
      stop("File not found: ", file_peploc)

    ans[[i]] <- if (is.null(df)) readRDS(file_peploc) else df
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
  saveRDS(out, file.path(out_path, "temp", "peploc.rds"))
  saveRDS(out_paths, file.path(out_path, "temp", "out_paths.rds"))
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


#' Combines the results of ion matches.
#' 
#' @param out_path A parent output path.
#' @param out_paths Sub output pathes.
#' @param type The type of data for combining.
combine_ion_matches <- function (out_path, out_paths, type = "ion_matches_") 
{
  out_path_temp <- create_dir(file.path(out_path, "temp"))
  out_paths_temp <- lapply(out_paths, function(x) file.path(x, "temp"))
  
  pat <- paste0(type, "[0-9]+\\.rds$")
  pat2 <- paste0(type, "([0-9]+)\\.rds$")
  
  files_mts <- local({
    xs <- list.files(out_paths_temp[[1]], pattern = pat)
    
    if (length(xs)) {
      idxes <- sort(as.integer(gsub(pat2, "\\1", xs)))
      files <- paste0(type, idxes, ".rds")
    }
    else {
      files <- NULL
    }
  })
  
  len_mts <- length(files_mts)
  
  if (!len_mts) {
    warning("Files not found: ", type)
    return(NULL)
  }
  
  ans_mts <- vector("list", len_mts)
  
  for (i in seq_along(ans_mts)) {
    ans_mts[[i]] <- lapply(out_paths_temp, function (path) {
      readRDS(file.path(path, files_mts[i]))
    }) %>% 
      dplyr::bind_rows()
    
    saveRDS(ans_mts[[i]], file.path(out_path_temp, files_mts[i]))
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
  combine_ion_matches(out_path, sub_paths, type = "ion_matches_rev_")
  
  invisible(NULL)
}


#' Checks the path of MGF files
#' 
#' @param error Character string; the level of error.
#' @inheritParams matchMS
#' @inheritParams matchMS_par_groups
checkMGF <- function (mgf_path = NULL, grp_args = NULL, 
                            error = c("stop", "warn")) 
{
  mgf_path <- find_dir(mgf_path)
  error <- match.arg(error)
  
  if (! error %in% c("warn", "stop"))
    stop("`error` needs to be one of \"error\" or \"stop\".")
  
  if (is.null(mgf_path)) 
    stop("`mgf_path` not found.", call. = FALSE)
  
  filelist <- list.files(path = mgf_path, pattern = "\\.mgf$")
  
  if (!length(filelist)) {
    if (error == "warn")
      warning("No `.mgf` files immediately under ", mgf_path, call. = FALSE)
    else
      stop("No `.mgf` files immediately under ", mgf_path, call. = FALSE)
  }

  rm(list = c("filelist"))
  
  invisible(mgf_path)
}


#' Maps raw_file and scan_title from indexes to real values.
#' 
#' @param df A data frame.
#' @inheritParams matchMS
map_raw_n_scan <- function (df, mgf_path) 
{
  file_raw <- file.path(mgf_path, "raw_indexes.rds")
  file_scan <- file.path(mgf_path, "scan_indexes.rds")
  
  if (file.exists(file_raw)) {
    raws <- readRDS(file_raw)
    raws2 <- names(raws)
    names(raws2) <- raws
    df$raw_file <- unname(raws2[df$raw_file])
  }
  else {
    stop("File not found: ", file_raw)
  }
  
  if (file.exists(file_scan)) {
    scans <- readRDS(file_scan)
    scans2 <- names(scans)
    names(scans2) <- scans
    df$scan_title <- unname(scans2[df$scan_title])
  }
  else {
    stop("File not found: ", file_scan)
  }
  
  invisible(df)
}


