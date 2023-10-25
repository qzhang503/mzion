#' An integrated facility for searches of mass spectrometry data.
#'
#' Database searches of MS/MS data (DDA).
#'
#' @section \code{Output columns}: \code{system.file("extdata",
#'   "column_keys.txt", package = "mzion")} \cr
#'
#' @section \code{Notes}: The annotation of protein attributes, including
#'   percent coverage, will be performed with \link[proteoQ]{normPSM} given that
#'   values will be affected with the combination of multiple PSM tables.
#'
#'   The search is a two-way match: (a) a forward matching of theoretical values
#'   to experiment ones and (b) a backward matching of the experimental values
#'   to the theoretical ones. This allows the establishment of one-to-one
#'   correspondences between experiments and theoreticals. The correspondences
#'   are made available to users in files of \code{ion_matches_1.rds...} (nested
#'   form) and \code{list_table_1.rds...} etc. (flat form). To open the files,
#'   use \code{qs::qread(...)}. A more self-contained output can be made
#'   available at the TRUE of \code{add_ms2theos}, \code{add_ms2theos2},
#'   \code{add_ms2moverzs} and \code{add_ms2ints}.
#'
#'   When there is no evidence to distinguish, e.g. distinct primary sequences
#'   of \code{P[EMPTY]EPTIDE} and \code{P[MTYPE]EPTIDE}, both will be reported
#'   by \code{mzion} and further kept by proteoQ. For peptides under the same
#'   primary sequence, the redundancy in the positions and/or neutral losses of
#'   \code{Anywhere} variable modifications are also kept in the outputs of
#'   \code{mzion} but removed with proteoQ.
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
#' @param rm_dup_term_anywhere Logical; if TRUE, removes combinations in
#'   variable modifications with site(s) in positions of both terminal and
#'   anywhere, e.g., "Gln->pyro-Glu (N-term = Q)" and "Deamidated (Q).
#' @param fixedlabs Character string(s) of fixed isotopic labels. See examples
#'   of SILAC for details. Can be but not typically used in standard alone
#'   searches of labeled residues.
#' @param varlabs Character string(s) of variable isotopic labels. See examples
#'   of SILAC for details. Can be but not typically used in standard alone
#'   searches of labeled residues.
#' @param locmods Among \code{varmods} for the consideration of localization
#'   probabilities; for instance, \code{locmods = NULL} for nothing,
#'   \code{locmods = c("Phospho (S)", "Phospho (T)", "Phospho (Y)")} for
#'   phosphopeptides, \code{locmods = "Acetyl (K)"} for lysine acetylation.
#'   \code{fixedmods} that were coerced to \code{varmods} will be added
#'   automatically to \code{locmods}.
#'
#'   For convenience, the default is set to look for applicable peptide
#'   phosphorylation (and may encounter warning messages if the data type is
#'   different to the default).
#' @param mod_motifs The motifs to restrict \code{Anywhere} variable
#'   modification. For example, provided the \code{Anywhere} variable
#'   modifications containing \code{c("Oxidation (M)", "Deamidated (N)")} and
#'   \code{mod_motifs = list(`Deamidated (N)` = c("NG", "NM"), `Oxidation (M)` =
#'   c("NM", "MP"))},
#'   variable modifications will only be considered at sites that satisfy the
#'   motifs.
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
#' @param nes_fdr_group A character string in one of \code{c("all",
#'   "all_cterm_tryptic", "all_cterm_nontryptic", "base", "base_cterm_tryptic",
#'   "base_cterm_nontryptic")}. All peptides will be used in the classifications
#'   of targets and decoys at \code{"all"}. Peptides with the chemistry of
#'   C-terminal K or R will be used at \code{"all_cterm_tryptic"} (peptides from
#'   protein C-terminals being excluded). Peptides without C-terminal K or R
#'   will be used at \code{"all_cterm_nontryptic"}. The same applied to
#'   \code{"base_cterm_tryptic"} and \code{"base_cterm_nontryptic"} with the
#'   difference of only peptides from the \code{base} group being used. See also
#'   parameter \code{fdr_group}.
#' @param noenzyme_maxn Non-negative integer; the maximum number of peptide
#'   lengths for sectional searches at \code{noenzyme} specificity. The argument
#'   may be used to guard against RAM exhaustion. At the zero default, The
#'   peptide lengths from \code{min_len} to \code{max_len} will be broken
#'   automatically into continuous sections. At value 1, searches will be
#'   performed against individual peptide lengths; at value 2, two adjacent
#'   lengths will be taken at a time, etc.
#' @param maxn_vmods_setscombi A non-negative integer; the maximum number of
#'   sets of combinatorial variable modifications. The default is 512.
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
#' @param maxn_fnl_per_seq A non-negative integer; the maximum number of
#'   permutative neutral losses per peptide sequence for fixed modifications. To
#'   bypass the combinatorial of neutral losses, set \code{maxn_fnl_per_seq = 1}
#'   or \code{0}.
#' @param maxn_vnl_per_seq A non-negative integer; the maximum number of
#'   permutative neutral losses per peptide sequence for variable modifications.
#'   To bypass the combinatorial of neutral losses, set \code{maxn_vnl_per_seq =
#'   1} or \code{0}.
#' @param maxn_vmods_sitescombi_per_pep A non-negative integer; the maximum
#'   number of combinatorial variable modifications per peptide sequence (per
#'   module). The ways include the permutations in neutral losses and
#'   modifications (e.g., \code{Acetyl (K) and TMT (K)}).
#' @param min_len A positive integer; the minimum length of peptide sequences
#'   for considerations. Shorter peptides will be excluded. The default is 7.
#' @param max_len A positive integer; the maximum length of peptide sequences
#'   for considerations. Longer peptides will be excluded. The default is 40.
#' @param max_miss A non-negative integer; the maximum number of mis-cleavages
#'   per peptide sequence for considerations. The default is 2.
#' @param min_mass A positive integer; the minimum precursor mass for
#'   interrogation. The default is an arbitrarily low value (the primary guard
#'   against low molecular-weight precursors is \code{min_len}).
#' @param max_mass A positive integer; the maximum precursor mass for
#'   interrogation.
#' @param min_ms2mass A positive integer; the minimum MS2 mass for
#'   interrogation. The default is 110.
#' @param max_ms2mass A positive integer; the maximum MS2 mass for
#'   interrogation.
#' @param n_13c Number(s) of 13C off-sets in precursor masses, for example, over
#'   the range of \code{-1:2}. The default is 0.
#' @param ms1_notches A numeric vector; notches (off-sets) in precursor masses,
#'   e.g., \code{c(-97.976896)} to account fo the loss of a phosphoric acid in
#'   precursor masses.
#' @param ms1_neulosses Character string(s) specifying the neutral losses of
#'   precursors. The nomenclature is the same as those in argument
#'   \code{varmods}, e.g., \code{c("Phospho (S)", "Phospho (T)", "Phospho
#'   (Y)")}.
#'
#'   The argument is a simplified (narrower) usage of argument
#'   \code{ms1_notches} with additional specificity in modifications and sites.
#' @param maxn_neulosses_fnl A positive integer used in conjunction with
#'   arguments \code{ms1_notches} and \code{ms1_neulosses}. The maximum number
#'   of fixed MS2 neutral losses to be considered when searching against peak
#'   lists with MS1 mass off-sets. The default is two (one MS2 neutral loss in
#'   addition to 0).
#' @param maxn_neulosses_vnl A positive integer used in conjunction with
#'   arguments \code{ms1_notches} and \code{ms1_neulosses}. The maximum number
#'   of variable MS2 neutral losses to be considered when searching against peak
#'   lists with MS1 mass off-sets. The default is two (one MS2 neutral loss in
#'   addition to 0).
#' @param par_groups A low-priority feature. Parameter(s) of \code{matchMS}
#'   multiplied by sets of values in groups. Multiple searches will be performed
#'   separately against the parameter groups. For instance with one set of
#'   samples in SILAC light and the other in SILAC heavy, the experimenters may
#'   specify two arguments for parameter \code{mgf_path} and two arguments for
#'   parameter \code{fixedmods} that link to the respective samples. In this
#'   way, there is no need to search against, e.g. heavy-isotope-labeled K8R10
#'   with the light samples and vice versa. Note that results will be combined
#'   at the end, with the group names indicated under column \code{pep_group}.
#'   The default is NULL without grouped searches. See the examples under SILAC
#'   and Group searches.
#' @param silac_mix A list of labels indicating SILAC groups in samples. The
#'   parameter is most relevant for SILAC experiments where peptides of heavy,
#'   light etc. were \emph{mixed} into one sample. The default is NULL
#'   indicating a none mixed-SILAC experiment. See also the examples under
#'   SILAC.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @param reproc_dda_ms1 A shortcut for \code{is_mdda = TRUE, maxn_mdda_precurs
#'   = 1L}. The net effect is to reprocess the default MS1 m-over-zs and charge
#'   states.
#' @param is_mdda Logical; if TRUE, consider the multiple (chimeric) precursors
#'   in DDA.
#' @param deisotope_ms2 Logical; if TRUE, de-isotope MS2 features.
#' @param max_ms2_charge Maximum charge states for consideration with MS2
#'   deisotoping.
#' @param use_defpeaks Logical; if TRUE, uses MS1 m-over-z's, intensities and
#'   charge states pre-calculated by other peak-picking algorithms.
#' @param maxn_dia_precurs Maximum number of precursors for consideration in a
#'   DIA scan.
#' @param maxn_mdda_precurs Maximum number of precursors for consideration in a
#'   multi-precursor DDA scan.
#' @param n_mdda_flanks The number of preceding and following MS1 scans for
#'   consideration when averaging isotope envelops of precursors.
#' @param ppm_ms1_deisotope Mass error tolerance in MS1 deisotoping.
#' @param ppm_ms2_deisotope Mass error tolerance in MS2 deisotoping.
#' @param grad_isotope Positive numeric; the gradient threshold between two
#'   adjacent peaks in an isotopic envelop. The smaller the value, the more
#'   stringent it is in calling an adjacent peak being a mono-isotopic
#'   precursor.
#' @param fct_iso2 A multiplication factor for the fuzzy discrimination of a
#'   secondary precursor in parallel to a primary precursor in an isotopic
#'   envelop. The smaller the value, the more stringent it is in calling an
#'   adjacent peak being a mono-isotopic precursor
#' @param topn_ms2ions A positive integer; the top-n species for uses in MS2 ion
#'   searches.
#' @param topn_ms2ion_cuts Advanced feature. Either \code{NA} or a named vector.
#'   For instance, at \code{topn_ms2ions = 100} and \code{topn_ms2ion_cuts =
#'   c(`1000` = 90, `1100` = 5, `4500` = 5)}, the maximum number of MS2 peaks
#'   that can be used is \eqn{90} at \eqn{m/z \le 1000}, \eqn{5} at \eqn{1000 <
#'   m/z < 1100} and \code{5} at \eqn{m/z > 1100}. The trailing \code{`4500` =
#'   5} can be skipped.
#'
#'   To exclude MS2 features such as at \eqn{m/z > 4500}: \code{topn_ms2ion_cuts
#'   = c(`4500` = 100)}.
#'
#'   It is also possible to make a zone of voids. For instance, features at
#'   \eqn{1200 < m/z < 1250} can be excluded at \code{topn_ms2ion_cuts =
#'   c(`1000` = 90, `1200` = 5, `1250` = 0)}.
#'
#'   The default is \code{NA} where \code{topn_ms2ions} are picked uniformly
#'   across the entire m/z range.
#' @param minn_ms2 A positive integer; the minimum number of matched MS2 ions
#'   for consideration as a hit. Counts of secondary ions, e.g. b0, b* etc., are
#'   not part of the threshold.
#' @param exclude_reporter_region Logical; if TRUE, excludes MS2 ions in the
#'   region of TMT reporter ions. The default is FALSE. The corresponding range
#'   of TMT reporter ions is informed by \code{tmt_reporter_lower} and
#'   \code{tmt_reporter_upper}. The argument affects only TMT data.
#' @param tmt_reporter_lower The lower bound of the region of TMT reporter ions.
#'   The default is \eqn{126.1}.
#' @param tmt_reporter_upper The upper bound of the region of TMT reporter ions.
#'   The default is \eqn{135.2}.
#' @param index_mgf_ms2 Depreciated. A low-priority feature. Logical; if TRUE,
#'   converts up-frontly MS2 m-over-z values from numeric to integers as opposed
#'   to \emph{on-the-fly} conversion during ion matches. The default is FALSE.
#'   The \code{index_mgf_ms2 = TRUE} might be useful for very large MS files by
#'   reducing RAM footprints.
#'
#'   At \code{index_mgf_ms2 = TRUE}, the resolution of mass deltas between
#'   theoretical and experimental MS2 m-over-z values is limited by the
#'   \code{bin_width}, which is the ceiling half of the \code{ppm_ms2}. For
#'   instance, the \code{bin_width} is 10 ppm at the default \code{ppm_ms2 =
#'   20}. Due to the low resolution in mass deltas at \code{index_mgf_ms2 = TRUE},
#'   the fields of \code{pep_ms2_deltas, pep_ms2_deltas2, pep_ms2_deltas_mean,
#'   pep_ms2_deltas_sd} are nullified in the outputs.
#'
#' @param min_ms1_charge A positive integer; the minimum MS1 charge state for
#'   considerations. The default is 2.
#' @param max_ms1_charge A positive integer; the maximum MS1 charge state for
#'   considerations. The default is 6.
#' @param min_scan_num Depreciated. A positive integer; the minimum scan number
#'   for considerations. The default is 1. The setting only applies to MGFs with
#'   numeric scan numbers. For example, it has no effects on Bruker's timsTOF
#'   data.
#' @param max_scan_num Depreciated. A positive integer; the maximum scan number
#'   for considerations. The default is the maximum machine integer. The setting
#'   only applies to MGFs with numeric scan numbers.
#' @param min_ret_time A non-negative numeric; the minimum retention time in
#'   seconds for considerations. The default is 0.
#' @param max_ret_time A non-negative numeric; the maximum retention time in
#'   seconds for considerations. The default is \code{Inf}.
#' @param ppm_ms1 A positive integer; the mass tolerance of MS1 species. The
#'   default is 20.
#' @param ppm_ms2 A positive integer; the mass tolerance of MS2 species. The
#'   default is 20.
#' @param calib_ms1mass Logical; if TRUE, calibrates precursor masses.
#' @param ppm_reporters A positive integer; the mass tolerance of MS2 reporter
#'   ions. The default is 10.
#' @param ppm_ms1calib A positive integer; the mass tolerance of MS1 species for
#'   precursor mass calibration. The argument has no effect at
#'   \code{calib_ms1mass = FALSE}.
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
#'   one of c("protein", "peptide", "psm"). The default is \code{protein}.
#'
#'   Note that \code{fdr_type = protein} is comparable to \code{fdr_type =
#'   peptide} with the additional filtration of data at \code{prot_tier == 1}.
#' @param fdr_group A character string; the modification group(s) for uses in
#'   peptide FDR controls. The value is in one of \code{c("all", "base")}. The
#'   \code{base} corresponds to the modification group with the largest number
#'   of matches.
#' @param max_pepscores_co A positive numeric; the upper limit in the cut-offs
#'   of peptide scores for discriminating significant and insignificant
#'   identities. Note that a probability p-value of, e.g., \eqn{1e-20} may not
#'   be interpreted as more probable than \eqn{1e-10} (beyond the precision of a
#'   probability test). Also note that the numeric limit when converting
#'   probability to scores (\eqn{-10 \times log10(p)}); e.g. \eqn{-10 \times
#'   log10(e-324)} yielded \code{Inf} instead of 3240. For reasons like these,
#'   a high value of \code{max_pepscores_co} is not recommended.
#' @param min_pepscores_co A non-negative numeric; the lower limit in the
#'   cut-offs of peptide scores for discriminating significant and insignificant
#'   identities.
#' @param max_protscores_co A positive numeric; the upper limit in the cut-offs
#'   of protein scores for discriminating significant and insignificant
#'   identities.  For higher quality and data-driven thresholds, choose the
#'   default \code{max_protscores_co = Inf}.
#' @param max_protnpep_co A positive integer; the maximum number of peptides
#'   under a protein (\code{prot_n_pep}) to warrant the protein significance.
#'   For instance, proteins with \code{prot_n_pep > max_protnpep_co} will have a
#'   protein significance score cutoff of zero and thus are significant. Choose
#'   \code{max_protnpep_co = Inf} to learn automatically the cut-off from data.
#'   Note that the the value of \code{prot_n_pep} includes the counts of shared
#'   peptides.
#' @param method_prot_es_co A low-priority setting. A character string; the
#'   method to calculate the cut-offs of protein enrichment scores. The value is
#'   in one of \code{"median", "mean", "max", "min"} with the default of
#'   \code{"median"}. For instance at the default, the median of
#'   \code{peptide_score - pep_score_cutoff} under a protein will be used to
#'   represent the threshold of a protein enrichment score. For more conserved
#'   thresholds, the statistics of \code{"max"} may be considered.
#' @param soft_secions Logical; if TRUE, collapses the intensities of secondary
#'   ions to primary ions even when the primaries are absent. The default is
#'   FALSE. For instance, the signal of \code{b5^*} will be ignored if its
#'   primary ion \code{b5} is not matched. The impacts of \code{soft_secions =
#'   TRUE} on search performance has not yet been assessed.
#' @param topn_seqs_per_query Positive integer; a threshold to discard peptide
#'   matches under the same MS query with scores beyond the top-n.
#'
#'   The same \code{MS query} refers to the identity in \code{MS scan number}
#'   and \code{MS raw file name}. Target and decoys matches are treated
#'   separately.
#' @param topn_mods_per_seq Positive integer; a threshold to discard variable
#'   modifications under the same peptide match with scores beyond the top-n.
#'
#'   The same \code{peptide match} refers to matches with identities in \code{MS
#'   scan number}, \code{MS raw file name} and \code{peptide sequence}. Target
#'   and decoys matches are treated separately.
#'
#'   For a variable modification with multiple neutral losses (NL), the
#'   best-scored NL will be used in the ranking.
#' @param combine_tier_three Logical; if TRUE, combines search results at tier-3
#'   to tier-1 to form the single output of \code{psmQ.txt}. The default is
#'   FALSE in that data will be segregated into the three quality tiers (shown
#'   below) by the choice of \code{fdr_type}. Note that the argument affects
#'   only at the \code{fdr_type} of \code{psm} or \code{peptide} where there are
#'   no tier-2 outputs. In general, the tier-3 results correspond to
#'   one-hit-wonders and setting \code{combine_tier_three = TRUE} is
#'   discouraged.
#'
#'   In subproteome analysis, such as phosphoproteome analysis, some proteins
#'   may be well established globally, but fail the significance assessment by
#'   protein FDR on the local scale. In situations like this, it may be suitable
#'   to apply \code{fdr_type = "peptide"} or \code{fdr_type = "psm"} other than
#'   incurring \code{combine_tier_three = TRUE}.
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
#' @param max_n_prots Softly depreciated. A positive integer to threshold the
#'   maximum number of protein entries before coercing \code{fdr_type} from
#'   \code{psm} or \code{peptide} to \code{protein}. The argument has no effect
#'   if \code{fdr_type} is already \code{protein}. In general, there is no need
#'   to change the default.
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
#'
#' @param svm_reproc Logical; if TRUE, reprocesses peptide data for significance
#'   thresholds with a support vector machine (SVM) approach analogous to
#'   \href{https://www.nature.com/articles/nmeth1113}{Percolator}.
#' @param svm_kernel The SVM kernel. See also \link[e1071]{svm}.
#' @param svm_feats Features used for SVM classifications.
#' @param svm_iters A positive integer; the number of iterations in
#'   \link[e1071]{svm}.
#' @param svm_cv Logical; if TRUE, performs cross validation for the
#'   regularization cost.
#' @param svm_k A positive integer; specifies the k-number of folds in cross
#'   validation.
#' @param svm_costs The cost constraints for k-fold cross validation.
#' @param svm_def_cost The default cost for SVM.
#' @param svm_iters The number of iteration in SVM learning.
#'
#' @param .path_cache The file path of cached search parameters. The parameter
#'   is for the users' awareness of the underlying structure of file folders and
#'   the use of default is suggested. Occasionally experimenters may remove the
#'   file folder for disk space or under infrequent events of modified framework
#'   incurred by the developer.
#'
#' @param .path_fasta The parent file path to the theoretical masses of MS1
#'   precursors. At the NULL default, the path is \code{gsub("(.*)\\.[^\\.]*$",
#'   "\\1", get("fasta", envir = environment())[1])}. The parameter is for the
#'   users' awareness of the structure of file folders and the use of default is
#'   suggested. Occasionally experimenters may remove the file folder for disk
#'   space or under infrequent events of modified framework incurred by the
#'   developer.
#' @param make_speclib Makes spectrum library from the search results of
#'   \code{psmQ.txt}.
#' @param by_modules Not used. Logical. At the TRUE default, searches MS data by
#'   individual modules of combinatorial fixed and variable modifications. If
#'   FALSE, search all modules together. The later would probably need more than
#'   32G RAM if the number of modules is over 96.
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
#' @section \code{mzTab}: \link{make_mztab} converts outputs from the mzion ->
#'   proteoQ pipeline to mzTab files.
#' @return A list of complete PSMs in \code{psmC.txt}; a list of quality PSMs in
#'   \code{psmQ.txt}.
#' @examples
#' \dontrun{
#' ## All examples are hypothetical
#' ## (Users are responsible for supplying FASTA and peak lists in MGF or mzML)
#'
#' # TMT-10plex
#' matchMS(
#'   fasta    = c("~/mzion/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                "~/mzion/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                "~/mzion/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   max_miss = 2,
#'   quant    = "tmt10",
#'   fdr_type = "protein",
#'   out_path = "~/mzion/examples",
#' )
#'
#' # TMT-16plex, phospho
#' matchMS(
#'   fixedmods = c("TMTpro (N-term)", "TMTpro (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)", "Phospho (S)", "Phospho (T)",
#'                 "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"),
#'   locmods   = c("Phospho (S)", "Phospho (T)", "Phospho (Y)"),
#'   quant     = "tmt16",
#'   fdr_type  = "psm",
#'   out_path  = "~/mzion/examples",
#' )
#'
#' # TMT-18plex
#' matchMS(
#'   fixedmods = c("TMTpro (N-term)", "TMTpro (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'                 "Deamidated (N)", "Deamidated (Q)",
#'                 "Gln->pyro-Glu (N-term = Q)"),
#'   quant     = "tmt18",
#'   out_path  = "~/mzion/examples",
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
#'   out_path  = "~/mzion/examples",
#' )
#'
#' # Wrapper of matchMS(enzyme = noenzyme, ...) without sectional searches
#' #   by ranges of peptide lengths
#' matchMS_NES(
#'   fasta    = c("~/mzion/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
#'                "~/mzion/dbs/fasta/refseq/refseq_mm_2013_07.fasta",
#'                "~/mzion/dbs/fasta/crap/crap.fasta"),
#'   acc_type = c("refseq_acc", "refseq_acc", "other"),
#'   quant    = "tmt10",
#'   fdr_type = "protein",
#'   out_path = "~/mzion/examples",
#' )
#'
#'
#' # Custom Unimod (Oxi+Carbamidomethyl)
#' # (see also calc_unimod_compmass)
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
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'                 "Oxi+Carbamidomethyl (M)"),
#'   quant    = "tmt10",
#' )
#'
#'
#' ## Stable isotope-labeled K and R (not SILAC mixtures)
#' # Anywhere K8, R10
#' matchMS(
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'                 "Label:13C(6)15N(2) (K)", "Label:13C(6)15N(4) (R)"),
#'   quant     = "none",
#' )
#'
#' # K8, R10 + TMT10
#' matchMS(
#'   fixedmods = c("TMT6plex (N-term)", "TMT10plex+K8 (K)", "Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'                 "R10 (R)"),
#'   quant     = "tmt10",
#' )
#'
#'
#' #######################################
#' # SILAC
#' #######################################
#'
#' ## 1. heavy and light mixed into one sample (classical SILAC)
#' # (i) K8R10
#' matchMS(
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'                 "Gln->pyro-Glu (N-term = Q)"),
#'   silac_mix = list(base  = c(fixedlabs = NULL, varlabs = NULL),
#'                    K8R10 = c(fixedlabs = c("Label:13C(6)15N(2) (K)",
#'                                            "Label:13C(6)15N(4) (R)"),
#'                              varlabs   = NULL)),
#'   ...
#'   )
#'
#'
#' # (ii) base: unlabeled; grpC: 13C; grpN: 15N
#' # (example Dong-Ecoli-QE, Nat. Biotech. 2018, 36, 1059-1061)
#' matchMS(
#'   fixedmods = c("Carbamidomethyl (C)"),
#'   varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'                 "Gln->pyro-Glu (N-term = Q)"),
#'
#'   silac_mix = list(base = c(fixedlabs = NULL, varlabs   = NULL),
#'
#'                    grpC = c(fixedlabs = c("Label:13C(3) (A)", "Label:13C(6) (R)",
#'                                           "Label:13C(4) (N)", "Label:13C(4) (D)",
#'                                           "Label:13C(3) (C)", "Label:13C(5) (E)",
#'                                           "Label:13C(5) (Q)", "Label:13C(2) (G)",
#'                                           "Label:13C(6) (H)", "Label:13C(6) (I)",
#'                                           "Label:13C(6) (L)", "Label:13C(6) (K)",
#'                                           "Label:13C(5) (M)", "Label:13C(9) (F)",
#'                                           "Label:13C(5) (P)", "Label:13C(3) (S)",
#'                                           "Label:13C(4) (T)", "Label:13C(11) (W)",
#'                                           "Label:13C(9) (Y)", "Label:13C(5) (V)"),
#'                             varlabs = c("Label:13C(2) (Protein N-term)")),
#'
#'                    grpN = c(fixedlabs = c("Label:15N(1) (A)", "Label:15N(4) (R)",
#'                                           "Label:15N(2) (N)", "Label:15N(1) (D)",
#'                                           "Label:15N(1) (C)", "Label:15N(1) (E)",
#'                                           "Label:15N(2) (Q)", "Label:15N(1) (G)",
#'                                           "Label:15N(3) (H)", "Label:15N(1) (I)",
#'                                           "Label:15N(1) (L)", "Label:15N(2) (K)",
#'                                           "Label:15N(1) (M)", "Label:15N(1) (F)",
#'                                           "Label:15N(1) (P)", "Label:15N(1) (S)",
#'                                           "Label:15N(1) (T)", "Label:15N(2) (W)",
#'                                           "Label:15N(1) (Y)", "Label:15N(1) (V)"),
#'                             varlabs = c("Label:15N(-1) (N)",
#'                                         "Label:15N(-1) (N-term = Q)"))),
#'   ...)
#'
#'
#' # (iii) labeled base
#' # (first to add exemplary K4 Unimod if not yet available)
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
#'   silac_mix = list(base   = c(fixedlabs = c("K4 (K)"), varlabs = NULL),
#'                    median = c(fixedlabs = c("K6 (K)", "R6 (R)"), varlabs = NULL),
#'                    heavy  = c(fixedlabs = c("K8 (K)", "R10 (R)"), varlabs = NULL)),
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
#'     light = list(mgf_path  = "~/mzion/my_project/mgf/light",
#'                  fixedmods = "Carbamidomethyl (C)"),
#'     heavy = list(mgf_path  = "~/mzion/my_project/mgf/heavy",
#'                  fixedmods = c("Carbamidomethyl (C)", "K8 (K)", "R10 (R)"))
#'   ),
#'   quant = "none",
#'   ...
#' )
#'
#' }
#' @export
matchMS <- function (out_path = "~/mzion/outs",
                     mgf_path = file.path(out_path, "mgf"),
                     fasta = c("~/mzion/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
                               "~/mzion/dbs/fasta/crap/crap.fasta"),
                     acc_type = c("uniprot_acc", "other"),
                     acc_pattern = NULL,
                     fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                   "Carbamidomethyl (C)"),
                     varmods = c("Acetyl (Protein N-term)",
                                 "Oxidation (M)", "Deamidated (N)",
                                 "Gln->pyro-Glu (N-term = Q)"),
                     rm_dup_term_anywhere = TRUE, 
                     
                     ms1_neulosses = NULL, 
                     maxn_neulosses_fnl = 2L, 
                     maxn_neulosses_vnl = 2L, 

                     fixedlabs = NULL, 
                     varlabs = NULL, 
                     locmods = c("Phospho (S)", "Phospho (T)", "Phospho (Y)"), 
                     mod_motifs = NULL, 
                     enzyme = c("Trypsin_P", "Trypsin", "LysC", "LysN", "ArgC", 
                                "LysC_P", "Chymotrypsin", "GluC", "GluN", 
                                "AspC", "AspN", "SemiTrypsin_P", "SemiTrypsin", 
                                "SemiLysC", "SemiLysN", "SemiArgC", 
                                "SemiLysC_P", "SemiChymotrypsin", "SemiGluC", 
                                "SemiGluN", "SemiAspC", "SemiAspN", "Noenzyme", 
                                "Nodigest"),
                     custom_enzyme = c(Cterm = NULL, Nterm = NULL), 
                     nes_fdr_group = c("base", "base_cterm_tryptic", 
                                       "base_cterm_nontryptic", 
                                       "all", "all_cterm_tryptic", 
                                       "all_cterm_nontryptic", 
                                       "top3", "top3_cterm_tryptic", 
                                       "top3_cterm_nontryptic"), 
                     noenzyme_maxn = 0L, 
                     maxn_vmods_setscombi = 512L,
                     maxn_vmods_per_pep = 5L,
                     maxn_sites_per_vmod = 3L,
                     maxn_fnl_per_seq = 3L, 
                     maxn_vnl_per_seq = 3L, 
                     maxn_vmods_sitescombi_per_pep = 64L,
                     min_len = 7L, max_len = 40L, max_miss = 2L, 
                     min_mass = 200L, max_mass = 4500L, 
                     ppm_ms1 = 20L, 
                     n_13c = 0L, 
                     ms1_notches = 0, 
                     
                     par_groups = NULL, 
                     silac_mix = NULL, 
                     
                     type_ms2ions = "by", 
                     min_ms2mass = 115L, 
                     max_ms2mass = 4500L, 
                     minn_ms2 = 6L, 
                     ppm_ms2 = 20L, 
                     tmt_reporter_lower = 126.1, 
                     tmt_reporter_upper = 135.2, 
                     exclude_reporter_region = FALSE, 
                     index_mgf_ms2 = FALSE, 
                     
                     ppm_reporters = 10L,
                     quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16", "tmt18"),
                     
                     target_fdr = 0.01,
                     fdr_type = c("protein", "peptide", "psm"),
                     fdr_group = c("base", "all", "top3"), 
                     max_pepscores_co = 50, min_pepscores_co = 0, 
                     max_protscores_co = Inf, 
                     max_protnpep_co = 10L, 
                     method_prot_es_co = c("median", "mean", "max", "min"), 
                     soft_secions = FALSE, 
                     
                     topn_mods_per_seq = 1L, 
                     topn_seqs_per_query = 1L, 
                     
                     combine_tier_three = FALSE,
                     max_n_prots = 60000L, 
                     use_ms1_cache = TRUE, 
                     .path_cache = "~/mzion/.MSearches (1.3.0.1)/Cache/Calls", 
                     .path_fasta = NULL,
                     
                     reproc_dda_ms1 = TRUE, is_mdda = FALSE, 
                     deisotope_ms2 = TRUE, max_ms2_charge = 3L, 
                     use_defpeaks = FALSE, maxn_dia_precurs = 300L, 
                     maxn_mdda_precurs = 5L, n_mdda_flanks = 6L, 
                     ppm_ms1_deisotope = 8L, ppm_ms2_deisotope = 8L, 
                     grad_isotope = 1.6, fct_iso2 = 3.0, 
                     
                     topn_ms2ions = 150L,
                     topn_ms2ion_cuts = NA, 
                     min_ms1_charge = 2L, max_ms1_charge = 4L, 
                     min_scan_num = 1L, max_scan_num = .Machine$integer.max, 
                     min_ret_time = 0, max_ret_time = Inf, 
                     calib_ms1mass = FALSE, 
                     ppm_ms1calib = 20L,

                     add_ms2theos = FALSE, add_ms2theos2 = FALSE, 
                     add_ms2moverzs = FALSE, add_ms2ints = FALSE,
                     
                     svm_reproc = FALSE,
                     svm_kernel = "radial",
                     svm_feats  = c("pep_score", "pep_ret_range", 
                                    "pep_delta", "pep_n_ms2", 
                                    "pep_expect", "pep_exp_mz", # "pep_exp_z", 
                                    "pep_exp_mr", "pep_tot_int", 
                                    "pep_n_matches2", "pep_ms2_deltas_mean"), 
                     svm_cv = TRUE, svm_k  = 3L, 
                     svm_costs = c(.1, .3, 1, 3, 10), svm_def_cost = 1, 
                     svm_iters  = 10L, 
                     
                     make_speclib = FALSE, 
                     
                     by_modules = TRUE, 
                     digits = 4L, ...) 
{
  options(digits = 9L)

  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) {
        # Don't: "fun = fun"; seem name collide of `fun` when called from Shiny
        tryCatch(save_call2(path = file.path(out_path, "Calls"), fun = "matchMS"), 
                 error = function(e) NA)
      }
    },
    add = TRUE
  )
  
  message("Started at: ", Sys.time())
  max_integer <- .Machine$integer.max
  
  # Shiny compatibles
  if (is.na(max_protscores_co)) max_protscores_co <- Inf
  if (is.na(max_scan_num)) max_scan_num <- max_integer
  if (is.na(max_ret_time)) max_ret_time <- max_integer
  if ((!is.na(topn_ms2ion_cuts)) && (topn_ms2ion_cuts == "")) topn_ms2ion_cuts <- NA
  if (!(is.numeric(ms1_notches) && length(ms1_notches))) ms1_notches <- 0
  if (is.null(noenzyme_maxn)) noenzyme_maxn <- 0L
  if ((!is.null(custom_enzyme)) && custom_enzyme == "")
    custom_enzyme = c(Cterm = NULL, Nterm = NULL)

  oks <- fasta != ""
  
  if (!all(is.null(acc_pattern)) && length(acc_pattern) == length(fasta))
    acc_pattern <- acc_pattern[oks]
  else
    acc_pattern <- NULL
  
  if (length(acc_type) == length(fasta))
    acc_type <- acc_type[oks]
  
  fasta <- fasta[oks]
  
  # fixedmods <- gsub("Protein N-term = N-term", "Protein N-term", fixedmods)
  # fixedmods <- gsub("Protein C-term = C-term", "Protein C-term", fixedmods)
  # varmods <- gsub("Protein N-term = N-term", "Protein N-term", varmods)
  # varmods <- gsub("Protein C-term = C-term", "Protein C-term", varmods)

  # Calls
  this_call <- match.call()
  fun <- as.character(this_call[1])
  this_fml <- formals()
  
  ## Match arguments
  method_prot_es_co <- match.arg(method_prot_es_co)

  ## Developer's dots
  dots <- as.list(substitute(...()))

  if (is.null(aa_masses <- dots$aa_masses)) {
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
  
  if (!is.null(fixedlabs))
    aa_masses <- add_fixedlab_masses(fixedlabs, aa_masses)

  suppressWarnings(
    rm(list = c(".path_cache", ".path_fasta", ".path_ms1masses", 
                ".time_stamp", ".time_bin", ".path_bin"), 
       envir = .GlobalEnv))

  ## Preparation
  # modifications
  fixedmods <- sort(fixedmods)
  varmods <- sort(varmods)
  locmods <- check_locmods(locmods, fixedmods, varmods, ms1_neulosses)
  
  # accession pattern
  db_ord <- order(fasta)
  fasta  <- fasta[db_ord]
  acc_type <- acc_type[db_ord]
  
  if ((!is.null(acc_pattern)) && all(acc_pattern == "")) 
    acc_pattern <- NULL
  
  if (!is.null(acc_pattern)) {
    if (length(acc_pattern) != length(acc_type))
      stop("The length of `acc_pattern` needs to be the same as `acc_type`.")
    else
      acc_pattern <- acc_pattern[db_ord]
  }
  rm(list = "db_ord")

  # logical types
   
  
  stopifnot(vapply(c(soft_secions, combine_tier_three, calib_ms1mass, 
                     use_ms1_cache, add_ms2theos, add_ms2theos2, add_ms2moverzs, 
                     add_ms2ints, exclude_reporter_region, index_mgf_ms2, 
                     svm_reproc, svm_cv, rm_dup_term_anywhere, reproc_dda_ms1 , 
                     make_speclib, is_mdda, deisotope_ms2, use_defpeaks), 
                   is.logical, logical(1L)))

  # numeric types 
  stopifnot(vapply(c(maxn_vmods_setscombi, maxn_vmods_per_pep, 
                     maxn_sites_per_vmod, maxn_fnl_per_seq, maxn_vnl_per_seq, 
                     ms1_notches, maxn_neulosses_fnl, maxn_neulosses_vnl, 
                     maxn_vmods_sitescombi_per_pep, 
                     min_len, max_len, max_miss, topn_ms2ions, minn_ms2, 
                     min_mass, max_mass, min_ms2mass, max_ms2mass, n_13c, 
                     ppm_ms1, ppm_ms2, ppm_reporters, max_n_prots, digits, 
                     target_fdr, max_pepscores_co, min_pepscores_co, 
                     max_protscores_co, max_protnpep_co, topn_mods_per_seq, 
                     topn_seqs_per_query, tmt_reporter_lower, tmt_reporter_upper, 
                     max_ms2_charge, maxn_dia_precurs, maxn_mdda_precurs, 
                     n_mdda_flanks, ppm_ms1_deisotope, ppm_ms2_deisotope, 
                     grad_isotope, fct_iso2), 
                   is.numeric, logical(1L)))

  # (a) integers casting for parameter matching when calling cached)
  if (is.infinite(max_len)) max_len <- max_integer
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

  maxn_vmods_setscombi <- as.integer(maxn_vmods_setscombi)
  maxn_vmods_per_pep <- as.integer(maxn_vmods_per_pep)
  maxn_sites_per_vmod <- as.integer(maxn_sites_per_vmod)
  maxn_vmods_sitescombi_per_pep <- as.integer(maxn_vmods_sitescombi_per_pep)
  maxn_fnl_per_seq <- as.integer(maxn_fnl_per_seq)
  maxn_vnl_per_seq <- as.integer(maxn_vnl_per_seq)
  maxn_neulosses_fnl <- as.integer(maxn_neulosses_fnl)
  maxn_neulosses_vnl <- as.integer(maxn_neulosses_vnl)
  min_len <- as.integer(min_len)
  max_len <- as.integer(max_len)
  max_miss <- as.integer(max_miss)
  topn_ms2ions <- as.integer(topn_ms2ions)
  minn_ms2 <- as.integer(minn_ms2)
  min_mass <- as.integer(min_mass)
  max_mass <- as.integer(max_mass)
  min_ms2mass <- as.integer(min_ms2mass)
  max_ms2mass <- as.integer(max_ms2mass)
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
  
  max_ms2_charge <- as.integer(max_ms2_charge)
  maxn_dia_precurs <- as.integer(maxn_dia_precurs)
  maxn_mdda_precurs <- as.integer(maxn_mdda_precurs)
  n_mdda_flanks <- as.integer(n_mdda_flanks)
  ppm_ms1_deisotope <- as.integer(ppm_ms1_deisotope)
  ppm_ms2_deisotope <- as.integer(ppm_ms2_deisotope)
  digits <- as.integer(digits)
  
  stopifnot(min_len >= 1L, max_len >= min_len, max_miss <= 10L, minn_ms2 >= 2L, 
            min_mass >= 1L, max_mass >= min_mass, 
            min_ms2mass >= 1L, max_ms2mass > min_ms2mass, 
            maxn_vmods_sitescombi_per_pep >= 2L, noenzyme_maxn >= 0L, 
            maxn_fnl_per_seq >= 0L, maxn_vnl_per_seq >= 0L, 
            maxn_neulosses_fnl >= 0L, maxn_neulosses_vnl >= 0L,
            maxn_vmods_per_pep >= maxn_sites_per_vmod, max_n_prots > 1000L, 
            min_ms1_charge >= 1L, max_ms1_charge >= min_ms1_charge, 
            min_scan_num >= 1L, max_scan_num >= min_scan_num, 
            topn_mods_per_seq >= 1L, topn_seqs_per_query >= 1L, 
            tmt_reporter_lower < tmt_reporter_upper, max_ms2_charge >= 1L, 
            maxn_dia_precurs >= 1L, maxn_mdda_precurs >= 1L, n_mdda_flanks >= 1L, 
            ppm_ms1_deisotope >= 1L, ppm_ms2_deisotope >= 1L)

  # (b) doubles
  target_fdr <- round(as.double(target_fdr), digits = 2L)
  
  if (target_fdr > .25) 
    stop("Choose a smaller `target_fdr`.")
  
  min_ret_time <- round(min_ret_time, digits = 2L)
  max_ret_time <- round(max_ret_time, digits = 2L)
  max_pepscores_co <- round(max_pepscores_co, digits = 2L)
  min_pepscores_co <- round(min_pepscores_co, digits = 2L)
  max_protscores_co <- round(max_protscores_co, digits = 2L)

  stopifnot(max_pepscores_co >= 0, min_pepscores_co >= 0, max_protscores_co >= 0, 
            min_ret_time >= 0, max_pepscores_co >= min_pepscores_co, 
            max_ret_time >= min_ret_time, max_protnpep_co >= 1L, 
            grad_isotope >= 1.0, grad_isotope <= 5.0, 
            fct_iso2 >= 1.0, fct_iso2 <= 6.0)
  
  # named vectors
  if (any(is.na(topn_ms2ion_cuts)))
    mgf_cutmzs <- mgf_cutpercs <- numeric()
  else {
    if (is.infinite(topn_ms2ions))
      stop("Choose a finite \"topn_ms2ions\" value to enable \"topn_ms2ion_cuts\".")
    
    mgf_cutmzs <- as.numeric(names(topn_ms2ion_cuts))
    len <- length(topn_ms2ion_cuts)
    
    if (!identical(mgf_cutmzs, sort(mgf_cutmzs)))
      stop("\"mgf_cutmzs\" is not in an ascending order.")
    
    if (anyDuplicated(mgf_cutmzs))
      warning("Duplicated m-over-z cutpoints in \"topn_ms2ion_cuts\".")
    
    s_topn <- sum(topn_ms2ion_cuts)
    
    if (s_topn > topn_ms2ions) {
      stop("\"sum(topn_ms2ion_cuts) = ", s_topn, "\" is greater than ", 
           "\"topn_ms2ions = ", topn_ms2ions, ".\"")
    }
    else if (s_topn < topn_ms2ions) {
      mgf_cutpercs <- c(unname(topn_ms2ion_cuts), topn_ms2ions - s_topn)
      mgf_cutmzs <- c(mgf_cutmzs, max_ms2mass)
    }
    else {
      mgf_cutpercs <- unname(topn_ms2ion_cuts)
    }
    
    rm(list = c("len", "s_topn"))
    
    if (mgf_cutpercs[length(mgf_cutpercs)] != 0) {
      mgf_cutpercs <- c(mgf_cutpercs, 0)
      mgf_cutmzs <- c(mgf_cutmzs, max_ms2mass)
    }
  }

  # enzyme
  if ((!is.null(custom_enzyme)) && custom_enzyme == "")
    custom_enzyme <- NULL

  if (is.null(custom_enzyme)) {
    enzyme <- tolower(match.arg(enzyme))
  }
  else {
    warning("Overrule `enzyme` with `custom_enzyme`.")
    enzyme <- NULL
  }

  if ((!is.null(enzyme)) && (enzyme == "noenzyme"))
    max_miss <- 0L
  
  if ((!is.null(enzyme)) && (enzyme == "nodigest")) {
    max_miss <- 0L
    max_len <- max_integer
  }

  # fdr_type
  fdr_type <- match.arg(fdr_type)
  fdr_group <- match.arg(fdr_group)
  nes_fdr_group <- match.arg(nes_fdr_group)
  
  # quant
  quant <- match.arg(quant)

  # TMT
  check_tmt_pars(fixedmods, varmods, quant)
  
  # MS1 off-sets
  if (!all(ms1_neulosses %in% varmods))
    stop("Not all `ms1_neulosses` found in `varmods`.")
  
  check_notches(ms1_notches = ms1_notches, ms1_neulosses = ms1_neulosses)
  ms1_offsets <- find_ms1_offsets(n_13c = n_13c, ms1_notches = ms1_notches)
  is_notched  <- length(unique(c(ms1_offsets, ms1_neulosses))) > 1L

  # system paths
  homedir <- find_dir("~")
  
  if (is.null(.path_cache))
    .path_cache <- "~/mzion/.MSearches/Cache/Calls/"

  if (is.null(.path_fasta))
    .path_fasta <- file.path(gsub("(.*)\\.[^\\.]*$", "\\1", fasta[1]))

  .path_cache <- create_dir(.path_cache)
  .path_fasta <- create_dir(.path_fasta)
  .path_ms1masses <- create_dir(file.path(.path_fasta, "ms1masses"))

  fasta <- lapply(fasta, function (x) if (grepl("~", x)) gsub("~", homedir, x) else x)
  fasta <- unlist(fasta, recursive = FALSE, use.names = FALSE)
  
  # Output path
  out_path <- create_dir(out_path)
  dir.create(file.path(out_path, "Calls"), showWarnings = FALSE, recursive = FALSE)
  dir.create(file.path(out_path, "temp"), showWarnings = FALSE, recursive = FALSE)

  # grouped searches 
  # (this step before checking mgf_path)
  if (length(par_groups)) {
    if ("out_path" %in% names(par_groups))
      stop("Do not include `out_path` in `par_groups`.\n", 
           "The same parent `out_path` is assumed.")
    
    if ("fasta" %in% names(par_groups))
      stop("Do not include `fasta` in `par_groups`.\n", 
           "The same set of `fasta` files is assumed.")
    
    grp_args <- local({
      nms <- lapply(par_groups, names)
      all_nms <- sort(unique(unlist(nms, use.names = FALSE, recursive = FALSE)))
      nms_1 <- sort(nms[[1]])
      
      if (!identical(nms_1, all_nms))
        stop("Not all names are identical to those in the first group: ", 
             paste(nms_1, collapse = ", "))
      
      fargs <- formalArgs(fun)
      bads <- nms_1[! nms_1 %in% fargs]
      
      if (length(bads)) 
        stop("Arguments in `par_groups` not defined in `", fun, "`:\n  ", 
             paste(bads, collapse = ", "))
      
      cargs <- names(this_call)
      cargs <- cargs[cargs != ""]
      dups <- nms_1[nms_1 %in% cargs]
      
      if (length(dups))
        stop("Arguments in `par_groups` already in the call", ":\n  ", 
             paste(dups, collapse = ", "))
      
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
  
  if (isTRUE(enzyme == "noenzyme") && exec_noenzyme) {
    matchMS_noenzyme(this_call = this_call, min_len = min_len, max_len = max_len, 
                     fasta = fasta, out_path = out_path, mgf_path = mgf_path, 
                     noenzyme_maxn = noenzyme_maxn, quant = quant, 
                     silac_noenzyme = if (!is.null(silac_mix)) TRUE else FALSE, 
                     groups_noenzyme = if (!is.null(par_groups)) TRUE else FALSE)

    return(NULL)
  }

  ## Mixed SILAC
  exec_silac_mix <- if (isTRUE(dots$bypass_silac_mix)) FALSE else TRUE
  
  if (length(silac_mix) && exec_silac_mix) {
    if (!is.null(fixedlabs)) {
      stop("Arguments \"fixedlabs\" and \"silac_mix\" both are non-NULL.\n", 
           "  Set up \"fixedlabs\" under \"silac_mix\" for SILAC;\n", 
           "  Use directly \"fixedlabs\" for direct searches with labels.\n", 
           "The same applies to \"varlabs\".")
    }
    
    matchMS_silac_mix(silac_mix = silac_mix, 
                      this_call = this_call, 
                      out_path = out_path, 
                      mgf_path = mgf_path, 
                      aa_masses = aa_masses)

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
    
    .savecall <- TRUE
    
    return(df)
  }
  
  ## Theoretical MS1 masses
  if (is.null(bypass_pepmasses <- dots$bypass_pepmasses)) 
    bypass_pepmasses <- FALSE

  ## temporary fix
  maxn_fasta_seqs <- 300000L
  mod_motifs <- NULL
  par_groups = NULL
  tmt_reporter_lower <- 126.1, 
  tmt_reporter_upper <- 135.2, 
  ppm_reporters <- 10L,
  ##
  
  if (!bypass_pepmasses)
    res <- calc_pepmasses2(
      aa_masses = aa_masses, 
      fasta = fasta,
      acc_type = acc_type,
      acc_pattern = acc_pattern,
      fixedmods = fixedmods,
      varmods = varmods,
      rm_dup_term_anywhere = rm_dup_term_anywhere, 
      fixedlabs = fixedlabs, 
      varlabs = varlabs, 
      mod_motifs = mod_motifs, 
      enzyme = enzyme,
      custom_enzyme = custom_enzyme, 
      noenzyme_maxn = noenzyme_maxn, 
      maxn_fasta_seqs = maxn_fasta_seqs,
      maxn_vmods_setscombi = maxn_vmods_setscombi,
      maxn_vmods_per_pep = maxn_vmods_per_pep,
      maxn_sites_per_vmod = maxn_sites_per_vmod,
      min_len = min_len,
      max_len = max_len,
      max_miss = max_miss,
      min_mass = min_mass, 
      max_mass = max_mass, 
      out_path = out_path,
      digits = digits,
      use_ms1_cache = use_ms1_cache, 
      .path_cache = .path_cache, 
      .path_fasta = .path_fasta, 
      .path_ms1masses = .path_ms1masses)

  ## Bin theoretical peptides
  if (is.null(bypass_bin_ms1 <- dots$bypass_bin_ms1)) 
    bypass_bin_ms1 <- FALSE
  
  reframe_mgfs <- calib_ms1mass && ppm_ms1calib != ppm_ms1

  if (!bypass_bin_ms1) {
    .path_bin <- 
      bin_ms1masses(res = res, 
                    min_mass = min_mass, 
                    max_mass = max_mass, 
                    min_len = min_len,
                    max_len = max_len,
                    ppm_ms1 = ppm_ms1, 
                    use_ms1_cache = use_ms1_cache, 
                    .path_cache = .path_cache, 
                    .path_ms1masses = .path_ms1masses, 
                    enzyme = enzyme, 
                    out_path = out_path)
    
    .path_bin_calib <- if (reframe_mgfs)
      bin_ms1masses(res = res, 
                    min_mass = min_mass, 
                    max_mass = max_mass, 
                    min_len = min_len,
                    max_len = max_len,
                    ppm_ms1 = ppm_ms1calib, 
                    use_ms1_cache = use_ms1_cache, 
                    .path_cache = .path_cache, 
                    .path_ms1masses = .path_ms1masses, 
                    enzyme = enzyme, 
                    out_path = out_path)
    else
      .path_bin

    if (exists("res"))
      rm(list = "res")
  }

  ## MGFs
  if (reproc_dda_ms1) {
    is_mdda <- TRUE
    
    if (maxn_mdda_precurs > 1L)
      warning("Coerce to `maxn_mdda_precurs = 1` at `reproc_dda_ms1 = TRUE`. ")
    
    maxn_mdda_precurs = 1L
  }
  
  if (is.null(bypass_mgf <- dots$bypass_mgf)) 
    bypass_mgf <- FALSE
  
  if (!bypass_mgf)
    load_mgfs(out_path = out_path, 
              mgf_path = mgf_path,
              min_mass = min_mass,
              max_mass = max_mass, 
              min_ms2mass = min_ms2mass,
              max_ms2mass = max_ms2mass,
              topn_ms2ions = topn_ms2ions,
              min_ms1_charge = min_ms1_charge, 
              max_ms1_charge = max_ms1_charge, 
              min_scan_num = min_scan_num, 
              max_scan_num = max_scan_num, 
              min_ret_time = min_ret_time,
              max_ret_time = max_ret_time, 
              ppm_ms1 = ppm_ms1, 
              ppm_ms2 = ppm_ms2,
              mgf_cutmzs = mgf_cutmzs, 
              mgf_cutpercs = mgf_cutpercs, 
              enzyme = enzyme, 
              exclude_reporter_region = exclude_reporter_region, 
              tmt_reporter_lower = tmt_reporter_lower, 
              tmt_reporter_upper = tmt_reporter_upper, 
              index_mgf_ms2 = index_mgf_ms2, 
              is_mdda = is_mdda, 
              deisotope_ms2 = deisotope_ms2, 
              max_ms2_charge = max_ms2_charge, 
              use_defpeaks = use_defpeaks, 
              maxn_dia_precurs = maxn_dia_precurs, 
              maxn_mdda_precurs = maxn_mdda_precurs, 
              n_mdda_flanks = n_mdda_flanks, 
              ppm_ms1_deisotope = ppm_ms1_deisotope, 
              ppm_ms2_deisotope = ppm_ms2_deisotope, 
              grad_isotope = grad_isotope, 
              fct_iso2 = fct_iso2, 
              quant = quant, 
              digits = digits)

  ## MSMS matches
  if (is.null(bypass_ms2match <- dots$bypass_ms2match)) 
    bypass_ms2match <- FALSE

  if (length(.time_stamp <- find_ms1_times(out_path)) == 1L) {
    path_time <- file.path(.path_ms1masses, .time_stamp)
    file_aams <- file.path(path_time, "aa_masses_all.rds")
    file_mods <- file.path(path_time, "mod_indexes.txt")
    aa_masses_all <- qs::qread(file_aams)
    mod_indexes <- find_mod_indexes(file_mods)
    
    file.copy(file_aams, file.path(out_path, "aa_masses_all.rds"), overwrite = TRUE)
    file.copy(file_mods, file.path(out_path, "mod_indexes.txt"), overwrite = TRUE)
    rm(list = c("path_time", "file_aams", "file_mods"))
  }
  else {
    # only with group searches (low priority)
    aa_masses_all <- NULL
    mod_indexes <- NULL
  }
  
  if (calib_ms1mass)
    calib_mgf(mgf_path = mgf_path, aa_masses_all = aa_masses_all[1], # base
              out_path = out_path, .path_bin = .path_bin_calib, 
              mod_indexes = mod_indexes[names(mod_indexes) %in% fixedmods], 
              type_ms2ions = type_ms2ions, 
              maxn_vmods_per_pep = maxn_vmods_per_pep,
              maxn_sites_per_vmod = maxn_sites_per_vmod, 
              maxn_fnl_per_seq = maxn_fnl_per_seq, 
              maxn_vnl_per_seq = maxn_vnl_per_seq, 
              maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep,
              minn_ms2 = minn_ms2, 
              ppm_ms1 = ppm_ms1calib, 
              reframe_mgfs = reframe_mgfs, 
              ppm_ms2 = ppm_ms2, min_mass = min_mass, max_mass = max_mass, 
              min_ms2mass = min_ms2mass, quant = quant, 
              ppm_reporters = ppm_reporters, index_mgf_ms2 = index_mgf_ms2, 
              by_modules = by_modules, fasta = fasta, acc_type = acc_type, 
              acc_pattern = acc_pattern, topn_ms2ions = topn_ms2ions, 
              fixedmods = fixedmods, varmods = NULL, # the first search
              enzyme = enzyme, maxn_fasta_seqs = maxn_fasta_seqs, 
              maxn_vmods_setscombi = maxn_vmods_setscombi,
              min_len = min_len, max_len = max_len, max_miss = max_miss)

  if (!bypass_ms2match) {
    if (min_ms2mass < 5L) 
      warning("Maybe out of RAM at \"min_ms2mass < 5L\".")
    
    ms2match(mgf_path = mgf_path,
             aa_masses_all = aa_masses_all,
             out_path = out_path,
             .path_bin = .path_bin, 
             mod_indexes = mod_indexes,
             type_ms2ions = type_ms2ions,
             maxn_vmods_per_pep = maxn_vmods_per_pep,
             maxn_sites_per_vmod = maxn_sites_per_vmod,
             maxn_fnl_per_seq = maxn_fnl_per_seq, 
             maxn_vnl_per_seq = maxn_vnl_per_seq, 
             maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep,
             minn_ms2 = minn_ms2,
             ppm_ms1 = ppm_ms1,
             reframe_mgfs = FALSE, 
             ppm_ms2 = ppm_ms2,
             min_mass = min_mass, 
             max_mass = max_mass, 
             min_ms2mass = min_ms2mass,
             quant = quant,
             ppm_reporters = ppm_reporters,
             index_mgf_ms2 = index_mgf_ms2, 
             by_modules = by_modules, 
             ms1_offsets = ms1_offsets, 
             ms1_neulosses = ms1_neulosses, 
             maxn_neulosses_fnl = maxn_neulosses_fnl, 
             maxn_neulosses_vnl = maxn_neulosses_vnl, 
             deisotope_ms2 = deisotope_ms2, 

             # dummy for argument matching
             fasta = fasta,
             acc_type = acc_type,
             acc_pattern = acc_pattern,
             topn_ms2ions = topn_ms2ions,
             fixedmods = fixedmods,
             varmods = varmods,
             enzyme = enzyme,
             maxn_fasta_seqs = maxn_fasta_seqs,
             maxn_vmods_setscombi = maxn_vmods_setscombi,
             min_len = min_len,
             max_len = max_len,
             max_miss = max_miss)
  }

  ## Peptide scores
  if (is.null(bypass_from_pepscores <- dots$bypass_from_pepscores)) 
    bypass_from_pepscores <- FALSE

  if (bypass_from_pepscores) 
    return(NULL)
  
  if (is.null(bypass_pepscores <- dots$bypass_pepscores)) 
    bypass_pepscores <- FALSE
  
  if (!bypass_pepscores) {
    if (is.null(tally_ms2ints <- dots$tally_ms2ints)) 
      tally_ms2ints <- TRUE
    
    calc_pepscores(topn_ms2ions = topn_ms2ions,
                   type_ms2ions = type_ms2ions,
                   target_fdr = target_fdr,
                   min_len = min_len,
                   max_len = max_len,
                   ppm_ms2 = ppm_ms2,
                   soft_secions = soft_secions, 
                   out_path = out_path,
                   min_ms2mass = min_ms2mass,
                   index_mgf_ms2 = index_mgf_ms2, 
                   tally_ms2ints = tally_ms2ints, 

                   # dummies
                   mgf_path = mgf_path,
                   maxn_vmods_per_pep = maxn_vmods_per_pep,
                   maxn_sites_per_vmod = maxn_sites_per_vmod,
                   maxn_vmods_sitescombi_per_pep = 
                     maxn_vmods_sitescombi_per_pep,
                   minn_ms2 = minn_ms2,
                   ppm_ms1 = ppm_ms1,
                   quant = quant,
                   ppm_reporters = ppm_reporters,
                   fasta = fasta,
                   acc_type = acc_type,
                   acc_pattern = acc_pattern,
                   fixedmods = fixedmods,
                   varmods = varmods,
                   enzyme = enzyme,
                   maxn_fasta_seqs = maxn_fasta_seqs,
                   maxn_vmods_setscombi = maxn_vmods_setscombi,
                   add_ms2theos = add_ms2theos, 
                   add_ms2theos2 = add_ms2theos2, 
                   add_ms2moverzs = add_ms2moverzs, 
                   add_ms2ints = add_ms2ints,
                   by_modules = by_modules, 
                   digits = digits)
  }
  
  if (is.null(bypass_primatches <- dots$bypass_primatches)) 
    bypass_primatches <- FALSE
  
  if (!bypass_primatches)
    hadd_primatches(out_path = out_path, 
                    is_notched = is_notched, 
                    add_ms2theos = add_ms2theos, 
                    add_ms2theos2 = add_ms2theos2, 
                    add_ms2moverzs = add_ms2moverzs, 
                    add_ms2ints = add_ms2ints, 
                    by_modules = by_modules, 
                    index_mgf_ms2 = index_mgf_ms2)

  ## Peptide FDR 
  if (is.null(bypass_pepfdr <- dots$bypass_pepfdr)) 
    bypass_pepfdr <- FALSE
  
  if (!bypass_pepfdr) {
    prob_cos <- calc_pepfdr(target_fdr = target_fdr, 
                            fdr_type = fdr_type, 
                            min_len = min_len, 
                            max_len = max_len, 
                            is_notched = is_notched, 
                            max_pepscores_co = max_pepscores_co, 
                            min_pepscores_co = min_pepscores_co, 
                            enzyme = enzyme, 
                            fdr_group = fdr_group, 
                            nes_fdr_group = nes_fdr_group, 
                            out_path = out_path)
    
    ans <- post_pepfdr(prob_cos, out_path)

    if (svm_reproc) {
      message("SVM reprocessing of peptide probabilities.")
      
      prob_cos <- perco_svm(out_path = out_path, df = ans, prob_cos = prob_cos, 
                            target_fdr = target_fdr, fdr_type = fdr_type, 
                            min_len = min_len, max_len = max_len, 
                            max_pepscores_co = max_pepscores_co, 
                            min_pepscores_co = min_pepscores_co, enzyme = enzyme, 
                            fdr_group = fdr_group, nes_fdr_group = nes_fdr_group, 
                            svm_kernel = svm_kernel, svm_feats = svm_feats, 
                            cross_valid = svm_cv, k  = svm_k, 
                            costs = svm_costs, 
                            def_cost = svm_def_cost, 
                            svm_iters = svm_iters)

      # post_pepfdr(prob_cos, out_path)
      message("Completed SVM reprocessing.")
    }

    rm(list = c("ans", "prob_cos"))
  }

  ## Peptide ranks and score deltas between `pep_ivmod`
  if (is.null(bypass_peploc <- dots$bypass_peploc)) 
    bypass_peploc <- FALSE
  
  if (!bypass_peploc) {
    calc_peploc(out_path = out_path, 
                mod_indexes = mod_indexes, 
                locmods = locmods, 
                is_notched = is_notched,
                topn_mods_per_seq = topn_mods_per_seq, 
                topn_seqs_per_query = topn_seqs_per_query)
  }

  ## Protein accessions
  if (is.null(bypass_from_protacc <- dots$bypass_from_protacc)) 
    bypass_from_protacc <- FALSE
  
  if (bypass_from_protacc) 
    return(NULL)

  if (is.null(bypass_protacc <- dots$bypass_protacc)) 
    bypass_protacc <- FALSE
  
  temp_dir <- file.path(out_path, "temp")
  file_protacc <- file.path(temp_dir, "df_protacc.rds")
  
  if (bypass_protacc && file.exists(file_protacc))
    df <- qs::qread(file_protacc)
  else {
    if (enzyme != "noenzyme" || isTRUE(dots[["direct_prot_acc"]]))
      df <- add_protacc(out_path = out_path, 
                        .path_cache = .path_cache, 
                        .path_fasta = .path_fasta)
    else {
      silac_noenzyme <- if (isTRUE(dots$silac_noenzyme)) TRUE else FALSE
      
      # see matchMS_noenzyme for nested silac under noenzyme
      df <- if (silac_noenzyme)
        add_protacc(out_path = out_path, 
                    .path_cache = .path_cache, 
                    .path_fasta = .path_fasta)
      else
        add_protacc2(out_path = out_path, 
                     .path_cache = .path_cache, 
                     .path_fasta = .path_fasta)
      
      rm(list = c("silac_noenzyme"))
    }
    
    qs::qsave(df, file_protacc, preset = "fast")
  }
  
  rm(list = "file_protacc")
  
  ## Protein FDR
  if (is.null(bypass_protfdr <- dots$bypass_protfdr)) 
    bypass_protfdr <- FALSE
  
  file_protfdr <- file.path(temp_dir, "df_protfdr.rds")
  
  if (bypass_protfdr && file.exists(file_protfdr)) {
    df <- qs::qread(file_protfdr)
  }
  else {
    df <- calc_protfdr(df = df, 
                       target_fdr = target_fdr, 
                       max_protscores_co = max_protscores_co, 
                       max_protnpep_co = max_protnpep_co, 
                       method_prot_es_co = method_prot_es_co, 
                       out_path = out_path)
    qs::qsave(df, file_protfdr, preset = "fast")
  }
  
  df <- add_rptrs(df, quant, out_path)

  ## Clean-ups
  # (raw_file etc. already mapped if `from_group_search`)
  if (!isTRUE(from_group_search <- dots$from_group_search)) 
    df <- map_raw_n_scan(df, mgf_path)
  
  df <- dplyr::mutate(df, pep_expect = 10^((pep_score_co - pep_score)/10) * target_fdr)
  df[["pep_score_co"]] <- NULL
  df$pep_delta <- df$pep_exp_mr - df$pep_calc_mr

  nms <- names(df)
  df  <- dplyr::bind_cols(
    df[grepl("^prot_", nms)],
    df[grepl("^pep_", nms)],
    df[grepl("^psm_", nms)],
    df[!grepl("^prot_|^pep_|^psm_", nms)], )
  rm(list = "nms")
  
  df <- reloc_col_after(df, "pep_exp_z", "pep_exp_mr")
  df <- reloc_col_after(df, "pep_calc_mr", "pep_exp_z")
  df <- reloc_col_after(df, "pep_delta", "pep_calc_mr")
  
  # e.g. realization with Acetyl (K) but no TMT (K)
  cols_tmt <- grepl("^I[0-9]{3}[Nc]{0,1}", names(df))
  rows_tmt <- grepl("TMT", df[["pep_fmod"]]) | grepl("TMT", df[["pep_vmod"]])
  df[!rows_tmt, cols_tmt] <- NA_real_
  rm(list = c("cols_tmt", "rows_tmt"))
  
  local({
    df$pep_exp_mz  <- round(df$pep_exp_mz, digits = 4L)
    df$pep_exp_mr  <- round(df$pep_exp_mr, digits = 4L)
    df$pep_calc_mr <- round(df$pep_calc_mr, digits = 4L)
    df$pep_delta   <- round(df$pep_delta, digits = 4L)
    df$pep_tot_int <- round(df$pep_tot_int, digits = 1L)
    df$pep_expect  <- format(df$pep_expect, digits = 3L)
    
    readr::write_tsv(df, file.path(out_path, "psmC.txt"))
    session_info <- sessionInfo()
    save(session_info, file = file.path(out_path, "Calls", "mzion.rda"))
  })

  ## psmC to psmQ
  df <- df[, c("prot_acc", "pep_seq", "pep_issig", "pep_isdecoy", 
               "prot_issig", "prot_n_pep")]
  
  df <- dplyr::filter(df, pep_issig, !pep_isdecoy, !grepl("^-", prot_acc))

  df <- try_psmC2Q(df, 
                   out_path = out_path,
                   fdr_type = fdr_type, # for workflow controls 
                   combine_tier_three = combine_tier_three, 
                   max_n_prots = max_n_prots)

  message("Completed at: ", Sys.time())
  
  .savecall <- TRUE

  invisible(df)
}


#' Helper of \link{psmC2Q}.
#' 
#' "n_peps" and "n_prots" including both targets and decoys: \cr
#' "n_prots" about 1:1 \cr
#' "n_peps" about 1.8:1
#' 
#' @inheritParams psmC2Q
#' @importFrom magrittr %>% %T>%
try_psmC2Q <- function (df = NULL, out_path = NULL, fdr_type = "protein",
                        combine_tier_three = FALSE, max_n_prots = 60000L) 
{
  n_peps <- length(unique(df$pep_seq))
  n_prots <- length(unique(df$prot_acc))
  
  if (n_prots == 1L) {
    message("No grouping with the number of of proteins = ", n_prots, ".\n",
            "Search completed successfully.")
    options(show.error.messages = FALSE)
    stop()
  }
  
  if (n_peps > 1000000L && n_prots > 100000L)
    df <- NA
  else
    df <- tryCatch(
      psmC2Q(df,
             out_path = out_path,
             fdr_type = fdr_type,
             combine_tier_three = combine_tier_three, 
             max_n_prots = max_n_prots),
      error = function(e) NA)

  if (length(df) == 1L && is.na(df)) {
    message("Retry with a new R session: \n\n",
            "Manual execution of the following codes if not start automatically.\n\n", 
            "mzion::reproc_psmC(\n",
            "  out_path = \"", out_path, "\",\n",
            "  fdr_type = \"", fdr_type, "\",\n",
            "  combine_tier_three  = ", combine_tier_three, ",\n",
            "  max_n_prots  = ", max_n_prots, "\n",
            ")\n")
    
    fileConn <- file(file.path("~/post_psmC.R"))
    
    lines <- c(
      "library(mzion)\n",
      "mzion::reproc_psmC(",
      paste0("  out_path = \"", out_path, "\","),
      paste0("  fdr_type = \"", fdr_type, "\","),
      paste0("  combine_tier_three = ", combine_tier_three, ","),
      paste0("  max_n_prots = ", max_n_prots),
      ")\n",
      "unlink(\"~/post_psmC.R\")"
    )
    
    writeLines(lines, fileConn)
    close(fileConn)
    
    rstudioapi::restartSession(command = 'source("~/post_psmC.R")')
  } 
  else {
    suppressWarnings(
      rm(list = c(".path_cache", ".path_ms1masses", ".time_stamp"), 
         envir = .GlobalEnv))

    message("Done.")
  }
  
  invisible(df)
}


#' Reprocessing of \code{psmC.txt}.
#'
#' Protein grouping from \code{psmC.txt} to \code{psmQ.txt}.
#'
#' May solve some memory shortage issues for large data sets by restarting An
#' Rstudio session.
#'
#' The score cut-offs are different among the \code{fdr_type} of "psm",
#' "peptide" and "protein". An experimenter need to match the value of
#' \code{fdr_type}.
#'
#' @param fct A factor for data splitting into chunks. May consider a greater
#'   value for a larger data set.
#' @inheritParams matchMS
#' @export
reproc_psmC <- function (out_path = NULL, fdr_type = "protein",
                         combine_tier_three = FALSE, max_n_prots = 60000L, 
                         fct = 4L) 
{
  if (is.null(out_path)) 
    stop("`out_path` cannot be NULL.", call. = FALSE)

  message("Leave the session open and wait for the `Search completed` message.")

  df <- suppressWarnings(
    readr::read_tsv(file.path(out_path, "psmC.txt"), 
                    col_types = get_mzion_coltypes()))

  df <- df[, c("pep_seq", "prot_acc", "prot_issig", "prot_n_pep",
               "pep_issig", "pep_isdecoy")]
  
  df <- dplyr::filter(df, pep_issig, !pep_isdecoy, !grepl("^-", prot_acc))
  gc()

  psmC2Q(df, out_path = out_path,
         fdr_type = fdr_type,
         combine_tier_three = combine_tier_three, 
         max_n_prots = max_n_prots, 
         fct = fct)

  message("Done.")
}


#' From \code{psmC.txt} to \code{psmQ.txt}.
#'
#' Non-significant and decoy peptides should have been removed from the input
#' \code{df}, as well as decoy proteins.
#'
#' @param df A result of \code{psmC.txt} with the removals of non-significant
#'   or decoy peptides, as well as decoy proteins.
#' @param fct A factor for data splitting into chunks. May consider a greater
#'   value for a larger data set.
#' @inheritParams matchMS
#' @importFrom fastmatch %fin%
psmC2Q <- function (df = NULL, out_path = NULL, fdr_type = "protein",
                    combine_tier_three = FALSE, max_n_prots = 60000L, 
                    fct = 4L) 
{
  options(warn = 1L)
  
  # if (!all(df[["pep_issig"]])) stop("Developer: filter data by \"pep_issig\" first.")
  # if (any(df[["pep_isdecoy"]])) stop("Developer: remove decoy peptide first.")
  # if (any(grepl("^-", df["prot_acc"]))) stop("Developer: remove decoy proteins first.")

  message("\n=================================\n",
          "prot_tier  prot_issig  prot_n_pep \n",
          "    1          [y]          \n",
          "    2          [n]          > 1\n",
          "    3          [n]          = 1\n",
          "=================================\n")
  
  # Set aside one-hit wonders
  df3 <- dplyr::filter(df, !prot_issig, prot_n_pep == 1L)
  df3 <- dplyr::mutate(df3, prot_tier = 3L)

  df <- dplyr::bind_rows(
    dplyr::filter(df, prot_issig),
    dplyr::filter(df, !prot_issig, prot_n_pep >= 2L))

  df <- dplyr::mutate(df, prot_tier = ifelse(prot_issig, 1L, 2L))

  # the same peptide can be present in all three protein tiers; 
  # steps up if pep_seq(s) in tier 3 also in tiers 1, 2
  if (FALSE) {
    rows <- df3$pep_seq %in% df$pep_seq
    df <- dplyr::bind_rows(df, df3[rows, ])
    df3 <- df3[!rows, ]
    rm(list = "rows")
  }
  
  # Protein groups
  message("Building protein-peptide maps.")
  
  len_prots <- length(unique(df$prot_acc))
  
  if (len_prots > max_n_prots && fdr_type != "protein") {
    warning("Large number of proteins at ", len_prots, ".\n", 
            "Coerce to `fdr_type = protein` ",
            "and save peptide results of tier-2 proteins in `psmT2.txt`.",
            call. = FALSE)
    
    fdr_type <- "protein"
    
    df2 <- dplyr::filter(df, prot_tier == 2L)
    df  <- dplyr::filter(df, prot_tier == 1L)
  } 
  else {
    if (fdr_type == "protein") {
      df2 <- dplyr::filter(df, prot_tier == 2L)
      df  <- dplyr::filter(df, prot_tier == 1L)
    } 
    else {
      message("No tier-2 outputs at `fdr_type = ", fdr_type, "`.")
      
      if (len_prots > max_n_prots) {
        warning("The number of proteins is ", len_prots, ".\n", 
                "Consider `fdr_type = protein`.",
                call. = FALSE)
      }
      
      df2 <- df[0, ]
      df <- df # prot_tiers: 1 + 2
    }
  }
  
  # df may have both prot_tier 1 and 2 if fdr_type != "protein"
  df_tier12 <- unique(df[, c("prot_acc", "prot_tier")])
  
  df  <- unique(df [, c("prot_acc", "pep_seq")])
  df2 <- unique(df2[, c("prot_acc", "pep_seq")])
  df3 <- unique(df3[, c("prot_acc", "pep_seq")])

  nms <- c("prot_acc", "pep_seq", "prot_isess", "prot_hit_num", 
           "prot_family_member", "pep_literal_unique", "pep_razor_unique")
  
  if (nrow(df)) {
    df <- groupProts(df, out_path = file.path(out_path, "temp1"), fct = fct)
    df <- dplyr::left_join(df, df_tier12, by = "prot_acc")
  }
  else
    df <- make_zero_df(nms)

  df2 <- if (nrow(df2)) 
    groupProts(df2, out_path = file.path(out_path, "temp2"), fct = fct)
  else 
    make_zero_df(nms)
  
  df3 <- if (nrow(df3)) 
    groupProts(df3, out_path = file.path(out_path, "temp3"), fct = fct)
  else 
    make_zero_df(nms)
  
  rm(list = c("nms", "df_tier12"))

  # Cleanup
  dfC <- suppressWarnings(
    read_tsv(file.path(out_path, "psmC.txt"), col_types = get_mzion_coltypes()))
  dfC <- dplyr::filter(dfC, pep_issig, !pep_isdecoy, !grepl("^-", prot_acc))
  dfC <- tidyr::unite(dfC, uniq_id, prot_acc, pep_seq, sep = ".", remove = FALSE)

  df  <- post_psmC2Q(df,  dfC, tier = NULL)
  df2 <- post_psmC2Q(df2, dfC, tier = 2L)
  df3 <- post_psmC2Q(df3, dfC, tier = 3L)
  
  rm(list = "dfC")

  # Three-tier combines
  nms_df <- names(df)
  df2 <- df2[, nms_df]
  df3 <- df3[, nms_df]
  rm(list = "nms_df")
  
  max <- max(df$prot_hit_num, na.rm = TRUE)
  
  if (fdr_type == "protein" && combine_tier_three) {
    warning("Coerce to `combine_tier_three = FALSE` at `fdr_type = protein`.",
            call. = FALSE)
    combine_tier_three <- FALSE
  }
  
  if (combine_tier_three) {
    df3 <- df3[!df3[["pep_seq"]] %fin% df[["pep_seq"]], ]
    df <- dplyr::bind_rows(list(df, df3)) # df2 should have no rows
    df <- dplyr::arrange(df, prot_acc, pep_seq)
    readr::write_tsv(df, file.path(out_path, "psmQ.txt"))

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
    df <- dplyr::arrange(df, prot_acc, pep_seq)
    readr::write_tsv(df, file.path(out_path, "psmQ.txt"))

    if (nrow(df2)) {
      df2 <- dplyr::mutate(df2[names(df)], prot_hit_num = prot_hit_num + max)
      readr::write_tsv(df2, file.path(out_path, "psmT2.txt"))
      max <- max(df2[["prot_hit_num"]], na.rm = TRUE)
    }
    
    if (nrow(df3)) {
      df3 <- dplyr::mutate(df3[names(df)], prot_hit_num = prot_hit_num + max)
      readr::write_tsv(df3, file.path(out_path, "psmT3.txt"))
    }
  }

  #  No pepQ.txt and prnQ.txt; use proteoQ for data mining

  invisible(df)
}


#' Post \link{psmC2Q}.
#'
#' @param df A data frame of protein-peptide map.
#' @param dfC A \code{psmQ} data with the removal of non-significant peptides
#'   etc.
#' @param tier The tier of proteins in \code{df}.
post_psmC2Q <- function (df, dfC, tier = NULL) 
{
  if (!is.null(tier))
    df <- dplyr::mutate(df, prot_tier = tier)

  df <- tidyr::unite(df, uniq_id, prot_acc, pep_seq, sep = ".", remove = TRUE)
  df <- dplyr::left_join(df, dfC, by = "uniq_id")
  df <- dplyr::select(df, -uniq_id)

  ord_prots <- c("prot_acc", "prot_issig")
  
  df <- dplyr::bind_cols(
    df[, ord_prots, drop = FALSE], 
    df[, !names(df) %in% ord_prots, drop = FALSE])

  ord_peps <- c("pep_seq", "pep_issig", "pep_literal_unique", 
                "pep_razor_unique", "pep_score", "pep_expect")
  
  df <- dplyr::bind_cols(
    df[, ord_peps, drop = FALSE], 
    df[, !names(df) %in% ord_peps, drop = FALSE])

  df <- dplyr::bind_cols(
    df[grepl("^prot_", names(df))],
    df[grepl("^pep_", names(df))],
    df[grepl("^psm_", names(df))],
    df[!grepl("^prot_|^pep_|^psm_", names(df))])
  
  
  df$pep_exp_mz  <- round(df$pep_exp_mz, digits = 4L)
  df$pep_exp_mr  <- round(df$pep_exp_mr, digits = 4L)
  df$pep_calc_mr <- round(df$pep_calc_mr, digits = 4L)
  df$pep_delta   <- round(df$pep_delta, digits = 4L)
  df$pep_tot_int <- round(df$pep_tot_int, digits = 1L)
  df$pep_expect  <- format(df$pep_expect, digits = 3L)

  df <- dplyr::select(df, -which(names(df) %in% c("prot_n_psm", "prot_n_pep")))
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
  
  if (!grepl("^tmt[0-9]+", quant))
    return(NULL)
  
  tmts <- fvmods[grepl("^TMT", fvmods)]
  
  if (quant == "tmt18") {
    ok <- all(grepl("TMTpro18.* |TMT18plex.* ", tmts))
    
    if (!ok) 
      warning("All TMT modifications need to be `TMTpro18` or `TMT18plex` at `", 
              quant, "`.\n", 
              tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3)
  } 
  else if (quant == "tmt16") {
    ok <- all(grepl("TMTpro.* |TMT16plex.* ", tmts))
    
    if (!ok) 
      warning("All TMT modifications need to be `TMTpro` or `TMT16plex` at `", 
              quant, "`.\n", 
              tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3)
  } 
  else {
    ok <- all(grepl("TMT6plex.* |TMT10plex.* |TMT11plex.* ", tmts))
    
    if (!ok) 
      warning("All TMT modifications need to be `TMT6plex`, `TMT10plex` or `TMT11plex` at `", 
              quant, "`.\n", 
              tmt_msg_1, "\n", tmt_msg_2, "\n", tmt_msg_3)
  }
  
  invisible(NULL)
}


#' Checks the path of MGF files
#' 
#' @param error Character string; the level of error.
#' @inheritParams matchMS
#' @inheritParams matchMS_par_groups
checkMGF <- function (mgf_path = NULL, grp_args = NULL, error = c("stop", "warn")) 
{
  mgf_path <- find_dir(mgf_path)
  error <- match.arg(error)
  
  if (! error %in% c("warn", "stop"))
    stop("`error` needs to be one of \"error\" or \"stop\".")
  
  if (is.null(mgf_path)) 
    stop("`mgf_path` not found.")
  
  fi_mgf <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mgf$")
  fi_mzml <- list.files(path = file.path(mgf_path), pattern = "^.*\\.mzML$")
  len_mgf <- length(fi_mgf)
  len_mzml <- length(fi_mzml)
  
  if (len_mgf && len_mzml)
    stop("Peak lists need to be in either MGF or mzML, but not both.")
  
  if (!(len_mgf || len_mzml)) {
    if (error == "warn")
      warning("No `.mgf` files immediately under ", mgf_path)
    else
      stop("No `.mgf` files immediately under ", mgf_path)
  }
  
  invisible(mgf_path)
}


#' Checks \code{locmods}
#' 
#' Coerced \code{fixedmods} not considered.
#' 
#' @inheritParams matchMS
check_locmods <- function (locmods, fixedmods, varmods, ms1_neulosses = NULL)
{
  if (!length(locmods))
    return(NULL)
  
  if (!is.null(ms1_neulosses)) {
    if (sum(bads <- !ms1_neulosses %in% locmods)) {
      warning("\nPLEASE READ: \n\n", 
              "\"ms1_neulosses\": ", paste(ms1_neulosses[bads], collapse = ", "), 
              " not found in \"locmods\": ", paste(locmods, collapse =, ""), 
              "\n!!! Consider matching some of the \"locmods\" setting to", 
              " \"ms1_neulosses\". !!!\n")
      
      if (FALSE) {
        nls <- lapply(ms1_neulosses, find_unimod)
        nlresids <- lapply(nls, `[[`, "position_site")
        nlresids <- unlist(nlresids, use.names = FALSE, recursive = FALSE)
        
        vs <- lapply(varmods, find_unimod)
        vresids <- lapply(vs, `[[`, "position_site")
        vresids <- unlist(vresids, use.names = FALSE, recursive = FALSE)
      }
    }
  }

  if (!all(oks <- locmods %in% c(fixedmods, varmods))) {
    warning("Ignore \"locmods\" not in \"varmods\" or \"fixedmods\": ", 
            paste(locmods[!oks], collapse = ", "))
    locmods <- locmods[oks]
  }

  if (!length(locmods))
    return(NULL)

  # locmods are only among fixedmods
  if (!length(vids <- which(varmods %in% locmods))) 
    stop("No \"varmods\" matched to \"locmods\": ", paste(locmods, collapse = ", "))
  
  vmods <- find_modps(varmods)
  vsites <- unlist(vmods[vids], recursive = FALSE, use.names = FALSE)
  fmods <- find_modps(fixedmods)
  fids <- which(fmods %in% vsites)
  
  if (length(fids))
    locmods <- c(fixedmods[fids], locmods)

  locmods
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
    raws <- qs::qread(file_raw)
    raws2 <- names(raws)
    names(raws2) <- raws
    df$raw_file <- unname(raws2[df$raw_file])
  }
  else {
    stop("File not found: ", file_raw)
  }
  
  if (file.exists(file_scan)) {
    scans <- qs::qread(file_scan)
    scans2 <- names(scans)
    names(scans2) <- scans
    df$pep_scan_title <- unname(scans2[df$pep_scan_title])
  }
  else {
    stop("File not found: ", file_scan)
  }
  
  invisible(df)
}


#' Checks the values of \code{fdr_group}
#' 
#' Not yet used. Takes values of integers or character strings.
#' 
#' @param oks A vector of allowed modification groups.
#' @inheritParams matchMS
check_fdr_group <- function (fdr_group = c("base", "all", "top3"), 
                             oks = c("base", "all"))
{
  is_trivial <- all(is.null(fdr_group)) || all(is.na(fdr_group)) || 
    all(fdr_group == "")
  
  if (is_trivial)
    return(oks[[1]])
  
  len  <- length(fdr_group <- unique(fdr_group))
  oks2 <- fdr_group %in% oks
  
  if (len > 1L)
    fdr_group <- if (all(oks2)) oks[1] else fdr_group[!oks2]

  as.character(fdr_group)
}


#' Checks the compatibility between ms1_notches and ms1_neulosses.
#' 
#' @inheritParams matchMS
check_notches <- function (ms1_notches, ms1_neulosses)
{
  n_notches   <- length(ms1_notches)
  n_neulosses <- length(ms1_neulosses)
  
  if (n_neulosses) {
    cdn_1 <- n_notches == 1L && ms1_notches != 0
    cdn_2 <- n_notches > 1L
    
    if (cdn_1 || cdn_2)
      stop("Not support simultaneous non-trival `ms1_notches` and `ms1_neulosses`.")
  }
}


