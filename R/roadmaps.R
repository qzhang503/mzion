## matchMS
#   calc_pepmasses2 (ms1_precursors.R)
#   bin_ms1masses (bin_masses.R)
#   load_mgfs (mgf.R)
#   ms2match (msmsmatches2.R)
#   calc_pepscores (scores.R)
#   calc_pepfdr
#   post_pepfdr
#   calc_peploc (scores)
#   add_protacc (quant2.R)
#   calc_protfdr (scores.R)
#   add_rptrs (quant2.R)
#   try_psmC2Q
#     psmC2Q
#       grp_prots (quant2.R)
# 
# ======================================
# bin_masses.R
#   - bin_ms1masses
#     - (i) binTheoSeqs (small dataset)
#       - bin_theoseqs (-> export)
#         - find_ms1_cutpoints (-> export)
#     - (ii) binTheoSeqs_i (larger dataset; -> export)
#       - binTheoSeqs2 (-> export)
#         - bin_theoseqs (-> export)
#           - find_ms1_cutpoints (-> export)
# 
# mgf.R
#   - readMGF
#    - read_mgf_chunks
#     // START parallel
#     - proc_mgf_chunks_i
#      - proc_mgf_chunks
#       - proc_mgfs
#         - which_topx2
#         - index_mz
#         - sub_mgftopn
#         - integerize_ms2ints
#     // END parallel
# 
# calc_tmtint (quant2.R)
#   find_reporter_ints
#     find_reporters_ppm
#   find_int_cols
# 
# vmod_ms1_labels.R
#   - make_ms1vmod_i
#     - make_ms1_vmodsets
#       - bacth_vmods_combi
#         - make_unique_sets
#           - find_unique_sets
#     - find_intercombi2
#       - expand_grid_rows (utils_engine.R)
# 
# vmod_ms2_labels.R
#     - combi_namesiteU
#       - find_vmodposU
#         - vec_to_list
#         - sim_combn
#     - combi_namesiteM
#       - find_vmodposM
#         - vec_to_list
#         - sim_combn
#       - match_aas_indexes


# ms1_precursors.R: 
#   - calc_pepmasses2
#    - find_aa_masses
#       - calc_aamasses
#         - add_fix_masses
#         - add_var_masses
#         - parse_aamasses
#     - split_fastaseqs
#       - load_fasta2 (fastas.R)
#       - chunksplit (msmsmatches.R)
#         - make_fastapeps0 (-> export)
#           - keep_n_misses (-> export; ms1_precursors.R)
#     - distri_fpeps (fixedmods)
#     - ms1masses_bare
#       - ms1masses_noterm
#         - calcms1mass_noterm (-> export)
#           - calcms1mass_noterm_byprot (-> export)
#             - calcms1mass_noterm_bypep (-> export)
#       - roll_sum (-> export; ms1_precursors.R)
#     - find_motif_pat
#     - distri_peps (varmods)
#       - subpeps_by_vmods (dispatch.R)
#         - find_nmodtree
#           ...
#         - find_cmodtree
#           ...
#       - rm_char_in_nfirst (ms1_precursors.R)
#       - rm_char_in_nlast (ms1_precursors.R)
#     - simple_prots_peps
#     - flat_pepseqs
#     - add_term_mass
#     - helpers below
# 
# helpers at sets of realized modifications: 
#   (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
#     - hms1_a0_vnl0_fnl1
#       - ms1_a0_vnl0_fnl1
#         - expand_grid_rows
#         - delta_ms1_a0_fnl1 (-> export)
#   
#  (7-8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
#    (9-10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
#    (11-12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
#    (13-14) "amods+ tmod- vnl+ fnl+", "amods+ tmod+ vnl+ fnl+"
#         [by nested combinatorial conditions; no explicit functions]
#         - hms1_a1_vnl0_fnl0
#           - ms1_a1_vnl0_fnl0
#             - match_mvmods (vmods_ms1_labels.R)
#             - expand_grid_rows (utils_engine.R)
#             - delta_ms1_a0_fnl1 (ms1_precursors.R)

## ms2match (msmsmatches2.R)
# 
# ms2base.R: (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-"
#   ms2match_base 
#     purge_search_space (utils_engine.R)
#       subset_theoframes (msmsmatches.R)
#     
#       frames_adv (frame-advancing)
#         gen_ms2ions_base (for specific pep_seq)
#           ms2ions_by_type (ion_ladder.R)
#             byions, czions, axions (ion_ladder.R)
#           calc_rev_ms2 (utils_engine.R)
#         search_mgf
#           find_ms2_bypep2
#             fuzzy_match_one
#             fuzzy_match_one2
#       post_frame_adv (utils_engine.R)
# 
# ms2_a0_vnl0_fnl1.R: (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
#   ms2match_a0_vnl0_fnl1 
#     purge_search_space (utils_engine.R)
#     
#       frames_adv (ms2_base.R)
#         gen_ms2ions_a0_vnl0_fnl1
#           // early return
#           gen_ms2ions_base (ms2base.R)
#             calc_rev_ms2 (utils_engine.R)
#             ms2ions_by_type (ion_ladder.R)
#               byions, czions, axions (ion_ladder.R)
#             calc_rev_ms2 (utils_engine.R)
#           // regular return
#           ms2ions_by_type (ion_ladder.R)
#             byions, czions, axions (ion_ladder.R)
#         search_mgf (ms2base.R)
#           find_ms2_bypep2 (ms2base.R)
#             fuzzy_match_one (ms2base.R)
#             fuzzy_match_one2 (ms2base.R)
#       post_frame_adv (utils_engine.R)
# 
# ms2_a1_vnl0_fnl0.R: (7, 8) "amods+ tmod+ vnl- fnl-", "amods+ tmod- vnl- fnl-"
#   "ms2match_a1_vnl0_fnl0"
#     "purge_search_space" (utils_engine.R)
#     // START parallel
#     
#       "frames_adv" (ms2_base.R)
#         "gen_ms2ions_a1_vnl0_fnl0" (ms2_a1_vnl0_fnl0.R)
#           "match_mvmods" (vmods_ms1_labels.R)
#             "expand_grid_rows" (utils_engine.R)
#             "combi_namesiteU" (vmods_ms2_labels.R)
#               "find_vmodposU" (vmods_ms2_labels.R)
#                 "vec_to_list" (utils_engine.R)
#                 "sim_combn" (vmod_ms2_labels.R)
#             "combi_namesiteM" (vmods_ms2_labels.R)
#               "find_vmodposM" (vmods_ms2_labels.R)
#                 "vec_to_list" (vmods_ms2_labels.R)
#                 "sim_combn" (vmod_ms2_labels.R)
#               "match_aas_indexes" (vmods_ms2_labels.R)
#           "check_ms1_mass_vmods2" (ms2_a1_vnl0_fnl0.R)
#           "calc_rev_ms2" (utils_engine.R)
#           "calc_ms2ions_a1_vnl0_fnl0" (ms2_a1_vnl0_fnl0.R)
#             "ms2ions_by_type" (ion_ladder.R)
#               "byions", "czions", "axions"
#                 "bions_base", "yions_base",
#                 "cions_base", "zions_base", 
#                 "aions_base", "xions_base", 
#           "add_hexcodes" (ms2_a1_vnl0_fnl0.R)
#         "search_mgf" (ms2base.R)
#           "find_ms2_bypep2" (ms2base.R)
#             "fuzzy_match_one" (ms2base.R)
#             "fuzzy_match_one2" (ms2base.R)
#       "post_frame_adv" (utils_engine.R)
#     // END parallel
#       (saveRDS ion_matches_[i].rds)
# 
# ms2_a1_vnl1_fnl0.R: (9, 10) "amods+ tmod+ vnl+ fnl-", "amods+ tmod- vnl+ fnl-"
#   "ms2match_a1_vnl1_fnl0" 
#     "purge_search_space" (utils_engine.R)
#     // START parallel
#     
#       "frames_adv" (ms2_base.R)
#         "gen_ms2ions_a1_vnl1_fnl0"
#           "match_mvmods" (vmods_ms1_labels.R)
#             "expand_grid_rows" (utils_engine.R)
#             "combi_namesiteU" (vmods_ms2_labels.R)
#               "find_vmodposU" (vmods_ms2_labels.R)
#                 "vec_to_list" (utils_engine.R)
#                 "sim_combn" (vmod_ms2_labels.R)
#             "combi_namesiteM" (vmods_ms2_labels.R)
#               "find_vmodposM" (vmods_ms2_labels.R)
#                 "vec_to_list" (utils_engine.R)
#                 "sim_combn" (vmod_ms2_labels.R)
#               "match_aas_indexes" (vmods_ms2_labels.R)
#           "check_ms1_mass_vmods2" (ms2_a1_vnl0_fnl0.R)
#           "expand_grid_rows" (utils_engine.R)
#           "calc_rev_ms2" (utils_engine.R)
#           "calc_ms2ions_a1_vnl1_fnl0"
#             "ms2ions_by_type" (ion_ladder.R)
#               "byions", "czions", "axions"
#                 "bions_base", "yions_base",
#                 "cions_base", "zions_base", 
#                 "aions_base", "xions_base", 
#           "add_hexcodes_vnl2"
#         "search_mgf" (ms2base.R)
#           "find_ms2_bypep2" (ms2base.R)
#             "fuzzy_match_one" (ms2base.R)
#             "fuzzy_match_one2" (ms2base.R)
#       "post_frame_adv" (utils_engine.R)
#     // END parallel
# 
# ms2_a1_vnl0_fnl1.R: (11, 12) "amods+ tmod+ vnl- fnl+", "amods+ tmod- vnl- fnl+"
#   ms2match_a1_vnl0_fnl1 
#     purge_search_space (utils_engine.R)
#     
#       frames_adv (ms2_base.R)
#         gen_ms2ions_a1_vnl0_fnl1
#           - match_mvmods (vmods_ms1_labels.R)
##             - combi_namesiteU
#               - find_vmodposU
#                 - "sim_combn" (vmod_ms2_labels.R)
#             - combi_namesiteM
#               - find_vmodposM
#                 - "sim_combn" (vmod_ms2_labels.R)
#               - match_aas_indexes
#           check_ms1_mass_vmods2 (ms2_a1_vnl0_fnl0.R)
#           calc_rev_ms2 (utils_engine.R)
#           calc_ms2ions_a1_vnl0_fnl1
#             ms2ions_by_type (ion_ladder.R)
#               byions, czions, axions
#           add_hexcodes_fnl2
#         search_mgf (ms2base.R)
#           find_ms2_bypep2 (ms2base.R)
#             fuzzy_match_one (ms2base.R)
#             fuzzy_match_one2 (ms2base.R)
#       post_frame_adv (utils_engine.R)
# 

## calc_pepscores (scores.R)
#   calcpepsc
#   calc_pepfdr
#     calc_pepprobs_i
#       scalc_pepprobs
#         calc_probi
#           calc_probi_bypep
#             calc_probi_byvmods
#               add_seions
#               find_ppm_outer_bycombi
#               match_ex2th2
#       add_primatches

## calc_protfdr (scores.R)
#   calc_protfdr_i
#   fit_protfdr

## calc_peploc (scores.R)
#   which_topx3
#   findLocFracsDF
#   concatFracs

## (quant2.R)
#   groupProts
#     map_pepprot
#       pcollapse_sortpeps
#         collapse_sortpeps
#         find_group_breaks
#     cut_proteinGroups
#       as_lgldist
#         proteoCpp::to_lgldistC
#     greedysetcover3
#     

#################################
# mapMS2ions.R
#################################
# mapMS2ions
# match_mgf_path
# match_raw_id
# add_raw_ids
# find_secion_types
# find_psm_rows
# find_psm_rows1
# find_psm_rows2
# find_theoexpt_pair
# find_mgf_query
# combine_prisec_matches
# check_existed_psms
# get_mzion_coltypes

# batch_ms2ions
#   hbatch_ms2ions
#     mgen_ms2ions
#     make_ms2frames
#       make_ms2frames_bypars
#         create_dir
#         find_ms1_cutpoints
#         hash_frame_nums



#################################
# utils_engine.R
#################################
# which_topx
# which_topx2
# get_topn_vals
# insVal
# topx
# find_ppm_error
# find_mass_error_range
# `%+%`
# post_frame_adv
# purge_search_space
# subset_theoframes
# subset_neuloss_peps
# find_nterm_mass
# find_cterm_mass
# detect_cores
# find_free_mem
# find_mod_indexes
# is_equal_sets
# purge_decoys
# expand_grid_rows
# count_elements
# vec_to_list
# split_vec
# accumulate_char
# combi_mat

#################################
# utils_os.R
#################################
# `names_pos<-`
# find_int_cols
# ins_cols_after
# add_cols_at
# replace_cols_at
# reloc_col_after
# reloc_col_after_last
# reloc_col_after_first
# reloc_col_before
# reloc_col_before_last
# reloc_col_before_first
# find_preceding_colnm
# recur_flatten
# chunksplit
# chunksplitLB
# find_dir
# create_dir
# save_call2
# find_callarg_vals
# match_calltime
# delete_files
# is_nulllist
# add_nulllist
#################################

#################################
# utils_ui.R (user interfaces)
#################################
# calc_monopeptide
# calc_monopep
# check_aaseq
# calc_ms2ionseries
# calc_ms2ions
# unique_mvmods
# vmods_elements
# find_intercombi
#################################



#################################
# fastas.R
#################################
# read_fasta
# write_fasta
# load_fasta
# load_fasta2
# find_acc_type
#################################



#################################
# unimods.R
#################################
# parse_unimod
# find_unimod
# table_unimods
#################################



#################################
# dispatch.R
#################################
# find_pos_site
# find_protntsite
# find_anyntsite
# ...
# subpeps_by_vmods
# find_nmodtree
# find_cmodtree
# subset_by_prps
#################################



#################################
# pasef.R
#################################
# hextract_pasef
#   extract_pasef_ms1 (helper)
#     extract_pasef_frame (+ collapse_pasef_xys of all MS1 slices)
#       centroid_pasefms
#   extract_pasef_ms2 (helper)
#     collapse_pasef_ms2scans
#       extract_pasef_frame (not yet same-slice data collapsion)
#         collapse_pasef_xys (of MS2 data within the same slice)
#     group_ms2pasef_by_precursors
#################################


#######################################################################
## MS2 permutations
#
# aas - "H" "Q" "G" "V" "M" "N" "V" "G" "M" "G" "Q" "K" "M" "N" "S"
# 
# Ma - Carbamidomethyl (M)
# Mb - Carbamyl (M)
# N - Deamidated (N)
# 
# Level 1 - sets of MS1 labels (lists of 6)
# 
# Level 2 - permutation of Level-1 with  the positions of residues in aas
#    (lists of 6, the number of total permutations in `n2_perm`; 
# 
#       Carbamidomethyl (M) Carbamyl (M) Deamidated (N) n1_perm n2_perm
# L1.1                   1            1              1      6      36
# L1.2                   2            1              1     12      24
# L1.3                   1            2              1     12      24
# L1.4                   1            1              2     12      36
# L1.5                   2            1              2     30      30
# L1.6                   1            2              2     30      30
#######################################################################



## (Tentative) same-site rules: no additive varmods
# 
# No additive terminal mods (fixed/fixed; var/var; fixed/var)!!!
# 
# (1) No more than one fixedmod on the same residue or N/C terminal.
#   [N] fixed "Oxidation (M)" + fixed "Carbamidomethyl (M)"
#   [N] fixed "TMT6plex (N-term)" + fixed "Acetyl (Protein N-term)"
# 
#   # otherwise additive
#   [Y] fixed "Oxidation (M)" + fixed "TMT6plex (N-term)" with M on the N-term
# (2) OK different variable mods to the same residue at different sites
#   [Y] variable "Oxidation (M) @3" + variable "Carbamidomethyl (M) @4"
#   [Y] variable "Oxidation (M)" + variable "TMT6plex (N-term)" with M on the N-term
# (3) no conflict between fixedmods and variable mods (additive)
#   [Y] fixed "Oxidation (M)" + variable "TMT6plex (N-term)" with M on the N-term
#   [N] fixed "Oxidation (M) @3" + variable "Carbamidomethyl (M) @3"
#   [N] fixed "TMT6plex (N-term)" + variable "Acetyl (Protein N-term)"



