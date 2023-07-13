# $bin_masses.R
# [1] "bin_ms1masses"      "binTheoSeqs_i"      "binTheoSeqs2"       "bin_theoseqs"       "binTheoSeqs"       
# [6] "find_ms1_cutpoints" "s_readRDS"          "set_bin_ncores"    
# 
# $dispatch.R
#  [1] "find_pos_site"       "contain_pos_site"    "contain_termpos_any" "subset_by_prps"      "subset_protntsite"  
#  [6] "subset_protntany"    "subset_anyntsite"    "subset_anyntany"     "subset_anysite"      "subset_protctsite"  
# [11] "subset_protctany"    "subset_anyctsite"    "subset_anyctany"     "find_nmodtree"       "find_cmodtree"      
# [16] "subpeps_by_vmods"   
# 
# $fastas.R
# [1] "read_fasta"       "write_fasta"      "load_fasta"       "load_fasta2"      "find_acc_pattern" "find_acc_type"   
# 
# $funs.R
# character(0)
# 
# $ion_ladder.R
#  [1] "ms2ions_by_type" "byions"          "czions"          "axions"          "bions_base"      "yions_base"     
#  [7] "cions_base"      "zions_base"      "c2ions"          "z2ions"          "aions_base"      "xions_base"     
# [13] "a2ions"          "astarions"       "astar2ions"      "a0ions"          "a02ions"         "x2ions"         
# 
# $mapMS2ions.R
#  [1] "mapMS2ions"             "match_mgf_path"         "match_raw_id"           "add_raw_ids"           
#  [5] "find_secion_types"      "find_psm_rows"          "find_psm_rowsQ"         "find_psm_rowsC"        
#  [9] "find_theoexpt_pair"     "find_mgf_query"         "combine_prisec_matches" "check_existed_psms"    
# [13] "get_mzion_coltypes"    
# 
# $mgfs.R
#  [1] "load_mgfs"          "readMGF"            "post_readmgf"       "readlineMGFs"       "  f"               
#  [6] "read_mgf_chunks"    "proc_mgf_chunks"    "proc_mgfs"          "sub_mgftopn"        "integerize_ms2ints"
# [11] "extract_mgf_rptrs"  "find_ms1_interval"  "index_mz"           "find_mgf_type"      "readmzML"          
# [16] "proc_mzml"          "read_mzml"          "prepBrukerMGF"      "mprepBrukerMGF"    
# 
# $ms1_precursors.R
#  [1] "calc_pepmasses2"           "find_aa_masses"            "find_motif_pat"            "simple_prots_peps"        
#  [5] "flat_pepseqs"              "find_aa_site"              "calc_aamasses"             "finalize_aamasses"        
#  [9] "save_mod_indexes"          "check_dupfvmods"           "coerce_fvmods"             "find_f_to_v"              
# [13] "check_mod_motifs"          "find_aamasses_vmodscombi"  "add_var_masses"            "add_fixed_masses"         
# [17] "find_except_sites"         "find_modps"                "extract_umods"             "check_resunimod"          
# [21] "check_fmods_pos_site"      "check_dup_term_any"        "add_aamasses_neulosses"    "add_aamasses_motifs"      
# [25] "parse_aamasses"            "split_fastaseqs"           "make_fastapeps0"           "split_fastaseqs_noenz"    
# [29] "mmake_noenzpeps"           "make_noenzpeps"            "hmake_noenzpeps"           "ms1masses_bare_noenz"     
# [33] "keep_n_misses"             "exclude_n_misses"          "str_exclude_count"         "rm_char_in_nfirst"        
# [37] "rm_char_in_nlast"          "adj_base_masses"           "adj_anywhere_masses"       "add_term_mass"            
# [41] "ms1masses_bare"            "add_ms1_13c"               "add_ms1_notches"           "ms1masses_noterm"         
# [45] "calcms1mass_noterm"        "calcms1mass_noterm_byprot" "calcms1mass_noterm_bypep"  "distri_peps"              
# [49] "ct_counts"                 "distri_fpeps"              "roll_sum"                  "hsemipeps_byprots"        
# [53] "semipeps_byprots"          "calc_semipepmasses"        "delta_ms1_a0_fnl1"         "hms1_a0_vnl0_fnl1"        
# [57] "ms1_a0_vnl0_fnl1"          "hms1_a1_vnl0_fnl0"         "ms1_a1_vnl0_fnl0"         
# 
# $ms2_gen.R
# [1] "gen_ms2ions_base"          "gen_ms2ions_a0_vnl0_fnl1"  "gen_ms2ions_a1_vnl0_fnl0"  "calc_ms2ions_a1_vnl0_fnl0"
# [5] "check_ms1_mass_vmods"      "gen_ms2ions_a1_vnl0_fnl1"  "calc_ms2ions_a1_vnl0_fnl1" "gen_ms2ions_a1_vnl1_fnl0" 
# [9] "calc_ms2ions_a1_vnl1_fnl0"
# 
# $ms2frames.R
#  [1] "pair_mgftheos"  "hpair_mgths"    "hms2match"      "ms2match_all"   "mframes_adv"    "find_ms2_bypep" "search_mgf"    
#  [8] "hms2match_one"  "ms2match_one"   "frames_adv"    
# 
# $msmsmatches.R
#  [1] "matchMS"         "try_psmC2Q"      "reproc_psmC"     "psmC2Q"          "post_psmC2Q"     "check_tmt_pars" 
#  [7] "checkMGF"        "check_locmods"   "map_raw_n_scan"  "check_fdr_group" "check_notches"  
# 
# $msmsmatches2.R
# [1] "ms2match"              "reverse_peps_in_frame" "reverse_seqs"          "calib_mgf"             "calib_ms1"            
# [6] "cv_ms1err"             "post_calib"            "find_ms1_offsets"      "comb_ms1_offsets"     
# 
# $mzion.R
# character(0)
# 
# $mztab.R
# [1] "make_mztab"
# 
# $percolator.R
# [1] "create_folds" "cv_svm"       "perco_svm"   
# 
# $quant2.R
#  [1] "hcalc_tmtint"       "calc_tmtint"        "add_rptrs"          "find_int_cols"      "find_reporter_ints"
#  [6] "find_reporters_ppm" "msub_protpep"       "sub_protpep"        "add_protacc2"       "add_protacc"       
# [11] "hannot_decoys"      "groupProts"         "map_pepprot"        "collapse_sortpeps"  "pcollapse_sortpeps"
# [16] "chunksplit_spmat"   "find_group_breaks"  "cut_proteinGroups"  "sparseD_fourquad"   "as_dist"           
# [21] "greedysetcover3"   
# 
# $roadmaps.R
# character(0)
# 
# $scores.R
#  [1] "add_seions"             "list_leftmatch"         "calc_probi_byvmods"     "calc_probi_bypep"      
#  [5] "calc_probi"             "scalc_pepprobs"         "calc_pepprobs_i"        "calc_pepscores"        
#  [9] "split_im"               "order_fracs"            "order_fracs3"           "combine_fracs"         
# [13] "move_scfiles"           "find_decoy"             "find_targets"           "calcpepsc"             
# [17] "hadd_primatches"        "add_primatches"         "collapse_vecs"          "post_pepscores"        
# [21] "find_pepscore_co1"      "find_pepscore_co2"      "probco_bypeplen"        "sub_td_byfdrtype"      
# [25] "find_optlens"           "find_probco_valley"     "prep_pepfdr_td"         "keep_pepfdr_best"      
# [29] "calc_pepfdr"            "fill_probco_nas"        "fill_probs"             "post_pepfdr"           
# [33] "calc_protfdr"           "aggr_prot_es"           "calc_protfdr_i"         "fit_protfdr"           
# [37] "  f"                    "find_ppm_outer_bycombi" "match_ex2th2"           "calc_peploc"           
# [41] "calcpeprank_1"          "calcpeprank_2"          "calcpeprank_3"          "find_bestnotch"        
# [45] "find_chunkbreaks"       "findLocFracsDF"         "concatFracs"            "na.interp"             
# [49] "is.constant"            "tsoutliers"            
# 
# $silac.R
# [1] "matchMS_silac_mix"   "matchMS_par_groups"  "add_fixedlab_masses" "matchMS_noenzyme"    "combine_ion_matches"
# [6] "comine_PSMsubs"      "matchMS_ms1calib"   
# 
# $unimods.R
#  [1] "parse_unimod"             "find_unimod"              "hfind_unimod"             "table_unimods"           
#  [5] "htable_unimods"           "add_unimod"               "add_modification"         "add_specificy"           
#  [9] "add_delta"                "add_neuloss"              "hadd_neuloss"             "add_comp_elements"       
# [13] "remove_unimod"            "standardize_unimod_ps"    "remove_unimod_title"      "calc_unimod_compmass"    
# [17] "parse_unimod_composition"
# 
# $utils_engine.R
#  [1] "which_topx"            "which_topx2"           "get_topn_vals"         "insVal"                "topx"                 
#  [6] "find_ppm_error"        "find_mass_error_range" "`%+%`"                 "`%+%`"                 "post_frame_adv"       
# [11] "subset_theoframes"     "subset_neuloss_peps"   "find_nterm_mass"       "find_cterm_mass"       "quick_rightjoin"      
# [16] "quick_leftjoin"        "detect_cores"          "find_free_mem"         "find_mod_indexes"      "is_equal_sets"        
# [21] "expand_grid_rows"      "expand_grid"           "expand_gr"             "expand_grid_rows0"     "count_elements"       
# [26] "vec_to_list"           "split_matrix"          "split_vec"             "fold_vec"              "rep_vec"              
# [31] "accumulate_char"       "combi_mat"             "make_zero_df"          "calc_threeframe_ppm"   "get_ms1charges"       
# [36] "finds_uniq_vec"        "my_dataframe"          "flatten_list"          "calc_rev_ms2"          "bind_dfs"             
# 
# $utils_os.R
#  [1] "`names_pos<-`"          "ins_cols_after"         "add_cols_at"            "replace_cols_at"       
#  [5] "reloc_col_after"        "reloc_col_after_last"   "reloc_col_after_first"  "reloc_col_before"      
#  [9] "reloc_col_before_last"  "reloc_col_before_first" "find_preceding_colnm"   "recur_flatten"         
# [13] "chunksplit"             "chunksplitLB"           "find_dir"               "create_dir"            
# [17] "save_call2"             "find_callarg_vals"      "match_calltime"         "delete_files"          
# [21] "find_ms1_times"         "is_nulllist"            "add_nulllist"          
# 
# $utils_ui.R
# [1] "calc_monopeptide"  "calc_monopep"      "check_aaseq"       "calc_ms2ionseries" "calc_ms2ions"      "unique_mvmods"    
# [7] "vmods_elements"    "find_intercombi"  
# 
# $vmod_ms1_labels.R
#  [1] "match_mvmods"      "make_ms1vmod_i"    "make_ms1_vmodsets" "bacth_vmods_combi" "make_unique_sets"  "find_unique_sets" 
#  [7] "gtools_combn"      "    sub"           "  else sub"        "find_intercombi2" 
# 
# $vmod_ms2_labels.R
# [1] "find_vmodposU"   "find_vmodposM"   "make_ms2vmods"   "find_ms2resids"  "find_perm_sets"  "add_one_permlab"
# [7] "add_one_label"   "ins_permlab"     "sim_combn"      
# 
# $wrappers.R
# [1] "my_dist"       "cos_sim"       "matchMS_NES"   "rematchMS_NES"
# 
# $zzz.R
# [1] ".onAttach"
# 
