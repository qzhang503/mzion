# $batchMS2.R
# [1] "batch_ms2ions"         "hbatch_ms2ions"        "mgen_ms2ions"          "make_ms2frames"        "make_ms2frames_bypars"
# [6] "check_ms2frames"      
# 
# $bin_masses.R
# [1] "bin_ms1masses"      "binTheoSeqs_i"      "binTheoSeqs2"       "bin_theoseqs"       "binTheoSeqs"       
# [6] "find_ms1_cutpoints" "s_readRDS"         
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
#  [7] "b2ions_base"     "bstarions"       "bstar2ions"      "b0ions"          "b02ions"         "y2ions"         
# [13] "ystarions"       "ystar2ions"      "y0ions"          "y02ions"         "cions_base"      "c2ions"         
# [19] "zions_base"      "z2ions"          "aions_base"      "a2ions"          "astarions"       "astar2ions"     
# [25] "a0ions"          "a02ions"         "xions_base"      "x2ions"         
# 
# $mapMS2ions.R
#  [1] "mapMS2ions"             "match_mgf_path"         "match_raw_id"           "add_raw_ids"           
#  [5] "find_secion_types"      "find_psm_rows"          "find_psm_rows1"         "find_psm_rows2"        
#  [9] "find_theoexpt_pair"     "find_mgf_query"         "combine_prisec_matches" "check_existed_psms"    
# [13] "get_proteoM_coltypes"  
# 
# $mgfs.R
#  [1] "load_mgfs"          "readMGF"            "post_readmgf"       "readlineMGFs"       "  f"               
#  [6] "read_mgf_chunks"    "proc_mgf_chunks"    "proc_mgfs"          "sub_mgftopn"        "integerize_ms2ints"
# [11] "extract_mgf_rptrs"  "find_ms1_interval"  "index_mz"           "find_mgf_type"      "readmzML"          
# [16] "proc_mzml"          "read_mzml"         
# 
# $ms1_precursors.R
#  [1] "calc_pepmasses2"           "find_aa_masses"            "find_motif_pat"            "tbl_prots_peps"           
#  [5] "flat_pepseqs"              "find_aa_site"              "calc_aamasses"             "finalize_aamasses"        
#  [9] "save_mod_indexes"          "check_dupfvmods"           "coerce_fvmods"             "find_f_to_v"              
# [13] "check_mod_motifs"          "find_aamasses_vmodscombi"  "add_var_masses"            "add_fixed_masses"         
# [17] "find_except_sites"         "find_modps"                "extract_umods"             "check_resunimod"          
# [21] "check_fmods_pos_site"      "add_aamasses_neulosses"    "add_aamasses_motifs"       "parse_aamasses"           
# [25] "split_fastaseqs"           "make_fastapeps0"           "split_fastaseqs_noenz"     "mmake_noenzpeps"          
# [29] "make_noenzpeps"            "hmake_noenzpeps"           "ms1masses_bare_noenz"      "keep_n_misses"            
# [33] "exclude_n_misses"          "str_exclude_count"         "rm_char_in_nfirst"         "rm_char_in_nlast"         
# [37] "adj_base_masses"           "adj_anywhere_masses"       "add_term_mass"             "ms1masses_bare"           
# [41] "add_ms1_13c"               "add_ms1_13c_orig"          "ms1masses_noterm"          "calcms1mass_noterm"       
# [45] "calcms1mass_noterm_byprot" "calcms1mass_noterm_bypep"  "distri_peps"               "ct_counts"                
# [49] "distri_fpeps"              "roll_sum"                  "hsemipeps_byprots"         "semipeps_byprots"         
# [53] "calc_semipepmasses"        "delta_ms1_a0_fnl1"         "hms1_a0_vnl0_fnl1"         "ms1_a0_vnl0_fnl1"         
# [57] "hms1_a1_vnl0_fnl0"         "ms1_a1_vnl0_fnl0"         
# 
# $ms2_a0_vnl0_fnl1.R
# [1] "ms2match_a0_vnl0_fnl1"    "gen_ms2ions_a0_vnl0_fnl1"
# 
# $ms2_a1_vnl0_fnl0.R
# [1] "ms2match_a1_vnl0_fnl0"     "gen_ms2ions_a1_vnl0_fnl0"  "calc_ms2ions_a1_vnl0_fnl0" "check_ms1_mass_vmods2"    
# [5] "add_hexcodes"             
# 
# $ms2_a1_vnl0_fnl1.R
# [1] "ms2match_a1_vnl0_fnl1"     "gen_ms2ions_a1_vnl0_fnl1"  "calc_ms2ions_a1_vnl0_fnl1" "add_hexcodes_fnl2"        
# 
# $ms2_a1_vnl1_fnl0.R
# [1] "ms2match_a1_vnl1_fnl0"     "gen_ms2ions_a1_vnl1_fnl0"  "calc_ms2ions_a1_vnl1_fnl0" "add_hexcodes_vnl2"        
# 
# $ms2_base.R
# [1] "ms2match_base"    "frames_adv"       "gen_ms2ions_base" "fuzzy_match_one"  "fuzzy_match_one2" "find_ms2_bypep"  
# [7] "search_mgf"      
# 
# $msmsmatches.R
# [1] "matchMS"        "try_psmC2Q"     "reproc_psmC"    "psmC2Q"         "post_psmC2Q"    "check_tmt_pars" "checkMGF"      
# [8] "check_locmods"  "map_raw_n_scan"
# 
# $msmsmatches2.R
# [1] "ms2match"              "hcalc_tmtint"          "reverse_peps_in_frame" "reverse_seqs"          "calib_ms1masses"      
# 
# $mztab.R
# [1] "make_mztab"
# 
# $proteoM.R
# character(0)
# 
# $quant2.R
#  [1] "calc_tmtint"        "add_rptrs"          "find_reporter_ints" "find_reporters_ppm" "add_prot_acc"      
#  [6] "hfwd_prps"          "hrev_prps"          "add_prot_acc2"      "hadd_prot_acc"      "groupProts"        
# [11] "map_pepprot"        "collapse_sortpeps"  "pcollapse_sortpeps" "chunksplit_spmat"   "find_group_breaks" 
# [16] "cut_proteinGroups"  "sparseD_fourquad"   "as_dist"            "as_lgldist"         "greedysetcover3"   
# 
# $roadmaps.R
# character(0)
# 
# $scores.R
#  [1] "add_seions"             "list_leftmatch"         "calc_probi_byvmods"     "calc_probi_bypep"      
#  [5] "calc_probi"             "scalc_pepprobs"         "calc_pepprobs_i"        "calc_pepscores"        
#  [9] "find_decoy"             "find_targets"           "calcpepsc"              "add_primatches"        
# [13] "collapse_vecs"          "post_pepscores"         "find_pepscore_co1"      "find_pepscore_co2"     
# [17] "probco_bypeplen"        "find_optlens"           "find_probco_valley"     "calc_pepfdr"           
# [21] "fill_probco_nas"        "fill_probs"             "post_pepfdr"            "calc_protfdr"          
# [25] "aggr_prot_es"           "calc_protfdr_i"         "fit_protfdr"            "  f"                   
# [29] "find_ppm_outer_bycombi" "match_ex2th2"           "calc_peploc"            "calcpeprank_1"         
# [33] "calcpeprank_2"          "calcpeprank_3"          "find_chunkbreaks"       "findLocFracsDF"        
# [37] "concatFracs"            "na.interp"              "is.constant"            "tsoutliers"            
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
#  [6] "find_ppm_error"        "find_mass_error_range" "`%+%`"                 "post_ms2match"         "post_frame_adv"       
# [11] "purge_search_space"    "subset_theoframes"     "subset_neuloss_peps"   "find_nterm_mass"       "find_cterm_mass"      
# [16] "quick_rightjoin"       "quick_leftjoin"        "detect_cores"          "find_free_mem"         "find_mod_indexes"     
# [21] "is_equal_sets"         "purge_decoys"          "expand_grid_rows"      "count_elements"        "vec_to_list"          
# [26] "split_vec"             "accumulate_char"       "combi_mat"             "make_zero_df"          "hash_frame_nums"      
# [31] "find_from_hash"        "calc_threeframe_ppm"   "check_ms1calib"        "save_ms1calib"         "get_ms1charges"       
# 
# $utils_os.R
#  [1] "`names_pos<-`"          "find_int_cols"          "ins_cols_after"         "add_cols_at"           
#  [5] "replace_cols_at"        "reloc_col_after"        "reloc_col_after_last"   "reloc_col_after_first" 
#  [9] "reloc_col_before"       "reloc_col_before_last"  "reloc_col_before_first" "find_preceding_colnm"  
# [13] "recur_flatten"          "chunksplit"             "chunksplitLB"           "find_dir"              
# [17] "create_dir"             "save_call2"             "find_callarg_vals"      "match_calltime"        
# [21] "delete_files"           "find_ms1_times"         "get_globalvar"          "load_cache_info"       
# [25] "is_nulllist"            "add_nulllist"          
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
#  [1] "find_vmodscombi"   "combi_namesiteU"   "find_vmodposU"     "combi_namesiteM"   "find_vmodposM"     "match_aas_indexes"
#  [7] "make_ms2vmods"     "find_ms2resids"    "find_perm_sets"    "add_one_permlab"   "add_one_label"     "ins_permlab"      
# [13] "sim_combn"        
# 
# $wrappers.R
# [1] "my_dist"     "cos_sim"     "matchMS_NES"
# 
# $zzz.R
# [1] ".onAttach"
# 
