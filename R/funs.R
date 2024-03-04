# $bin_masses.R
# [1] "bin_ms1masses"      "binTheoSeqs_i"      "binTheoSeqs2"       "bin_theoseqs"       "binTheoSeqs"       
# [6] "find_ms1_cutpoints" "s_readRDS"          "set_bin_ncores"    
# 
# $deisotope.R
# [1] "find_ms1stat"      "find_charge_state" "check_chduo"       "is_true"           "find_dbl_z"       
# [6] "find_lcpeaks"      "sep_ms1rts"        "bin_ms1rts"        "gen_isoenvlope"   
# 
# $dispatch.R
#  [1] "find_pos_site"       "contain_pos_site"    "contain_termpos_any" "subset_by_prps"     
#  [5] "subset_protntsite"   "subset_protntany"    "subset_anyntsite"    "subset_anyntany"    
#  [9] "subset_anysite"      "subset_protctsite"   "subset_protctany"    "subset_anyctsite"   
# [13] "subset_anyctany"     "find_nmodtree"       "find_cmodtree"       "subpeps_by_vmods"   
# 
# $fastas.R
# [1] "read_fasta"       "write_fasta"      "load_fasta"       "load_fasta2"      "find_acc_pattern"
# [6] "find_acc_type"   
# 
# $funs.R
# character(0)
# 
# $ion_ladder.R
#  [1] "ms2ions_by_type" "byions"          "czions"          "axions"          "bions_base"      "yions_base"     
#  [7] "cions_base"      "zions_base"      "c2ions"          "z2ions"          "aions_base"      "xions_base"     
# [13] "a2ions"          "astarions"       "astar2ions"      "a0ions"          "a02ions"         "x2ions"         
# 
# $lfq.R
# [1] "subMSfull"    "prep_traceXY" "htraceXY"     "traceXY"      "updateMS1Int" "traceMS1"     "getMS1Int"   
# 
# $mapMS2ions.R
#  [1] "mapMS2ions"         "plotMS2ions"        "match_mgf_path"     "match_raw_id"       "add_raw_ids"       
#  [6] "find_secion_types"  "find_mgf_query"     "make_speclib"       "get_mzion_coltypes" "check_ggname"      
# 
# $mgfs.R
#  [1] "load_mgfs"          "readMGF"            "post_readmgf"       "readlineMGFs"       "  f"               
#  [6] "read_mgf_chunks"    "proc_mgf_chunks"    "proc_mgfs"          "reset_rettimes"     "sub_mgftopn"       
# [11] "integerize_ms2ints" "extract_mgf_rptrs"  "index_mz"           "find_mgf_type"      "prepBrukerMGF"     
# [16] "mprepBrukerMGF"    
# 
# $ms1_precursors.R
#  [1] "calc_pepmasses2"           "find_aa_masses"            "find_motif_pat"           
#  [4] "simple_prots_peps"         "flat_pepseqs"              "find_aa_site"             
#  [7] "calc_aamasses"             "finalize_aamasses"         "save_mod_indexes"         
# [10] "check_dupfvmods"           "coerce_fvmods"             "find_f_to_v"              
# [13] "check_mod_motifs"          "find_aamasses_vmodscombi"  "add_var_masses"           
# [16] "add_fixed_masses"          "find_except_sites"         "find_modps"               
# [19] "extract_umods"             "check_resunimod"           "check_fmods_pos_site"     
# [22] "check_dup_term_any"        "add_aamasses_neulosses"    "add_aamasses_motifs"      
# [25] "parse_aamasses"            "split_fastaseqs"           "make_fastapeps0"          
# [28] "split_fastaseqs_noenz"     "mmake_noenzpeps"           "make_noenzpeps"           
# [31] "hmake_noenzpeps"           "ms1masses_bare_noenz"      "keep_n_misses"            
# [34] "exclude_n_misses"          "str_exclude_count"         "rm_char_in_nfirst"        
# [37] "rm_char_in_nlast"          "adj_base_masses"           "adj_anywhere_masses"      
# [40] "add_term_mass"             "ms1masses_bare"            "add_ms1_13c"              
# [43] "add_ms1_notches"           "ms1masses_noterm"          "calcms1mass_noterm"       
# [46] "calcms1mass_noterm_byprot" "calcms1mass_noterm_bypep"  "distri_peps"              
# [49] "ct_counts"                 "distri_fpeps"              "roll_sum"                 
# [52] "hsemipeps_byprots"         "semipeps_byprots"          "calc_semipepmasses"       
# [55] "delta_ms1_a0_fnl1"         "hms1_a0_vnl0_fnl1"         "ms1_a0_vnl0_fnl1"         
# [58] "hms1_a1_vnl0_fnl0"         "ms1_a1_vnl0_fnl0"         
# 
# $ms2_gen.R
# [1] "gen_ms2ions_base"          "gen_ms2ions_a0_vnl0_fnl1"  "gen_ms2ions_a1_vnl0_fnl0" 
# [4] "calc_ms2ions_a1_vnl0_fnl0" "check_ms1_mass_vmods"      "gen_ms2ions_a1_vnl0_fnl1" 
# [7] "calc_ms2ions_a1_vnl0_fnl1" "gen_ms2ions_a1_vnl1_fnl0"  "calc_ms2ions_a1_vnl1_fnl0"
# 
# $ms2frames.R
#  [1] "pair_mgftheos"  "hpair_mgths"    "make_dia_mgfs"  "hms2match"      "ms2match_all"   "mframes_adv"   
#  [7] "find_ms2_bypep" "search_mgf"     "hms2match_one"  "ms2match_one"   "frames_adv"    
# 
# $msmsmatches.R
#  [1] "matchMS"            "try_psmC2Q"         "reproc_psmC"        "psmC2Q"             "post_psmC2Q"       
#  [6] "check_tmt_pars"     "checkMGF"           "check_locmods"      "map_raw_n_scan"     "map_raw_n_scan_old"
# [11] "check_fdr_group"    "check_notches"     
# 
# $msmsmatches2.R
# [1] "ms2match"              "reverse_peps_in_frame" "reverse_seqs"          "calib_mgf"            
# [5] "calib_ms1"             "cv_ms1err"             "post_calib"            "find_ms1_offsets"     
# [9] "comb_ms1_offsets"     
# 
# $mzion.R
# character(0)
# 
# $mzml.R
#  [1] "readmzML"          "hloadMZML"         "loadMZML"          "extrDDA"           "hdeisoDDA"        
#  [6] "deisoDDA"          "getMSrowIndexes"   "getMS1xyz"         "getMS2xyz"         "extrDIA"          
# [11] "hdeisoDIA"         "deisoDIA"          "hsubDIAMS1"        "subDIAMS1"         "htraceDIA"        
# [16] "traceDIA"          "flattenMSxyz"      "spreadMSohw"       "spreadMS_v1"       "comb_mstraces"    
# [21] "find_gates"        "find_gate_edges"   "traceLCMS"         "collapse_xyz"      "mapcoll_xyz"      
# [26] "find_lc_gates"     "fill_lc_gaps"      "collapse_mms1ints" "calc_ms1xys"       "find_mdda_mms1s"  
# [31] "find_ms1byms2"    
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
# [16] "find_group_breaks"  "cut_proteinGroups"  "sparseD_fourquad"   "as_dist"            "greedysetcover3"   
# 
# $roadmaps.R
# character(0)
# 
# $scores.R
#  [1] "add_seions"             "list_leftmatch"         "calc_probi_byvmods"     "calc_probi_bypep"      
#  [5] "calc_probi"             "scalc_pepprobs"         "calc_pepprobs_i"        "calc_pepscores"        
#  [9] "split_im"               "order_fracs"            "order_fracs3"           "combine_fracs"         
# [13] "move_scfiles"           "find_decoy"             "find_targets"           "calcpepsc"             
# [17] "find_iexunv"            "addChim"                "hadd_primatches"        "add_primatches"        
# [21] "collapse_vecs"          "post_pepscores"         "find_pepscore_co1"      "find_pepscore_co2"     
# [25] "probco_bypeplen"        "sub_td_byfdrtype"       "find_optlens"           "find_probco_valley"    
# [29] "prep_pepfdr_td"         "keep_pepfdr_best"       "calc_pepfdr"            "fill_probco_nas"       
# [33] "find_fdr_fits"          "fill_probs"             "post_pepfdr"            "calc_protfdr"          
# [37] "aggr_prot_es"           "calc_protfdr_i"         "fit_protfdr"            "  f"                   
# [41] "find_ppm_outer_bycombi" "match_ex2th2"           "calc_peploc"            "calcpeprank_1"         
# [45] "calcpeprank_2"          "calcpeprank_3"          "find_bestnotch"         "find_chunkbreaks"      
# [49] "findLocFracsDF"         "concatFracs"            "na.interp"              "is.constant"           
# [53] "tsoutliers"             "rm_dup13c"             
# 
# $silac.R
# [1] "matchMS_silac_mix"   "matchMS_par_groups"  "add_fixedlab_masses" "matchMS_noenzyme"   
# [5] "combine_ion_matches" "comine_PSMsubs"      "matchMS_ms1calib"   
# 
# $unimods.R
#  [1] "parse_unimod"             "find_unimod"              "hfind_unimod"             "table_unimods"           
#  [5] "htable_unimods"           "add_unimod"               "add_modification"         "add_specificy"           
#  [9] "add_delta"                "add_neuloss"              "hadd_neuloss"             "add_comp_elements"       
# [13] "remove_unimod"            "standardize_unimod_ps"    "remove_unimod_title"      "calc_unimod_compmass"    
# [17] "parse_unimod_composition"
# 
# $utils_engine.R
#  [1] "which_topx"            "which_topx2"           "get_topn_vals"         "insVal"               
#  [5] "topx"                  "find_ppm_error"        "find_mass_error_range" "`%+%`"                
#  [9] "`%+%`"                 "post_frame_adv"        "subset_theoframes"     "subset_neuloss_peps"  
# [13] "find_nterm_mass"       "find_cterm_mass"       "quick_rightjoin"       "quick_leftjoin"       
# [17] "detect_cores"          "find_free_mem"         "find_mod_indexes"      "is_equal_sets"        
# [21] "expand_grid_rows"      "expand_grid"           "expand_gr"             "expand_grid_rows0"    
# [25] "count_elements"        "vec_to_list"           "split_matrix"          "split_vec"            
# [29] "fold_vec"              "fold_vec2"             "sep_vec"               "rep_vec"              
# [33] "accumulate_char"       "combi_mat"             "make_zero_df"          "calc_threeframe_ppm"  
# [37] "get_ms1charges"        "finds_uniq_vec"        "my_dataframe"          "flatten_list"         
# [41] "calc_rev_ms2"          "bind_dfs"              "find_min_ncores"      
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
# [1] "calc_monopeptide"  "calc_monopep"      "check_aaseq"       "calc_ms2ionseries" "calc_ms2ions"     
# [6] "unique_mvmods"     "vmods_elements"    "find_intercombi"  
# 
# $vmod_ms1_labels.R
#  [1] "match_mvmods"      "make_ms1vmod_i"    "make_ms1_vmodsets" "bacth_vmods_combi" "make_unique_sets" 
#  [6] "find_unique_sets"  "gtools_combn"      "    sub"           "  else sub"        "find_intercombi2" 
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
