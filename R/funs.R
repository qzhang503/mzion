# $bin_masses.R
# [1] "bin_ms1masses"      "binTheoSeqs_i"      "binTheoSeqs2"       "bin_theoseqs"       "binTheoSeqs"        "find_ms1_cutpoints"
# [7] "s_readRDS"          "set_bin_ncores"    
# 
# $deisotope.R
# [1] "find_ms1stat"      "find_charge_state" "check_chduo"       "is_true"           "find_dbl_z"        "find_lcpeaks"     
# [7] "sep_ms1rts"        "bin_ms1rts"        "gen_isoenvlope"   
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
#  [1] "ms2ions_by_type" "byions"          "czions"          "axions"          "bions_base"      "yions_base"      "cions_base"     
#  [8] "zions_base"      "c2ions"          "z2ions"          "aions_base"      "xions_base"      "a2ions"          "astarions"      
# [15] "astar2ions"      "a0ions"          "a02ions"         "x2ions"         
# 
# $lfq.R
# [1] "subMSfull"     "pretraceXY"    "htraceXY"      "traceXY"       "updateMS1Int"  "updateMS1Int2" "traceMS1"      "getMS1Int"    
# 
# $mapMS2ions.R
#  [1] "mapMS2ions"         "plotMS2ions"        "match_mgf_path"     "match_raw_id"       "add_raw_ids"        "find_secion_types" 
#  [7] "find_mgf_query"     "make_speclib"       "get_mzion_coltypes" "check_ggname"      
# 
# $mgfs.R
#  [1] "load_mgfs"          "readMGF"            "post_readmgf"       "readlineMGFs"       "  f"                "read_mgf_chunks"   
#  [7] "proc_mgf_chunks"    "proc_mgfs"          "reset_rettimes"     "sub_mgftopn"        "integerize_ms2ints" "extract_mgf_rptrs" 
# [13] "index_mz"           "find_mgf_type"      "prepBrukerMGF"      "mprepBrukerMGF"    
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
#  [1] "pair_mgftheos"   "thin_mgf"        "hpair_mgths"     "make_dia_mgfs"   "clean_flat_mgfs" "hms2match"       "find_ms2_bypep" 
#  [8] "search_mgf"      "hms2match_one"   "ms2match_one"    "frames_adv"     
# 
# $msfilereader.R
# [1] "readRAW"                   "proc_raws"                 "exeReadRAW"                "acceptMSFileReaderLicense"
# 
# $msmsmatches.R
#  [1] "matchMS"            "try_psmC2Q"         "reproc_psmC"        "psmC2Q"             "post_psmC2Q"        "check_tmt_pars"    
#  [7] "checkMGF"           "check_locmods"      "map_raw_n_scan"     "map_raw_n_scan_old" "check_fdr_group"    "check_notches"     
# 
# $msmsmatches2.R
#  [1] "ms2match"              "reverse_peps_in_frame" "reverse_seqs"          "calib_mgf"             "calib_ms1"            
#  [6] "substract_ms1mass"     "cv_ms1err"             "post_calib"            "find_ms1_offsets"      "comb_ms1_offsets"     
# 
# $mzion.R
# character(0)
# 
# $mzml.R
#  [1] "readmzML"          "hloadMZML"         "loadMZML"          "extrDDA"           "hdeisoDDA"         "deisoDDA"         
#  [7] "get_ms1xs_space"   "estimate_rtgap"    "predeisoDDA"       "deisoDDAMS2"       "pasefMS1xyz"       "getMSrowIndexes"  
# [13] "find_ms2ends"      "getMS1xyz"         "getMS2xyz"         "extrDIA"           "hdeisoDIA"         "deisoDIA"         
# [19] "hsubDIAMS1"        "subDIAMS1"         "htraceDIA"         "traceDIA"          "flattenMSxyz"      "spreadMSohw"      
# [25] "spreadMS_v1"       "comb_mstraces"     "find_gates"        "find_gate_edges"   "traceLCMS"         "collapse_xyz"     
# [31] "mapcoll_xyz"       "find_lc_gates"     "fill_lc_gaps"      "collapse_mms1ints" "calc_ms1xys"       "find_mdda_mms1s"  
# [37] "find_ms1byms2"    
# 
# $mztab.R
# [1] "make_mztab"
# 
# $pasefreader.R
#  [1] "readPASEF"                    "proc_pasefs"                  "add_pasef_precursors"         "sep_pasef_ms2info"           
#  [5] "hextract_pasef"               "extract_pasef_ms1"            "extract_pasef_ms2"            "centroid_pasefms"            
#  [9] "sum_pasef_ms1"                "collapse_rawtims_xys"         "group_ms2pasef_by_precursors" "add_pasef_ms2iso"            
# [13] "extract_pasef_frame"          "collapse_pasef_xys"           "exeReadPASEF"                 "acceptBrukerLicense"         
# 
# $percolator.R
# [1] "create_folds" "cv_svm"       "perco_svm"   
# 
# $quant2.R
#  [1] "hcalc_tmtint"       "calc_tmtint"        "add_rptrs"          "find_int_cols"      "find_reporter_ints" "find_reporters_ppm"
#  [7] "msub_protpep"       "sub_protpep"        "add_protacc2"       "add_protacc"        "hannot_decoys"      "groupProts"        
# [13] "map_pepprot"        "collapse_sortpeps"  "pcollapse_sortpeps" "find_group_breaks"  "cut_proteinGroups"  "sparseD_fourquad"  
# [19] "as_dist"            "greedysetcover3"   
# 
# $roadmaps.R
# character(0)
# 
# $scores.R
#  [1] "add_seions"             "list_leftmatch"         "calc_probi_byvmods"     "calc_probi_bypep"       "calc_probi"            
#  [6] "scalc_pepprobs"         "calc_pepprobs_i"        "calc_pepscores"         "split_im"               "order_fracs"           
# [11] "order_fracs3"           "combine_fracs"          "move_scfiles"           "find_decoy"             "find_targets"          
# [16] "calcpepsc"              "find_iexunv"            "addChim"                "hadd_primatches"        "add_primatches"        
# [21] "collapse_vecs"          "post_pepscores"         "find_pepscore_co1"      "find_pepscore_co2"      "probco_bypeplen"       
# [26] "sub_td_byfdrtype"       "find_optlens"           "find_probco_valley"     "prep_pepfdr_td"         "keep_pepfdr_best"      
# [31] "calc_pepfdr"            "fill_probco_nas"        "find_fdr_fits"          "fill_probs"             "post_pepfdr"           
# [36] "calc_protfdr"           "aggr_prot_es"           "calc_protfdr_i"         "fit_protfdr"            "  f"                   
# [41] "find_ppm_outer_bycombi" "match_ex2th2"           "calc_peploc"            "calcpeprank_1"          "calcpeprank_2"         
# [46] "calcpeprank_3"          "find_bestnotch"         "find_chunkbreaks"       "findLocFracsDF"         "concatFracs"           
# [51] "na.interp"              "is.constant"            "tsoutliers"             "rm_dup13c"             
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
# [11] "subset_theoframes"     "subset_neuloss_peps"   "find_nterm_mass"       "find_cterm_mass"       "quick_join"           
# [16] "detect_cores"          "find_free_mem"         "find_mod_indexes"      "is_equal_sets"         "expand_grid_rows"     
# [21] "expand_grid"           "expand_gr"             "expand_grid_rows0"     "count_elements"        "vec_to_list"          
# [26] "split_matrix"          "split_vec"             "fold_vec"              "fold_vec2"             "sep_vec"              
# [31] "rep_vec"               "accumulate_char"       "combi_mat"             "make_zero_df"          "calc_threeframe_ppm"  
# [36] "get_ms1charges"        "finds_uniq_vec"        "my_dataframe"          "flatten_list"          "calc_rev_ms2"         
# [41] "bind_dfs"              "find_min_ncores"      
# 
# $utils_os.R
#  [1] "`names_pos<-`"          "ins_cols_after"         "add_cols_at"            "replace_cols_at"        "reloc_col_after"       
#  [6] "reloc_col_after_last"   "reloc_col_after_first"  "reloc_col_before"       "reloc_col_before_last"  "reloc_col_before_first"
# [11] "find_preceding_colnm"   "recur_flatten"          "chunksplit"             "chunksplitLB"           "find_dir"              
# [16] "create_dir"             "save_call2"             "find_callarg_vals"      "match_calltime"         "delete_files"          
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
# [1] "find_vmodposU"   "find_vmodposM"   "make_ms2vmods"   "find_ms2resids"  "find_perm_sets"  "add_one_permlab" "add_one_label"  
# [8] "ins_permlab"     "sim_combn"      
# 
# $wrappers.R
# [1] "my_dist"       "cos_sim"       "matchMS_NES"   "rematchMS_NES"
# 
# $zzz.R
# [1] ".onAttach"
# 
