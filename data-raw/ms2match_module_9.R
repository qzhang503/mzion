# ms2match_a1_vnl1_fnl0 (Khare_001)
i=100
x <- frames_adv(mgf_frames[[i]][39:39], theopeps[[i]][54:56], 
                aa_masses = aa_masses, 
                ms1vmods = ms1vmods, 
                ms2vmods = ms2vmods, 
                ntmod = ntmod, 
                ctmod = ctmod, 
                ntmass = ntmass, 
                ctmass = ctmass, 
                amods = amods, 
                vmods_nl = vmods_nl, 
                fmods_nl = NULL, 
                mod_indexes = mod_indexes, 
                type_ms2ions = type_ms2ions, 
                maxn_vmods_per_pep = 
                  maxn_vmods_per_pep, 
                maxn_sites_per_vmod = 
                  maxn_sites_per_vmod, 
                maxn_vmods_sitescombi_per_pep = 
                  maxn_vmods_sitescombi_per_pep, 
                minn_ms2 = minn_ms2, 
                ppm_ms1 = ppm_ms1, 
                ppm_ms2 = ppm_ms2, 
                min_ms2mass = min_ms2mass, 
                digits = digits, 
                FUN = gen_ms2ions_a1_vnl1_fnl0)

## Module 1 (base)
# [1] "PGACITNSAR" "DETNLDLAR"  "IDGSVPSSER" "QLTEEDGVR"  "SLAAEEEAAR"
i=217
x <- mgf_frames[[i]]; y <- theopeps[[i]]
a <- frames_adv(x[11:11], y[14:16], 
                aa_masses = aa_masses, 
                ms1vmods = ms1vmods, 
                ms2vmods = ms2vmods, 
                ntmod = NULL, 
                ctmod = NULL, 
                ntmass = ntmass, 
                ctmass = ctmass, 
                amods = NULL, 
                vmods_nl = NULL, 
                fmods_nl = NULL, 
                mod_indexes = mod_indexes, 
                type_ms2ions = type_ms2ions, 
                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                maxn_sites_per_vmod = maxn_sites_per_vmod, 
                maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
                minn_ms2 = minn_ms2, 
                ppm_ms1 = ppm_ms1, 
                ppm_ms2 = ppm_ms2, 
                min_ms2mass = min_ms2mass, 
                digits = digits, 
                FUN = gen_ms2ions_base)

# Advanced to search_mgf
x <- mapply(
  search_mgf, 
  expt_mass_ms1 = exptmasses_ms1[1:1], 
  expt_moverz_ms2 = exptmoverzs_ms2[1:1], 
  # ex = exptimoverzs_ms2[1:1], 
  MoreArgs = list(
    theomasses_bf_ms1 = theomasses_bf_ms1, 
    theomasses_cr_ms1 = theomasses_cr_ms1, 
    theomasses_af_ms1 = theomasses_af_ms1, 
    theos_bf_ms2 = theos_bf_ms2, 
    theos_cr_ms2 = theos_cr_ms2, 
    theos_af_ms2 = theos_af_ms2, 
    minn_ms2 = minn_ms2, 
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
    min_ms2mass = min_ms2mass
  ), 
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
)

search_mgf(expt_mass_ms1, expt_moverz_ms2, # ex, 
                        theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                        theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                        minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                        min_ms2mass = 115L) 

# Advanced to theos_cr_ms2 
x <- find_ms2_bypep(theos_cr_ms2[[30]], 
                    expts = expt_moverz_ms2, 
                    ex = ex, 
                    d = d2, 
                    ppm_ms2 = ppm_ms2, 
                    min_ms2mass = min_ms2mass, 
                    minn_ms2 = minn_ms2)