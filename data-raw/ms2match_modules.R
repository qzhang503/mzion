x = mgf_frames[[1]]
y = theopeps[[1]]
from = 133
to = 133
sta = as.integer(names(x[from])) - 1L
end = as.integer(names(x[to])) + 1L

df <- frames_adv(
  x[from:to], y[which(names(y) == sta):which(names(y) == end)], 
  aa_masses = aa_masses, 
  ms1vmods = ms1vmods, 
  ms2vmods = ms2vmods, 
  ntmod = ntmod, 
  ctmod = ctmod, 
  ntmass = ntmass, 
  ctmass = ctmass, 
  amods = amods, 
  vmods_nl = vmods_nl, 
  fmods_nl = fmods_nl, 
  pep_mod_group = pep_mod_group, 
  mod_indexes = mod_indexes, 
  type_ms2ions = type_ms2ions, 
  maxn_vmods_per_pep = maxn_vmods_per_pep, 
  maxn_sites_per_vmod = maxn_sites_per_vmod, 
  maxn_fnl_per_seq = maxn_fnl_per_seq, 
  maxn_vnl_per_seq = maxn_vnl_per_seq, 
  maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
  minn_ms2 = minn_ms2, 
  ppm_ms1 = ppm_ms1, 
  ppm_ms2 = ppm_ms2, 
  min_ms2mass = min_ms2mass, 
  index_mgf_ms2 = index_mgf_ms2, 
  ms1_offsets = ms1_offsets, 
  FUN = FUN
)

## Module 9
# ms2match_a1_vnl1_fnl0 (Kh_001)
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

## Another module 1: ill-form match
i=217
x <- mgf_frames[[i]]; y <- theopeps[[i]]
z = x[[31]]
z = z[5:7, ]
x[[31]] <- z
a <- frames_adv(x[31:31], y[37:39], 
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

which(theopeps_cr_ms1 == "NCIDITGVR") # 23

mapply(
  search_mgf, 
  expt_mass_ms1 = exptmasses_ms1[1], 
  expt_moverz_ms2 = exptmoverzs_ms2[1], 
  ###
  # ex = exptimoverzs_ms2, 
  ###
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

bf_allowed <- which(theomasses_bf_ms1 > max(theomasses_bf_ms1)) # should be empty
af_allowed <- which(theomasses_af_ms1 < min(af_allowed))

ans_cr <- lapply(theos_cr_ms2[23:23], find_ms2_bypep, 
                 expts = expt_moverz_ms2, 
                 ex = ex, 
                 d = d2, 
                 ppm_ms2 = ppm_ms2, 
                 min_ms2mass = min_ms2mass, 
                 minn_ms2 = minn_ms2)

## Module 5
# frame 104634
# "VFPPDEMEQVSNK"    "TFVVQGFGNVGLHSMR"
i=125
theopeps[[217]] <- list(NULL)
out <- frames_adv(mgf_frames[[i]][167], theopeps[[i]][217:219], 
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

row <- which(theos_cr_ms1$pep_seq == "VFPPDEMEQVSNK") # 18
# List of three: (1) NULL; (2) and (3) results
x <- mapply(
  search_mgf, 
  expt_mass_ms1 = exptmasses_ms1, 
  expt_moverz_ms2 = exptmoverzs_ms2, 
  MoreArgs = list(
    theomasses_bf_ms1 = NULL, 
    theomasses_cr_ms1 = theomasses_cr_ms1[row], 
    theomasses_af_ms1 = NULL, 
    theos_bf_ms2 = list(NULL), 
    theos_cr_ms2 = theos_cr_ms2[row], 
    theos_af_ms2 = list(NULL), 
    minn_ms2 = minn_ms2, 
    ppm_ms1 = ppm_ms1, 
    ppm_ms2 = ppm_ms2, 
    min_ms2mass = min_ms2mass
  ), 
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
)


