#' Calculates the mono-isotopic mass of a peptide sequence.
#'
#' Only for direct uses from an R console (with trade-offs in speed).
#'
#' @inheritParams calc_monopep
#' @inheritParams calc_aamasses
#'
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = NULL,
#'                       varmods = NULL)
#'
#' stopifnot((unlist(x$mass) - 1594.5369) < 1e-4)
#'
#' # (2-a)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = "Oxidation (M)",
#'                       varmods = NULL)
#'
#' m <- unlist(x$mass)
#'
#' stopifnot((m[1] - 1626.5267) < 1e-4)
#'
#' # (2-b) combinatorial NL for fixed modifications
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = "Oxidation (M)",
#'                       varmods = NULL)
#'
#' m <- unlist(x$mass)
#'
#' stopifnot(length(m) == 2L)
#' stopifnot((m[1] - m[2] - 127.9965) < 1e-4)
#'
#' ## With variable modifications
#' # (3-a)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = NULL,
#'                       varmods = "Oxidation (M)")
#'
#' m <- unlist(x$mass)
#'
#' stopifnot((m[1] - 1594.5369) < 1e-4)
#' stopifnot((m[2] - m[1]) == (m[3] - m[2]))
#'
#' # x$vmods_ps
#'
#' # (3-b)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = NULL,
#'                       varmods = "Oxidation (M)")
#'
#' x$mass
#'
#' # (4-a)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       c("TMT6plex (N-term)",
#'                         "TMT6plex (K)",
#'                         "Carbamidomethyl (C)"),
#'                       c("Acetyl (N-term)",
#'                         "Gln->pyro-Glu (N-term = Q)",
#'                         "Oxidation (M)"))
#'
#' x$mass
#'
#' # The N-term M realizes with acetylation
#' x$vmods_ps[[1]]
#'
#' # The N-term M realizes with TMT
#' x$vmods_ps[[2]]
#'
#'
#' # (4-b)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       fixedmods = c("TMT6plex (N-term)",
#'                                     "TMT6plex (K)",
#'                                     "Carbamidomethyl (C)"),
#'                       varmods = c("Acetyl (N-term)",
#'                                   "Gln->pyro-Glu (N-term = Q)",
#'                                   "Oxidation (M)"))
#'
#' x$mass
#' }
#' @export
calc_monopeptide <- function (aa_seq, fixedmods, varmods,
                              maxn_vmods_setscombi = 64L,
                              maxn_vmods_per_pep = Inf,
                              maxn_sites_per_vmod = Inf,
                              min_mass = 200L, 
                              max_mass = 4500L, 
                              digits = 4L) 
{
  options(digits = 9L)
  
  aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                                 varmods = varmods,
                                 maxn_vmods_setscombi = maxn_vmods_setscombi)
  
  peps <- check_aaseq(aa_seq, aa_masses_all, fixedmods, varmods)
  oks <- purrr::map_lgl(peps, function (x) !purrr::is_empty(x))
  peps <- peps[oks]
  aa_masses_all <- aa_masses_all[oks]
  
  ms <- purrr::map2(peps, aa_masses_all, ~ {
    calc_monopep(.x, .y,
                 maxn_vmods_per_pep = maxn_vmods_per_pep,
                 maxn_sites_per_vmod = maxn_sites_per_vmod,
                 min_mass = min_mass, 
                 max_mass = max_mass, 
                 digits = digits)
  })
  
  attrs <- purrr::map(aa_masses_all, attributes)
  vmods_ps <- map(attrs, `[[`, "vmods_ps")
  
  list(mass = ms, vmods_ps = vmods_ps)
}


#' Calculates the mono-isotopic mass of a peptide sequence.
#'
#' Only used for calc_monopeptide at a user's interface. Typically coupled to
#' \link{subpeps_by_vmods} for automatic dispatching of peptide sequences by
#' sets of fixed and variable modifications. For manual calculations, uses
#' \link{calc_monopeptide}.
#'
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @inheritParams add_var_masses
#' @inheritParams calc_pepmasses2
#' @import purrr
#' @importFrom stringr str_split
calc_monopep <- function (aa_seq, aa_masses,
                          maxn_vmods_per_pep = 5L,
                          maxn_sites_per_vmod = 3L,
                          min_mass = 200L, 
                          max_mass = 4500L, 
                          digits = 5L) 
{
  if (is.na(aa_seq)) return(NULL)
  
  aas <- stringr::str_split(aa_seq, "", simplify = TRUE)
  type <- attr(aa_masses, "type", exact = TRUE)
  
  # bare
  mass <- aas %>%
    aa_masses[.] %>%
    sum() %>%
    `+`(aa_masses["N-term"]) %>%
    `+`(aa_masses["C-term"]) %>%
    setNames(aa_seq) %>%
    round(digits = digits)
  
  if (type == "amods- tmod- vnl- fnl-") {
    return(mass)
  }
  
  # adds terminal mass
  if (grepl("tmod+", type, fixed = TRUE)) {
    mass <- add_term_mass(mass, aa_masses, min_mass, max_mass)
  }
  
  # --- Mass of variable mods and/or NLs ---
  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)
  fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
  vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
  amods <- attr(aa_masses, "amods", exact = TRUE)
  tmod <- attr(aa_masses, "tmod", exact = TRUE)
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
  if (TRUE) {
    if (type %in% c("amods- tmod- vnl- fnl+", "amods- tmod+ vnl- fnl+")) {
      fnl_combi <- expand_grid_rows(fmods_nl)
      deltas <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
      masses <- round(mass - deltas, digits = digits)
    }
    else {
      masses <- mass
    }
  } 
  else {
    masses <- mass
  }
  
  # `amods+`; (9-14) nested under (7-8)
  #
  # (7-8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #   (9-10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #   (11-12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  #   (13-14) "amods+ tmod- vnl+ fnl+", "amods+ tmod+ vnl+ fnl+"
  
  ok <- type %in% c("amods+ tmod- vnl- fnl-",
                    "amods+ tmod+ vnl- fnl-",
                    "amods+ tmod- vnl+ fnl-",
                    "amods+ tmod+ vnl+ fnl-",
                    "amods+ tmod- vnl- fnl+",
                    "amods+ tmod+ vnl- fnl+",
                    "amods+ tmod- vnl+ fnl+",
                    "amods+ tmod+ vnl+ fnl+")
  
  if (ok) {
    vmods_combi <- unique_mvmods(amods = amods, ntmod = NULL, ctmod = NULL,
                                 aa_masses = aa_masses, aas = aas,
                                 maxn_vmods_per_pep = maxn_vmods_per_pep,
                                 maxn_sites_per_vmod = maxn_sites_per_vmod,
                                 digits = digits) %>% 
      find_intercombi()

    deltas <- lapply(vmods_combi, function (x) sum(aa_masses[x]))
    
    masses <- 
      sapply(deltas, function (x) round(unlist(masses) + x, digits = digits))
  }
  
  masses
}


#' Checks the validity of a peptide sequence and dispatched it by fixedmods and
#' varmods.
#'
#' A sequence may be invalid at a given set of fixedmods and varmods. For
#' example, "MAKEMASSPECFUN" cannot have a mod of "Gln->pyro-Glu (N-term = Q)".
#' 
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param aa_masses_all All the amino acid lookup tables.
#' @inheritParams matchMS
check_aaseq <- function (aa_seq, aa_masses_all, fixedmods, varmods) 
{
  if (any(grepl("Protein N-term", fixedmods))) 
    stop("Need to change fixed 'Protein N-term' modification to variable.", 
         call. = FALSE)

  if (any(grepl("Protein C-term", fixedmods))) 
    stop("Need to change fixed 'Protein C-term' modification to variable.", 
         call. = FALSE)

  # ---
  has_prot_nt <- any(grepl("Protein N-term", varmods))
  
  if (has_prot_nt) {
    aa_seq <- paste0("-", aa_seq)
  }
  
  has_prot_ct <- any(grepl("Protein C-term", varmods))
  
  if (has_prot_ct) 
    aa_seq <- paste0(aa_seq, "-")
  
  peps <- lapply(aa_masses_all, subpeps_by_vmods, aa_seq) %>%
    purrr::flatten()
  
  if (has_prot_nt) 
    peps <- purrr::map(peps, ~ gsub("^-", "", .x))

  if (has_prot_ct) 
    peps <- purrr::map(peps, ~ gsub("-$", "", .x))

  invisible(peps)
}


#' Calculates the mono-isotopic mass of a MS2 ions.
#'
#' For direct uses from an R console (with trade-offs in speed).
#'
#' @inheritParams calc_ms2ions
#' @inheritParams calc_monopeptide
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' x <- calc_ms2ionseries("MAKEMASSPECFUN", 
#'                        fixedmods = NULL, 
#'                        varmods = NULL)
#' 
#' x$mass
#' 
#' # (2) no combinatorial NL for fixed modifications
#' x <- calc_ms2ionseries("MAKEMASSPECFUN", 
#'                        fixedmods = "Oxidation (M)", 
#'                        varmods = NULL)
#'                        
#' x$mass
#' 
#' ## With variable modifications
#' # (3) combinatorial sites and NL available
#' aa_seq <- "MAKEMASSPECFUN"
#' fixedmods <- NULL
#' varmods <- "Oxidation (M)"
#' ms <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = ms$mass[[2]][1])
#' 
#' x$mass
#' # x$vmods_ps
#' 
#' # (4)
#' aa_seq <- "MAKEMASSPECFUN"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)")
#' ms <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = ms$mass[[4]][[1]])
#'                       
#' x$mass
#' 
#' # (5) Neutral losses for occurrences of both fixed 
#' aa_seq <- "MAKEMASSPECFUN"
#' fixedmods <- c("TMT6plex (N-term)", "Oxidation (M)", "Deamidated (N)")
#' varmods <- c("dHex (S)")
#' ms <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' # and variable modifications ignored
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = 2002.7316)
#'                       
#' stopifnot(is.null(x$mass[[2]]))
#' 
#' # Change from fixed to variable for full combinatorials
#' fixedmods <- c("TMT6plex (N-term)", "Deamidated (N)")
#' varmods <- c("Acetyl (Protein N-term)", "Oxidation (M)", "dHex (S)")
#' 
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = 2002.7316)
#'                       
#' x$mass[[8]]
#' x$vmods_ps[[8]]
#' 
#' # (6) A lot of S
#' aa_seq <- "MAKEMASSSSSSPECFUNSS"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Acetyl (N-term)", "Oxidation (M)", "Deamidated (N)", 
#'              "Phospho (S)", "Phospho (T)", "Phospho (Y)", 
#'              "Gln->pyro-Glu (N-term = Q)")
#' ms <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = 2649.0553)
#' 
#' # (7) A lot of S and Y
#' aa_seq <- "MAKEMASSSSSSPECFUNSSYYYYYYY"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Acetyl (N-term)", "Oxidation (M)", "Deamidated (N)", 
#'              "Phospho (S)", "Phospho (T)", "Phospho (Y)", 
#'              "Gln->pyro-Glu (N-term = Q)")
#' ms <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' x <- calc_ms2ionseries(aa_seq, fixedmods = fixedmods, varmods = varmods, 
#'                        ms1_mass = ms$mass[[21]][[3]])
#' }
#' 
#' @export
calc_ms2ionseries <- function (aa_seq, fixedmods, varmods, 
                               type_ms2ions = "by", ms1_mass = NULL, 
                               maxn_vmods_setscombi = 512L,
                               maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               digits = 5L) 
{
  options(digits = 9L)
  
  aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
                                 varmods = varmods,
                                 maxn_vmods_setscombi = maxn_vmods_setscombi)
  
  peps <- check_aaseq(aa_seq, aa_masses_all, fixedmods, varmods)
  oks <- purrr::map_lgl(peps, ~ !purrr::is_empty(.x))
  peps <- peps[oks]
  aa_masses_all <- aa_masses_all[oks]
  
  mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
    as.hexmode() %>% 
    `names<-`(c(fixedmods, varmods))
  
  ms <- purrr::map2(peps, aa_masses_all, function (x, y) {
    pri <- calc_ms2ions(x, ms1_mass, y, mod_indexes, type_ms2ions, 
                        maxn_vmods_per_pep, maxn_sites_per_vmod, 
                        maxn_vmods_sitescombi_per_pep, digits)
    
    sec <- lapply(pri, add_seions, type_ms2ions = type_ms2ions, digits = digits)
    
    list(pri = pri, sec = sec)
  })
  
  attrs <- lapply(aa_masses_all, attributes)
  vmods_ps <- lapply(attrs, `[[`, "vmods_ps")
  
  list(mass = lapply(ms, `[[`, "pri"), 
       sec_mass = lapply(ms, `[[`, "sec"), 
       vmods_ps = vmods_ps)
}


#' Calculates the masses of MS2 ion series.
#'
#' For a given type of fragmentation. Minimal error handling for speeds.
#'
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @inheritParams calc_monopep
#' @inheritParams calc_aamasses
#' @inheritParams ms2match
#' @import purrr
#'
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' library(magrittr)
#'
#' fixedmods = NULL
#' varmods = NULL
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' x <- calc_ms2ions("MAKEMASSPECFUN", NULL, aa_masses_all[[1]], mod_indexes)
#'
#' }
calc_ms2ions <- function (aa_seq, ms1_mass = NULL, aa_masses, mod_indexes = NULL, 
                          type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                          maxn_sites_per_vmod = 3L, 
                          maxn_vmods_sitescombi_per_pep = 64L, digits = 5L) 
{
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # carbon <- 12.0
  # oxygen <- 15.99491462
  # nitrogen <- 14.003074
  # h3o_p <- 19.0178415
  # electron <- 0.000549
  
  # moverz <- (mass + h2o + z*proton)/z
  # moverz*z 
  
  if (is.na(aa_seq) || is.null(aa_seq)) 
    return(NULL)
  
  aas <- stringr::str_split(aa_seq, "", simplify = TRUE)
  type <- attr(aa_masses, "type", exact = TRUE)
  
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  if (type %in% c("amods- tmod- vnl- fnl-", 
                  "amods- tmod+ vnl- fnl-")) {
    
    ans <- gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                            aa_masses = aa_masses, 
                            ms1vmods = NULL, ms2vmods = NULL, 
                            ntmod = NULL, ctmod = NULL, 
                            ntmass = NULL, ctmass = NULL, 
                            amods = NULL, vmods_nl = NULL, fmods_nl = NULL, 
                            mod_indexes = mod_indexes, 
                            type_ms2ions = type_ms2ions, 
                            maxn_vmods_per_pep = maxn_vmods_per_pep, 
                            maxn_sites_per_vmod = maxn_sites_per_vmod, 
                            maxn_vmods_sitescombi_per_pep = 
                              maxn_vmods_sitescombi_per_pep, 
                            digits = digits)
    
    return(ans)
  }
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  #        (mutual exclusive btw. (1, 2) and (5, 6)
  #         "ANY" fmod has neuloss -> 5, 6;
  #         "ALL" fmods have no neuloss -> 1, 2)
  
  if (type %in% c("amods- tmod- vnl- fnl+", 
                  "amods- tmod+ vnl- fnl+")) {
    
    ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
    
    ntmass <- if (length(ntmod)) 
      aa_masses[names(ntmod)] + 1.00727647
    else 
      aa_masses["N-term"] - 0.000549

    ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
    
    ctmass <- if (length(ctmod)) 
      aa_masses[names(ctmod)] + 2.01510147
    else 
      aa_masses["C-term"] + 2.01510147

    fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
    
    ans <- gen_ms2ions_a0_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, ntmass = ntmass, 
                                    ctmass = ctmass, fmods_nl = fmods_nl, 
                                    mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    
    return(ans)
  }
  
  ms1vmods <- make_ms1vmod_i(aa_masses = aa_masses, 
                             maxn_vmods_per_pep = maxn_vmods_per_pep,
                             maxn_sites_per_vmod = maxn_sites_per_vmod)
  
  ms2vmods <- lapply(ms1vmods, make_ms2vmods)
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #        (ALL amods are vnl-)

  if (type %in% c("amods+ tmod- vnl- fnl-", 
                  "amods+ tmod+ vnl- fnl-")) {
    
    ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
    
    ntmass <- if (length(ntmod)) 
      aa_masses[names(ntmod)] + 1.00727647
    else 
      aa_masses["N-term"] - 0.000549
    
    ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
    
    ctmass <- if (length(ctmod)) 
      aa_masses[names(ctmod)] + 2.01510147
    else 
      aa_masses["C-term"] + 2.01510147
    
    amods <- attr(aa_masses, "amods", exact = TRUE)
    
    ans <- gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ms1vmods = ms1vmods, ms2vmods = ms2vmods,
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, 
                                    vmods_nl = NULL, fmods_nl = NULL,
                                    mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    
    return(ans)
  }
  
  # (9, 10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #         (ANY amod is vnl+)

  if (type %in% c("amods+ tmod- vnl+ fnl-", 
                  "amods+ tmod+ vnl+ fnl-")) {
    
    ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
    
    ntmass <- if (length(ntmod)) 
      aa_masses[names(ntmod)] + 1.00727647
    else 
      aa_masses["N-term"] - 0.000549
    
    ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
    
    ctmass <- if (length(ctmod)) 
      aa_masses[names(ctmod)] + 2.01510147
    else 
      aa_masses["C-term"] + 2.01510147

    amods <- attr(aa_masses, "amods", exact = TRUE)
    vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
    
    ans <- gen_ms2ions_a1_vnl1_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ms1vmods = ms1vmods, ms2vmods = ms2vmods,
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, vmods_nl = vmods_nl, 
                                    mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    
    return(ans)
  }
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  #          (mutual exclusive btw. (11, 12) and (7, 8);
  #           logicial ANY versus ALL)

  if (type %in% c("amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+")) {
    
    ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
    
    ntmass <- if (length(ntmod)) 
      aa_masses[names(ntmod)] + 1.00727647
    else 
      aa_masses["N-term"] - 0.000549

    ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
    
    ctmass <- if (length(ctmod)) 
      aa_masses[names(ctmod)] + 2.01510147
    else 
      aa_masses["C-term"] + 2.01510147

    amods <- attr(aa_masses, "amods", exact = TRUE)
    fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
    
    ans <- gen_ms2ions_a1_vnl0_fnl1(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ms1vmods = ms1vmods, ms2vmods = ms2vmods,
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, fmods_nl = fmods_nl, 
                                    mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    
    return(ans)
  }
  
  ans <- NULL
}


#' The unique combinations of variable modifications.
#'
#' The same residue, e.g. M, at different modifications, c("Carbamyl (M",
#' "Oxidation (M)")).
#'
#' Goes over all the \code{Anywhere} modifications specified in \code{amods} for
#' a given \code{aa_masses}.
#'
#' @param amods Anywhere modifications.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @param .ms1_vmodsets Not used.
#' @param .base_ent Not used.
#' @inheritParams matchMS
#' @inheritParams add_var_masses
#' @import purrr
#' @return Lists by residues in \code{amods}.
#' @seealso \link{ms1_a1_vnl0_fnl0} for examples.
#'
#' @examples
#' \donttest{
#' ## M
#' fixedmods = c("TMT6plex (K)", "dHex (S)")
#' varmods = c("Carbamidomethyl (M)", "Carbamyl (M)", "Acetyl (Protein N-term)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses <- aa_masses_all[[8]]
#'
#' amods <- list(`Carbamidomethyl (M)` = c(Anywhere = "M"),
#'               `Carbamyl (M)` = c(Anywhere = "M"))
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#' ans <- unique_mvmods(amods = amods, ntmod = NULL, ctmod = NULL,
#'                      aa_masses = aa_masses, aas = aas)
#'
#' stopifnot(length(ans) == 1L,
#'           length(ans[[1]]) == 3L)
#'
#' ## M and N
#' fixedmods = c("TMT6plex (K)", "dHex (S)")
#' varmods = c("Carbamidomethyl (M)", "Carbamyl (M)",
#'             "Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#'
#' aa_masses <- aa_masses_all[[16]]
#'
#' amods <- list(`Carbamidomethyl (M)` = c(Anywhere = "M"),
#'               `Carbamyl (M)` = c(Anywhere = "M"),
#'               `Deamidated (N)` = c(Anywhere = "N"))
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#' ans <- unique_mvmods(amods = amods, ntmod = NULL, ctmod = NULL,
#'                      aa_masses = aa_masses, aas = aas)
#'
#' stopifnot(length(ans) == 2L,
#'           length(ans[[1]]) == 3L,
#'           length(ans[[2]]) == 2L)
#' }
unique_mvmods <- function (amods, ntmod, ctmod, aa_masses, aas,
                           maxn_vmods_per_pep = 5L,
                           maxn_sites_per_vmod = 3L,
                           .ms1_vmodsets = NULL, 
                           .base_ent = NULL, 
                           digits = 5L) 
{
  # (6) "amods- tmod- vnl- fnl+"
  if (!length(amods)) 
    return(NULL)
  
  residue_mods <- .Internal(unlist(amods, recursive = FALSE, use.names = FALSE))
  names(residue_mods) <- names(amods)
  residue_mods <- split_vec(residue_mods)
  
  lapply(residue_mods, function (x) {
    vmods_elements(aas = aas, residue_mods = x, 
                   ntmod = ntmod, ctmod = ctmod,
                   maxn_vmods_per_pep = maxn_vmods_per_pep,
                   maxn_sites_per_vmod = maxn_sites_per_vmod,
                   .ms1_vmodsets = .ms1_vmodsets, 
                   .base_ent = .base_ent, 
                   digits = digits)
  })
}


#' Find the sets of variable modifications.
#'
#' The same residue, e.g. M, at different modifications, c("Carbamyl (M",
#' "Oxidation (M)")). 
#' 
#' Excluding position differences, i.e., \code{A, B} and \code{B, A} is the
#' same set.
#'
#' @param residue_mods Amino-acid residues with Unimod names. For example
#'   rownames of \code{Carbamidomethyl (M)} and \code{Oxidation (M)} and a
#'   column residues of \code{M, M}.
#' @inheritParams unique_mvmods
#' @import purrr
#' @examples 
#' \donttest{
#' ntmod <- list(`Acetyl (Protein N-term)` = c(`Protein N-term` = "N-term"))
#' 
#' ctmod <- list()
#' names(ctmod) <- character()
#' 
#' aas <- unlist(strsplit("HQGVMNVGMGQKSMNS", ""))
#' residue_mods <- c(`Carbamidomethyl (M)` = "M", `Carbamyl (M)` = "M")
#' 
#' x <- vmods_elements(aas, residue_mods, ntmod, ctmod)
#' }
vmods_elements <- function (aas,
                            residue_mods,
                            ntmod,
                            ctmod,
                            maxn_vmods_per_pep = 5L,
                            maxn_sites_per_vmod = 3L,
                            .ms1_vmodsets = NULL, 
                            .base_ent = NULL, 
                            digits = 5L) 
{
  residue <- residue_mods[[1]]
  
  ns <- names(residue_mods)
  len_n <- length(ns)
  
  # the exact positions not needed
  len_p <- sum(aas == residue)
  
  # i.e., btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  # MFGMFNVSMR cannot have three `Oxidation (M)` and `Acetyl (N-term)`
  
  len_nt <- length(ntmod)
  len_ct <- length(ctmod)
  
  if (len_nt && len_ct) {
    len_aas <- length(aas)
    aas_1 <- aas[1]
    aas_n <- aas[len_aas]
    if (aas_1 == residue && aas_n == residue) {
      len_p <- len_p - 2
    } else if ((aas_1 == residue) || (aas_n == residue)) {
      len_p <- len_p - 1
    }
  } else if (len_nt) {
    aas_1 <- aas[1]
    if (aas_1 == residue) {
      len_p <- len_p - 1
    }
  } else if (len_ct) {
    aas_n <- aas[len_aas]
    if (aas_n == residue) {
      len_p <- len_p - 1
    }
  }
  
  if (len_p <= 0) 
    return(list())
  
  len_p <- min(len_p, maxn_vmods_per_pep)
  
  if (is.null(.ms1_vmodsets) || is.null(.base_ent)) {
    if (len_p > len_n) {
      x <- lapply((len_n + 1):len_p, function (x) find_unique_sets(x, ns))
      x <- .Internal(unlist(x, recursive = FALSE, use.names = FALSE))
      x <- c(list(ns), x)
    } else {
      x <- list(ns)
    }
    
    maxn_vmod <- lapply(x, count_elements)
    maxn_vmod <- lapply(maxn_vmod, max)
    rows <- (maxn_vmod <= maxn_sites_per_vmod)
    
    x <- x[rows]
  } else {
    # x <- extract_vmodsets(.ms1_vmodsets, .base_ent, len_p, ns)
  }
  
  invisible(x)
}


#' Finds the combinations across residues.
#'
#' For uses with MS1 precursors. For multiple residues (each residue one to
#' multiple modifications).
#'
#' @param intra_combis The results from \link{unique_mvmods}.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' C <- list(c("Carbamidomethyl (C)"),
#'           rep("Carbamidomethyl (C)", 2))
#'
#' N <- list(c("Deamidated (N)"),
#'           rep("Deamidated (N)", 2))
#'
#' intra_combis <- list(C = C, N = N)
#'
#' ans <- find_intercombi(intra_combis)
#' 
#' # three large lists
#' S <- list(c("Carbamidomethyl (S)", "Phospho (S)", "Phospho (S)"),
#'            c("Carbamidomethyl (S)", "Phospho (S)", "Phospho (S)", "Phospho (S)"))
#' 
#' M <- list(c("Oxidation (M)", "Carbamidomethyl (M)"), 
#'           c("Oxidation (M)", "Carbamidomethyl (M)", "Carbamyl (M)"))
#' 
#' N <- list(c("Deamidated (N)"),
#'           rep("Deamidated (N)", 2))
#' 
#' ans <- find_intercombi(list(S = S, M = M, N = N))
#' }
find_intercombi <- function (intra_combis, maxn_vmods_per_pep = 5L) 
{
  len <- length(intra_combis)
  
  if (!len) { # scalar
    v_out <- list()
  } else if (any(.Internal(unlist(lapply(intra_combis, purrr::is_empty), 
                                  recursive = FALSE, use.names = FALSE)))) { # list
    v_out <- list()
  } else if (len > 1L) {
    v_out <- expand_grid_rows(intra_combis, use.names = FALSE)
    
    lens <- lapply(v_out, length)
    lens <- .Internal(unlist(lens, recursive = FALSE, use.names = FALSE))
    oks <- (lens <= maxn_vmods_per_pep)
    v_out <- v_out[oks]
    
    ## DONT: 
    ## each aa_masses is a realization of a set of combinatorial amods;
    ## if length > maxn_vmods_per_pep, the combination should be dropped;
    ## duplicated entries if subset by 1:min(length(x), maxn_vmods_per_pep)
    
    # v_out <- lapply(v_out, function (x) x[1:min(length(x), maxn_vmods_per_pep)])
  } else {
    v_out <- purrr::flatten(intra_combis)
  }
  
  invisible(v_out)
}

