#' Matches to the pre-calculated labels of combinatorial MS1 variable
#' modifications.
#'
#' No checking of conflicting [NC]-term and Anywhere sites. Even so, it only
#' sets a overall limit on the counts. Still need to check the MS2 permutation
#' with indexes in aas for exact exclusion of certain permutations.
#'
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @param ms1vmods The i-th result from lapply(aa_masses_all, make_ms1vmod_i).
#' @param amods \code{Anywhere} variable modifications.
#' @examples
#' \donttest{
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)",
#'                "Carbamidomethyl (C)")
#'
#' varmods <- c("Acetyl (Protein N-term)", "Oxidation (M)",
#'              "Deamidated (N)",
#'              "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' i <- 11L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#' ms1vmods <- ms1vmods_all[[i]]
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#' vmods_combi <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
#'
#' ## M and N
#' fixedmods = c("TMT6plex (K)", "dHex (S)")
#' varmods = c("Carbamidomethyl (M)", "Carbamyl (M)",
#'             "Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' i <- 16L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#' ntmod <- attr(aa_masses, "ntmod")
#' ctmod <- attr(aa_masses, "ctmod")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#' vmods_combi <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)
#'
#' stopifnot(sapply(vmods_combi$ms1, function (x) all(names(amods) %in% x)),
#'           length(vmods_combi$ms1) == 6L)
#'
#' }
match_mvmods <- function (aas = NULL, ms1vmods = NULL, amods = NULL) {
  
  ## stronger check by each residues
  resids <- unique(amods)
  len_r <- length(resids)
  max_rs <- vector("integer", len_r)
  
  for (i in seq_len(len_r)) {
    max_rs[i] <- sum(aas == resids[[i]])
  }
  
  ## `make_ms1_vmodsets` obtained from the same `amods`
  ##   -> no mess up in the order of `amods` -> no name sorting
  
  len <- length(ms1vmods)
  rows <- vector("logical", len)
  
  for (i in seq_len(len)) {
    ps <- attr(ms1vmods[[i]], "ps")
    rows[i] <- all(ps <= max_rs) # && all(labs <= max_ls)
  }
  
  rows <- which(rows)
  
  list(ms1 = ms1vmods[rows], inds = rows)
}


#' Makes the sets of MS1 labels for a given \code{aa_masses}.
#'
#' For \code{Anywhere} variable modifications (\code{amods}) across all
#' residues.
#'
#' For the universe of labels, loops through \code{aa_masses_all}.
#'
#' By the design of \code{aa_masses}, the \code{amods} in a \code{aa_masses} are
#' all realized. Therefore, each list in the resulted labels should contain at
#' least one of the \code{amods} residues from the \code{aa_masses}. For
#' example, if \code{amods} contain M, N and S, each list in the result should
#' contains at least one of the residues.
#'
#' Currently, variable terminal modifications are exempted from the restriction
#' by \code{maxn_vmods_per_pep}. If taken into account (e.g. variable N-term),
#' the combination of \code{Anywhere} can be further reduced.
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' ## Simple
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)",
#'                "Carbamidomethyl (C)")
#'
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)",
#'             "Deamidated (N)",
#'             "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' stopifnot(length(ms1vmods_all[[1]]) == 0L,
#'           length(ms1vmods_all[[2]]) == 0L,
#'           length(ms1vmods_all[[3]]) == 0L)
#'
#' # (M, N)
#' len_ps <- lapply(ms1vmods_all[[12]], function (x) attr(x, "ps"))
#'
#' # (M)
#' len_ps <- lapply(ms1vmods_all[[4]], function (x) attr(x, "ps"))
#'
#' ## More complex
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")
#'
#' varmods <- c("Acetyl (Protein N-term)", "Deamidated (N)",
#'              "Oxidation (M)", "Carbamidomethyl (M)", "Carbamyl (M)",
#'              "Carbamidomethyl (S)", "Phospho (S)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' # No duplication within each aa_masses
#' any_dups <- lapply(ms1vmods_all,
#'                    function (x) anyDuplicated(x))
#'
#' stopifnot(all(unlist(any_dups) == 0L))
#'
#' # Can have identical sets of labels at different aa_masses
#' identical(ms1vmods_all[[3]], ms1vmods_all[[9]])
#'
#' attr(aa_masses_all[[3]], "tmod")
#' attr(aa_masses_all[[9]], "tmod")
#'
#' ms1vmods <- ms1vmods_all[[64]]
#' n <- length(ms1vmods[[1]])
#' 
#' 
#' ## "Oxidation (M)" on the N-term
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#'
#' varmods = c("Acetyl (Protein N-term)", "Oxidation (M)", "Deamidated (N)",
#'             "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#' 
#' i <- 6L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
#' 
#' aas <- unlist(strsplit("MHQGVMNVGMGQKNS", ""))
#' 
#' # Subset from ms1vmods by aas
#' oks <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)$inds
#' ms2vmods <- ms2vmods[oks]
#'
#' vmods_combi <- find_vmodscombi(aas, ms2vmods)
#' 
#' # "Oxidation (M)" on N-term kept in the label space
#' lapply(vmods_combi, function (x) any(names(x) == 1))
#' 
#' # Otherwise check if any attr(aa_masses, "amods") in the N-term
#' amods <- attr(aa_masses, "amods")
#' # may do it during `combi_byvmodsM`...
#' }
make_ms1vmod_i <- function (aa_masses = NULL, maxn_vmods_per_pep = 5L, 
                            maxn_sites_per_vmod = 3L) {
  
  # stopifnot(maxn_vmods_per_pep >= 2L)
  
  vmodsets <- make_ms1_vmodsets(aa_masses_all = list(aa_masses), 
                                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                maxn_sites_per_vmod = maxn_sites_per_vmod)
  
  find_intercombi2(vmodsets = vmodsets, maxn_vmods_per_pep = maxn_vmods_per_pep)
}


#' Makes the sets of labels of variable modifications.
#' 
#' No position permutation (for MS1 masses).
#' 
#' @param aa_masses_all All the amino acid lookup tables.
#' @inheritParams matchMS
make_ms1_vmodsets <- function (aa_masses_all = NULL, maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L) {
  
  # stopifnot(maxn_vmods_per_pep >= 2L)
  
  if (!is.list(aa_masses_all)) aa_masses_all <- list(aa_masses_all)
  
  amods_all <- lapply(aa_masses_all, attr, "amods", exact = TRUE)
  
  len <- length(amods_all)
  resmods_all <- vector("list", len)
  
  for (i in seq_len(len)) {
    aa_masses <- aa_masses_all[[i]]
    amods <- attr(aa_masses, "amods", exact = TRUE)
    
    # split into lists by individual residues
    if (length(amods)) {
      x <- .Internal(unlist(amods, recursive = FALSE, use.names = FALSE))
      names(x) <- names(amods)
      resmods_all[[i]] <- split_vec(x)
    } else {
      resmods_all[[i]] <- NULL
    }
  }
  
  # unique sets of modifications by individual residues
  resmods_all <- unlist(resmods_all, recursive = FALSE, use.names = FALSE)
  resmods_all <- unique(resmods_all)
  
  out <- lapply(resmods_all, bacth_vmods_combi, 
                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                maxn_sites_per_vmod = maxn_sites_per_vmod)
  
  resids <- lapply(resmods_all, `[`, 1)
  resids <- unlist(resids, recursive = FALSE, use.names = FALSE)
  names(out) <- resids
  
  invisible(out)
}


#' Helper in cycling through \code{2:maxn_vmods_per_pep} numbers of positions.
#'
#' For the same residue (M) at different names of modifications.
#'
#' @param resmods A vector of amino-acid residues with Unimod in names. For
#'   example, c(`Oxidation (M)` = "M", `Carbamyl` = "M").The residues are the
#'   same and the names are different.
#' @inheritParams matchMS
#' @examples 
#' \donttest{
#' resmods <- c(`Oxidation (M)` = "M", `Carbamidomethyl (M)` = "M")
#' ans <- bacth_vmods_combi(resmods)
#' }
bacth_vmods_combi <- function (resmods = NULL, maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L) {
  
  # stopifnot(maxn_vmods_per_pep >= 2L)
  
  resid <- unname(resmods[[1]])

  ans <- make_unique_sets(p = maxn_vmods_per_pep, 
                          n = length(resmods), 
                          labs = names(resmods), 
                          # resid = resmods, 
                          maxn_vmods_per_pep = maxn_vmods_per_pep, 
                          maxn_sites_per_vmod = maxn_sites_per_vmod)
  
  lapply(ans, function (x) {
    names(x) <- rep(resid, length(x))
    x
  })
}


#' Makes the unique sets of modifications by individual AA residues. 
#' 
#' @param p The number of open positions for filling.
#' @param n The number of labels.
#' @param labs The labels of balls.
#' @inheritParams matchMS
#' @examples 
#' \donttest{
#' # multiple outs
#' p <- 5
#' labs <- c("Oxidation (M)", "Carbamidomethyl (M)", "Carbamyl (M)")
#' n <- length(labs)
#' stopifnot(n == length(labs))
#' 
#' ans5 <- make_unique_sets(5, n, labs)
#' ans7 <- make_unique_sets(7, n, labs)
#' 
#' stopifnot(identical(ans5, ans7))
#' 
#' # single out
#' p <- 3
#' labs <- c("Oxidation (M)", "Carbamidomethyl (M)", "Carbamyl (M)")
#' n <- length(labs)
#' stopifnot(n == length(labs))
#' 
#' ans <- make_unique_sets(p, n, labs)
#' }
make_unique_sets <- function (p = 5L, n = 2L, labs = c("X", "Y"), 
                              maxn_vmods_per_pep = 5L, 
                              maxn_sites_per_vmod = 3L) {
  
  # stopifnot(n == length(labs))
  
  if (p < n) return(NULL)
  
  p <- min(p, maxn_vmods_per_pep)
  ps <- n:p
  
  combs <- lapply(ps, function (i) {
    ans <- find_unique_sets(i, labs)
    
    n_vmod <- lapply(ans, count_elements)
    maxn_vmod <- lapply(n_vmod, max)
    rows <- (maxn_vmod <= maxn_sites_per_vmod)
    
    ans <- ans[rows]
  })
  
  lens <- lapply(combs, length)
  lens <- unlist(lens, recursive = FALSE, use.names = FALSE)
  nms_p <- rep(ps, lens)
  nms_n <- rep(n, length(nms_p))
  
  combs <- .Internal(unlist(combs, recursive = FALSE, use.names = FALSE))
  names(combs) <- nms_p
  
  invisible(combs)
}


#' Finds the sets of unique labels for \code{n} balls in \code{p-positions}.
#'
#' At least one occupancy for each balls in \code{ns}.
#'
#' @param p The number of positions.
#' @param labs The names to be filled into the \code{p}-number of positions.
#' @importFrom gtools combinations
#' @examples 
#' \donttest{
#' find_unique_sets(5, c("Oxidation (M)", 
#'                       "Carbamidomethyl (M)", 
#'                       "Carbamyl (M)"))
#' }
find_unique_sets <- function (p = 5L, labs = c("A", "B", "C")) {
  
  ps <- seq_len(p)
  n <- length(labs)
  r <- p - n
  
  if (r == 0L) return(list(labs))
  if (r < 0L) return(NULL)
  
  x <- gtools::combinations(n, r, labs, repeats = TRUE)
  
  n_row <- nrow(x)
  out <- vector("list", n_row)
  
  for (i in seq_len(n_row)) {
    out[[i]] <- c(labs, x[i, ])
  }
  
  out
}


#' Finds the combinations across residues (with attributes).
#'
#' Used after \link{make_ms1_vmodsets} across residues.
#'
#' @param vmodsets The results from \link{make_ms1_vmodsets}.
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' ## Simple
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")
#'
#' varmods <- c("Acetyl (Protein N-term)",
#'              "Oxidation (M)", "Deamidated (N)", 
#'              "Gln->pyro-Glu (N-term = Q)")
#' 
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' # M and N
#' i <- 10L
#' aa_masses <- aa_masses_all[[i]]
#'
#' vmodsets <- make_ms1_vmodsets(aa_masses_all = aa_masses,
#'                               maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                               maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' ms1vmods <- find_intercombi2(vmodsets)
#' 
#' ## More complex
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")
#'
#' varmods <- c("Acetyl (Protein N-term)",
#'              "Oxidation (M)", "Carbamidomethyl (M)",
#'              "Deamidated (N)", "Carbamyl (M)",
#'              "Gln->pyro-Glu (N-term = Q)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' # Only `Oxidation (M)` 
#' i <- 8L
#' aa_masses <- aa_masses_all[[i]]
#'
#' vmodsets <- make_ms1_vmodsets(aa_masses_all = aa_masses,
#'                               maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                               maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' ms1vmods <- find_intercombi2(vmodsets)
#' 
#' ## `Oxidation (M)`, "Carbamidomethyl (M)"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")
#'
#' varmods <- c("Acetyl (Protein N-term)",
#'              "Oxidation (M)", "Carbamidomethyl (M)",
#'              "Deamidated (N)")
#'
#' aa_masses_all <- calc_aamasses(fixedmods = fixedmods,
#'                                varmods = varmods,
#'                                maxn_vmods_setscombi = 64,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE,
#'                                exclude_phospho_nl = TRUE,
#'                                out_path = NULL)
#' 
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#' 
#' i <- 16L
#' aa_masses <- aa_masses_all[[i]]
#'
#' vmodsets <- make_ms1_vmodsets(aa_masses_all = aa_masses,
#'                               maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                               maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' ms1vmods <- find_intercombi2(vmodsets)
#' 
#' }
find_intercombi2 <- function (vmodsets = NULL, maxn_vmods_per_pep = 5L) {
  
  len <- length(vmodsets)
  
  # scalar
  if (!len) return(list())
  
  # empty if any is empty
  if (any(.Internal(unlist(lapply(vmodsets, purrr::is_empty), 
                           recursive = FALSE, use.names = FALSE)))) 
    return(list())
  
  # all lists are non-empty
  vmod_nms <- lapply(vmodsets, function (x) lapply(x, names))
  len_ps <- lapply(vmodsets, function (x) lapply(x, length))
  
  v_out <- expand_grid_rows(vmodsets, use.names = FALSE)
  nm_out <- expand_grid_rows(vmod_nms, use.names = FALSE)

  v_out <- mapply(function (x, y) {
    names(x) <- y
    x
  }, v_out, nm_out, 
  SIMPLIFY = FALSE, USE.NAMES = FALSE)

  ps <- expand_grid_rows(len_ps, use.names = TRUE)
  
  ## Characteristics of expand grids: 
  # (residues will be in sections, e.g. M, M, M, N)
  #  M                      M                  M                   N       
  # "Carbamidomethyl (M)", "Carbamyl (M)", "Carbamidomethyl (M)", "Deamidated (N)"
  # 
  # The same order for ps (M > N)
  # M N 
  # 3 1
  # 
  # Therefore can safely tell from ps that the first three are M(s)

  lens <- lapply(v_out, length)
  lens <- .Internal(unlist(lens, recursive = FALSE, use.names = FALSE))
  oks <- (lens <= maxn_vmods_per_pep)
  
  v_out <- v_out[oks]
  ps <- ps[oks]
  labs <- lapply(v_out, count_elements)
  
  mapply(function (x, y, z) {
    attr(x, "ps") <- y
    attr(x, "labs") <- z
    x
  }, v_out, ps, labs, 
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


