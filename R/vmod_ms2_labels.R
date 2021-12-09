#' Finds the combinatorial MS2 variable modifications.
#'
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @param ms2vmods The i-th result from
#'
#'   lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods)).
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' ## One-to-one correspondance between Names and Sites
#' #  (no need to permutate MS1 labels)
#'
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
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#'
#' i <- 11L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#' # Subset from ms1vmods by aas
#' oks <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)$inds
#' ms2vmods <- ms2vmods[oks]
#'
#' vmods_combi <- find_vmodscombi(aas, ms2vmods)
#'
#'
#' ## 'Carbamidomethyl (M)',  'Carbamyl (M)' and N
#' #  (need permutation of MS1 labels)
#'
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
#' ms2vmods_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#'
#' i <- 16L
#' aa_masses <- aa_masses_all[[i]]
#' amods <- attr(aa_masses, "amods")
#' ntmod <- attr(aa_masses, "ntmod")
#' ctmod <- attr(aa_masses, "ctmod")
#'
#' ms1vmods <- ms1vmods_all[[i]]
#' ms2vmods <- ms2vmods_all[[i]]
#'
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#'
#'
#' # Subset from ms1vmods by aas
#' oks <- match_mvmods(aas = aas, ms1vmods = ms1vmods, amods = amods)$inds
#' ms2vmods <- ms2vmods[oks]
#'
#' vmods_combi <- find_vmodscombi(aas, ms2vmods)
#' 
#' n_pos <- lapply(vmods_combi, function (x) names(x[x == "Deamidated (N)"]))
#' stopifnot(all(sapply(n_pos, function (x) all(x %in% c("6", "14")))))
#' }
find_vmodscombi <- function (aas = NULL, ms2vmods = NULL, 
                             maxn_vmods_sitescombi_per_pep = 64L) 
{
  # Starts from sets of combinatorial MS1 labels
  # 
  #   if multiple names for the same site 
  #     (e.g., Carbamidomethyl (M)", "Oxidation (M)")
  #     -> permutation of the names
  #   else one-to-one correspondence between sites and names 
  #     (e.g., M <-> "Oxidation (M)")
  #     -> single-row matrix
  # 
  #   -> permutation (by indexes matched to `aas`)
  #        (5, 9 | 5, 13 | 9, 13...)
  
  len <- length(ms2vmods)
  pos <- vector("list", len)
  tot <- 0L
  
  for (i in 1:len) {
    M <- ms2vmods[[i]]
    nrows <- nrow(M)
    
    if (nrows == 1L) {
      ans <- combi_namesiteU(M = M, aas = aas)
    } else {
      ans <- combi_namesiteM(M = M, aas = aas, nrows = nrows)
      ans <- .Internal(unlist(ans, recursive = FALSE, use.names = FALSE))
    }

    len_a <- length(ans)
    tot <- tot + len_a
    
    if (tot > maxn_vmods_sitescombi_per_pep) {
      lags <- len_a + maxn_vmods_sitescombi_per_pep - tot
      pos[[i]] <- ans[seq_len(lags)]
      
      break
    } else {
      pos[[i]] <- ans
    }
  }
  
  .Internal(unlist(pos[1:i], recursive = FALSE, use.names = FALSE))

  ## Level 1 (permutated by MS1 labels): 
  # 
  # x = out[[1]] 
  # 
  # Belongs to the same MS1 label set, e.g. 
  #   "Carbamidomethyl (M)", "Carbamyl (M)", "Deamidated (N)"
  # at (six) different permutations of 
  #   "Carbamidomethyl (M)", "Carbamyl (M)", "Deamidated (N)"
  #   "Carbamyl (M)", "Carbamidomethyl (M)", "Deamidated (N)"
  #   ...                                          
  #                                                ^
  ## Levle 2 (permutated by `aas` indexes):        ^
  #                                                ^
  # x1 = x[[1]]                                    ^
  #                                                ^
  # Belongs to the same permutated MS1 label, e.g.
  #   "Carbamyl (M)", "Carbamidomethyl (M)", "Deamidated (N)"
  # 
  # with sub lists at different `aas` indexes 
  #         5                 9                    14
  #   "Carbamyl (M)", "Carbamidomethyl (M)", "Deamidated (N)"
  #         5                 13                    14
  #   "Carbamyl (M)", "Carbamidomethyl (M)", "Deamidated (N)"
}



#' Helper of \link{combi_vmodsMat} (by each rows of labels).
#'
#' One-to-one correspondence between Names and Sites. Finds the positions of
#' residues (sites) from a given amino acid sequence (aas).
#'
#' @param M A matrix of labels of modification names with permutation. Each
#'   permutation in a row. The corresponding matrix of site permutations in the
#'   attribute of \code{resids}.
#'
#'   Note that M is a matrix other than lists of vectors, which allows the
#'   application of one copy of attributes to all rows.
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
combi_namesiteU <- function (M, aas) 
{
  m <- attr(M, "resids")
  ps <- attr(M, "ps")
  
  ans <- find_vmodposU(m, ps, aas)
  combi <- ans$combi
  vpos <- ans$vpos
  
  len_out <- nrow(combi)
  out <- rep(list(M[1, ]), len_out)
  
  for (i in seq_along(vpos)) { # by residue
    ansi <- combi[[i]] # list of six: 5, 9; 9, 13 etc.
    pi <- vpos[[i]] # 1, 3
    
    for (j in seq_len(len_out)) 
      names(out[[j]])[pi] <- ansi[[j]] # by combi
  }
  
  out
}


#' Helper of \link{combi_namesiteU}.
#'
#' One-to-one correspondence between Names and Sites.
#' 
#' Custom functions: vec_to_list.
#'
#' @param vec A vector of labels.
#' @param ps Named vector; counts for each site. Sites in names and counts in
#'   values.
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
find_vmodposU <- function (vec, ps, aas) 
{
  nres <- length(ps)
  M <- vpos <- vector("list", nres)
  
  for (i in seq_len(nres)) {
    resid <- names(ps)[i] # M
    aapos <- which(aas == resid) # M:5, 9, 13; N: 6, 14
    
    ct <- ps[[i]] # M: 2; N: 1
    vpos[[i]] <- which(vec == resid) # M: 1, 2; N: 3
    
    if (ct == 1L) 
      M[[i]] <- vec_to_list(aapos) 
    else 
      M[[i]] <- sim_combn(aapos, ct)
  }
  
  list(
    combi = expand.grid(M, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), 
    vpos = vpos
  )
}


#' Helper of \link{combi_vmodsMat} (by each rows of labels).
#'
#' Multiple Names to the same Site. Finds the positions of
#' residues (sites) from a given amino acid sequence (aas).
#' 
#' @param nrows The number of rows of M.
#' @inheritParams combi_namesiteU
combi_namesiteM <- function (M, aas, nrows) 
{
  ps <- attr(M, "ps")
  m <- attr(M, "resids")
  
  # convert to vectors for speed
  mv <- vector("list", nrows)
  for (i in 1:nrows) mv[[i]] <- m[i, ]
  uniqs <- !duplicated.default(mv)
  umv <- mv[uniqs]
  
  len <- length(mv)
  cache <- ans <- vector("list", len)

  for (i in 1:len) {
    Vec <- M[i, ]
    vec <- mv[[i]]
    is_new <- uniqs[i]
    
    if (is_new) {
      ans[[i]] <- find_vmodposM(Vec = Vec, vec = vec, ps = ps, aas = aas)
      cache[[i]] <- vec
    } else {
      # must have a preceding match by the way of `duplicated`
      for (j in 1:len) {
        cj <- cache[[j]]
        
        if ((!is.null(cj)) && identical(vec, cj)) {
          ans[[i]] <- match_aas_indexes(ans[[j]], Vec)
          break
        }
      }
    }
  }
  
  ans
}


#' Helper of \link{combi_namesiteM} (by each rows of labels).
#'
#' Multiple Names to the same Site.
#' 
#' @param Vec A vector of names (lower-case vec for sites).
#' @inheritParams find_vmodposU
find_vmodposM <- function (Vec, vec, ps, aas) 
{
  nres <- length(ps)
  M <- vpos <- vector("list", nres)
  
  for (i in seq_len(nres)) { # by residues
    resid <- names(ps)[i] # M
    aapos <- which(aas == resid) # M:5, 9, 13; N: 6, 14
    
    ct <- ps[[i]] # M: 2; N: 1
    vpos[[i]] <- which(vec == resid) # M: 1, 2; N: 3
    
    if (ct == 1L) 
      M[[i]] <- vec_to_list(aapos) 
    else 
      M[[i]] <- sim_combn(aapos, ct)
  }
  
  ans <- expand.grid(M, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  len_out <- nrow(ans)
  out <- rep(list(Vec), len_out)
  
  for (i in seq_along(vpos)) { # by residue
    ansi <- ans[[i]] # list of six: 5, 9; 9, 13 etc.
    pi <- vpos[[i]] # 1, 3
    
    for (j in seq_len(len_out)) 
      names(out[[j]])[pi] <- ansi[[j]] # by combi
  }
  
  out
}


#' Matches the indexes of amino-acid residues to cached results.
#' 
#' @param X Lists of cached results.
#' @param Vec A vector of names (lower-case vec for sites).
match_aas_indexes <- function (X, Vec) 
{
  len <- length(X)
  out <- rep(list(Vec), len)
  
  for (i in 1:len) 
    names(out[[i]]) <- names(X[[i]])

  out
}


#' Makes the sets of MS2 labels (with permutations) for an \code{aa_masses}.
#'
#' For \code{Anywhere} variable modifications (\code{amods}) across all
#' residues.
#'
#' For the universe of labels, loops through \code{aa_masses_all} and the
#' indexes of lists are in one-to-one correspondence to \code{aa_masses_all}.
#'
#' By the design of \code{aa_masses}, the \code{amods} in a \code{aa_masses} are
#' all realized. Therefore, each list in the resulted labels should contain at
#' least one of the \code{amods} residues from the \code{aa_masses}. For
#' example, if \code{amods} contain M, N and S, each list in the result should
#' contains at least one of the residues.
#'
#' @param vec A named vector of labels.
#' @examples
#' \donttest{
#' ## with a bare vector
#' # One-to-one correspondance between Names and Sites
#' vec <- c(M = "Oxidation (M)", N = "Deamidated (N)")
#' attr(vec, "ps") <- c(M = 3L, N = 2L)
#' attr(vec, "labs") <- c(`Oxidation (M)` = 3L, 
#'                        `Deamidated (N)` = 2L)
#' 
#' ans <- make_ms2vmods(vec)
#' 
#' # Multiple Names to the same Site (S)
#' vec <- c(M = "Oxidation (M)", 
#'          S = "Carbamidomethyl (S)", 
#'          S = "Phospho (S)")
#' attr(vec, "ps") <- c(M = 2L, S = 3L)
#' attr(vec, "labs") <- c(`Oxidation (M)` = 2L, 
#'                        `Carbamidomethyl (S)` = 2L, 
#'                        `Phospho (S)` = 1L)
#' 
#' ans <- make_ms2vmods(vec)
#' stopifnot(nrow(ans) == 6L, ncol(ans) == length(vec))
#' 
#' # Another one-to-multiple
#' vec <- c(M = "Oxidation (M)", M = "Carbamyl (M)", 
#'          S = "Carbamidomethyl (S)", S = "Phospho (S)")
#' attr(vec, "ps") <- c(M = 3L, S = 2L)
#' attr(vec, "labs") <- c(`Oxidation (M)` = 2L, 
#'                        `Carbamyl (M)` = 1L, 
#'                        `Carbamidomethyl (S)` = 1L, 
#'                        `Phospho (S)` = 1L)
#' 
#' ans <- make_ms2vmods(vec)
#' stopifnot(nrow(ans) == 24L, ncol(ans) == length(vec), 
#'           nrow(ans) == nrow(unique(ans)))
#' 
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
#' # M, N (up to 3 M's or N's and 5 in total)
#' len_ps <- lapply(ms1vmods_all[[12]], function (x) attr(x, "ps"))
#' 
#' # M (up to 3 M's)
#' len_ps <- lapply(ms1vmods_all[[4]], function (x) attr(x, "ps"))
#' 
#' x <- ms1vmods_all[[12]] # list of 8
#' # x <- ms1vmods_all[[1]] # empty list
#' 
#' ans <- lapply(x, make_ms2vmods) # list of 8 matrices
#'
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
#' x <- ms1vmods_all[[64]] # list of 10 vectors
#' 
#' # [[1]]
#' #  M                     S                     S 
#' # "Carbamyl (M)"  "Carbamidomethyl (S)"  "Phospho (S)" 
#' # [[2]]
#' # M                     M                     S               S 
#' # "Carbamyl (M)"  "Carbamyl (M)"  "Carbamidomethyl (S)"  "Phospho (S)" 
#' # [[3]]
#' # ...
#' 
#' ans <- lapply(x, make_ms2vmods) # list of 10 matrices
#' 
#' ans_all <- lapply(ms1vmods_all, function (x) lapply(x, make_ms2vmods))
#' }
make_ms2vmods <- function (vec = NULL) 
{
  # stopifnot(maxn_vmods_per_pep >= 2L)
  
  n_nms <- length(unique(names(vec)))
  n_vals <- length(unique(vec))
  
  if (n_nms == n_vals)
    ans <- matrix(vec, nrow = 1L)
  else 
    ans <- find_perm_sets(vec)
  
  attr(ans, "ps") <- attr(vec, "ps")
  attr(ans, "labs") <- attr(vec, "labs")
  attr(ans, "resids") <- find_ms2resids(ans, vec)
  
  ans
}


#' Helper of \link{make_ms2vmods}.
#' 
#' Makes the corresponding residue matrix from \link{make_ms2vmods}.
#' 
#' @param M A modification-name matrix from \link{make_ms2vmods}.
#' @param vec A vector of labels.
find_ms2resids <- function (M, vec) 
{
  vecinv <- names(vec)
  names(vecinv) <- vec
  vecinv <- vecinv[unique(names(vecinv))]
  
  for (i in 1:nrow(M)) M[i, ] <- vecinv[M[i, ]]

  invisible(M)
}




##
# labels: n
# positions: p
# 
# n = c("A", "A", "B"); n1 = 2, n2 = = 1
# p = c(1, 3, 5, 9)
# n = c("Carbamidomethyl (M)",  "Carbamidomethyl (M)", "Oxidation (M)"); 
# p = c(1, 3, 5, 9)
# 
# choose(p, n1 + n2) * choose(n1 + n2, n2)
# = choose(4, 3) * choose(3, 2)
# 
# n = c("A", "A", "B", "B", "C"); n1 = 2, n2 = = 2, n3 = 1
# p = c(1, 3, 5, 7, 8, 9)
# 
# p!/(n1!n2!n3!)/(p-n1-n2-n3)!
# = choose(p, n1+n2+n3) * (n1+n2+n3)!/(n1!n2!n3!)
# = choose(p, n1+n2+n3) * choose(n1+n2+n3, n1+n2) * choose(n1+n2, n1)
##



#' Helper of \link{find_perm_sets}.
#' 
#' Adds one more labels to a permutation table.
#' 
#' @param M A permutation table with unique sets by rows.
#' @param x A label to be added to M for additional permutations.
add_one_label <- function (M, x) 
{
  ncols <- ncol(M)
  nrows <- nrow(M)
  
  out <- matrix(nrow = ncols * nrows, ncol = ncols + 1L)
  
  k = 1L
  
  for (i in 1:nrows) {
    mi <- M[i, ]
    
    if (ncols > 1L) {
      for (j in 1:(ncols-1L)) {
        out[k, ] <- c(mi[1:j], x, mi[(j+1):ncols])
        k <- k + 1L
      }
    }
    
    out[k, ] <- c(mi, x)
    k <- k + 1L
  }
  
  out[!duplicated.matrix(out), ]
}


#' A (faster alternative) to \link[gtools]{permutations} with duplicated labels.
#'
#' @param labs A vector of labels.
#'
#' @examples
#' \donttest{
#' labs <- c("A", "A", "A", "B", "B", "C")
#' x <- find_perm_sets(labs)
#'
#' # Compares to gtools::permutations
#' len <- length(labs)
#' y <- gtools::permutations(len, len, labs, set = FALSE)
#' y <- unique(y)
#'
#' stopifnot(dim(x) == dim(y))
#'
#' # Compares to the theoretical number of entries
#' len_labs <- length(labs)
#' ns <- count_elements(labs)
#' theo_cts <- factorial(len_labs)/prod(factorial(ns))
#' stopifnot(nrow(x) == theo_cts)
#'
#' ## Others
#' x <- find_perm_sets(rep("A", 3))
#' }
find_perm_sets <- function (labs = c("A", "A", "A", "B", "B", "C")) 
{
  grps <- split_vec(labs)
  len_grps <- length(grps)
  
  ns <- count_elements(labs)
  unilabs <- names(ns)
  
  # extra counts with duplicated labels
  ns2 <- ns - 1L
  ns2 <- ns2[ns2 > 0L]
  
  # base perm without duplicated labels
  M <- gtools::permutations(len_grps, len_grps, unilabs)
  
  for (i in seq_along(ns2)) {
    ct <- ns2[[i]]
    x <- names(ns2)[i]
    
    for (j in seq_len(ct)) M <- add_one_label(M, x)
  }
  
  invisible(M)
}



#' From \link[utils]{combn}.
#' 
#' To avoid the conversion, e.g., \code{3} to \code{1:3}.
#' 
#' @param x Integer; the vector source for combinations.
#' @param m Integer; the number of elements to choose.
sim_combn <- function (x, m) 
{
  if (length(m) > 1L) 
    stop("length(m) > 1", domain = NA)
  
  if (m < 0L) 
    stop("m < 0", domain = NA)
  
  m <- as.integer(m)
  
  n <- length(x)
  if (n < m) stop("n < m", domain = NA)
  
  x0 <- x
  e <- 0
  h <- m
  a <- seq_len(m)
  
  len.r <- length(r <- x[a])
  count <- as.integer(choose(n, m))
  
  out <- vector("list", count)
  out[[1L]] <- r
  
  if (m > 0) {
    i <- 2L
    nmmp1 <- n - m + 1L
    
    while (a[1L] != nmmp1) {
      if (e < n - h) {
        h <- 1L
        e <- a[m]
        j <- 1L
      } else {
        e <- a[m - h]
        h <- h + 1L
        j <- 1L:h
      }
      
      a[m - h + j] <- e + j
      r <- x[a]
      out[[i]] <- r
      
      i <- i + 1L
    }
  }
  
  out
}


