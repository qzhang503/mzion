#' Combinatorial vmods (permutation of positions but not labels).
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
#' @param nmax The maximum number of combinations.
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' aa_seq <- "MHQGVMNVNMGQKMNS"
#' aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
#' aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
#' ms <- c("M", "M", "N", "N")
#' labs <- ps <- c(2, 2)
#' names(ps) <- c("M", "N")
#' names(labs) <- c("Oxidation (M)", "Deamidated (N)")
#' M <- c("Oxidation (M)", "Oxidation (M)", "Deamidated (N)", "Deamidated (N)")
#' attr(M, "ps") <- ps
#' attr(M, "resids") <- ms
#' 
#' mzion:::find_vmodposU(M, aas)
#' mzion:::find_vmodposU(M, aas, nmax = 1L)
#' 
#' # one residue
#' aa_seq <- "MHQGVMNVNMGQKMSS"
#' aas <- .Internal(strsplit(aa_seq, "", fixed = TRUE, perl = FALSE, useBytes = FALSE))
#' aas <- .Internal(unlist(aas, recursive = FALSE, use.names = FALSE))
#' ms <- c("N")
#' labs <- ps <- c(1)
#' names(ps) <- c("N")
#' names(labs) <- c("Deamidated (N)")
#' M <- c("Deamidated (N)")
#' attr(M, "ps") <- ps
#' attr(M, "resids") <- ms
#' 
#' mzion:::find_vmodposU(M, aas)
#' mzion:::find_vmodposU(M, aas, nmax = 1L)
#' }
find_vmodposU <- function (M, aas, nmax = 64L) 
{
  rs <- attr(M, "resids")
  nr <- length(rs)
  ps <- attr(M, "ps")
  ss <- names(ps)

  if (nmax == 1L) {
    ns <- length(ss)
    combi <- character(nr)
    ends <- cumsum(ps)
    r1 <- 1L

    for (i in 1:ns) {
      si <- ss[[i]]
      pi <- ps[[i]]
      ei <- ends[[i]]
      combi[r1:ei] <- .Internal(which(aas == si))[1:pi]
      r1 <- pi + 1L
    }
    
    # one-row matrix
    combi <- .Internal(matrix(combi, nrow = 1L, ncol = 1L, byrow = FALSE, 
                              dimnames = NULL, FALSE, TRUE))
    # remove attributes
    attr(combi, "mods") <- M[1:nr]
    
    return(combi)
  }
  
  if (nr == 1L) {
    combi <- .Internal(which(aas == ss))
    combi <- .Internal(matrix(combi, nrow = 1L, ncol = 1L, byrow = FALSE, 
                              dimnames = NULL, TRUE, FALSE))
  }
  else {
    ns <- length(ss)
    X  <- vector("list", ns)
    
    for (i in 1:ns) {
      si <- .Internal(which(aas == names(ps)[i]))
      ni <- ps[[i]]
      
      X[[i]] <- if (ni == 1L) 
        matrix(si) 
      else 
        arrangements::combinations(si, ni, nitem = nmax, layout = "row")
    }
    
    combi <- expand_gr(X, nmax = nmax)
  }
  
  attr(combi, "mods") <- M[1:nr]
  
  combi
}


#' Combinatorial vmods (permutation of both positions and labels).
#'
#' Multiple Names to the same Site.
#' 
#' @param M A vector of names (lower-case vec for sites).
#' @param nmax The maximum number of combinations.
#' @inheritParams find_vmodposU
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' Vec <- c("Carbamidomethyl (M)", "Deamidated (N)", "Carbamyl (M)")
#' vec <- c("M", "N", "M")
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#' ps <- c(2, 1)
#' names(ps) <- c("M", "N")
#' attr(Vec, "ps") <- ps
#' attr(Vec, "resids") <- vec
#' rs <- unique(vec)
#' inds <- lapply(rs, function (x) which(vec == x))
#' attr(Vec, "inds") <- inds
#' mzion:::find_vmodposM(Vec, aas)
#' mzion:::find_vmodposM(Vec, aas, nmax = 1L)
#' 
#' Vec <- c("Carbamidomethyl (M)", "Carbamidomethyl (M)", "Carbamyl (M)", "Deamidated (N)")
#' vec <- c("M", "M", "M", "N")
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNS", ""))
#' ps <- c(3, 1)
#' names(ps) <- c("M", "N")
#' attr(Vec, "ps") <- ps
#' attr(Vec, "resids") <- vec
#' rs <- unique(vec)
#' inds <- lapply(rs, function (x) which(vec == x))
#' attr(Vec, "inds") <- inds
#' mzion:::find_vmodposM(Vec, aas)
#' mzion:::find_vmodposM(Vec, aas, nmax = 1L)
#' 
#' Vec <- c("Carbamidomethyl (M)", "Carbamidomethyl (M)", "Carbamyl (M)", "Deamidated (N)", 
#'          "Carbamidomethyl (S)", "Phospho (S)")
#' vec <- c("M", "M", "M", "N", "S", "S")
#' aas <- unlist(strsplit("HQGVMNVGMGQKMNSSS", ""))
#' ps <- c(3, 1, 2)
#' names(ps) <- c("M", "N", "S")
#' attr(Vec, "ps") <- ps
#' attr(Vec, "resids") <- vec
#' rs <- unique(vec)
#' inds <- lapply(rs, function (x) which(vec == x))
#' attr(Vec, "inds") <- inds
#' mzion:::find_vmodposM(Vec, aas)
#' mzion:::find_vmodposM(Vec, aas, nmax = 1L)
#' }
find_vmodposM <- function (M, aas, nmax = 64L) 
{
  rs <- attr(M, "resids", exact = TRUE)
  ps <- attr(M, "ps", exact = TRUE)
  inds <- attr(M, "inds", exact = TRUE)
  
  np <- length(ps)
  A <- P <- vector("list", np)
  
  for (i in 1:np) {
    ri <- names(ps)[i]
    si <- .Internal(which(aas == ri))
    ni <- ps[[i]]
    P[[i]] <- arrangements::combinations(si, ni, nitem = nmax, layout = "row")

    pi <- inds[[i]]
    ai <- unique(arrangements::permutations(M[pi], nitem = nmax, layout = "list"))
    A[[i]] <- do.call(rbind, ai)
  }
  
  combP <- expand_gr(P, nmax = nmax)
  combA <- expand_gr(A, nmax = nmax)
  attr(combP, "mods") <- combA

  combP
}


#' Makes the sets of MS2 labels (with permutations) at an \code{aa_masses}.
#'
#' For \code{Anywhere} variable modifications (\code{amods}) across all
#' residues.
#'
#' The universe of labels. It loops through \code{aa_masses_all}. The indexes of
#' lists are in one-to-one correspondence to \code{aa_masses_all}.
#'
#' By the design of \code{aa_masses}, the \code{amods} in an \code{aa_masses}
#' are all realized. Therefore, each list in the resulted labels should contain
#' at least one of the \code{amods} residues from the \code{aa_masses}. For
#' example, if \code{amods} contain M, N and S, each list in the result should
#' contains at least one of the residues.
#'
#' @param vec A named vector of labels. Site in name, modification in value.
#'
#' @examples
#' \donttest{
#' library(mzion)
#' 
#' ## with a bare vector
#' # One-to-one correspondence between Names and Sites
#' vec <- c(M = "Oxidation (M)", N = "Deamidated (N)")
#' attr(vec, "ps")   <- c(M = 3L, N = 2L)
#' attr(vec, "labs") <- c(`Oxidation (M)` = 3L, `Deamidated (N)` = 2L)
#'
#' ans <- mzion:::make_ms2vmods(vec)
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
#' ans <- mzion:::make_ms2vmods(vec)
#'
#' # Another one-to-multiple
#' vec <- c(M = "Oxidation (M)", M = "Carbamyl (M)",
#'          S = "Carbamidomethyl (S)", S = "Phospho (S)")
#' attr(vec, "ps") <- c(M = 2L, S = 2L)
#' attr(vec, "labs") <- c(`Oxidation (M)` = 1L,
#'                        `Carbamyl (M)` = 1L,
#'                        `Carbamidomethyl (S)` = 1L,
#'                        `Phospho (S)` = 1L)
#'
#' ans <- mzion:::make_ms2vmods(vec)
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
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
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
#' ans <- lapply(x, mzion:::make_ms2vmods) # list of 8 matrices
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
#'                                out_path = NULL)
#'
#' maxn_vmods_per_pep <- 5L
#' maxn_sites_per_vmod <- 3L
#'
#' ms1vmods_all <- lapply(aa_masses_all, mzion:::make_ms1vmod_i,
#'                        maxn_vmods_per_pep = maxn_vmods_per_pep,
#'                        maxn_sites_per_vmod = maxn_sites_per_vmod)
#'
#' # No duplication within each aa_masses
#' any_dups <- lapply(ms1vmods_all, function (x) anyDuplicated(x))
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
#' ans <- lapply(x, mzion:::make_ms2vmods) # list of 10 matrices
#'
#' ans_all <- lapply(ms1vmods_all, function (x) lapply(x, mzion:::make_ms2vmods))
#' }
make_ms2vmods <- function (vec = NULL) 
{
  n_nms  <- length(unique(names(vec)))
  n_vals <- length(unique(vec))
  oks <- n_nms == n_vals

  ans <- vec
  attr(ans, "ps") <- ps <- attr(vec, "ps")
  attr(ans, "labs") <- attr(vec, "labs")
  attr(ans, "single") <- oks

  rs <- find_ms2resids(ans, vec)
  attr(rs, "single") <- attr(rs, "ps") <- attr(rs, "labs") <- NULL
  attr(ans, "resids") <- rs

  ss <- names(ps)
  inds <- lapply(ss, function (p) which(rs == p))
  inds <- lapply(inds, unname)
  names(inds) <- names(ps)
  attr(ans, "inds") <- inds

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
  
  for (i in seq_along(M))
    M[[i]] <- vecinv[M[[i]]]

  invisible(M)
}




##
# labels: n
# positions: p
# 
# n <- c("A", "A", "B"); n1 <- 2; n2 <- 1
# p <- c(1, 3, 5, 9)
# n <- c("Carbamidomethyl (M)",  "Carbamidomethyl (M)", "Oxidation (M)")
# p <- c(1, 3, 5, 9)
# 
# choose(p, n1 + n2) * choose(n1 + n2, n2)
# = choose(4, 3) * choose(3, 2)
# 
# n <- c("A", "A", "B", "B", "C"); n1 <- 2; n2 <- 2; n3 = 1
# p <- c(1, 3, 5, 7, 8, 9)
# 
# p!/(n1!n2!n3!)/(p-n1-n2-n3)!
# = choose(p, n1+n2+n3) * (n1+n2+n3)!/(n1!n2!n3!)
# = choose(p, n1+n2+n3) * choose(n1+n2+n3, n1+n2) * choose(n1+n2, n1)
##


#' A (faster alternative) to \link[gtools]{permutations} with duplicated labels.
#'
#' @param labs A vector of labels.
#'
#' @examples
#' \dontrun{
#' library(gtools)
#' library(mzion)
#' library(dplyr)
#' 
#' # four positions for four balls in three colors
#' labs  <- c("A", "A", "B", "C")
#' len   <- length(labs)
#' g <- unique(permutations(len, len, labs, set = FALSE, repeats.allowed = FALSE))
#' m <- mzion:::find_perm_sets(labs)
#' 
#' cns <- paste0("p", 1:len)
#' g <- data.frame(g)
#' m <- data.frame(m)
#' colnames(m) <- colnames(g) <- cns
#' 
#' stopifnot(identical(dplyr::arrange_all(g), dplyr::arrange_all(m)))
#' 
#' # example-2
#' labs <- c("A", "A", "A", "B", "B", "C")
#' m <- mzion:::find_perm_sets(labs)
#'
#' # Compares to gtools::permutations
#' len <- length(labs)
#' g <- unique(permutations(len, len, labs, set = FALSE))
#'
#' cns <- paste0("p", 1:len)
#' g <- data.frame(g)
#' m <- data.frame(m)
#' colnames(m) <- colnames(g) <- cns
#' 
#' stopifnot(identical(dplyr::arrange_all(g), dplyr::arrange_all(m)))
#'
#' # Compares to the theoretical number of entries
#' len_labs <- length(labs)
#' ns <- mzion:::count_elements(labs)
#' theo_cts <- factorial(len_labs)/prod(factorial(ns))
#' stopifnot(nrow(m) == theo_cts)
#'
#' ## Others
#' x <- mzion:::find_perm_sets(rep("A", 3))
#' }
find_perm_sets <- function (labs = c("A", "A", "A", "B", "B", "C")) 
{
  grps <- split_vec(labs)
  lens <- count_elements(labs)
  
  ## base perm with unique labels
  ulabs <- names(lens)
  M  <- matrix(ulabs[1])
  
  if (length(ulabs) > 1L) {
    for (ul in ulabs[-1]) 
      M <- add_one_permlab(M, ul)
  }
  
  ## extended perm with duplicated labels
  dlens <- lens - 1L
  dlens <- dlens[dlens > 0L]
  lend  <- length(dlens)

  if (lend) {
    dlabs <- names(dlens)
    
    for (i in 1:lend) {
      ct <- dlens[[i]]
      dl <- dlabs[i]
      
      for (j in seq_len(ct)) 
        M <- add_one_label(M, dl)
    }
  }

  M
}


#' Helper in making permuation table.
#' 
#' @param M A permutation matrix with non-redundant labels.
#' @param x A new label not in M for permutation.
#' 
#' @examples
#' \dontrun{
#' library(mzion)
#' library(gtools)
#' library(dplyr)
#' 
#' unilabs <- c("A", "B", "C")
#' nu <- length(unilabs)
#' 
#' G <- gtools::permutations(nu, nu, unilabs)
#' 
#' M  <- matrix(unilabs[1])
#' for (ulab in unilabs[-1]) M <- mzion:::add_one_permlab(M, ulab)
#' 
#' cns <- paste0("p", 1:nu)
#' M <- data.frame(M)
#' G <- data.frame(G)
#' colnames(G) <- colnames(M) <- cns
#' 
#' stopifnot(identical(dplyr::arrange_all(G), dplyr::arrange_all(M)))
#' }
add_one_permlab <- function (M, x) 
{
  ncol  <- ncol(M)
  nrow  <- nrow(M)
  ncol2 <- ncol + 1L
  nrow2 <- ncol2 * nrow
  
  ins_permlab(M, x, nrow, ncol, nrow2, ncol2, rm_dup = FALSE)
}


#' Helper of \link{find_perm_sets}.
#' 
#' Adds one more labels to a permutation table.
#' 
#' @param M A permutation table with unique sets by rows.
#' @param x A duplicated label that can be found in M.
add_one_label <- function (M, x) 
{
  ncol  <- ncol(M)
  nrow  <- nrow(M)
  ncol2 <- ncol + 1L
  nrow2 <- ncol * nrow
  
  ins_permlab(M, x, nrow, ncol, nrow2, ncol2, rm_dup = TRUE)
}


#' Helper of \code{add_one_label} and \code{add_one_permlab}.
#' 
#' No duplicated labels for base permutation and \code{rm_dup = FALSE}. 
#' 
#' @inheritParams add_one_label
#' @param nrow The number of rows of M.
#' @param ncol The number of columns of M.
#' @param nrow2 The number of rows of the output.
#' @param ncol2 The number of columns of the output.
#' @param rm_dup Logical; remove duplicated rows or not. 
ins_permlab <- function (M, x, nrow, ncol, nrow2, ncol2, rm_dup = FALSE)
{
  out <- matrix(nrow = nrow2, ncol = ncol2)
  
  k <- 1L
  
  for (i in 1:nrow) {
    mi <- M[i, ]
    
    if (ncol > 1L) {
      for (j in 1:(ncol-1L)) {
        out[k, ] <- c(mi[1:j], x, mi[(j+1):ncol])
        k <- k + 1L
      }
    }
    
    out[k, ] <- c(mi, x)
    k <- k + 1L
  }
  
  if (rm_dup)
    out <- out[!duplicated.matrix(out), , drop = FALSE]
  else
    out[k:nrow2, ] <- cbind(x, M)
  
  out
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
  
  if (n < m) 
    stop("n < m", domain = NA)
  
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


