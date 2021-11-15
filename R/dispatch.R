#' Function factories for finding the amino-acid residue at a position.
#'
#' If the the value in a \code{vmods} is one of the LETTERS -> a site is
#' specified in the corresponding variable modification. Otherwise, a site is in
#' one of \code{c("C-term", "N-term")}. For example, \code{`Oxidation (M)` =
#' "M"} for the former and \code{`Protein N-term` = "N-term"} for the later. In
#' other words, \code{pos} can be a simple upper case letter or a string of
#' "\code{[NC]-term}".
#'
#' @param pos The UniMod position of a variable modification, which must be one
#'   of \code{c("Protein N-term", "Protein C-term", "Any N-term", "Any C-term",
#'   "Anywhere")}. No need to call stopifnot for the way the function is called
#'   (by the developer).
#' @return A function for finding the residue at the position specified by the
#'   argument \code{pos}. For each function, it takes a list of variable
#'   modifications specified by argument \code{vmods} and the corresponding
#'   "positions" (for vectorization) as inputs.
find_pos_site <- function (pos) 
{
  force(pos)
  
  # "^" for regex, not the indicator of an "N-term"
  p <- paste0("^", pos)
  
  function (vmods, posns) {
    idxes <- .Internal(grepl(p, posns, ignore.case = FALSE, 
                             FALSE, perl = FALSE, fixed = FALSE, 
                             useBytes = FALSE, FALSE)) & 
      vmods %in% LETTERS
    
    vmods[idxes]
  }
}


#' Finds the sites of amino-acid residues at the position of \code{Protein
#' N-term}.
#'
#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @param vmods A named list of variable modifications.
#' @inheritParams find_nmodtree
#' @seealso contain_protntany
#' @examples
#' \donttest{
#' ## `Protein N-term = P`
#' sites <- list(`Dimethyl (Protein N-term = P)` = "P",
#'               `Acetyl (Protein N-term)` = "N-term", 
#'               `Oxidation (M)` = "M",
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein N-term", "Protein N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#'
#' # (pretend `vmods` in a real workflow)
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' # stopifnot(identical(positions, posns))
#'
#' stopifnot(contain_protntsite(vmods, names(vmods), length(vmods)))
#' 
#' ans <- find_protntsite(vmods, posns)
#' stopifnot(identical(ans, vmods[1]))
#' }
find_protntsite <- find_pos_site("Protein N-term")


#' Finds the sites of amino-acid residues at the position of \code{Any N-term}.
#'
#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#'
#' @rdname  find_protntsite
#' @examples
#' \donttest{
#' ## Gln->pyro Glu (N-term = Q)
#' sites <- list(`Gln->pyro Glu (N-term = Q)` = "Q",
#'               `Acetyl (N-term)` = "N-term",
#'               `Oxidation (M)` = "M",
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#'
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#'
#' stopifnot(contain_anyntsite(vmods, posns, length(vmods)))
#'
#' ans <- find_anyntsite(vmods, posns)
#' stopifnot(identical(ans, vmods[1])) # M and N
#' }
find_anyntsite <- find_pos_site("Any N-term")


#' Finds amino-acid sites for a given set of variable modifications at the
#' position of \code{Anywhere}.
#'
#' Flowchart (5): `Oxidation (M)`
#'
#' In the less common cases of multiple modifications to the same
#' \code{Anywhere} site, e.g., "Oxidation (M)" and "Carbamyl (M)", an amino-acid
#' sequence should contain at least two "M"s. This is due to the design that
#' \code{vmods} in an \code{aa_masses} are realized. The additional subsetting
#' by the counts of residues is applied at this stage of the peptide
#' distributions. A downstream step in the calculations of MS2 ions will also
#' guard against this.
#'
#' In the case of, e.g., "Oxidation (M)" and "Protein (N-term)", and the peptide
#' sequence has only one "M" on the N-term, the entry is still kept with
#' the "additive" effect of modifications. 
#'
#' @rdname  find_protntsite
#'
#' @examples
#' \donttest{
#' ## `Oxidation (M)` and `Deamidated (N)`
#' sites <- list(`Acetyl (N-term)` = "N-term",
#'               `Oxidation (M)` = "M",
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#'
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#'
#' stopifnot(contain_anysite(vmods, posns, length(vmods)))
#'
#' ans <- find_anysite(vmods, posns)
#' stopifnot(length(ans) == 2L)
#'
#'
#' ## `Oxidation (M)` and `Carbamyl (M)`
#' sites <- list(`Dimethyl (Protein N-term = P)` = "P",
#'               `Oxidation (M)` = "M", `Carbamyl (M)` = "M",
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#'
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' vmods <- vmods[!duplicated.default(vmods)]
#'
#' posns <- names(vmods)
#' stopifnot(length(posns) < length(positions))
#'
#' ans <- find_anysite(vmods, posns)
#' stopifnot(length(ans) == 2L)
#' }
find_anysite <- find_pos_site("Anywhere")


#' Finds amino-acid sites at the position of \code{Protein C-term}.
#' 
#' Flowchart f(1-ct): `Dehydrated (Protein C-term = N)`
#' 
#' @rdname  find_protntsite
#' 
#' @examples 
#' \donttest{
#' ## `Dehydrated (Protein C-term = N)`
#' sites <- list(`Dehydrated (Protein C-term = N)` = "N", 
#'               `Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein C-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' stopifnot(contain_protctsite(vmods, posns, length(vmods)))
#' 
#' ans <- find_protctsite(vmods, posns)
#' stopifnot(identical(ans, vmods[1]))
#' }
find_protctsite <- find_pos_site("Protein C-term")


#' Finds amino-acid sites at the position of \code{Any C-term}.
#' 
#' Flowchart f(3-ct): `Oxidation (C-term = G)`
#' 
#' @rdname  find_protntsite
#' 
#' @examples 
#' \donttest{
#' ## `Oxidation (C-term = G)`
#' sites <- list(`Oxidation (C-term = G)` = "G", 
#'               `Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any C-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' stopifnot(contain_anyctsite(vmods, posns, length(vmods)))
#' 
#' ans <- find_anyctsite(vmods, posns)
#' stopifnot(identical(ans, vmods[1]))
#' }
find_anyctsite <- find_pos_site("Any C-term")


#' Function factories for checking the existence of an amino-acid residue at a
#' position.
#'
#' @inheritParams find_pos_site
#' @return A function for checking the existence of a residue at the position
#'   specified by the argument \code{pos}. For each function, it takes a list of
#'   variable modifications specified by argument \code{vmods}, and the
#'   corresponding postitions and counts as inputs.
contain_pos_site <- function (pos) 
{
  force(pos)
  
  # "^" for regex, not the indicator of an "N-term"
  p <- paste0("^", pos)

  function (vmods, posns, len) {
    if (!len) 
      return(FALSE)
    if (len == 1L && vmods == "") 
      return(FALSE)
    
    idxes <- .Internal(grepl(p, posns, ignore.case = FALSE, 
                             FALSE, perl = FALSE, fixed = FALSE, 
                             useBytes = FALSE, FALSE)) & 
      vmods %in% LETTERS

    any(idxes)
  }
}


#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @rdname find_protntsite
contain_protntsite <- contain_pos_site("Protein N-term")


#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#' 
#' @rdname  find_protntsite
contain_anyntsite <- contain_pos_site("Any N-term")


#' Flowchart (5): `Oxidation (M)`
#' 
#' @rdname find_protntsite
contain_anysite <- contain_pos_site("Anywhere")


#' Flowchart (1-ct): 'Dehydrated (Protein C-term = N)'
#' 
#' @rdname find_protntsite
contain_protctsite <- contain_pos_site("Protein C-term")


#' Flowchart (3-ct): 'Oxidation (C-term = G)'
#' 
#' @rdname find_protntsite
contain_anyctsite <- contain_pos_site("Any C-term")


#' Function factories for checking the existence of an amino-acid residue at a
#' \code{terminal} position.
#'
#' Note to the developer: when manufacturing, a \code{pos} must be terminal in
#' one of \code{c("Protein N-term", "Any N-term", "Protein C-term", "Any
#' C-term")}. \code{"Anywhere"} is not one of them.
#'
#' @param pos The position. It must be terminal in \code{c("Protein N-term",
#'   "Any N-term", "Protein C-term", "Any C-term")}.
#' @return A function for checking the existence of a residue at the
#'   \code{terminal} position specified by the argument \code{pos}. For each
#'   function, it takes a list of variable modifications specified by argument
#'   \code{vmods} as inputs.
contain_termpos_any <- function (pos) 
{
  force(pos)
  
  # "^" for regex, not the indicator of an "N-term"
  p <- paste0("^", pos)

  function (vmods, posns, len) {
    if (!len) 
      return(FALSE)
    if (len == 1L && vmods == "") 
      return(FALSE)
    
    idxes <- .Internal(grepl(p, posns, ignore.case = FALSE, 
                             FALSE, perl = FALSE, fixed = FALSE, 
                             useBytes = FALSE, FALSE))

    any(idxes)
  }
}


#' Checks if variable modifications are on \code{Protein N-term}.
#'
#' Flowchart (2-nt): `Acetyl (Protein N-term)`
#'
#' @inheritParams find_nmodtree
#'
#' @examples
#' \donttest{
#' ## `Acetyl (Protein N-term)`
#' sites <- list(`Acetyl (Protein N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' contain_protntany(vmods, posns, length(vmods))
#' }
contain_protntany <- contain_termpos_any("Protein N-term")


#' Checks if variable modifications are on \code{Any N-term}.
#' 
#' Flowchart (4-nt): `Acetyl (N-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' \donttest{
#' ## `Acetyl (N-term)`
#' sites <- list(`Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' contain_anyntany(vmods, posns, length(vmods))
#' }
contain_anyntany <- contain_termpos_any("Any N-term")


#' Checks if variable modifications are on \code{Protein C-term}.
#' 
#' Flowchart (2-ct): `Amidated (Protein C-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' \donttest{
#' ## `Amidated (Protein C-term)`
#' sites <- list(`Amidated (Protein C-term)` = "C-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein C-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' contain_protctany(vmods, posns, length(vmods))
#' }
contain_protctany <- contain_termpos_any("Protein C-term")


#' Checks if variable modifications are on \code{C-term}.
#' 
#' Flowchart (4-ct): `Amidated (C-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' \donttest{
#' ## `Amidated (C-term)`
#' sites <- list(`Amidated (C-term)` = "C-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any C-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' contain_anyctany(vmods, posns, length(vmods))
#' }
contain_anyctany <- contain_termpos_any("Any C-term")


#' Subsets proteins by variable modifications.
#'
#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @param site A site. No need to specify "any terminal" sites.
#' @inheritParams find_nmodtree
#' @examples
#' \donttest{
#' ## Protein N-term (Site)
#' sites <- list(`Dimethyl (Protein N-term = P)` = "P", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' prps <- list(PROT_A = c("-MAKEMASSPECFUN", "-PAKEKASSPECFUN"), 
#'              PROT_B = c("PAKEKASSPECFUN", "NKAKEKASSPECFU", 
#'                         "-NKAKEKASSPECFU"))
#' 
#' subset_protntsite(prps, find_protntsite(vmods, posns))
#' } 
subset_protntsite <- function (prps, site) 
{
  s <- paste0("^-", site)
  lapply(prps, function (ps) ps[grepl(s, ps)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (2-nt): `Acetyl (Protein N-term)`
#' 
#' @rdname subset_protntsite
subset_protntany <- function (prps) lapply(prps, function (ps) ps[grepl("^-", ps)])


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#' @rdname subset_protntsite
#' @examples 
#' \donttest{
#' ## Any N-term (Site)
#' sites <- list(`Dimethyl (N-term = P)` = "P", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' posns <- names(vmods)
#' 
#' prps <- list(PROT_A = c("-MAKEMASSPECFUN", "-PAKEKASSPECFUN"), 
#'              PROT_B = c("PAKEKASSPECFUN", "NKAKEKASSPECFU", 
#'                         "-NKAKEKASSPECFU"))
#' 
#' # "-PAKEKASSPECFUN" went with `subset_protntsite` (see flow charts)
#' subset_anyntsite(prps, find_anyntsite(vmods, posns))
#' }
subset_anyntsite <- function (prps, site = "Q") 
{
  s <- paste0("^", site)
  lapply(prps, function (ps) ps[grepl(s, ps)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart f(4-nt): `Acetyl (N-term)`
#' @rdname subset_protntsite
subset_anyntany <- function (peps) peps


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (5): `Oxidation (M)`
#' @rdname subset_protntsite
#' @inheritParams find_nmodtree
#' @examples 
#' \donttest{
#' ## Anywhere
#' min_n_res <- c(P = 1, M = 1, N = 1)
#' 
#' sites <- list(`Dimethyl (N-term = P)` = "P", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' vmods <- unname(vmods)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' 
#' prps <- list(PROT_A = c("-MAKEMASSPECFUN", "-PAKEKASSPECFUN"), 
#'              PROT_B = c("PAKEKASSPECFUN", "NKAKEKASSPECFU", 
#'                         "-NKAKEKASSPECFU"))
#' 
#' # should contain both M and N
#' subset_anysite(prps, sites, min_n_res)
#' 
#' ## Multiple mods to the same site
#' # (mimic from aa_masses)
#' min_n_res <- c(P = 1, M = 2, N = 1)
#' 
#' sites <- list(`Dimethyl (N-term = P)` = "P", 
#'               `Oxidation (M)` = "M", `Carbamyl (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' 
#' # No need of modification names
#' vmods <- unname(vmods)
#' sites <- unname(sites)
#' vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
#' 
#' is_same <- any(length(min_n_res) > 1L)
#' if (is_same) {
#'   ok <- !duplicated.default(vmods)
#'   vmods <- vmods[ok]
#'   sites <- sites[ok]
#' }
#' 
#' prps <- list(PROT_A = c("-MAKEMASSPECFUN", "-PAKEKASSPECFUN"), 
#'              PROT_B = c("PAKEKASSPECFUN", "NKAKEKASSPECFU", 
#'                         "-NKAKEKASSPECFU"))
#' 
#' ans <- subset_anysite(prps, sites, min_n_res)
#' }
subset_anysite <- function (prps, sites, min_n_res) 
{
  lapply(prps, function (ps) {
    oks <- lapply(ps, function (p) {
      a <- .Call(stringi:::C_stri_count_fixed, 
                 str = p, pattern = sites, opts_fixed = NULL)
      all(a >= min_n_res)
    })
    oks <- .Internal(unlist(oks, recursive = FALSE, use.names = FALSE))

    ps[oks]
  })
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (1-ct): `Dehydrated (Protein C-term = N)`
#' @rdname subset_protntsite
subset_protctsite <- function (prps, site = "V") 
{
  s <- paste0(site, "-$")
  lapply(prps, function (ps) ps[grepl(s, ps)])
}


# f(2-ct): 
subset_protctany <- function (prps) lapply(prps, function (ps) ps[grepl("-$", ps)])


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (3-ct): `Oxidation (C-term = G)`
#' @rdname subset_protntsite
subset_anyctsite <- function (prps, site = "Q") 
{
  s <- paste0(site, "$")
  lapply(prps, function (ps) ps[grepl(s, ps)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (4-ct): `Amidated (C-term)`
#' @rdname subset_protntsite
subset_anyctany <- function (peps) peps


#' Find and subset peptides.
#'
#' (1-nt) Protein N-term + site -> (2-nt) Any protein N-term -> (3-nt) Any
#' N-term + site -> (4-nt) Any N-term -> (5) Anywhere + site.
#'
#' @param vmods A named list of variable modifications. See also
#'   \link{find_protntsite} for examples of \code{vmods}.
#' @param min_n_res The minimum numbers of residues for a given set of sites.
#' @param posns The position (e.g., \code{Protein N-term}, \code{Anywhere},
#'   etc.) of \code{vmods}. The argument can be obtained from \code{vmods} but
#'   passed as a parameter for vectorization.
#' @param len The count of \code{vmods} (passed as a parameter for
#'   vectorization).
#' @inheritParams distri_peps
find_nmodtree <- function (prps, min_n_res, vmods, posns, len) 
{
  if (contain_protntsite(vmods, posns, len)) { # level_1: Protein N-term + Site
    prps <- subset_protntsite(prps, find_protntsite(vmods, posns))
    
    if (contain_anysite(vmods, posns, len)) {
      # (1) -|* .. |
      prps <- subset_anysite(prps, find_anysite(vmods, posns), min_n_res)
    } else {
      # (2) -|*    |
      prps <- prps
    }
  } else {
    if (contain_protntany(vmods, posns, len)) { # level_2: Protein N-term 
      prps <- subset_protntany(prps)
      
      if (contain_anysite(vmods, posns, len)) {
        # (3) -|o .. |
        prps <- subset_anysite(prps, find_anysite(vmods, posns), min_n_res)
      } else {
        # (4) -|o    |
        prps <- prps
      }
    } else { 
      if (contain_anyntsite(vmods, posns, len)) { # level_3: Any N-term + Site
        prps <- subset_anyntsite(prps, find_anyntsite(vmods, posns))
        
        if (contain_anysite(vmods, posns, len)) {
          # (5) |* .. |
          prps <- subset_anysite(prps, find_anysite(vmods, posns), min_n_res)
        } else {
          # (6) |*    |
          prps <- prps
        }
      } else { 
        if (contain_anyntany(vmods, posns, len)) { # level_4: Any N-term 
          prps <- prps 
          
          if (contain_anysite(vmods, posns, len)) {
            # (7) |o .. |
            prps <- subset_anysite(prps, find_anysite(vmods, posns), min_n_res)
          } else {
            # (8) |o    |
            prps <- prps
          }
        } else { 
          if (contain_anysite(vmods, posns, len)) { # level_5: Anywhere
            # (9) |  .. |
            prps <- subset_anysite(prps, find_anysite(vmods, posns), min_n_res)
          } else {
            # (10) |     |
            prps <- prps
          }
        }
      }
    }
  }
}


#' Find and subset peptides.
#'
#' (1-ct) Protein C-term + site -> (2-ct) Any protein C-term -> (3-ct) Any
#' C-term + site -> (4-ct) Any C-term.
#'
#' @inheritParams find_nmodtree
find_cmodtree <- function (prps, min_n_res, vmods, posns, len) 
{
  if (contain_protctsite(vmods, posns, len)) { # level_1: Protein C-term + Site
    # (1) -|* .. *|-, (2) -|*    *|-, (3) -|o .. *|-, (4) -|o    *|-, (5) |* .. *|-, 
    # (6) |*    *|-, (7) |o .. *|-, (8) |o    *|-, (9) |  .. *|-, (10) |     *|-
    prps <- subset_protctsite(prps, find_protctsite(vmods, posns))
  } else {
    if (contain_protctany(vmods, posns, len)) { # level_2: Protein C-term 
      # (1) -|* .. o|-, (2) -|*    o|-, (3) -|o .. o|-, # (4) -|o    o|-, # (5) |* .. o|-, 
      # (6) |*    o|-, (7) |o .. o|-, (8) |o    o|-, (9) |  .. o|-, (10) |     o|-
      prps <- subset_protctany(prps)
    } else {
      if (contain_anyctsite(vmods, posns, len)) { # level_3: Any C-term + Site
        # (1) -|* .. *|, (2) -|*    *|, (3) -|o .. *|, # (4) -|o    *|, # (5) |* .. *|, 
        # (6) |*    *|, (7) |o .. *|, (8) |o    *|, (9) |  .. *|, (10) |     *|
        prps <- subset_anyctsite(prps, find_anyctsite(vmods, posns))
      } else {
        if (contain_anyctany(vmods, posns, len)) { # level_4: Any C-term
          # (1) -|* .. o|, (2) -|*    o|, (3) -|o .. o|, # (4) -|o    o|, # (5) |* .. o|, 
          # (6) |*    o|, (7) |o .. o|, (8) |o    o|, (9) |  .. o|, (10) |     o|
          prps <- prps 
        } # else {
          # No use
        # }
      }
    }
  }
  
  invisible(prps)
}


#' Subsets peptides by variable modifications.
#'
#' From N-term to C-term.
#'
#' @inheritParams add_fixvar_masses
#' @inheritParams distri_peps
subpeps_by_vmods <- function(aa_masses, prps) 
{
  vmods <- attr(aa_masses, "vmods_ps", exact = TRUE) 
  min_n_res <- attr(aa_masses, "min_n_res", exact = TRUE)
  is_same <- attr(aa_masses, "is_same", exact = TRUE)
  
  if (is.list(vmods)) {
    vmods <- unname(vmods)
    vmods <- unlist(vmods, recursive = FALSE, use.names = TRUE)
  }
  
  if (is_same) 
    vmods <- vmods[!duplicated.default(vmods)]

  posns <- names(vmods)
  len <- length(vmods)

  prps <- find_nmodtree(prps, min_n_res, vmods, posns, len)
  prps <- find_cmodtree(prps, min_n_res, vmods, posns, len)
}

