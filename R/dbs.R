### Also in proteoQ

#' Reads a file in fasta format
#'
#' Reads a file in fasta format by line.
#'
#' @param file A character string to the name of a protein fasta file.
#' @param acc_pattern A regular expression describing the pattern to separate
#'   the header lines of fasta entries. The default is to separate a header and
#'   keep the character string before the first space where the so kept will be
#'   used as the name of an entry. The character ">" at the beginning of the
#'   header will also be removed.
#' @param comment_char Character: a character or an empty string. Use "" to turn
#'   off the interpretation of comment lines.
#' @examples
#' \donttest{
#' # assume the file and location of "uniprot_hs_2020_05.fasta"
#' fasta <- read_fasta("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' head(names(fasta))
#'
#' # use the first fifty characters
#' fasta <- read_fasta("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                     ">(.{50}).*")
#' head(names(fasta))
#'
#' # uniprot_acc
#' fasta <- read_fasta("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                     ">..\\|([^\\|]+)\\|.*")
#' head(names(fasta))
#'
#' # use all characters in the header
#' fasta <- read_fasta("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                     ">(.*)")
#' head(names(fasta))
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @seealso \code{\link{write_fasta}}
#' @export
read_fasta <- function (file = NULL, acc_pattern = ">([^ ]+?) .*", 
                        comment_char = "") {
  
  lines <- readLines(file)
  
  # removes empty lines
  empties <- grep("^\\s*$", lines)
  
  if (length(empties)) {
    lines <- lines[-empties]
  }
  
  rm(list = c("empties"))
  
  # removes comment lines
  if (nchar(comment_char)) {
    lines <- lines[!grepl(paste0("^", comment_char), lines)]
  }

  # begins and ends
  headers <- grep(">", lines)
  begins <- headers + 1L
  ends <- c(headers[-1L] - 1L, length(lines))

  seqs <- mapply(function (x, y) {
    Reduce(paste0, lines[x : y])
  }, begins, ends, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  hdrs <- lines[headers]

  db <- mapply(function (x, y) {
    attr(x, "header") <- y
    return(x)
  }, seqs, hdrs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  names(db) <- gsub(acc_pattern, "\\1", hdrs)
  
  invisible(db)
}


#' Writes fasta
#'
#' Writes a fasta file (Not yet used).
#'
#' @param fasta_db A list of protein entries from \code{\link{read_fasta}}.
#' @inheritParams read_fasta
#' @examples
#' \donttest{
#' fasta_db <- read_fasta(file = "~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' write_fasta(fasta_db, "~/proteoM/examples/my.fasta")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
write_fasta <- function (fasta_db, file) {
  
  filepath <- gsub("(^.*/).*$", "\\1", file)
  dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
  
  res <- lapply(fasta_db, function (x) paste(attr(x, "header"), x, sep = "\n"))
  res <- unlist(res)
  writeLines(res, file)
}


#' Loads fasta
#' 
#' Not used in proteoM.
#'
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. There is no default and the experimenters need to
#'   supply the files.
#' @examples
#' \donttest{
#' fasta_db <- load_fasta("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' }
load_fasta <- function (fasta = NULL) {
  
  if (is.null(fasta)) {
    stop("FASTA file(s) are required.", call. = FALSE)
  }

  if (!all(file.exists(fasta))) {
    stop("Missing FASTA file(s): \n",
         purrr::reduce(fasta %>% .[!file.exists(.)], paste, sep = "\n"),
         call. = FALSE)
  }

  lapply(fasta, function (x) read_fasta(x)) %>%
    do.call(`c`, .) %>%
    `names<-`(gsub(">", "", names(.))) %>%
    .[!duplicated(names(.))]
}

### End of also in proteoQ




#' Loads fasta (with parsing rule).
#'
#' The length of \code{acc_type} needs to match the length of \code{fasta};
#' otherwise, the first value will be used for all \code{fasta} files.
#' 
#' @param acc_type Character string(s); the types of protein accessions in one
#'   of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
#'   the corresponding regular expressions need to be supplied via argument
#'   \code{acc_pattern}.
#' @param acc_pattern Regular expression(s) describing the patterns to separate
#'   the header lines of fasta entries. At the \code{NULL} default, the pattern
#'   will be automated when \code{acc_type} are among c("uniprot_acc",
#'   "uniprot_id", "refseq_acc", "other").
#' @inheritParams matchMS
#' @examples
#' \donttest{
#' fasta_db <- load_fasta2(
#'               c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                 "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'               c("uniprot_acc", "other")
#' )
#'
#' # Need `acc_pattern` as "crap" is not one of the default acc_type
#' load_fasta2(
#'    c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'      "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'    c("uniprot_acc", "crap")
#' )
#'
#' # ok
#' fasta_db2 <- load_fasta2(
#'                c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                  "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'                c("uniprot_acc", "crap"),
#'                c("^>..\\|([^\\|]+)\\|[^\\|]+", "(.*)")
#' )
#'
#' fasta_db3 <- load_fasta2(
#'                c("~/proteoM/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                  "~/proteoM/dbs/fasta/crap/crap.fasta"),
#'                c("my_acc", "crap"),
#'                c("^>..\\|([^\\|]+)\\|[^\\|]+", "(.*)")
#' )
#'
#' stopifnot(identical(fasta_db, fasta_db2),
#'           identical(fasta_db, fasta_db3))
#' }
#' @export
load_fasta2 <- function (fasta = NULL, acc_type = NULL, acc_pattern = NULL) {
  
  if (is.null(fasta)) {
    stop("FASTA file(s) are required.", call. = FALSE)
  }

  if (!all(file.exists(fasta))) {
    stop("Missing FASTA file(s): \n",
         paste(fasta %>% .[!file.exists(.)], collapse = "\n"),
         call. = FALSE)
  }

  len_f <- length(fasta)
  len_a <- length(acc_type)
  len_p <- length(acc_pattern)

  if (len_f < len_a) {
    stop("More accession types than fasta files.",
         call. = FALSE)
  }

  if (len_f < len_p) {
    stop("More acc_pattern types than fasta files.",
         call. = FALSE)
  }

  if (len_a && len_a < len_f) {
    warning("More fasta files than accession types; ",
            "the first accession type will be used for all fastas.",
            call. = FALSE)
    acc_type <- rep(acc_type[1], len_f)
  }

  if (len_p && len_p < len_f) {
    warning("More fasta files than acc_pattern expressions; ",
            "the first acc_pattern expression will be used for all fastas.",
            call. = FALSE)
    acc_pattern <- rep(acc_pattern[1], len_f)
  }

  if (! (is.null(acc_type) || is.null(acc_pattern))) {
    acc_type <- acc_type
    acc_pattern <- acc_pattern
  } else if (is.null(acc_type) && is.null(acc_pattern)) {
    acc_type <- rep("other", len_f)
    acc_pattern <- rep("(.*)", len_f)
  } else if (!is.null(acc_type)) {
    acc_pattern <- map_chr(acc_type, find_acc_pattern)
  } else {
    acc_type <- map_chr(acc_pattern, find_acc_type)
  }

  stopifnot(length(acc_pattern) == len_f)
  
  # Not to USE.NAMES; otherwise fasta names prefix to accession names
  # this is different to map2 where names are NULL for each fasta_db
  
  mapply(function (x, y) read_fasta(x, y), fasta, acc_pattern, 
         SIMPLIFY = FALSE, USE.NAMES = FALSE) %>%
    do.call(`c`, .) %>%
    `names<-`(gsub(">", "", names(.))) %>%
    .[!duplicated(names(.))]
}


#' Helper for \link{load_fasta2}.
#'
#' Not used for custom acc_type, i.e. acc_type = "my_acctype".
#'
#' @inheritParams load_fasta2
find_acc_pattern <- function (acc_type) {
  
  stopifnot(length(acc_type) == 1L)
  stopifnot(acc_type %in% c("uniprot_acc", "uniprot_id", "refseq_acc", "other"))

  if (acc_type == "uniprot_acc") {
    acc_pattern <- "^>..\\|([^\\|]+)\\|[^\\|]+"
  } else if (acc_type == "uniprot_id") {
    acc_pattern <- "^>..\\|[^\\|]+\\|([^ ]+) .*"
  } else if (acc_type == "refseq_acc") {
    acc_pattern <- "^>([^ ]+?) .*"
  } else if (acc_type == "other") {
    acc_type <- "other"
    acc_pattern <- "(.*)"
  } else {
    stop("Unknown `acc_type`.",
         call. = FALSE)
  }

  invisible(acc_pattern)
}


#' Helper for \link{load_fasta2}.
#'
#' Not used for custom acc_pattern, i.e. acc_pattern = "...".
#'
#' @inheritParams load_fasta2
find_acc_type <- function (acc_pattern) {
  
  stopifnot(length(acc_pattern) == 1L)

  pat_upacc <- "^>..\\|([^\\|]+)\\|[^ ]+?"
  pat_upid <- "^>..\\|[^\\|]+\\|([^ ]+?)"
  pat_rsacc <- "^>([^ ]+?) "
  pat_other <- "(.*)"

  stopifnot(acc_pattern %in% c("pat_upacc", "pat_upid", "pat_rsacc", "pat_other"))

  if (acc_pattern == pat_upacc) {
    acc_type <- "uniprot_acc"
  } else if (acc_pattern == pat_upid) {
    acc_type <- "uniprot_id"
  } else if (acc_pattern == pat_rsacc) {
    acc_type <- "refseq_acc"
  } else if (acc_pattern == pat_other) {
    acc_type <- "other"
  } else {
    stop("Unknown `acc_pattern`.",
         call. = FALSE)
  }

  invisible(acc_type)
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
#' @inheritParams matchMS
#' @inheritParams add_fixvar_masses
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
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
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
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
                           digits = 5L) {
  
  # (6) "amods- tmod- vnl- fnl+"
  if (!length(amods)) return(NULL)
  
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
                            digits = 5L) {

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

  if (len_p <= 0) return(list())
  
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
    x <- extract_vmodsets(.ms1_vmodsets, .base_ent, len_p, ns)
  }
  
  invisible(x)
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
#' find_unique_sets(5, c("Oxidation (M)", "Carbamidomethyl (M)", "Carbamyl (M)"))
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


###

#' Extracts sets of combinatorial vmods labels.
#'
#' For faster calculations of precursor masses.
#'
#' @param p The number of open positions for filling.
#' @param labs The labels of balls.
#' @param .ms1_vmodsets The complete set of combinatorial variable modifications
#'   for MS1.
#' @param .base_ent The base entries of .ms1_vmodsets. 
#' @seealso make_ms1_vmodsets
extract_vmodsets <- function (.ms1_vmodsets, .base_ent, p = 4L, 
                              labs = c("Carbamyl (M)")) {
  
  n <- length(labs)
  
  idx <- lapply(.base_ent, is_equal_sets, labs)
  idx <- which(unlist(idx, recursive = FALSE, use.names = FALSE))
  
  vmodsets <- .ms1_vmodsets[[idx]]
  vmodsets[names(vmodsets) %in% n:p]
}


#' Makes the sets of labels of variable modifications.
#' 
#' No position permutation (for MS1 masses).
#' 
#' @param aa_masses_all All the amino acid lookup tables.
#' @inheritParams matchMS
#' @examples 
#' \donttest{
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)")
#' 
#' varmods <- c("Acetyl (Protein N-term)", "Oxidation (M)", "Carbamidomethyl (M)", 
#'             "Deamidated (N)", "Carbamyl (M)", 
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
#' .ms1_vmodsets <- make_ms1_vmodsets(aa_masses_all)
#' .base_ent <- lapply(.ms1_vmodsets, `[[`, 1)
#' 
#' extract_vmodsets(.ms1_vmodsets, .base_ent, p = 4, 
#'                  labs = c("Carbamidomethyl (M)", "Carbamyl (M)"))
#' }
make_ms1_vmodsets <- function (aa_masses_all = NULL, maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L) {
  
  # stopifnot(maxn_vmods_per_pep >= 2L)
  
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
  
  out <- lapply(resmods_all, bacth_vmods_combi, maxn_vmods_per_pep, 
                maxn_sites_per_vmod)
  
  resids <- lapply(resmods_all, `[`, 1)
  resids <- unlist(resids, recursive = FALSE, use.names = FALSE)
  names(out) <- resids
  
  # assign(".ms1_vmodsets", out, envir = .GlobalEnv)
  # .base_ent <- lapply(.ms1_vmodsets, `[[`, 1)
  # assign(".base_ent", .base_ent, envir = .GlobalEnv)

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
  
  make_unique_sets(p = maxn_vmods_per_pep, 
                   n = length(resmods), 
                   labs = names(resmods), 
                   maxn_vmods_per_pep = maxn_vmods_per_pep, 
                   maxn_sites_per_vmod = maxn_sites_per_vmod)
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
  # for (i in seq_along(combs)) attr(combs[[i]], "p") <- nms_p[i]
  names(combs) <- nms_p

  invisible(combs)
}



###




#' Finds the combinations across residues.
#'
#' For uses with MS1 precursors. For multiple residues (each residue one to
#' multiple modifications).
#'
#' @param intra_combis The results from \link{unique_mvmods}.
#' @seealso \link{find_intercombi_p2} for MS2 ions.
#' @examples
#' C <- list(c("Carbamidomethyl (C)"),
#'           rep("Carbamidomethyl (C)", 2))
#'
#' N <- list(c("Deamidated (N)"),
#'           rep("Deamidated (N)", 2))
#'
#' intra_combis <- list(C = C, N = N)
#'
#' ans <- find_intercombi(intra_combis)
find_intercombi <- function (intra_combis) {
  
  len <- length(intra_combis)

  if (!len) { # scalar
    v_out <- list()
  } else if (any(.Internal(unlist(lapply(intra_combis, purrr::is_empty), 
                                  recursive = FALSE, use.names = FALSE)))) { # list
    v_out <- list()
  } else if (len > 1L) {
    v_out <- expand_grid_rows(intra_combis, use.names = FALSE)
  } else {
    v_out <- purrr::flatten(intra_combis)
  }

  invisible(v_out)
}


#' Concatenates adjacent peptides in a list (with mass).
#' 
#' @param peps A list of peptide sequences with a one-letter representation of
#'   amino acid residues.
#' @param n The number of mis-cleavages for consideration.
#' @param include_cts Logical; the list, \code{peps}, includes the protein
#'   C-terminal sequence or not. At the default of TRUE, mis-cleaved peptides at
#'   the end of the protein C-terms will be added as they should. The arguments
#'   would be typically at FALSE, for example, when used for generating
#'   mis-cleaved peptides from the N-terminal of peptides with the removal of a
#'   starting residue \code{M}.
#' @examples
#' \donttest{
#' peps <- 1:26
#' names(peps) <- LETTERS
#' res <- roll_sum(peps, 2)
#' 
#' # length shorter than n
#' peps <- c(a = 1)
#' res <- roll_sum(peps, 2)
#' 
#' peps <- c(a = 1, b = 2, c = 3)
#' res <- roll_sum(peps, 4)
#' }
roll_sum <- function (peps = NULL, n = 2L, include_cts = TRUE) {
  
  len <- length(peps)
  
  if (!len) return(NULL)
  
  if (n >= len) n <- len - 1L
  
  res <- lapply(seq_len((len - n)), function (x) {
    ranges <- x:(x + n)
    
    psub <- peps[ranges]
    nms <- accumulate_char(names(psub), paste0)

    vals <- cumsum(psub)
    names(vals) <- nms
    
    vals
  }) 
  
  res <- .Internal(unlist(res, recursive = FALSE, use.names = TRUE))
  
  if (include_cts && n >= 1L) {
    ends <- peps[(len - n + 1L):len]
    
    res2 <- lapply(n:1L, function (x) {
      y <- tail(ends, x)
      nms <- accumulate_char(names(y), paste0)
      
      vals <- cumsum(y)
      names(vals)  <- nms
      
      vals
    })
    
    res2 <- .Internal(unlist(res2, recursive = FALSE, use.names = TRUE))
  } else {
    res2 <- NULL
  }
  
  c(res, res2)
}


#' Parse the name of a Unimod
#'
#' The general format: \code{parse_unimod("title (position = site)")}.
#'
#' @param unimod The name of a \href{https://www.unimod.org/}{Unimod} modification.
#' @seealso \link{table_unimods}, \link{find_unimod}.
#' @examples
#' \donttest{
#' # "dot" for anywhere (either position or site)
#' x1 <- parse_unimod("Carbamidomethyl (. = C)")
#' x2 <- parse_unimod("Carbamidomethyl (Anywhere = C)")
#' x3 <- parse_unimod("Carbamidomethyl (C)")
#'
#' identical(x1, x2); identical(x2, x3)
#'
#' # Any residue on protein N-term
#' x1 <- parse_unimod("Acetyl (Protein N-term = .)")
#' x2 <- parse_unimod("Acetyl (Protein N-term)")
#'
#' identical(x1, x2)
#'
#' # Any N-term residue
#' x1 <- parse_unimod("Acetyl (N-term)")
#' x2 <- parse_unimod("Acetyl (N-term = .)")
#'
#' identical(x1, x2)
#'
#' # N-term Q
#' x1 <- parse_unimod("Gln->pryo-Glu (N-term = Q)")
#' x2 <- parse_unimod("Gln->pryo-Glu (N-term Q)")
#'
#' identical(x1, x2)
#'
#' # ok with parenthesis in the 'title'
#' x <- parse_unimod("Hex(5)HexNAc(2) (N)")
#'
#' # ok with spaces in the 'title'
#' x <- parse_unimod("Met-loss (Protein N-term = M)")
#' }
#'
#' \dontrun{
#' # No modification of anywhere and anything
#' # (to every position and site)
#' x <- parse_unimod("Carbamidomethyl (Anywhere = .)")
#' x <- parse_unimod("Carbamidomethyl")
#' x <- parse_unimod("Carbamidomethyl (. = .)")
#'
#' # Prefer an "=" sign between 'N-term' and 'Q'
#' x <- parse_unimod("Gln->pyro-Glu (N-term Q)")
#' }
#' @export
parse_unimod <- function (unimod = "Carbamyl (M)") {
  
  # unimod = "Carbamidomethyl (Protein N-term = C)" # --> pos_site = "Protein N-term = C"
  # unimod = "Carbamidomethyl (Any N-term = C)" # --> pos_site = "Any N-term = C"
  # unimod = "Carbamidomethyl (N-term = C)" # --> pos_site = "N-term = C"
  # unimod = "Carbamidomethyl (. = C)" # --> pos_site = ". = C"
  # unimod = "Carbamidomethyl (C)" # --> pos_site = "C"
  # unimod = "Carbamidomethyl ()" # --> pos_site = ""
  # unimod = "Carbamidomethyl" # --> pos_site = ""
  # unimod = "" # --> pos_site = ""

  ## any N-term residue
  # unimod = "Carbamidomethyl (Protein N-term = .)"

  # unimod = "Hex(5)HexNAc(2) (N)"

  ## dual parentheses
  # unimod = "Carbamidomethyl ((. = C))" # --> pos_site = ". = C"

  if (grepl("([NC]{1}-term|Anywhere) [A-Z]{1}", unimod)) {
    unimod <- unimod %>%
      gsub("^(.*[NC]{1}-term|.*Anywhere)\\s*([A-Z]{1})", "\\1 = \\2", .)
  }

  # (assumed) no space in `title`
  # title <- gsub("(.*)\\s\\([^\\(]*\\)$", "\\1", unimod)
  title <- gsub("^([^ ]+?) .*", "\\1", unimod)

  pos_site <- unimod %>%
    gsub("^[^ ]+", "", .) %>%
    gsub("^[^\\(]+[\\(]*([^\\)]*)[\\)]*$", "\\1", .)

  if (grepl("=", pos_site)) {
    pos <- pos_site %>%
      gsub("^([^=]+?)[=].*", "\\1", .) %>%
      gsub("^[ ]*", "\\1", .) %>%
      gsub(" *$", "", .)

    site <- pos_site %>%
      gsub("^[^=]+?[=](.*)", "\\1", .) %>%
      gsub("^[ ]*", "\\1", .) %>%
      gsub(" *$", "", .)
  } else {
    pos <- "."
    site <- pos_site
  }

  if (site == "") {
    site = "."
  }

  if (site %in% c("Protein N-term", "Protein C-term",
                  "Anywhere N-term", "Anywhere C-term",
                  "N-term", "C-term")) {
    pos <- site
    site <- site %>% gsub("^(Protein|Anywhere) ", "", .)
  }

  # standardize `position`
  pos <- pos %>%
    gsub("^([NC]){1}-term", "Any \\1-term", .)

  if (pos %in% c(".", "")) {
    pos <- "Anywhere"
  }

  pos_allowed <- c("Anywhere", "Protein N-term", "Protein C-term",
                   "Any N-term", "Any C-term")

  if (! pos %in% pos_allowed) {
    stop("`pos` needs to be one of ", 
         paste0("\n  '", pos_allowed, collapse = "'"), 
         "'",  
         call. = FALSE)
  }

  # standardize terminal sites
  if (site == ".") {
    if (pos %in% c("Protein N-term", "Any N-term")) {
      site <- "N-term"
    } else if (pos %in% c("Protein C-term", "Any C-term")) {
      site <- "C-term"
    }
  }

  if (pos == "Anywhere" && site == ".") {
    stop("'position' or 'site' cannot be both 'Anywhere'.",
         call. = FALSE)
  }

  invisible(list(title = title, position = pos, site = site))
}


#' Find a Unimod.
#'
#' Find the mono-isotopic mass, position, site and neutral losses of a
#' modification.
#'
#' In the field of \code{position_site}, \code{position} is the name and
#' \code{site} is the value.
#'
#' @inheritParams parse_unimod
#' @seealso \link{table_unimods}, \link{parse_unimod}.
#' @examples
#' \donttest{
#' x <- find_unimod("Carbamidomethyl (C)")
#' x <- find_unimod("Carbamidomethyl (M)")
#' x <- find_unimod("Acetyl (Protein N-term)")
#' x <- find_unimod("Gln->pyro-Glu (N-term = Q)")
#' x <- find_unimod("Hex(5)HexNAc(2) (N)")
#' }
#'
#' \dontrun{
#' # Prefer an "=" sign between 'N-term' and 'Q'
#' x <- find_unimod("Gln->pyro-Glu (N-term Q)")
#' }
#' @export
find_unimod <- function (unimod = "Carbamidomethyl (C)") {
  
  options(digits = 9L)

  res <- parse_unimod(unimod)
  title <- res$title
  position <- res$position
  site <- res$site
  rm(res)

  parent <- system.file("extdata", "master.xml", package = "proteoM") %>%
    xml2::read_xml()

  # <umod:elements>
  # <umod:modifications>
  # <umod:amino_acids>
  # <umod:mod_bricks>

  children <- xml2::xml_children(parent)
  contents <- parent %>% xml2::xml_contents()

  elements <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "elements")]])
  modifications <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "modifications")]])
  amino_acids <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "amino_acids")]])
  mod_bricks <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "mod_bricks")]])

  this_mod <- local({
    idx <- which(xml2::xml_attr(modifications, "title") == title)

    if (purrr::is_empty(idx)) {
      stop("Modification not found: '", title, "'.\n",
           "For example, use 'Acetyl' (title) instead of 'Acetylation' (full_name).",
           call. = FALSE)
    }

    this_mod <- modifications[[idx]]
  })

  stopifnot(xml2::xml_attrs(this_mod) %>% .["title"] == title)

  # --- children of `this_mod` ---
  modch <- xml2::xml_children(this_mod)

  # --- find mass ---
  monomass <- xml2::xml_attrs(modch) %>%
    purrr::map(`[`, "mono_mass") %>%
    `[`(!is.na(.)) %>%
    unlist() %>%
    as.numeric()

  stopifnot(length(monomass) == 1)

  # --- find sites and positions ---
  positions_sites <- xml2::xml_attrs(modch) %>%
    purrr::map(`[`, c("site", "position")) %>%
    purrr::map(~ invisible(setNames(.x[1], .x[2])))

  sites <- xml2::xml_attrs(modch) %>%
    purrr::map(`[`, c("site")) %>%
    unlist()

  positions <- xml2::xml_attrs(modch) %>%
    purrr::map(`[`, c("position")) %>%
    unlist()

  # --- neutral loss ---
  idx_nl <- grep("NeutralLoss", modch)

  if (length(idx_nl) > 0L) {
    nls <- purrr::map(idx_nl, ~ {
      modnl <- xml2::xml_children(modch[.x])

      xml2::xml_attrs(modnl) %>%
        purrr::map(`[`, "mono_mass") %>%
        `[`(!is.na(.)) %>%
        as.numeric() %>%
        sort() %>% # ensures the first is 0
        setNames(paste("nl", 1:length(.), sep = "."))
    }) %>%
      setNames(rep("nl", length(.)))

    positions_sites[idx_nl] <-
      purrr::map2(idx_nl, nls, ~ {c(positions_sites[[.x]], .y)})
  }

  positions_sites <- positions_sites[sites == site & positions == position] %>%
    .[purrr::map_lgl(., function (x) !is.null(x))]

  if (purrr::is_empty(positions_sites)) {
    stop("'", unimod, "' not found.", call. = FALSE)
  }

  neulosses <- positions_sites[[1]][-1]
  if (purrr::is_empty(neulosses)) {
    neulosses <- 0
  } else {
    neulosses <- as.numeric(neulosses)
  }

  invisible(list(monomass = monomass,
                 position_site = positions_sites[[1]][1],
                 nl = neulosses))
}


#' Parses \href{https://www.unimod.org/}{Unimod} entries.
#'
#' For convenience findings of the \code{title}, \code{site} and
#' \code{position}.
#'
#' @param file A file path to a Unimod ".xml".
#' @param out_nm A name to outputs.
#' @seealso \link{find_unimod}, \link{parse_unimod}.
#' @export
#' @examples
#' \donttest{
#' ans <- table_unimods()
#' 
#' ans[with(ans, title == "TMTpro"), ]
#' this_mod <- parse_unimod("TMTpro (Anywhere = K)")
#' 
#' ans[with(ans, grepl("^Gln->pyro", title)), ]
#' this_mod <- parse_unimod("Gln->pryo-Glu (N-term = Q)")
#' }
table_unimods <- function (file = system.file("extdata", "master.xml", 
                                              package = "proteoM"), 
                           out_nm = "~/proteoM/unimods.txt") {
  
  parent <- xml2::read_xml(file)
  
  # <umod:elements>
  # <umod:modifications>
  # <umod:amino_acids>
  # <umod:mod_bricks>
  
  children <- xml2::xml_children(parent)
  contents <- xml2::xml_contents(parent)
  
  elements <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "elements")]])
  modifications <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "modifications")]])
  amino_acids <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "amino_acids")]])
  mod_bricks <-
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "mod_bricks")]])
  
  titles <- xml2::xml_attr(modifications, "title")
  
  lapply(titles, function (title) {
    idx <- which(xml2::xml_attr(modifications, "title") == title)
    this_mod <- modifications[[idx]]
    
    # --- children of `this_mod` ---
    modch <- xml2::xml_children(this_mod)
    mod_attrs <-  xml2::xml_attrs(modch)
    
    sites <- lapply(mod_attrs, `[`, c("site"))
    sites <- unlist(sites)
    
    positions <- lapply(mod_attrs, `[`, c("position"))
    positions <- unlist(positions)
    
    rows_s <- unlist(lapply(sites, is.na))
    rows_p <- unlist(lapply(positions, is.na))
    rows <- rows_s | rows_p
    
    data.frame(title = title, 
               position = positions[!rows], 
               site = sites[!rows])
  }) %>% 
    dplyr::bind_rows() %>% 
    `rownames<-`(NULL) %T>% 
    readr::write_tsv(out_nm)
}


#' Find mis-cleavages in a vector.
#'
#' A convenience utility may be used to extract the first \eqn{n+1} peptides
#' from 0 to n mis-cleavages. It also assumes that the data were already sorted
#' in a desirable way.
#'
#' @param x A vector of data.
#' @param n Integer. The number of mis-cleavages.
keep_n_misses <- function (x, n) {
  
  len <- length(x)
  
  if (n < 0L) {
    stop("`n` cannot be nagative integers: ", n)
  }
  
  if (!len) {
    stop("Length of `x` cannot be zero.")
  }

  x[1:min(n + 1, len)]
}


#' Exclude mis-cleavages in a vector.
#'
#' @inheritParams keep_n_misses
#' @seealso keep_n_misses
exclude_n_misses <- function (x, n) {
  
  len <- length(x)

  if (n < 0L) {
    stop("`n` cannot be nagative integers: ", n)
  }
  
  if (!len) {
    stop("Length of `x` cannot be zero.")
  }
  
  x[-(1:min(n + 1, len))]
}


#' Excludes a character in string counting
#'
#' @param x A character string
#' @param char A character to be excluded for counting.
#' @importFrom stringi stri_length stri_count_fixed
str_exclude_count <- function (x, char = "-") {
  stringi::stri_length(x) - stringi::stri_count_fixed(x, char)
}


#' Remove a starting character from the first \code{n} entries.
#'
#' @param x A list of character strings.
#' @param char A starting character to be removed.
#' @param n The number of beginning entries to be considered.
rm_char_in_nfirst2 <- function (x, char = "^-", n = (max_miss + 1L) * 2L) {
  
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- seq_len(n)
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms
  
  x
}


#' Remove a trailing character from the last \code{n} entries.
#'
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst2
rm_char_in_nlast2 <- function (x, char = "-$", n = (max_miss + 1L) * 2L) {
  
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- (len - n + 1L):len
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms
  
  x
}


