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
                        comment_char = "") 
{
  lines <- readLines(file)
  
  # removes empty lines
  empties <- grep("^\\s*$", lines)
  
  if (length(empties)) 
    lines <- lines[-empties]
  
  rm(list = c("empties"))
  
  # removes comment lines
  if (nchar(comment_char)) 
    lines <- lines[!grepl(paste0("^", comment_char), lines)]

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
    x
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
write_fasta <- function (fasta_db, file) 
{
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
load_fasta <- function (fasta = NULL) 
{
  if (is.null(fasta)) 
    stop("FASTA file(s) are required.")
  
  oks <- file.exists(fasta)

  if (!all(oks)) {
    bads <- fasta[!oks]
    stop("Missing FASTA file(s): \n", paste(bads, collapse = "\n"))
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
load_fasta2 <- function (fasta = NULL, acc_type = NULL, acc_pattern = NULL) 
{
  if (is.null(fasta)) 
    stop("FASTA file(s) are required.")

  oks <- file.exists(fasta)
  
  if (!all(oks)) {
    bads <- fasta[!oks]
    stop("Missing FASTA file(s): \n", paste(bads, collapse = "\n"))
  }

  len_f <- length(fasta)
  len_a <- length(acc_type)
  len_p <- length(acc_pattern)

  if (len_f < len_a) 
    stop("More accession types than fasta files.")

  if (len_f < len_p) 
    stop("More acc_pattern types than fasta files.",
         call. = FALSE)

  if (len_a && (len_a < len_f)) {
    warning("More fasta files than accession types; ",
            "the first accession type will be used for all fastas.")
    acc_type <- rep(acc_type[1], len_f)
  }

  if (len_p && (len_p < len_f)) {
    warning("More fasta files than acc_pattern expressions; ",
            "the first acc_pattern expression will be used for all fastas.")
    acc_pattern <- rep(acc_pattern[1], len_f)
  }

  if (! (is.null(acc_type) || is.null(acc_pattern))) {
    acc_type <- acc_type
    acc_pattern <- acc_pattern
  } 
  else if (is.null(acc_type) && is.null(acc_pattern)) {
    acc_type <- rep("other", len_f)
    acc_pattern <- rep("(.*)", len_f)
  } 
  else if (!is.null(acc_type)) {
    acc_pattern <- purrr::map_chr(acc_type, find_acc_pattern)
  } 
  else {
    acc_type <- purrr::map_chr(acc_pattern, find_acc_type)
  }

  if (length(acc_pattern) != len_f)
    stop("Unequal length between `acc_pattern` and `fasta`.")
  
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
find_acc_pattern <- function (acc_type) 
{
  if (length(acc_type) != 1L)
    stop("The length of `acc_type` is not one.")
  
  oks <- c("uniprot_acc", "uniprot_id", "refseq_acc", "other")
  
  if (!acc_type %in% oks)
    stop("`acc_type` is not one of ", paste(oks, collapse = ", "))

  acc_pattern <- if (acc_type == "uniprot_acc") {
    "^>..\\|([^\\|]+)\\|[^\\|]+"
  } 
  else if (acc_type == "uniprot_id") {
    "^>..\\|[^\\|]+\\|([^ ]+) .*"
  } 
  else if (acc_type == "refseq_acc") {
    "^>([^ ]+?) .*"
  } 
  else if (acc_type == "other") {
    "(.*)"
  } 
  else {
    stop("Unknown `acc_type`.")
  }

  invisible(acc_pattern)
}


#' Helper for \link{load_fasta2}.
#'
#' Not used for custom acc_pattern, i.e. acc_pattern = "...".
#'
#' @inheritParams load_fasta2
find_acc_type <- function (acc_pattern) 
{
  if (length(acc_pattern) != 1L)
    stop("The length of `acc_pattern` is not one.")
  
  oks <- c("pat_upacc", "pat_upid", "pat_rsacc", "pat_other")
  
  if (!acc_pattern %in% oks)
    stop("`acc_pattern` is not one of ", paste(oks, collapse = ", "))

  pat_upacc <- "^>..\\|([^\\|]+)\\|[^ ]+?"
  pat_upid <- "^>..\\|[^\\|]+\\|([^ ]+?)"
  pat_rsacc <- "^>([^ ]+?) "
  pat_other <- "(.*)"

  acc_type <- if (acc_pattern == pat_upacc) {
    "uniprot_acc"
  } 
  else if (acc_pattern == pat_upid) {
    "uniprot_id"
  } 
  else if (acc_pattern == pat_rsacc) {
    "refseq_acc"
  } 
  else if (acc_pattern == pat_other) {
    "other"
  } 
  else {
    stop("Unknown `acc_pattern`.")
  }

  invisible(acc_type)
}

