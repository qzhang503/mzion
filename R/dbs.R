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
    stop("FASTA file(s) are required.", call. = FALSE)

  if (!all(file.exists(fasta))) 
    stop("Missing FASTA file(s): \n",
         purrr::reduce(fasta %>% .[!file.exists(.)], paste, sep = "\n"),
         call. = FALSE)

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
    stop("FASTA file(s) are required.", call. = FALSE)

  if (!all(file.exists(fasta))) 
    stop("Missing FASTA file(s): \n",
         paste(fasta %>% .[!file.exists(.)], collapse = "\n"),
         call. = FALSE)

  len_f <- length(fasta)
  len_a <- length(acc_type)
  len_p <- length(acc_pattern)

  if (len_f < len_a) 
    stop("More accession types than fasta files.",
         call. = FALSE)

  if (len_f < len_p) 
    stop("More acc_pattern types than fasta files.",
         call. = FALSE)

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
find_acc_pattern <- function (acc_type) 
{
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
find_acc_type <- function (acc_pattern) 
{
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
parse_unimod <- function (unimod = "Carbamyl (M)") 
{
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

  if (grepl("([NC]{1}-term|Anywhere) [A-Z]{1}", unimod)) 
    unimod <- 
      gsub("^(.*[NC]{1}-term|.*Anywhere)\\s*([A-Z]{1})", "\\1 = \\2", unimod)

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

  if (site == "") site = "."

  if (site %in% c("Protein N-term", "Protein C-term",
                  "Anywhere N-term", "Anywhere C-term",
                  "N-term", "C-term")) {
    pos <- site
    site <- site %>% gsub("^(Protein|Anywhere) ", "", .)
  }

  # standardize `position`
  pos <- gsub("^([NC]){1}-term", "Any \\1-term", pos)

  if (pos %in% c(".", "")) pos <- "Anywhere"

  pos_allowed <- c("Anywhere", "Protein N-term", "Protein C-term",
                   "Any N-term", "Any C-term")

  if (! pos %in% pos_allowed) 
    stop("`pos` needs to be one of ", 
         paste0("\n  '", pos_allowed, collapse = "'"), 
         "'",  
         call. = FALSE)

  # standardize terminal sites
  if (site == ".") {
    if (pos %in% c("Protein N-term", "Any N-term")) {
      site <- "N-term"
    } else if (pos %in% c("Protein C-term", "Any C-term")) {
      site <- "C-term"
    }
  }

  if (pos == "Anywhere" && site == ".") 
    stop("'position' or 'site' cannot be both 'Anywhere'.",
         call. = FALSE)

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
find_unimod <- function (unimod = "Carbamidomethyl (C)") 
{
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

    if (purrr::is_empty(idx)) 
      stop("Modification not found: '", title, "'.\n",
           "For example, use 'Acetyl' (title) instead of 'Acetylation' (full_name).",
           call. = FALSE)

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

  if (length(idx_nl)) {
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

  if (purrr::is_empty(positions_sites)) 
    stop("'", unimod, "' not found.", call. = FALSE)

  neulosses <- positions_sites[[1]][-1]
  
  neulosses <- if (length(neulosses))
    as.numeric(neulosses)
  else 
    0
  
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
#' @examples
#' \donttest{
#' ans <- table_unimods()
#' 
#' ## TMT-6, -10 and -11 plexes 
#' # share the same Unimod entry at title "TMT6plex"
#' # (the same chemistry at tag mass 229.162932 Da)
#' ans[with(ans, title == "TMT6plex"), ]
#' this_mod1 <- parse_unimod("TMT6plex (Anywhere = K)")
#' 
#' # Convenience title, "TMT10plex", alias to "TMT6plex"
#' this_mod2 <- parse_unimod("TMT10plex (Anywhere = K)")
#' 
#' # Title "TMT11plex" alias to "TMT6plex"
#' this_mod3 <- parse_unimod("TMT11plex (Anywhere = K)")
#' 
#' ## TMT-16
#' ans[with(ans, title == "TMTpro"), ]
#' this_mod1 <- parse_unimod("TMTpro (Anywhere = K)")
#' 
#' # Both "TMTpro16" and "TMT16plex" alias to "TMTpro"
#' this_mod2 <- parse_unimod("TMTpro16 (Anywhere = K)")
#' this_mod3 <- parse_unimod("TMT16plex (Anywhere = K)")
#' 
#' ## TMT-18
#' ans[with(ans, title == "TMTpro18"), ]
#' this_mod1 <- parse_unimod("TMTpro18 (Anywhere = K)")
#' this_mod2 <- parse_unimod("TMTpro18 (Anywhere = K)")
#' this_mod3 <- parse_unimod("TMT18plex (Anywhere = K)")
#' 
#' ## Summary of TMT entries and alias
#' ans[with(ans, grepl("^TMT", title)), ]
#' 
#' 
#' ## Special characters in the title (e.g., "->")
#' ans[with(ans, grepl("^Gln->pyro", title)), ]
#' this_mod <- parse_unimod("Gln->pryo-Glu (N-term = Q)")
#' 
#' }
#' @export
table_unimods <- function (file = system.file("extdata", "master.xml", 
                                              package = "proteoM"), 
                           out_nm = "~/proteoM/unimods.txt") 
{
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

