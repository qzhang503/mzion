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


#' Finds a Unimod.
#'
#' Finds the mono-isotopic mass, position, site and neutral losses of a
#' modification.
#'
#' In the field of \code{position_site}, \code{position} is the name and
#' \code{site} is the value.
#'
#' @param xml_files Name(s) of Unimod ".xml" files. The file path is a system
#'   setting of \code{system.file("extdata", xml_file, package = "proteoM")}.
#' @inheritParams parse_unimod
#' @seealso \link{table_unimods}, \link{parse_unimod}.
#' @examples
#' \donttest{
#' x1 <- find_unimod("Carbamidomethyl (C)")
#' x2 <- find_unimod("Carbamidomethyl (M)")
#' x3 <- find_unimod("Acetyl (Protein N-term)")
#' x4 <- find_unimod("Gln->pyro-Glu (N-term = Q)")
#' x5 <- find_unimod("Hex(5)HexNAc(2) (N)")
#' }
#'
#' \dontrun{
#' # Prefer an "=" sign between 'N-term' and 'Q'
#' x <- find_unimod("Gln->pyro-Glu (N-term Q)")
#' }
#' @export
find_unimod <- function (unimod = "Carbamidomethyl (C)", 
                         xml_files = c("master.xml", "custom.xml")) 
{
  options(digits = 9L)
  
  res <- parse_unimod(unimod)
  title <- res$title
  position <- res$position
  site <- res$site
  rm(list = c("res"))
  
  this_mod <- hfind_unimod(xml_files = xml_files, unimod)
  
  # mass
  node_delta <- xml2::xml_find_all(this_mod, "umod:delta")
  monomass <- as.numeric(xml2::xml_attr(node_delta, "mono_mass"))
  
  if (!(length(monomass) == 1L))
    stop("The length of `mono_mass` is not one.", call. = FALSE)
  
  # sites and positions
  node_specs <- xml2::xml_find_all(this_mod, "umod:specificity")
  attrs_specs <- xml2::xml_attrs(node_specs)
  
  sites <- attrs_specs %>%
    lapply(`[`, c("site")) %>%
    unlist()
  
  positions <- attrs_specs %>%
    lapply(`[`, c("position")) %>%
    unlist()
  
  names(sites) <- positions
  
  idx_ps <- which(sites == site & positions == position)
  
  sites <- sites[idx_ps] %>%
    .[purrr::map_lgl(., function (x) !is.null(x))]
  
  local({
    len_sites <- length(sites)
    
    if (len_sites > 1L)
      stop("Multiple matches in site and position for '", unimod, "'.", 
           call. = FALSE)
    else if (len_sites == 0L)
      stop("No matches in site and position for '", unimod, "'.", 
           call. = FALSE)
  })

  # neutral loss
  idxes_nl <- grep("NeutralLoss", node_specs)
  idxes_nl <- idxes_nl[idxes_nl == idx_ps]
  
  if (length(idxes_nl)) {
    nodes_nl <- xml2::xml_children(node_specs[idxes_nl])
    
    neulosses <- xml2::xml_attr(nodes_nl, "mono_mass") %>%
      as.numeric() %>%
      sort() # ensures the first is 0
  }
  else 
    neulosses <- 0
  
  invisible(list(title = title, 
                 monomass = monomass,
                 position_site = sites,
                 nl = neulosses))
}


#' Helper of \link{find_unimod}.
#' 
#' @param xml_files A list of xml file names.
#' @inheritParams find_unimod
hfind_unimod <- function (xml_files = c("master.xml", "custom.xml"), unimod) 
{
  for (xml_file in xml_files) {
    title <- gsub("^([^ ]+?) .*", "\\1", unimod)
    
    xml_root <- xml2::read_xml(system.file("extdata", xml_file, package = "proteoM"))
    nodes_lev1_four <- xml2::xml_children(xml_root)
    node_modif <- xml2::xml_find_all(nodes_lev1_four, "//umod:modifications")
    modifications <- xml2::xml_children(node_modif)
    
    idx <- which(xml2::xml_attr(modifications, "title") == title)
    
    if (length(idx)) {
      this_mod <- modifications[[idx]]
      break
    } else
      this_mod <- NULL
  }
  
  if (is.null(this_mod))
    stop("Modification not found: '", title, "'.\n",
         "For example, use 'Acetyl' (title) instead of 'Acetylation' (full_name).",
         call. = FALSE)
  
  invisible(this_mod)
}


#' Tabulates \href{https://www.unimod.org/}{Unimod} entries.
#'
#' For convenience summary of the \code{title}, \code{site} and
#' \code{position}.
#'
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
table_unimods <- function (out_nm = "~/proteoM/unimods.txt") 
{
  files <- c("master.xml", "custom.xml")
  
  lapply(files, htable_unimods) %>% 
    do.call(rbind, .) %>% 
    `rownames<-`(NULL) %T>% 
    readr::write_tsv(out_nm)
  
}


#' Helper of \link{table_unimods}.
#' 
#' @param file A file path to a Unimod ".xml".
htable_unimods <- function (file) 
{
  xml_root <- xml2::read_xml(system.file("extdata", file, package = "proteoM"))
  nodes_lev1_four <- xml2::xml_children(xml_root)
  node_modif <- xml2::xml_find_all(nodes_lev1_four, "//umod:modifications")
  modifications <- xml2::xml_children(node_modif)
  
  titles <- xml2::xml_attr(modifications, "title")

  nodes_specs <- lapply(modifications, function (this_mod) 
    xml2::xml_find_all(this_mod, "umod:specificity"))

  sites <- lapply(nodes_specs, xml2::xml_attr, "site")
  positions <- lapply(nodes_specs, xml2::xml_attr, "position")
  
  nodes_delta <- lapply(modifications, function (this_mod) 
    xml2::xml_find_all(this_mod, "umod:delta"))
  
  mono_masses <- lapply(nodes_delta, xml2::xml_attr, "mono_mass")
  
  lens_specs <- lapply(nodes_specs, length)
  titles <- mapply(function (x, y) rep(x, y), titles, lens_specs)
  mono_masses <- mapply(function (x, y) rep(x, y), mono_masses, lens_specs)
  
  if (length(titles))
    data.frame(title = unname(unlist(titles)), 
               site = unlist(sites), 
               position = unlist(positions), 
               mono_masses = unlist(mono_masses))
  else
    NULL
}


#' Adds or modifies a Unimod entry.
#'
#' Changes are made to system files.
#'
#' @param header The header of a \href{https://www.unimod.org/}{Unimod}. The
#'   fields of \code{title} and \code{full_name} are required. No spaces are
#'   allowed in \code{title}.
#' @param specificity The specificity of \code{site} and \code{position} of a
#'   modification. The value of \code{site} is an upper-case, one-letter
#'   representation of an amino-acid residue, or "N-term" or "C-term". The value
#'   of \code{position} is one of "Anywhere", "Protein N-term", "Protein
#'   C-term", "Any N-term" or "Any C-term".
#'
#'   See also \link{remove_unimod} for a wildcard approach of \code{site = "."}
#'   and \code{position = "."} to remove all sites and positions under a title.
#' @param delta The mass and composition difference of a modification:
#'   \code{mono_mass}, the difference in mono-isotopic mass; \code{avge_mass},
#'   the difference in average mass; \code{composition}, the difference in
#'   chemical composition.
#'
#'   The precision of \code{mono_mass} is typically set with a decimal place of
#'   \eqn{\ge 5}.
#' @param neuloss The mass and composition \emph{loss} of a neutral species atop
#'   of the \code{delta}. Like many other search engines, \code{proteoM} adapts
#'   the common convention where the losses of masses and compositions are
#'   expressed in \emph{positive} forms for both the masses and the
#'   compositions.
#'
#'   See also \link{remove_unimod} for a wildcard approach of \code{neuloss =
#'   c(mono_mass = ".", ...)} to remove all neutral losses under a site and a
#'   position for a given title.
#' @import xml2
#' @seealso \link{table_unimods}, \link{parse_unimod}, \link{find_unimod}.
#' @return An xml object.
#' @examples
#' \dontrun{
#' library(proteoM)
#'
#' # To avoid unsound chemistries, proteoM prohibits
#' #   additive modifications to the same site.
#' # To enable cumulative effects, the solution is to
#' #   devise a "merged" modification.
#'
#' (system.file("extdata", "master.xml", package = "proteoM"))
#' (system.file("extdata", "custom.xml", package = "proteoM"))
#'
#' # Additive N-terminal modifications (though not chemically sound)
#' masses <- calc_unimod_compmass("H(17) C(8) 13C(4) 15N O(2)")
#' mono_mass <- masses$mono_mass
#' avge_mass <- masses$avge_mass
#' 
#' x <- add_unimod(header      = c(title       = "TMT10plexNterm+Gln->pyro-Glu",
#'                                 full_name   = "Additive N-term TMT10plex and Gln->pyro-Glu"),
#'
#'                 # think about this: site = "Q", not "N-term"
#'                 specificity = c(site        = "Q",
#'                                 position    = "Any N-term"),
#'                 delta       = c(mono_mass   = "212.136383",
#'                                 avge_mass   = "212.2329",
#'                                 composition = "H(17) C(8) 13C(4) 15N O(2)"),
#'                 neuloss     = c(mono_mass   = "0",
#'                                 avge_mass   = "0",
#'                                 composition = "0"))
#'
#' (ans <- table_unimods())
#' (ans[with(ans, title == "TMT10plexNterm+Gln->pyro-Glu"), ])
#'
#' # site `C`: Oxiation + Carbamidomethyl
#' # (without neutral losses)
#' x <- add_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'                                 full_name   = "Oxidation and iodoacetamide derivative"),
#'                 specificity = c(site        = "C",
#'                                 position    = "Anywhere"),
#'                 delta       = c(mono_mass   = "73.016379",
#'                                 avge_mass   = "73.0507",
#'                                 composition = "H(3) C(2) N O(2)"),
#'                 neuloss     = c(mono_mass   = "0",
#'                                 avge_mass   = "0",
#'                                 composition = "0"))
#'
#' # site `M`: Oxiation + Carbamidomethyl
#' # (with neutral losses)
#' x <- add_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'                                 full_name   = "Oxidation and iodoacetamide derivative"),
#'                 specificity = c(site        = "M",
#'                                 position    = "Anywhere"),
#'                 delta       = c(mono_mass   = "73.016379",
#'                                 avge_mass   = "73.0507",
#'                                 composition = "H(3) C(2) N O(2)"),
#'                 neuloss     = c(mono_mass   = "63.998285",
#'                                 avge_mass   = "64.1069",
#'                                 composition = "H(4) C O S"))
#' }
#' 
#' ## Heavy isotopes
#' # Lysine
#' K8 <- calc_unimod_compmass("13C(6) C(-6) 15N(2) N(-2)")
#' 
#' mono_mass <- K8$mono_mass
#' avge_mass <- K8$avge_mass
#' 
#' x <- add_unimod(header      = c(title       = "K8",
#'                                 full_name   = "Heavy lysine 13C(6) 15N(2)"),
#'                 specificity = c(site        = "K",
#'                                 position    = "Anywhere"),
#'                 delta       = c(mono_mass   = "8.0142",
#'                                 avge_mass   = "7.94272",
#'                                 composition = "13C(6) C(-6) 15N(2) N(-2)"),
#'                 neuloss     = c(mono_mass   = "0",
#'                                 avge_mass   = "0",
#'                                 composition = "0"))
#' 
#' # Arginine
#' R10 <- calc_unimod_compmass("13C(6) C(-6) 15N(4) N(-4)")
#' 
#' mono_mass <- R10$mono_mass
#' avge_mass <- R10$avge_mass
#' 
#' x <- add_unimod(header      = c(title       = "R10",
#'                                 full_name   = "Heavy arginine 13C(6) 15N(4)"),
#'                 specificity = c(site        = "R",
#'                                 position    = "Anywhere"),
#'                 delta       = c(mono_mass   = "10.00827",
#'                                 avge_mass   = "9.92954",
#'                                 composition = "13C(6) C(-6) 15N(4) N(-4)"),
#'                 neuloss     = c(mono_mass   = "0",
#'                                 avge_mass   = "0",
#'                                 composition = "0"))
#' 
#' 
#' # TMT reporter ions
#' electron <- 0.000549
#' 
#' tmt11_126 <- calc_unimod_compmass("C(8) N(1) H(16)")
#' tmt11_126 <- lapply(tmt11_126, `-`, electron)
#' 
#' tmt11_131n <- calc_unimod_compmass("13C(4) C(4) 15N(1) H(16)")
#' tmt11_131n <- lapply(tmt11_131n, `-`, electron)
#' 
#' tmt11_131c <- calc_unimod_compmass("13C(5) C(3) N(1) H(16)")
#' tmt11_131c <- lapply(tmt11_131c, `-`, electron)
#' 
#' @export
add_unimod <- function (header = c(title = "Foo", full_name = "Foo bar"), 
                        specificity = c(site  = "C", position = "Anywhere"), 
                        delta = c(mono_mass = "42.010565", 
                                  avge_mass = "42.0367", 
                                  composition = "H(2) C(2) O"), 
                        neuloss = c(mono_mass = "0", 
                                    avge_mass = "0", 
                                    composition = "0")) 
{
  options(digits = 9L)
  
  title <- header[["title"]]
  full_name <- header[["full_name"]]
  position <- specificity[["position"]]
  site <- specificity[["site"]]
  mod_mono_mass <- delta[["mono_mass"]]
  mod_avge_mass <- delta[["avge_mass"]]
  mod_composition <- delta[["composition"]]
  neuloss_mono_mass <- neuloss[["mono_mass"]]
  neuloss_avge_mass <- neuloss[["avge_mass"]]
  neuloss_composition <- neuloss[["composition"]]

  xml_file <- system.file("extdata", "custom.xml", package = "proteoM")
  xml_root <- xml2::read_xml(xml_file)
  nodes_lev1_four <- xml2::xml_children(xml_root)
  node_modif <- xml2::xml_find_all(nodes_lev1_four, "//umod:modifications")

  if (!length(node_modif)) {
    xml2::xml_add_child(xml_root, "umod:modifications")
    node_modif <- xml2::xml_find_all(xml_root, "//umod:modifications")
  }

  local({
    if (!(length(site) == 1L))
      stop("The length of `site` is not one.", call. = FALSE)
    
    if (!(length(position) == 1L))
      stop("The length of `position` is not one.", call. = FALSE)
    
    ok_sites <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
                  "P", "Q", "R", "S", "T", "U", "V", "W", "Y", 
                  "N-term", "C-term")
    
    ok_positions <- c("Anywhere", "Protein N-term", "Protein C-term", 
                      "Any N-term", "Any C-term")
    
    if (!site %in% c(LETTERS, "N-term", "C-term"))
      stop("`site` is not one of ", paste(ok_sites, collapse = ", "), ".", 
           call. = FALSE)

    if (!position %in% ok_positions)
      stop("`position` is not one of ", paste(ok_positions, collapse = ", "), ".", 
           call. = FALSE)
  })
  
  if (is.numeric(mod_mono_mass)) 
    mod_mono_mass <- as.character(mod_mono_mass)
  
  if (is.numeric(mod_avge_mass)) 
    mod_avge_mass <- as.character(mod_avge_mass)
  
  if (is.numeric(neuloss_mono_mass)) 
    neuloss_mono_mass <- as.character(neuloss_mono_mass)
  
  if (is.numeric(neuloss_avge_mass)) 
    neuloss_avge_mass <- as.character(neuloss_avge_mass)

  # individual nodes of modifications
  modifications <- xml2::xml_children(node_modif)
  idx_title <- which(xml2::xml_attr(modifications, "title") == title)
  len_title <- length(idx_title)
  
  if (!len_title) {
    # `title` not found -> new modifiation
    rec_ids <- as.integer(xml2::xml_attr(modifications, "record_id"))
    
    if (!length(rec_ids))
      rec_ids <- 0L
    
    record_id <- as.character(max(rec_ids) + 1L)
    
    add_modification(node_modif, title = title, full_name = full_name, 
                     site  = site, position = position, 
                     mod_mono_mass = mod_mono_mass, mod_avge_mass = mod_avge_mass, 
                     mod_composition = mod_composition, 
                     neuloss_mono_mass = neuloss_mono_mass, 
                     neuloss_avge_mass = neuloss_avge_mass, 
                     neuloss_composition = neuloss_composition, 
                     record_id = record_id)
  }
  else if (len_title == 1L) {
    # `title` found -> checks `site` and `position`
    this_mod <- modifications[[idx_title]]
    nodes_children <- xml2::xml_children(this_mod)
    attrs_children <- xml2::xml_attrs(nodes_children)
    
    local({
      this_full_name <- xml2::xml_attr(this_mod, "full_name")
      
      if (this_full_name != full_name)
        warning("The original `full_name = ", this_full_name, "` \n", 
                "replaced by `full_name = ", full_name, "`.", 
                call. = FALSE)
    })

    local({
      this_delta <- xml2::xml_find_all(this_mod, "umod:delta")

      stopifnot(length(this_delta) == 1L)
      
      this_mono_mass <- xml2::xml_attr(this_delta, "mono_mass")
      this_avge_mass <- xml2::xml_attr(this_delta, "avge_mass")
      this_composition <- xml2::xml_attr(this_delta, "composition")
      
      if (this_mono_mass != mod_mono_mass)
        stop("The original `mono_mass = ", mod_mono_mass, "` ", 
             "is different to the current `mono_mass = ", this_mono_mass, "`.\n", 
             "If you believe the new entry is correct: \n", 
             "  try `remove_unimod(title = ", title, "`, ", 
             "then repeat the current call.", 
             call. = FALSE)

      if (this_avge_mass != mod_avge_mass)
        stop("The original `avge_mass = ", mod_avge_mass, "` ", 
             "is different to the current `avge_mass = ", this_avge_mass, "`.\n", 
             "If you believe the new entry is correct: \n", 
             "  try `remove_unimod(title = ", title, "`, ", 
             "then repeat the current call.", 
             call. = FALSE)
      
      if (this_composition != mod_composition)
        stop("The original `composition = ", mod_composition, "` ", 
             "is different to the current `composition = ", this_composition, "`.\n", 
             "If you believe the new entry is correct: \n", 
             "  try `remove_unimod(title = ", title, "`, ", 
             "then repeat the current call.", 
             call. = FALSE)
    })
    
    sites <- unlist(lapply(attrs_children, `[`, c("site")))
    positions <- unlist(lapply(attrs_children, `[`, c("position")))
    ok_sitepos <- which((sites == site) & (positions == position))
    
    len_sitepos <- length(ok_sitepos)
    
    if (len_sitepos > 1L)
      stop("Multiple matches to `site = ", site, "` and ", 
           "`position = ", position, "` at ", 
           "`title = ", title, "`.\n", 
           "Fix the redundancy from ", xml_file, ".", 
           call. = FALSE)
    else if (len_sitepos) 
      # `site` and `position` found -> adds/checks neulosses
      add_neuloss(nodes_children[[ok_sitepos]], 
                  neuloss_mono_mass = neuloss_mono_mass, 
                  neuloss_avge_mass = neuloss_avge_mass, 
                  neuloss_composition = neuloss_composition)
    else 
      # new `site` and `position` with optional neuloss
      add_specificy(this_mod, site = site, position = position, 
                    neuloss_mono_mass = neuloss_mono_mass, 
                    neuloss_avge_mass = neuloss_avge_mass, 
                    neuloss_composition = neuloss_composition)
  }
  else if (len_title > 1L)
    stop("Multiple matches to `", title, "`.\n", 
         "Fix the redundancy from ", xml_file, ".", 
         call. = FALSE)
  

  xml2::write_xml(xml_root, xml_file)
  
  invisible(xml_root)
}


#' Adds a new node of modification.
#'
#' @param node A node of top-level \code{modification} (one of the top-4).
#' @param title The title of a modification.
#' @param full_name The full-name description of a modification.
#' @param position The position of a modification.
#' @param site The site of a modification.
#' @param mod_mono_mass The mono-isotopic mass delta of a modification.
#' @param mod_avge_mass The average mass delta of a modification.
#' @param mod_composition The chemical composition of a modification.
#' @param neuloss_mono_mass The mono-isotopic mass delta of neutral loss in a
#'   positive value.
#' @param neuloss_avge_mass The average mass delta of neutral loss in a positive
#'   value.
#' @param neuloss_composition The chemical composition of a neutral loss.
#' @param record_id The index of a modification.
add_modification <- function (node = NULL, title = "Acetyl", full_name = "", 
                              site  = "N-term", position = "Protein N-term", 
                              mod_mono_mass = "42.010565", mod_avge_mass = "42.0367", 
                              mod_composition = "H(2) C(2) O", 
                              neuloss_mono_mass = "0", neuloss_avge_mass = "0", 
                              neuloss_composition = "0", record_id = "0") 
{
  # (0) node: <umod:modifications>
  #   (1) <umod:mod title=... >
  #     (2.1) specificity
  #       (2.1.1) neuloss
  #         (2.1.1.1) element
  #     (2.2) delta
  #       (2.2.1) element

  ## (1) adds a "umod:mod" node
  xml2::xml_add_child(node, "umod:mod")
  modifications <- xml2::xml_find_all(node, "umod:mod")
  this_mod <- modifications[[length(modifications)]]

  xml2::xml_set_attrs(
    this_mod, 
    c(title = title, 
      full_name = full_name, 
      username_of_poster = "unimod", 
      group_of_poster = "user", 
      date_time_posted = as.character(Sys.time()), 
      date_time_modified = as.character(Sys.time()), 
      approved = "0", 
      record_id = record_id))
  
  # (2.1) specificity
  add_specificy(node = this_mod, site = site, position = position, 
                neuloss_mono_mass = neuloss_mono_mass, 
                neuloss_avge_mass = neuloss_avge_mass, 
                neuloss_composition = neuloss_composition) 

  # (2.2) delta
  add_delta(node = this_mod, mod_mono_mass = mod_mono_mass, 
            mod_avge_mass = mod_avge_mass, 
            mod_composition = mod_composition)
  
  invisible(node)
}


#' Adds a node of \code{specificity}.
#' 
#' @param node A node of \code{umod:mod}.
#' @inheritParams add_modification
add_specificy <- function (node = NULL, site = "C", position = "Anywhere", 
                           neuloss_mono_mass = "0", neuloss_avge_mass = "0", 
                           neuloss_composition = "0") 
{
  # (1) node: <umod:mod title=... >
  #   (2.1) specificity
  #     (2.1.1) neuloss
  #       (2.1.1.1) element
  #   (2.2) delta
  #     (2.2.1) element
  
  # (2.1) specificity
  xml2::xml_add_child(node, "umod:specificity")
  
  # DON't: "node_specs <- xml2::xml_children(node)"
  #   because of siblings other than `specificity`
  # 
  # [10] <umod:specificity hidden="1" site="U" position="Anywhere"
  # [11] <umod:specificity hidden="1" site="M" position="Anywhere" 
  # [12] <umod:delta mono_mass=
  # [13] <umod:alt_name>
  # [14] <umod:xref>\n 
  
  node_specs <- xml2::xml_find_all(node, "umod:specificity")

  len <- length(node_specs)
  this_spec <- node_specs[[len]]
  
  xml2::xml_set_attrs(
    this_spec, 
    c(hidden = "0", 
      site  = site, 
      position = position, 
      classification = "other", 
      spec_group = as.character(len)))
  
  # (2.1.1) neuloss
  add_neuloss(node = this_spec, neuloss_mono_mass = neuloss_mono_mass, 
              neuloss_avge_mass = neuloss_avge_mass, 
              neuloss_composition = neuloss_composition)
}


#' Adds a node of \code{delta}.
#' 
#' @param node A node of \code{umod:mod}.
#' @inheritParams add_modification
add_delta <- function (node = NULL, mod_mono_mass = "0", mod_avge_mass = "0", 
                       mod_composition = "0") 
{
  # (1) node: <umod:mod title=... >
  #   (2.1) specificity
  #     (2.1.1) neuloss
  #       (2.1.1.1) element
  #   (2.2) delta
  #     (2.2.1) element
  
  # (2.2) delta
  xml2::xml_add_child(node, "umod:delta")
  node_delta <- xml2::xml_find_all(node, "umod:delta")
  
  len <- length(node_delta)
  this_delta <- node_delta[[len]]
  
  xml2::xml_set_attrs(
    this_delta, 
    c(mono_mass = mod_mono_mass, 
      avge_mass  = mod_avge_mass, 
      composition = mod_composition))
  
  # (2.2.1) adds element children
  add_comp_elements(node_delta, mod_composition)
  
  invisible(node)
}


#' Adds a node of \code{NeutralLoss}.
#' 
#' Including the \code{0} NeutralLoss.
#' 
#' @param node A node of \code{umod:specificity}.
#' @inheritParams add_modification
add_neuloss <- function (node = NULL, neuloss_mono_mass = "0", 
                         neuloss_avge_mass = "0", neuloss_composition = "0") 
{
  #   (2.1) node: specificity
  #     (2.1.1) neuloss
  #       (2.1.1.1) element
  #   (2.2) delta
  #     (2.2.1) element
  
  if (neuloss_mono_mass != "0") {
    nodes_neuloss <- xml2::xml_find_all(node, "umod:NeutralLoss")
    attrs_neuloss <- xml2::xml_attrs(nodes_neuloss)
    masses <- unlist(lapply(attrs_neuloss, `[`, c("mono_mass")))
    
    if (any(masses == neuloss_mono_mass))
      warning("Pre-existed `NeutralLoss`; do nothing.")
    else {
      if (!any(masses == "0")) {
        this_neuloss <- hadd_neuloss(node = node, 
                                     neuloss_mono_mass = "0", 
                                     neuloss_avge_mass  = "0", 
                                     neuloss_composition = "0")
      }

      this_neuloss <- hadd_neuloss(node = node, 
                                   neuloss_mono_mass = neuloss_mono_mass, 
                                   neuloss_avge_mass  = neuloss_avge_mass, 
                                   neuloss_composition = neuloss_composition)
      
      # (2.1.1.1) neuloss elements
      add_comp_elements(this_neuloss, neuloss_composition)
    }
  }
  else {
    message("No neutral loss.")
  }
}


#' Helper of \link{add_neuloss}.
#' 
#' @inheritParams add_neuloss
hadd_neuloss <- function (node = NULL, neuloss_mono_mass = "0", 
                          neuloss_avge_mass = "0", neuloss_composition = "0") 
{
  xml2::xml_add_child(node, "umod:NeutralLoss")
  nodes_neuloss <- xml2::xml_find_all(node, "umod:NeutralLoss")
  this_neuloss <- nodes_neuloss[[length(nodes_neuloss)]]
  
  xml2::xml_set_attrs(
    this_neuloss,
    c(mono_mass = neuloss_mono_mass, 
      avge_mass  = neuloss_avge_mass, 
      flag = "false", 
      composition = neuloss_composition))
  
  invisible(this_neuloss)
}


#' Adds a series of nodes of \code{umod:element}.
#' 
#' @param node A node of either \code{NeutralLoss} or \code{delta}.
#' @param composition The chemical composition of a modification.
add_comp_elements <- function (node = NULL, composition = "0") 
{
  df <- parse_unimod_composition(composition)

  old_elems <- xml2::xml_find_all(node, "umod:element")
  xml2::xml_remove(old_elems)
  rm(list = c("old_elems"))
  
  seqs <- seq_len(nrow(df))
  
  for (i in seqs)
    xml2::xml_add_child(node, "umod:element")
  
  mod_elems <- xml2::xml_find_all(node, "umod:element")
  
  for (i in seqs) 
    xml2::xml_set_attrs(mod_elems[[i]], c(symbol = df[i, 1], number = df[i, 2]))

  invisible(node)
}


#' Removes an existing modification.
#'
#' @inheritParams add_unimod
#' @seealso remove_unimod_title
#' @examples 
#' \donttest{
#' # site `C`: Oxiation + Carbamidomethyl
#' # (without neutral losses)
#' x <- remove_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'                                    full_name   = "Oxidation and iodoacetamide derivative"),
#'                    specificity = c(site        = "C",
#'                                    position    = "Anywhere"),
#'                    delta       = c(mono_mass   = "73.016379",
#'                                    avge_mass   = "73.0507",
#'                                    composition = "H(3) C(2) N O(2)"),
#'                    neuloss     = c(mono_mass   = "0",
#'                                    avge_mass   = "0",
#'                                    composition = "0"))
#' 
#' # site `M`: Oxiation + Carbamidomethyl
#' # (with neutral losses)
#' x <- remove_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'                                    full_name   = "Oxidation and iodoacetamide derivative"),
#'                    specificity = c(site        = "M",
#'                                    position    = "Anywhere"),
#'                    delta       = c(mono_mass   = "73.016379",
#'                                    avge_mass   = "73.0507",
#'                                    composition = "H(3) C(2) N O(2)"),
#'                    neuloss     = c(mono_mass   = "63.998285",
#'                                    avge_mass   = "64.1069",
#'                                    composition = "H(4) C O S"))
#' 
#' x <- remove_unimod(header      = c(title       = "Oxi+Carbamidomethyl",
#'                                    full_name   = "Oxidation and iodoacetamide derivative"),
#'                    specificity = c(site        = "M",
#'                                    position    = "Anywhere"),
#'                    delta       = c(mono_mass   = "73.016379",
#'                                    avge_mass   = "73.0507",
#'                                    composition = "H(3) C(2) N O(2)"),
#'                    neuloss     = c(mono_mass   = ".",
#'                                    avge_mass   = ".",
#'                                    composition = "."))
#' 
#' x <- remove_unimod_title("TMT10plexNterm+Gln->pyro-Glu")
#' x <- remove_unimod_title("Oxi+Carbamidomethyl")
#' }
#' @export
remove_unimod <- function (header = c(title = "Foo", full_name = "Foo bar"), 
                           specificity = c(site  = "C", position = "Anywhere"), 
                           delta = c(mono_mass = "42.010565", 
                                     avge_mass = "42.0367", 
                                     composition = "H(2) C(2) O"), 
                           neuloss = c(mono_mass = "0", 
                                       avge_mass = "0", 
                                       composition = "0")) 
{
  options(digits = 9L)
  
  title <- header[["title"]]
  full_name <- header[["full_name"]]
  site <- specificity[["site"]]
  position <- specificity[["position"]]
  mod_mono_mass <- delta[["mono_mass"]]
  mod_avge_mass <- delta[["avge_mass"]]
  mod_composition <- delta[["composition"]]
  neuloss_mono_mass <- neuloss[["mono_mass"]]
  neuloss_avge_mass <- neuloss[["avge_mass"]]
  neuloss_composition <- neuloss[["composition"]]
  
  xml_file <- system.file("extdata", "custom.xml", package = "proteoM")
  xml_root <- xml2::read_xml(xml_file)
  nodes_lev1_four <- xml2::xml_children(xml_root)
  node_modif <- xml2::xml_find_all(nodes_lev1_four, "//umod:modifications")
  
  if (!length(node_modif)) 
    stop("Node `umod:modifications` not found and nothing to remove from.", 
         call. = FALSE)
  
  site <- standardize_unimod_ps(x = c(site = site))
  position <- standardize_unimod_ps(x = c(position = position))
  
  if (site == "." && position == ".")
    return(remove_unimod_title(title = title))
  
  if (site == "." || position == ".")
    stop("Specify both `site` and `position`.", call. = FALSE)
  
  if (is.numeric(mod_mono_mass)) 
    mod_mono_mass <- as.character(mod_mono_mass)
  
  if (is.numeric(mod_avge_mass)) 
    mod_avge_mass <- as.character(mod_avge_mass)
  
  if (is.numeric(neuloss_mono_mass)) 
    neuloss_mono_mass <- as.character(neuloss_mono_mass)
  
  if (is.numeric(neuloss_avge_mass)) 
    neuloss_avge_mass <- as.character(neuloss_avge_mass)
  
  # individual nodes of modifications
  modifications <- xml2::xml_children(node_modif)
  
  # title
  idx_title <- which(xml2::xml_attr(modifications, "title") == title)
  
  local({
    len_title <- length(idx_title)
    
    if (len_title > 1L)
      stop("Multiple matches to `", title, "`.\n", 
           "Fix the redundancy from ", xml_file, ".")
    else if (!len_title) 
      stop("Entry `", title, "` not found.")
  })
  
  # `site` and `position`
  this_mod <- modifications[[idx_title]]
  
  local({
    this_full_name <- xml2::xml_attr(this_mod, "full_name")
    
    if (this_full_name != full_name)
      warning("Ignored the difference between the original `full_name = ", 
              this_full_name, "` \nand the current ", 
              "`full_name = ", full_name, "`.", 
              call. = FALSE)
  })
  
  # (children can be `specification`, `delta` etc.)
  nodes_mod_ch <- xml2::xml_children(this_mod)
  attrs_mod_ch <- xml2::xml_attrs(nodes_mod_ch)
  
  sites <- unlist(lapply(attrs_mod_ch, `[`, c("site")))
  positions <- unlist(lapply(attrs_mod_ch, `[`, c("position")))
  ok_sitepos <- which((sites == site) & (positions == position))
  
  local({
    len_sitepos <- length(ok_sitepos)
    
    if (!len_sitepos)
      stop("No matches to `site = ", site, "` and ", 
           "`position = ", position, "` at ", 
           "`title = ", title, "`.\n", 
           call. = FALSE)
    else if (len_sitepos > 1L)
      stop("Multiple matches to `site = ", site, "` and ", 
           "`position = ", position, "` at ", 
           "`title = ", title, "`.\n", 
           "Fix the redundancy from ", xml_file, ".", 
           call. = FALSE)
  })

  this_spec <- nodes_mod_ch[ok_sitepos]
  
  # After "complete" removals of `specificity`, 
  # a modification may have `delta` without `specificity`.
  # This seems innocuous for now or simply `remove_unimod_title`.

  if (neuloss_mono_mass == "0") {
    xml2::xml_remove(this_spec)
    nodes_mod_ch <- xml2::xml_children(this_mod)
    xml2::write_xml(xml_root, xml_file)
    
    message("Entry `site = ", site, "` and `position = ", position, 
            "` removed from `title = ", title, "`.")
  }
  else if (neuloss_mono_mass == ".") {
    nodes_neuloss <- xml2::xml_find_all(this_spec, "umod:NeutralLoss")
    xml2::xml_remove(nodes_neuloss[seq_along(nodes_neuloss)])
    nodes_neuloss <- xml2::xml_children(this_spec)
    
    if (length(nodes_neuloss))
      stop("No NeutralLoss nodes expected; contact the developer for bugs.",
           call. = FALSE)
    
    xml2::write_xml(xml_root, xml_file)
    
    message("All NeutralLoss under ", 
            "`site = ", site, "` and `position = ", position, 
            "` removed from `title = ", title, "`.")
  }
  else {
    nodes_neuloss <- xml2::xml_find_all(this_spec, "umod:NeutralLoss")
    attrs_neuloss <- xml2::xml_attrs(nodes_neuloss)
    
    neuloss_mono_masses <- unlist(lapply(attrs_neuloss, `[`, c("mono_mass")))
    idx_neuloss <- which((neuloss_mono_masses == neuloss_mono_mass))
    
    local({
      len_neuloss <- length(idx_neuloss)
      
      if (len_neuloss > 1L)
        stop("Multiple matches to the `mono_mass` in `neuloss`.", call. = FALSE)
      else if (!len_neuloss) 
        stop("No matches to the `mono_mass` in `neuloss`.", call. = FALSE)
    })
    
    this_neuloss <- nodes_neuloss[idx_neuloss]
    xml2::xml_remove(this_neuloss)
    nodes_neuloss <- xml2::xml_children(this_spec)
    
    # only the "0" NeutralLoss
    if (length(nodes_neuloss) == 1L) {
      this_neuloss_0 <- nodes_neuloss[1]
      xml2::xml_remove(this_neuloss_0)
      nodes_neuloss <- xml2::xml_children(this_spec)
      
      if (length(nodes_neuloss))
        stop("No NeutralLoss nodes expected; contact the developer for bugs.",
             call. = FALSE)
    }
    
    xml2::write_xml(xml_root, xml_file)
    
    message("NeutralLoss at `mono_mass = ", neuloss_mono_mass, "`, ", 
            "`site = ", site, "` and `position = ", position, 
            "` removed from `title = ", title, "`.")
  }
  
  invisible(xml_root)
}


#' Standardizes the site and position of a modification.
#' 
#' @param x A name string. Note that \code{x} is either a site or a position.
standardize_unimod_ps <- function (x) 
{
  nm <- names(x)
  x <- unname(x)
  
  if (isFALSE(x)) 
    x <- "."
  
  if (length(x) != 1L) 
    stop("The length of `", nm, "` is not exactly one.", call. = FALSE)
  
  if (nchar(x) == 0L) 
    x <- "."
  
  if (nm == "site") {
    ok_sites <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
                  "P", "Q", "R", "S", "T", "U", "V", "W", "Y", 
                  "N-term", "C-term", ".")
    
    # site may later include "X", "Z" etc.
    if (! x %in% c(LETTERS, "N-term", "C-term", "."))
      stop("Invalid site = ", x, call. = FALSE)
  }
  
  if (nm == "position") {
    ok_positions <- c("Anywhere", "Protein N-term", "Protein C-term", 
                      "Any N-term", "Any C-term", ".")
    
    if (! x %in% ok_positions)
      stop("Invalid position = ", x, call. = FALSE)
  }
  
  invisible(x)
}


#' Removes an existing modification.
#'
#' @param title The title of a modification.
#' @examples 
#' \donttest{
#' x <- remove_unimod_title("Oxi+Carbamidomethyl")
#' }
#' @export
remove_unimod_title <- function (title = NULL) 
{
  if (isFALSE(title) || nchar(title) == 0L)
    stop("Provide a `title`.", call. = FALSE)
  
  xml_file <- system.file("extdata", "custom.xml", package = "proteoM")
  xml_root <- xml2::read_xml(xml_file)
  
  nodes_lev1_four <- xml2::xml_children(xml_root)
  node_modif <- xml2::xml_find_all(nodes_lev1_four, "//umod:modifications")
  modifications <- xml2::xml_children(node_modif)
  
  idx <- which(xml2::xml_attr(modifications, "title") == title)
  
  len <- length(idx)
  
  if (len > 1L)
    stop("Multiple matches to `", title, "`.\n", 
         "Fix the redundancy from ", xml_file, ".")
  else if (!len)
    stop("Entry `", title, "` not found.")
  
  xml2::xml_remove(modifications[idx])
  modifications <- xml2::xml_children(node_modif)
  xml2::write_xml(xml_root, xml_file)
  message("Entry `", title, "` removed from ", xml_file, ".")
  
  invisible(xml_root)
}


#' Calculates the masses of a chemical formula.
#'
#' @param composition A chemical composition.
#' @param digits A non-negative integer; the number of decimal places to be
#'   used.
#' @export
calc_unimod_compmass <- function (composition = "H(4) C O S", digits = 6L) 
{
  options(digits = 9L)
  
  nm <- system.file("extdata", "elem_masses.txt", package = "proteoM")
  
  if (file.exists(nm))
    lookup <- read.delim(file = nm, sep = "\t")
  else
    stop("Not found: ", nm, call. = FALSE)
  
  df <- parse_unimod_composition(composition)
  df$number <- as.numeric(df$number)
  
  df <- dplyr::left_join(df, lookup, by = "symbol")
  
  rows <- is.na(df$mono_mass)
  
  if (any(rows)) 
    stop("Unknown element(s): ", paste(df$symbol[rows], collapse = ", "), 
         call. = FALSE)
  
  nums <- df$number
  avges <- df$avge_mass
  monos <- df$mono_mass
  
  avge_mass <- round(sum(avges * nums), digits)
  mono_mass <- round(sum(monos * nums), digits)
  
  list(mono_mass = mono_mass, avge_mass = avge_mass)
}


#' Parses A unimod position.
#' 
#' @param composition A chemical composition.
parse_unimod_composition <- function (composition = "H(4) C O S") 
{
  options(digits = 9L)
  
  df <- composition %>% 
    stringr::str_replace_all("([:alnum:]+)$", paste0("\\1", "(1)")) %>% 
    stringr::str_replace_all("([:alnum:]+) ", paste0("\\1", "(1) ")) %>% 
    gsub("\\s+", "", .) %>% 
    stringr::str_match_all("(.*?)\\((-*[0-9]+)\\)") %>% 
    `[[`(1)
  
  df <- df[, -1]
  colnames(df) <- c("symbol", "number")
  
  data.frame(df)
}
