# combine all .R files into one
# foo_combine_codes(filepath = file.path("C:/Results/R/proteoQ/inst/extdata/examples"))
foo_combine_codes <- function (filepath = file.path("~/Github/mzion/R")) 
{
  filepath <- mzion:::find_dir(filepath)
  filenames <- dir(filepath, pattern = ".R$")

  dir.create(file.path(filepath, "temp"), showWarnings = FALSE)

  ans <- lapply(file.path(filepath, filenames), readLines)
  ans <- purrr::reduce(ans, `c`, init = NULL)
  writeLines(ans, file.path(filepath, "temp/all - mzion.R"))
}


foo_list_func <- function (filepath = file.path(file.path("~/Github/mzion/R"))) 
{
  filepath <- mzion:::find_dir(filepath)
  filenames <- dir(filepath, pattern = ".R$")
  
  dir.create(file.path(filepath, "temp"), showWarnings = FALSE)
  
  lines <- lapply(file.path(filepath, filenames), readLines)
  
  fns_all <- lapply(lines, function (x) {
    fn_lines <- x[grepl("<-\\s*function\\s*\\(", x)]
    fns <- gsub("^(.*)\\s*<- function\\s*\\(.*", "\\1", fn_lines)
    gsub("\\s*$", "", fns)
  })
  
  names(fns_all) <- filenames
  
  sink(file.path(filepath, "temp/funs.txt"))
  fns_all
  sink()
  
  ans <- readLines(file.path(filepath, "temp/funs.txt"))
  ans <- paste0("# ", ans)
  writeLines(ans, file.path(filepath, "temp/funs.R"))
}


foo_find_fmlmass <- function () 
{
  options(digits = 9L)
  
  xml_file <- system.file("extdata", "master.xml", package = "mzion")
  xml_root <- xml2::read_xml(xml_file)
  nodes_lev1_four <- xml2::xml_children(xml_root)
  node_elem <- xml2::xml_find_all(nodes_lev1_four, "//umod:elements")
  elements <- xml2::xml_find_all(node_elem, "umod:elem")
  
  attrs_elem <- xml2::xml_attrs(elements)
  
  ans <- data.frame(do.call(rbind, attrs_elem))
  ans$mono_mass <- round(as.numeric(ans$mono_mass), digits = 6L)
  ans$avge_mass <- round(as.numeric(ans$avge_mass), digits = 5L)
  
  col_title <- which(colnames(ans) == "title")
  
  stopifnot(length(col_title) == 1L)
  
  colnames(ans)[col_title] <- "symbol"
  
  saveRDS(ans, "~/mzion/elem_masses.rds")
  
  write.table(ans, file = "~/mzion/elem_masses.txt", sep = "\t", 
              row.names = FALSE)
}