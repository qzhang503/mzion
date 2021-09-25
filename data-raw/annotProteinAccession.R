# combine all .R files into one
# foo_combine_codes(filepath = file.path("C:/Results/R/proteoQ/inst/extdata/examples"))
foo_combine_codes <- function (filepath = file.path("E:/R/proteoM/R")) {
  filenames <- dir(filepath, pattern = ".R$")

  dir.create(file.path(filepath, "temp"))

  ans <- lapply(file.path(filepath, filenames), readLines)
  ans <- purrr::reduce(ans, `c`, init = NULL)
  writeLines(ans, file.path(filepath, "temp/all - proteoM.R"))
}
