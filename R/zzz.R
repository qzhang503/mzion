
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoM!\n\n",
                        "===============================================================================\n",
                        "NEW features (v0.1.0.2):\n",
                        "[x] Database searches separated from package proteoQ.\n",
                        "\n",
                        
                        "Notes:\n",
                        "[x] New handling of conflicting fixed/variable modifications and \n",
                        ",   requires the removals of older caches ",
                        "under `.MSearch` and `dbs/pepmasses`.\n",
                        "[x] Reprocessing of `matchMS` with a new R session supported.\n", 
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n",
                        "===============================================================================\n")
}
