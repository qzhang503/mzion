
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoM!\n\n",
                        "===============================================================================\n",
                        # "NEW features (v1.0.0.0):\n",
                        "[x] A search engine for mass spectrometry-based proteomics data.\n",
                        "[x] See also package `proteoQ` for downstream data QA and informatics.\n", 
                        "\n",
                        
                        # "Notes:\n",
                        "[x] Minimum system requirement: 32GB RAM and 8 (dual)-cores.\n",
                        "[x] Not yet tested under Linux or MacOS.\n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n",
                        "===============================================================================\n")
}
