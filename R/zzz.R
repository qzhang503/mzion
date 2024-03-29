
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Mzion.\n\n",
                        "============================================================================================\n",
                        # "NEW features (v1.3.3.4):\n",
                        # "[x] MS1, MS2 de-isotoping and chimeric peptide searches.\n\n",
                        "[x] For documents, enter \"?matchMS\".\n", 

                        # "Notes:\n",
                        "[x] Suggested configuration for large datasets: 32GB RAM and 8-cores.\n",
                        # "[x] May need to remove previously cached results (or use a new .path_cache and .path_fasta).\n",
                        # "    (OOps an accidental relocation of amino-acid lookups to a parent directory) \n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n",
                        "============================================================================================\n")
}
