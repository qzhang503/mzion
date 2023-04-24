
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to mzion.\n\n",
                        "============================================================================================\n",
                        # "NEW features (v1.2.4):\n",
                        # "[x] Incompatible with cached results from previous versions.\n\n",
                        "[x] For examples, enter \"?matchMS\".\n", 
                        # "[x] Please delete cached \"\temp\pep_score.rds\" for reprocessing wither older versions.\n",

                        # "[x] Added Percolator utility.\n",
                        # "[x] See also package `proteoQ` for downstream data QA and informatics.\n", 
                        # "\n",
                        
                        # "Notes:\n",
                        "[x] Suggested configuration for large datasets: 32GB RAM and 8-cores.\n",
                        # "[x] May need to remove previously cached results (or use a new .path_cache and .path_fasta).\n",
                        # "    (OOps an accidental relocation of amino-acid lookups to a parent directory) \n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n",
                        "============================================================================================\n")
}
