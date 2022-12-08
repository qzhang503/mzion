proteoM
================
true
2022-12-07

- <a href="#installation" id="toc-installation">Installation</a>
- <a href="#fastas-and-mgfs" id="toc-fastas-and-mgfs">FASTAs and MGFs</a>
- <a href="#database-searches" id="toc-database-searches">Database
  searches</a>
- <a href="#next-steps" id="toc-next-steps">Next steps</a>
- <a href="#other-utilities" id="toc-other-utilities">Other utilities</a>

## Installation

To install this package, start R and enter:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/proteoM")
```

## FASTAs and MGFs

``` r
# (feel free to skip, only for exemplary FASTAs and MGFs)
devtools::install_github("qzhang503/proteoQDA")

library(proteoQDA)

copy_refseq_hs("~/proteoM/dbs/fasta/refseq")
copy_refseq_mm("~/proteoM/dbs/fasta/refseq")
copy_uniprot_hsmm("~/proteoM/dbs/fasta/uniprot")
copy_crap("~/proteoM/dbs/fasta/crap")

# MGF (by MSConvert, Proteome Discoverer or Bruker's DataAnalysis)
copy_pd_mgf("~/proteoM/examples/mgfs")
copy_msconv_mgf("~/proteoM/examples_p/mgfs")
```

## Database searches

``` r
## Global, TMT-10plex

# If not using examples from proteoQDA
#   (1) link users' MGF or mzML files to 'mgf_path'
#   (2) link FASTA files to 'fasta'
#   (3) make addtional adjustments accordingly
library(proteoM)

matchMS(
  out_path  = "~/proteoM/examples", 
  mgf_path  = "~/proteoM/examples/mgfs",
  fasta     = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta", 
                "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta", 
                "~/proteoM/dbs/fasta/crap/crap.fasta"), 
  
  # (see also ?load_fasta2)
  acc_type  = c("refseq_acc", "refseq_acc", "other"), 
  
  # "TMT6plex" at mass 229.162932 Da for TMT-6, -10 and -11 
  # (see also ?table_unimods)
  fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"),
  varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)",
               "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)"),
  max_miss  = 4, 
  quant     = "tmt10", 
  fdr_type  = "protein", 
)


## Phospho, TMT
library(proteoM)

matchMS(
  out_path  = "~/proteoM/examples_p", 
  mgf_path  = "~/proteoM/examples_p/mgfs",
  
  fasta     = c("~/proteoM/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta", 
                "~/proteoM/dbs/fasta/crap/crap.fasta"), 
  acc_type  = c("uniprot_acc", "other"), 
  fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"), 
  varmods   = c("Acetyl (Protein N-term)", "Oxidation (M)", 
                "Deamidated (N)", "Phospho (S)", "Phospho (T)", 
                "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"), 
  locmods   = c("Phospho (S)", "Phospho (T)", "Phospho (Y)"), 
  max_miss  = 4, 
  quant     = "tmt10", 
  
  fdr_type  = "protein",
)

## See also ?matchMS for SILAC, acetylome etc.
```

## Next steps

- Search against full-length MGF or mzML files

  - Thermo’s MS
    - [x] `MSConvert`

      \[+\] Binary encode 64-bit

      \[+\] Write index

      \[+\] TPP compatibility

      \[-\] (*uncheck*) Use zlib compression

    - [x] `Proteome Discoverer`
  - Bruker’s MS
    - [x] `DataAnalysis`

- Data QC and mining with
  [proteoQ](https://github.com/qzhang503/proteoQ/)

## Other utilities

- `mapMS2ions`: visualizations of matched MS2 ions
- `table_unimods` summary of Unimod entries
- `add_unimod`: addition of a Unimod entry
