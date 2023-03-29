Mzion
================
true
2023-03-25

- <a href="#installation" id="toc-installation">Installation</a>
- <a href="#peaklist-formats" id="toc-peaklist-formats">Peaklist
  formats</a>
- <a href="#help-documents" id="toc-help-documents">Help documents</a>
- <a href="#database-searches" id="toc-database-searches">Database
  searches</a>
- <a href="#data-qc-and-mining" id="toc-data-qc-and-mining">Data QC and
  mining</a>
- <a href="#other-utilities" id="toc-other-utilities">Other utilities</a>

## Installation

To install this package, start R and enter:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/mzion")
```

## Peaklist formats

- Thermo’s MS
  - [x] `MSConvert mgf`

  - [x] `MSConvert mzML`

    \[+\] Binary encode 64-bit

    \[+\] Write index

    \[+\] TPP compatibility

    \[-\] Use zlib compression (do not compress)

  - [x] `Proteome Discoverer mgf`
- Bruker’s MS
  - [x] `DataAnalysis mgf`

## Help documents

Enter `?mzion::matchMS` from an R console.

## Database searches

``` r
## Global, TMT-10plex
library(mzion)

matchMS(
  out_path  = "~/mzion/examples", 
  mgf_path  = "~/mzion/examples/mgfs",
  fasta     = c("~/mzion/dbs/fasta/refseq/refseq_hs_2013_07.fasta", 
                "~/mzion/dbs/fasta/refseq/refseq_mm_2013_07.fasta", 
                "~/mzion/dbs/fasta/crap/crap.fasta"), 
  
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
matchMS(
  out_path  = "~/mzion/examples_p", 
  mgf_path  = "~/mzion/examples_p/mgfs",
  
  fasta     = c("~/mzion/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta", 
                "~/mzion/dbs/fasta/crap/crap.fasta"), 
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

## See also ?matchMS for SILAC, acetylome workflows and more.
```

## Data QC and mining

- [proteoQ](https://github.com/qzhang503/proteoQ/)

## Other utilities

- `mapMS2ions`: visualizes MS2 spectrum matches
- `table_unimods` summarizes Unimod entries
- `find_unimod`: finds a Unimod entry
- `parse_unimod`: parses a Unimod entry
- `calc_unimod_compmass`: calculates the composition masses of a Unimod
- `add_unimod`: adds a Unimod entry
- `remove_unimod`: removes a Unimod entry
- `remove_unimod_title`: removes a Unimod entry by title
- `make_mztab`: prepares a mzTab file from the search results
