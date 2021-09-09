proteoM
================
true
2021-09-09

-   [Installation](#installation)
-   [Database searches](#database-searches)
-   [Possible Next steps](#possible-next-steps)

## Installation

To install this package, start R (latest version) and enter:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/proteoM")
```

## Database searches

``` r
## An example of global TMT
library(proteoQDA)

# Fasta databases 
copy_refseq_hs("~/proteoQ/dbs/fasta/refseq")
copy_refseq_mm("~/proteoQ/dbs/fasta/refseq")
copy_crap("~/proteoQ/dbs/fasta/crap")

# MGF (e.g. by Proteome Discoverer)
copy_pd_mgf("~/proteoQ/examples/mgfs")

# Ion searches
library(proteoM)

matchMS(
  out_path = "~/proteoQ/examples", 
  mgf_path = "~/proteoQ/examples/mgfs",
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta", 
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta", 
            "~/proteoQ/dbs/fasta/crap/crap.fasta"), 
  acc_type = c("refseq_acc", "refseq_acc", "other"), 
  fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"),
  varmods = c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
             "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)"),
  max_miss = 4, 
  quant = "tmt10", 
  fdr_type = "protein", 
)


## An example: phospho TMT
library(proteoQDA)

# Try a different fasta database
copy_uniprot_hsmm("~/proteoQ/dbs/fasta/uniprot")

# Try MGFs, e.g., by MSConvert
copy_msconv_mgf("~/proteoQ/examples_p/mgfs")

# Ion searches
library(proteoM)

matchMS(
  out_path = "~/proteoQ/examples_p", 
  mgf_path = "~/proteoQ/examples_p/mgfs",
  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta", 
            "~/proteoQ/dbs/fasta/crap/crap.fasta"), 
  acc_type = c("uniprot_acc", "other"), 
  fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"), 
  varmods = c("TMT6plex (N-term)", "Acetyl (N-term)", "Oxidation (M)", 
              "Deamidated (N)", "Phospho (S)", "Phospho (T)", 
              "Phospho (Y)", "Gln->pyro-Glu (N-term = Q)"), 
  max_miss = 4, 
  quant = "tmt10", 
  fdr_type = "protein", 
)
```

## Possible Next steps

-   Search against real MGFs.
-   Data QC and informatics with `library(proteoQ)`.
