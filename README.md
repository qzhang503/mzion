proteoM
================
true
2021-09-15

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
copy_refseq_hs("~/proteoM/dbs/fasta/refseq")
copy_refseq_mm("~/proteoM/dbs/fasta/refseq")
copy_crap("~/proteoM/dbs/fasta/crap")

# MGF (e.g. by Proteome Discoverer)
copy_pd_mgf("~/proteoM/examples/mgfs")

# Ion searches
library(proteoM)

matchMS(
  out_path = "~/proteoM/examples", 
  mgf_path = "~/proteoM/examples/mgfs",
  fasta = c("~/proteoM/dbs/fasta/refseq/refseq_hs_2013_07.fasta", 
            "~/proteoM/dbs/fasta/refseq/refseq_mm_2013_07.fasta", 
            "~/proteoM/dbs/fasta/crap/crap.fasta"), 
  acc_type = c("refseq_acc", "refseq_acc", "other"), 
  fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"),
  varmods = c("Acetyl (Protein N-term)", "Oxidation (M)",
             "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)"),
  max_miss = 4, 
  quant = "tmt10", 
  fdr_type = "protein", 
)


## An example of phospho TMT
library(proteoQDA)

# Try a different fasta database
copy_uniprot_hsmm("~/proteoM/dbs/fasta/uniprot")

# Try MGFs, e.g., by MSConvert
copy_msconv_mgf("~/proteoM/examples_p/mgfs")

# Ion searches
library(proteoM)

matchMS(
  out_path = "~/proteoM/examples_p", 
  mgf_path = "~/proteoM/examples_p/mgfs",
  fasta = c("~/proteoM/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta", 
            "~/proteoM/dbs/fasta/crap/crap.fasta"), 
  acc_type = c("uniprot_acc", "other"), 
  fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"), 
  varmods = c("Acetyl (N-term)", "Oxidation (M)", 
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
