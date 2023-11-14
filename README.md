Mzion
================
true
2023-11-12

- [Installation](#installation)
- [Peaklist formats](#peaklist-formats)
- [Database searches (via an app)](#database-searches-via-an-app)
- [Database searches (via R scripts)](#database-searches-via-r-scripts)
- [Help documents](#help-documents)
- [Specifications of fixed and variable
  modifications](#specifications-of-fixed-and-variable-modifications)
- [Data QC and mining](#data-qc-and-mining)
- [Other utilities](#other-utilities)
- [Cite Mzion](#cite-mzion)

## Installation

To install this package, start R and enter:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/mzion")

# or install both Mzion and its Shiny App:
devtools::install_github("qzhang503/mzionShiny")
```

## Peaklist formats

- Thermo’s MS
  - [x] `MSConvert mgf`

  - [x] `MSConvert mzML`

    \[+\] Binary encode 64-bit

    \[+\] Write index

    \[+\] TPP compatibility

    \[-\] Use zlib compression (*do not compress*)

    \[+\] Filters

    1)  peakPicking: vendor msLevel = 1-

    2)  (optional) zeroSamples: removeExtra 1-
- Bruker’s MS
  - [x] `DataAnalysis mgf`

## Database searches (via an app)

``` r
library(mzionShiny)
run_app()
```

## Database searches (via R scripts)

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
  fdr_type  = "psm",
)

## See also ?matchMS for SILAC, acetylome workflows etc.
```

## Help documents

Enter `?mzion::matchMS` from an RStudio section.

## Specifications of fixed and variable modifications

The Unimod definition of positions and sites were adopted by Mzion for
specifying fixed and variable modifications. The value of a position is
in one of “Anywhere”, “Protein N-term”, “Protein C-term”, “Any N-term”
or “Any C-term”. The last two position labels can be shorthanded as
“N-term” and “C-term”. A site is a one-letter representation of the
twenty amino-acid residues, as well as the terminal sites of “N-term”
and “C-term”. The general format in specifying a fixed or variable
modification is `title (position = site)` where title is a unique
character string without space. At a position of “Anywhere”, the
modification can be shorthanded as `title (site)`, for example,
`TMT10plex (K)`. For a terminal modification at any site, it can be
abbreviated as `title (position)`, for examples,
`Acetyl (Protein N-term)` and `TMT10plex (N-term)`. There are
circumstances that both position and site are needed for specifying a
modification, for instance, `Gln->pyro-Glu (N-term = Q)`. More examples
are available in the help document of Mzion utility of `parse_unimod`.

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

## Cite Mzion

Zhang, Q. Mzion enables deep and precise identification of peptides in
data-dependent acquisition proteomics. Sci. Rep. 13, 7056 (2023).
