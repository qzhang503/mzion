Mzion
================
true
2023-04-30

- <a href="#installation" id="toc-installation">Installation</a>
- <a href="#peaklist-formats" id="toc-peaklist-formats">Peaklist
  formats</a>
- <a href="#help-documents" id="toc-help-documents">Help documents</a>
- <a href="#specifications-of-fixed-and-variable-modifications"
  id="toc-specifications-of-fixed-and-variable-modifications">Specifications
  of fixed and variable modifications</a>
- <a href="#database-searches" id="toc-database-searches">Database
  searches</a>
- <a href="#data-qc-and-mining" id="toc-data-qc-and-mining">Data QC and
  mining</a>
- <a href="#other-utilities" id="toc-other-utilities">Other utilities</a>
- <a href="#cite-mzion" id="toc-cite-mzion">Cite Mzion</a>

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
abbreviated as `title (position)`, for example,
`Acetyl (Protein N-term)` and `TMT10plex (N-term)`. There are
circumstances that both position and site are needed for specifying a
modification, for instance, `Gln->pyro-Glu (N-term = Q)`. More examples
are available in the help document of Mzion utility of `parse_unimod`.

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
  fdr_type  = "psm",
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

## Cite Mzion

Zhang, Q. Mzion enables deep and precise identification of peptides in
data-dependent acquisition proteomics. Sci. Rep. 13, 7056 (2023).
