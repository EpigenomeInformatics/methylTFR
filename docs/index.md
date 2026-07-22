# *methylTFR* : Quantification of DNA Methylation Patterns in TFBS

`methylTFR` is an R-package to analyze DNA methylation signatures in
transcription factor binding sites in each individual cells or samples.

## Installation instructions

Get the latest release `methylTFR` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("methylTFR")
```

And the development version from
[GitHub](https://github.com/EpigenomeInformatics/methylTFR) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("EpigenomeInformatics/methylTFR")
```

## Documentation

Full documentation and vignettes are hosted at
[epigenomeinformatics.github.io/methylTFR](https://epigenomeinformatics.github.io/methylTFR/):

- [Get
  started](https://epigenomeinformatics.github.io/methylTFR/articles/methylTFR.html)
  — reading data, computing deviations, and footprints.
- [Case study: memory vs. naive T
  cells](https://epigenomeinformatics.github.io/methylTFR/articles/memTcells.html)
  — differential TF activity on bundled example data.

## Quick Start

This is a basic example which shows you how to run `methylTFR` :

``` r
library(GenomicRanges)
library(dplyr)
library(methylTFRAnnotationHg38) # annotation package for hg38
library(methylTFR)

gcfreqs <- getGCfreq(motifSet = "jaspar2020")
gc_dist <- getGenomeGC("hg38")
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

sample_dir <- file.path("samples_dir")
sample_ann <- "samples.tsv" # should contain column name bedFile

# deviation score matrix
deviations <- run_methyltfr(sample_ann, # sample annotation file
  sample_dir, # where the EPP files are
  threads = 8, # number of threads
  chunkSize = 10, # number of chunks to process
  sampleColName = "bedFile", # column name for EPP file paths in sample_ann
  tf_bindsites = tf_bindsites, # TF binding sites
  gcfreqs = gcfreqs, # GC frequency
  gc_dist = gc_dist, # GC distribution
  filetype = "EPP" # file type
)
```

## Citation

Below is the citation output from using `citation('methylTFR')` in R.
Please run this yourself to check for any updates on how to cite
**methylTFR**.

``` r
print(citation("methylTFR"), bibtex = TRUE)
#> 
#> To cite package 'methylTFR' in publications use:
#> 
#>   Irem B. Gündüz, Sarath Kumar Murugan and Fabian Muller (2026).
#>   methylTFR: Quantification of DNA Methylation Signatures in TFBS. R
#>   package version 0.99.0.
#>   https://github.com/EpigenomeInformatics/methylTFR
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {methylTFR: Quantification of DNA Methylation Signatures in TFBS},
#>     author = {Irem B. Gündüz and Sarath Kumar Murugan and Fabian Muller},
#>     year = {2026},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/EpigenomeInformatics/methylTFR},
#>   }
```
