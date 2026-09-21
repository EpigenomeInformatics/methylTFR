<!-- README.md is generated from README.Rmd. Please edit that file -->

# *methylTFR* : Quantification of DNA Methylation Patterns in TFBS

<img src="man/figures/logo.png" align="right" height="139" alt="methylTFR logo" />

<!-- badges: start -->
[![Test R-universe](https://github.com/EpigenomeInformatics/methylTFR/actions/workflows/r-universe.yml/badge.svg)](https://github.com/EpigenomeInformatics/methylTFR/actions/workflows/r-universe.yml)
[![GitHub issues](https://img.shields.io/github/issues/EpigenomeInformatics/methylTFR)](https://github.com/EpigenomeInformatics/methylTFR/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/EpigenomeInformatics/methylTFR)](https://github.com/EpigenomeInformatics/methylTFR/pulls)
<!-- badges: end -->

`methylTFR` is an R-package to analyze DNA methylation signatures in transcription factor binding sites in each individual cells or samples.

## Installation instructions

Get the latest release `methylTFR` from [Bioconductor](http://bioconductor.org/) using the following code:


```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("methylTFR")
```

And the development version from [GitHub](https://github.com/EpigenomeInformatics/methylTFR) with:


```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("EpigenomeInformatics/methylTFR")
```
## Documentation

Full documentation and vignettes are hosted at
[epigenomeinformatics.github.io/methylTFR](https://epigenomeinformatics.github.io/methylTFR/):

- [Get started](https://epigenomeinformatics.github.io/methylTFR/articles/methylTFR.html) — reading data, computing deviations, and footprints.
- [Case study: memory vs. naive T cells](https://epigenomeinformatics.github.io/methylTFR/articles/memTcells.html) — differential TF activity on bundled example data.

## Quick Start

This is a basic example which shows you how to run `methylTFR` :


```r
library(GenomicRanges)
library(dplyr)
library(methylTFRAnnotationHg38) # annotation package for hg38
library(methylTFR)

gcfreqs <- getGCfreq(motifSet = "jaspar2020")
gc_dist <- getGenomeGC()
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

If you use `methylTFR` in your work, please cite:

> Gunduz IB, Murugan SK, Mueller F (2026). methylTFR: Quantification of DNA methylation signatures in TFBS. R package version 0.99.9. https://github.com/EpigenomeInformatics/methylTFR

```bibtex
@Manual{methylTFR,
  title = {methylTFR: Quantification of DNA methylation signatures in TFBS},
  author = {Irem B. Gunduz and Sarath Kumar Murugan and Fabian Mueller},
  year = {2026},
  note = {R package version 0.99.9},
  url = {https://github.com/EpigenomeInformatics/methylTFR},
  doi = {10.18129/B9.bioc.methylTFR},
}
```

Run `citation("methylTFR")` in R to get the entry for the version you have installed.
