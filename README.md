# methylTFR

`methylTFR` is an R-package to analyze the methylation data from whole genome bisulfite sequencing. It takes inputs from RnBeads processed bed file - EPP format. This tool aims to indentify the methylation patterns on transcription factor binding sites in each individual cells or samples. This tool is still under development and testing stage.

## Requirements

- R (>= 3.5.0)
- GenomicRanges
- data.table
- dplyr
- stringr
- parallel
- logger

```R

# To install dependencies 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages(c("data.table", "dplyr", "stringr"))

```

## Installation

```sh

git clone git@github.com:EpigenomeInformatics/methylTFR.git
cd methylTFR
Rscript -e "library(devtools); devtools::install('.')"

```

## QuickStart

The samples directory should contain methylation bed files and a tab delimited sample annotation file. The bedFile column in the sample annotation file should be used to hold bedfile names. 

```R

library(GenomicRanges)
library(dplyr)
library(methylTFRann)
library(methylTFR)

sample_dir <- file.path("samples_dir")
sample_ann <- "samples.tsv" # should contain column name bedFile

# deviation score matrix
deviations <- run_methyltfr(sample_ann,
                            sample_dir,
                            threads = 8)

```
