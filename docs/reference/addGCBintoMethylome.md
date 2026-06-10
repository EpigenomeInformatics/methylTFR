# addGCBintoMethylome

This function is used to add GC bin to the methylation sites and
calculate the average methylation for each bin.

## Usage

``` r
addGCBintoMethylome(msites, gcdist, ignoreStrand = TRUE)
```

## Arguments

- msites:

  Imported methylation sites using `read_methylome` function

- gcdist:

  a `GRanges` object contains Genome wide GC distribution

- ignoreStrand:

  - if TRUE, it ignores strand info from annotation

## Value

a `matrix` object with GC bin with corresponding avg methylation

## Author

Irem Gunduz

## Examples

``` r
library(methylTFR)

# Load the data
load(system.file("extdata",
    "gcdist_subset.rda",
    package = "methylTFR"
))
load(system.file("extdata",
    "example_data.rda",
    package = "methylTFR"
))

# Add GC bin
bin_meth <- addGCBintoMethylome(msites, gcdist)
```
