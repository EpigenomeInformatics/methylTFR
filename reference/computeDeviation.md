# computeDeviation

computeDeviation is a function to calculate the deviation in
transcription factor binding sites for a given motif

## Usage

``` r
computeDeviation(
  motif,
  msites,
  tf_bindsites,
  gcfreqs,
  enhancer = NULL,
  ignoreStrand = TRUE,
  binMsites
)
```

## Arguments

- motif:

  A character vector containing motif name

- msites:

  imported methylation sites

- tf_bindsites:

  a `GRangesList` object contains tf binding sites positions

- gcfreqs:

  a `list` of GC bin frequency tables (matrices for multiple motif)

- enhancer:

  a `GRanges` object specifying regions such as distal motif (optional)

- ignoreStrand:

  if TRUE, it ignores strand info from annotation

- binMsites:

  A matrix object with GC bin with corresponding avg methylation

## Value

a `numeric` deviation score for a given motif

## Examples

``` r
library(methylTFR)

# Load the data
load(system.file("extdata",
    "BATF_tf_bindsites.rda",
    package = "methylTFR"
))
load(system.file("extdata",
    "example_data.rda",
    package = "methylTFR"
))
load(system.file("extdata",
    "BATF_gcfreqs.rda",
    package = "methylTFR"
))
load(system.file("extdata",
    "gcdist_subset.rda",
    package = "methylTFR"
))

# Compute binMsites
bin_meth <- addGCBintoMethylome(msites, gcdist, TRUE)

# Compute the deviation
devs <- computeDeviation("BATF",
    msites,
    tf_bindsites,
    gcfreqs,
    enhancer = NULL,
    ignoreStrand = TRUE,
    binMsites = bin_meth
)
```
