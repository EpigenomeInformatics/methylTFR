# plotMotifFootprint

Creates a footprint plot of bias corrected methylation for a given motif
and sample

## Usage

``` r
plotMotifFootprint(
  motif,
  tf_bindsites,
  msites,
  sample_name = NULL,
  gc_dist,
  gcfreqs,
  enhancer = NULL,
  method = "division",
  flankNorm = 50
)
```

## Arguments

- motif:

  Motif name as a character string

- tf_bindsites:

  Transcript Factor binding sites from the annotation package

- msites:

  Methylation sites as a data frame

- sample_name:

  Optional sample name as a character string, defaults to the file name
  if not provided

- gc_dist:

  a `GRanges` object contains Genome wide GC distribution

- gcfreqs:

  GC frequency of the genome used to compute expected methylation

- enhancer:

  Specific region such as distal motif, proximal motif

- method:

  Method to calculate the difference, either "substraction" or
  "division"

- flankNorm:

  Numeric value indicating the flank normalization distance, default is
  50

## Value

A ggplot object of TF footprint difference plot for a single sample

## Examples

``` r
library(ggplot2)
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

# Plot the expected footprint
p <- plotMotifFootprint(
    motif = "BATF",
    tf_bindsites = tf_bindsites,
    msites = msites,
    sample_name = "Sample1",
    gc_dist = gcdist,
    gcfreqs = gcfreqs,
    enhancer = NULL,
    method = "division"
)
```
