# computeExpectedFootprint

Compute expected footprint returns a data.table required to create a
transcription factor footprint with label.

## Usage

``` r
computeExpectedFootprint(motif, gcfreqs, gc_dist, enhancer = NULL, msites)
```

## Arguments

- motif:

  chracter vector containing motif name eg: GATA string

- gcfreqs:

  GC frequency of the genome used to compute expected methylation

- gc_dist:

  a `GRanges` object contains Genome wide GC distribution

- enhancer:

  Specific regions such as distal motif, proximal motif

- msites:

  Methylation data

## Value

a `data.table` object to containing TF footprint
