# computeFootprint

Observed methylation profile around the binding sites of one motif.

## Usage

``` r
computeFootprint(
  motif_name,
  tf_bindsites,
  msites,
  enhancer = NULL,
  extendBy = 130
)
```

## Arguments

- motif_name:

  Name of the motif.

- tf_bindsites:

  A `GRangesList` named by motif. Each element is a `GRanges` of the
  occurrences of that motif. The names must include every name in
  `names(gcfreqs)`. Get it with
  `getTFbindsites(motifSet = "jaspar2020")` from the
  methylTFRAnnotationHg38 (or methylTFRAnnotationMm10) package, or build
  a custom set with methylTFRAnnotationBuilder.

- msites:

  A `GRanges` object with one range per CpG (width 1) and the metadata
  columns `score` (numeric methylation fraction between 0 and 1) and
  `coverage` (read count). Usually produced by
  [`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- extendBy:

  Number of base pairs added to the width of each binding site before
  overlapping it with the methylation sites, split evenly between both
  sides. Default 130.

## Value

a `data.table` with one row per distance to the site centre.
