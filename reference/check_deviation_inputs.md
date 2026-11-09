# check_deviation_inputs

Validate the inputs of `computeDeviation`.

## Usage

``` r
check_deviation_inputs(
  motif,
  msites,
  tf_bindsites,
  gcfreqs,
  binMsites,
  enhancer = NULL,
  extendBy = 130,
  windowBreaks = c(-250, -200, -25, 25, 200, 250)
)
```

## Arguments

- motif:

  Motif name as a character string.

- msites:

  A `GRanges` object with one range per CpG (width 1) and the metadata
  columns `score` (numeric methylation fraction between 0 and 1) and
  `coverage` (read count). Usually produced by
  [`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

- tf_bindsites:

  A `GRangesList` named by motif. Each element is a `GRanges` of the
  occurrences of that motif. The names must include every name in
  `names(gcfreqs)`. Get it with
  `getTFbindsites(motifSet = "jaspar2020")` from the
  methylTFRAnnotationHg38 (or methylTFRAnnotationMm10) package, or build
  a custom set with methylTFRAnnotationBuilder.

- gcfreqs:

  A named `list` with one element per motif, named as in `tf_bindsites`.
  Each element is a numeric matrix with the five GC bins in rows and the
  positions across the binding-site window in columns. Get it with
  `getGCfreq(motifSet = ...)` from the annotation package.

- binMsites:

  The matrix returned by
  [`addGCBinToMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/addGCBinToMethylome.md),
  with the columns `gcbin` and `avg_mscore`.

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- extendBy:

  Number of base pairs added to each binding site.

- windowBreaks:

  Break points of the distance bins.

## Value

Invisible `NULL`. Called for the errors it raises.
