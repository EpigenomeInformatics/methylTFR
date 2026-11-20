# plotExpectedFootprint

Plot the observed methylation profile around the binding sites of one
motif together with the profile expected from GC content.

## Usage

``` r
plotExpectedFootprint(
  motif,
  tf_bindsites,
  msites,
  sample_name = NULL,
  gc_dist,
  gcfreqs,
  enhancer = NULL,
  returnPlotData = FALSE
)
```

## Arguments

- motif:

  Name of the motif, one of `names(tf_bindsites)`.

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

- sample_name:

  Optional sample label used in the plot title. Default `"sample"`.

- gc_dist:

  A `GRanges` object of 30-bp genomic windows with an integer metadata
  column `GC_bin` (1 to 5). Get it with `getGenomeGC()` from the
  annotation package.

- gcfreqs:

  A named `list` with one element per motif, named as in `tf_bindsites`.
  Each element is a numeric matrix with the five GC bins in rows and the
  positions across the binding-site window in columns. Get it with
  `getGCfreq(motifSet = ...)` from the annotation package.

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- returnPlotData:

  Logical. If `TRUE`, return the plot data as well as the plot.

## Value

A `ggplot` object, or, if `returnPlotData` is `TRUE`, a list with the
elements `plot` and `plotDF`.

## See also

[`plotMotifFootprint`](https://epigenomeinformatics.github.io/methylTFR/reference/plotMotifFootprint.md)
for the bias-corrected profile,
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md)
for the score it summarises.

## Examples

``` r
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

p <- plotExpectedFootprint(
    motif = "BATF",
    tf_bindsites = tf_bindsites,
    msites = msites,
    sample_name = "Sample1",
    gc_dist = gcdist,
    gcfreqs = gcfreqs,
    enhancer = NULL,
    returnPlotData = FALSE
)
```
