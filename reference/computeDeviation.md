# computeDeviation

Compute the bias-corrected deviation score of one motif in one sample.

## Usage

``` r
computeDeviation(
  motif,
  msites,
  tf_bindsites,
  gcfreqs,
  enhancer = NULL,
  ignoreStrand = TRUE,
  binMsites,
  extendBy = 130,
  windowBreaks = c(-250, -200, -25, 25, 200, 250)
)
```

## Arguments

- motif:

  Name of the motif, one of `names(tf_bindsites)`.

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

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

- binMsites:

  The matrix returned by
  [`addGCBinToMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/addGCBinToMethylome.md),
  with the columns `gcbin` and `avg_mscore`.

- extendBy:

  Number of base pairs added to the width of each binding site before
  overlapping it with the methylation sites, split evenly between both
  sides. Default 130.

- windowBreaks:

  Six strictly increasing distances from the site centre, in base pairs,
  that define the five distance bins. Bins 1 and 5 are the outer flanks
  and bin 3 is the centre. Default `c(-250, -200, -25, 25, 200, 250)`.

## Value

A `data.table` with one row and the columns `dev` (bias-corrected
deviation) and `exp_dev` (expected deviation).

## Details

Each binding site is widened by `extendBy` base pairs and every CpG
inside it is placed by its distance to the site centre. The distances
are cut into five bins at `windowBreaks`: an outer left flank, an inner
left region, the centre, an inner right region and an outer right flank.
The observed score is the mean methylation in the centre bin divided by
the mean of the two outer flank bins. The same ratio is computed for the
methylation expected from GC content (see
[`addGCBinToMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/addGCBinToMethylome.md)),
and the deviation is the observed minus the expected ratio. If the
centre or either outer flank bin holds no CpG, the score is `NA`.

Every CpG counts once, whatever its read coverage: the bin means are
unweighted means of the per-CpG methylation fractions. Coverage only
enters through the `cov_threshold` filter applied when the data are
read.

## See also

[`runMethylTFR`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFR.md)
to score all motifs and samples at once,
[`plotMotifFootprint`](https://epigenomeinformatics.github.io/methylTFR/reference/plotMotifFootprint.md)
to see the profile behind the score.

## Examples

``` r
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

bin_meth <- addGCBinToMethylome(msites, gcdist, TRUE)

computeDeviation(
    motif = "BATF",
    msites = msites,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    enhancer = NULL,
    ignoreStrand = TRUE,
    binMsites = bin_meth
)
#>           dev   exp_dev
#>         <num>     <num>
#> 1: -0.1167365 0.9992794
```
