# computeExpectedFootprint

Methylation profile expected from GC content around the binding sites of
one motif.

## Usage

``` r
computeExpectedFootprint(motif, gcfreqs, gc_dist, enhancer = NULL, msites)
```

## Arguments

- motif:

  Name of the motif, one of `names(tf_bindsites)`.

- gcfreqs:

  A named `list` with one element per motif, named as in `tf_bindsites`.
  Each element is a numeric matrix with the five GC bins in rows and the
  positions across the binding-site window in columns. Get it with
  `getGCfreq(motifSet = ...)` from the annotation package.

- gc_dist:

  A `GRanges` object of 30-bp genomic windows with an integer metadata
  column `GC_bin` (1 to 5). Get it with `getGenomeGC()` from the
  annotation package.

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- msites:

  A `GRanges` object with one range per CpG (width 1) and the metadata
  columns `score` (numeric methylation fraction between 0 and 1) and
  `coverage` (read count). Usually produced by
  [`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

## Value

a `data.table` with the columns `x` and `avg_methyl`.
