# methyltfr_core

Internal engine shared by
[`runMethylTFR`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFR.md)
and
[`runMethylTFRFromRnBeads`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFRFromRnBeads.md).
It validates the motif set, allocates the on-disk sinks, iterates over
samples, computes per-motif deviations in chunks and assembles the
resulting `methylTFRdeviations` object.

The only difference between the two public entry points is where the
per-sample methylation calls come from. That difference is isolated in
the `msites_fun` argument, so both entry points share identical
numerical behaviour.

## Usage

``` r
methyltfr_core(
  sample_ids,
  msites_fun,
  samples,
  tf_bindsites,
  gcfreqs,
  gc_dist,
  chunkSize = 20,
  threads = 1,
  enhancer = NULL,
  ignoreStrand = TRUE,
  tempDir = file.path(tempdir(), "methylTFR")
)
```

## Arguments

- sample_ids:

  A character vector of sample identifiers. Used for the column names of
  the resulting object and to size the sinks.

- msites_fun:

  A function of a single integer `i` returning the methylation sites of
  sample `i` as a `GRanges` with the metadata columns `score` and
  `coverage`.

- samples:

  A `data.frame` of sample annotation with one row per entry of
  `sample_ids`, used as `colData`.

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

- gc_dist:

  A `GRanges` object of 30-bp genomic windows with an integer metadata
  column `GC_bin` (1 to 5). Get it with `getGenomeGC()` from the
  annotation package.

- chunkSize:

  Number of motifs processed together in one block. Default 20.

- threads:

  Number of workers used to process the motifs of a chunk in parallel.
  Default 1 (serial).

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

- tempDir:

  Directory for the temporary HDF5 files written during a run. Default
  is a folder under [`tempdir()`](https://rdrr.io/r/base/tempfile.html).
  The files are deleted when the run finishes.

## Value

a `methylTFRdeviations` object with bias-corrected deviations, row-wise
Z-scores and expected deviations.
