# process_core_sample

Compute and write the deviations of one sample, one motif chunk at a
time.

## Usage

``` r
process_core_sample(
  index,
  sample_ids,
  msites_fun,
  motif_chunks,
  tf_bindsites,
  gcfreqs,
  gc_dist,
  dev_grid,
  exp_grid,
  dev_sink,
  exp_sink,
  BPPARAM,
  enhancer,
  ignoreStrand
)
```

## Arguments

- index:

  Integer index of the sample within `sample_ids`.

- sample_ids:

  A character vector of sample identifiers.

- msites_fun:

  A function of a single integer sample index.

- motif_chunks:

  A `list` of character vectors of motif names.

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

- dev_grid, exp_grid:

  The grids the blocks are written on.

- dev_sink, exp_sink:

  The sinks the blocks are written to.

- BPPARAM:

  A `BiocParallelParam` object.

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

## Value

Invisible `NULL`. Called for its effect on the sinks.
