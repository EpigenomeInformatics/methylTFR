# methyltfr_core

Internal engine shared by
[`run_methyltfr`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methyltfr.md)
and
[`run_methylTFR_RnBeads`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methylTFR_RnBeads.md).
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
  ignoreStrand = TRUE
)
```

## Arguments

- sample_ids:

  A character vector of sample identifiers. Used for the column names of
  the resulting object and to size the sinks.

- msites_fun:

  A function of a single integer `i` returning a `GRanges` object of
  methylation calls for sample `i`, with a numeric `score` metadata
  column holding methylation levels in `[0, 1]`.

- samples:

  A `data.frame` of sample annotation with one row per entry of
  `sample_ids`, used as `colData`.

- tf_bindsites:

  a `GRangesList` object containing TF binding site positions.

- gcfreqs:

  a `list` of GC bin frequency tables.

- gc_dist:

  a `GRanges` object containing the genome-wide GC distribution.

- chunkSize:

  Chunk size for parallel processing of motifs.

- threads:

  Thread count for parallel processing.

- enhancer:

  a `GRanges` object restricting the analysis to a set of regions such
  as distal regulatory elements (optional).

- ignoreStrand:

  if TRUE, strand information is ignored.

## Value

a `methylTFRdeviations` object with bias-corrected deviations, row-wise
Z-scores and expected deviations.
