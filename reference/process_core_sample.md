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

  a `GRangesList` of TF binding site positions.

- gcfreqs:

  a `list` of GC bin frequency tables.

- gc_dist:

  a `GRanges` of the genome-wide GC distribution.

- dev_grid, exp_grid:

  The grids the blocks are written on.

- dev_sink, exp_sink:

  The sinks the blocks are written to.

- BPPARAM:

  A `BiocParallelParam` object.

- enhancer:

  a `GRanges` restricting the analysis (optional).

- ignoreStrand:

  if TRUE, strand information is ignored.

## Value

Invisible `NULL`. Called for its effect on the sinks.
