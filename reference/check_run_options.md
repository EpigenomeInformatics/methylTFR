# check_run_options

Validate the run-time options shared by both methylTFR entry points,
substituting documented defaults where a value is unusable.

## Usage

``` r
check_run_options(
  chunkSize = 20,
  threads = 1,
  ignoreStrand = TRUE,
  cov_threshold = 1
)
```

## Arguments

- chunkSize:

  Number of motifs processed together in one block. Default 20.

- threads:

  Number of workers used to process the motifs of a chunk in parallel.
  Default 1 (serial).

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

- cov_threshold:

  Numeric. Sites with fewer reads than this are dropped. Default 1.

## Value

A named list with the validated values.
