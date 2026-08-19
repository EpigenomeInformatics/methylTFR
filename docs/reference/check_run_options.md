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

  Chunk size for parallel processing of motifs.

- threads:

  Thread count for parallel processing.

- ignoreStrand:

  if TRUE, strand information is ignored.

- cov_threshold:

  numeric coverage threshold.

## Value

A named list with the validated values.
