# check_core_inputs

Validate the arguments shared by both methylTFR entry points.

## Usage

``` r
check_core_inputs(sample_ids, msites_fun, samples)
```

## Arguments

- sample_ids:

  A character vector of sample identifiers.

- msites_fun:

  A function of a single integer sample index.

- samples:

  A `data.frame` with one row per sample.

## Value

Invisible `NULL`. Called for the errors it raises.
