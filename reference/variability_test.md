# variability_test

Summarise the calibrated Z-scores into a variability per motif and test
each against the chi-squared null.

## Usage

``` r
variability_test(z, motif_names, padjMethod)
```

## Arguments

- z:

  The calibrated Z-score matrix.

- motif_names:

  Character vector of motif names.

- padjMethod:

  Multiple testing correction passed to
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A `data.frame` with one row per motif.
