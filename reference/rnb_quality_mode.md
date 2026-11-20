# rnb_quality_mode

Decide which per-site quality filter applies to an RnBeads object and
report it.

## Usage

``` r
rnb_quality_mode(rnb_set, dpval_threshold)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- dpval_threshold:

  numeric detection p-value threshold.

## Value

A list with the logical flags `has_covg` and `has_dpval`.
