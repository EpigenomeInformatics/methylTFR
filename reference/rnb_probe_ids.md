# rnb_probe_ids

Read the probe identifiers of an RnBeads object.

## Usage

``` r
rnb_probe_ids(rnb_set)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

## Value

A character vector of probe identifiers, or `NULL`.

## Details

The `sites` slot of an `RnBSet` is an index matrix whose row names carry
the probe identifiers.
