# rnb_sites_to_granges

Build a `GRanges` object of the site annotation of an RnBeads object.
The order of the ranges matches the row order of the methylation matrix
returned by
[`RnBeads::meth`](https://rdrr.io/pkg/RnBeads/man/meth-methods.html).

## Usage

``` r
rnb_sites_to_granges(rnb_set, ignoreStrand = TRUE)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- ignoreStrand:

  if TRUE, all ranges are returned with strand `"*"`.

## Value

A `GRanges` object with one range per site.
