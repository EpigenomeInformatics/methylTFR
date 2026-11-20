# rnb_sites_to_granges

Build a `GRanges` object of the site or probe annotation of an RnBeads
object.

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

A `GRanges` object with one range per site or probe.
