# rnb_sample_ids

Determine the sample identifiers of an RnBeads object.

## Usage

``` r
rnb_sample_ids(rnb_set)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

## Value

A character vector of sample identifiers.

## Details

RnBeads exports `samples` with `exportMethods` rather than `export`, so
the generic is not reachable as
[`RnBeads::samples`](https://rdrr.io/pkg/RnBeads/man/samples-methods.html)
and referring to it that way fails `R CMD check`. The identifiers are
therefore taken from the column names of the methylation matrix, falling
back to the row names of the sample annotation.
