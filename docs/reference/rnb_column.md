# rnb_column

Extract a single sample column from an RnBeads accessor, falling back to
full extraction on RnBeads versions that do not support column
subsetting.

## Usage

``` r
rnb_column(accessor, rnb_set, index)
```

## Arguments

- accessor:

  An RnBeads accessor function, either
  [`RnBeads::meth`](https://rdrr.io/pkg/RnBeads/man/meth-methods.html)
  or
  [`RnBeads::covg`](https://rdrr.io/pkg/RnBeads/man/covg-methods.html).

- rnb_set:

  An `RnBSet` object.

- index:

  Integer index of the sample to extract.

## Value

A numeric vector with one value per site.
