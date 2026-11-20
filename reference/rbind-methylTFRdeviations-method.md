# rbind

Combine methylTFRdeviations objects by row.

## Usage

``` r
# S4 method for class 'methylTFRdeviations'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  methylTFRdeviations objects.

- deparse.level:

  passed to the underlying method.

## Value

A methylTFRdeviations object.

## Details

See [`cbind`](https://rdrr.io/pkg/data.table/man/cbindlist.html) for why
this is defined explicitly rather than inherited.

## Examples

``` r
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
devs <- rbind(tc_mem[seq_len(2), ], tc_mem[3:4, ])
```
