# cbind

Combine methylTFRdeviations objects by column.

## Usage

``` r
# S4 method for class 'methylTFRdeviations'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  methylTFRdeviations objects.

- deparse.level:

  passed to the underlying method.

## Value

A methylTFRdeviations object.

## Details

Defined explicitly rather than left to inheritance. Although
`methylTFRdeviations` extends `SummarizedExperiment`, whether that
class's `cbind` method is found for a subclass defined in another
package depends on the search path and on the Bioconductor version:
dispatch can fall through to the `S4Vectors` method, which then fails
with "unable to find an inherited method for function 'bindCOLS'".

The objects are therefore coerced to `SummarizedExperiment` explicitly,
combined there, and re-wrapped. This is defined with `setMethod` on the
existing generic rather than `setGeneric`, so it does not shadow the
base function and any number of objects can be combined.

## Examples

``` r
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
devs <- cbind(tc_mem, tc_naive)
```
