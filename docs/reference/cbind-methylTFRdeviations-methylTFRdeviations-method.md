# cbind

Combine two methylTFRdeviations objects by column.

## Usage

``` r
# S4 method for class 'methylTFRdeviations,methylTFRdeviations'
cbind(x, y)
```

## Arguments

- x:

  methylTFRdeviations object

- y:

  methylTFRdeviations object

## Value

A methylTFRdeviations object.

## Examples

``` r
# Load example data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
# Combine the two objects
devs <- cbind(tc_mem, tc_naive)
```
