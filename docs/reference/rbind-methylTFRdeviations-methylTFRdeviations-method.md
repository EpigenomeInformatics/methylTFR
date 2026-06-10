# rbind

Combine two methylTFRdeviations objects by row.

## Usage

``` r
# S4 method for class 'methylTFRdeviations,methylTFRdeviations'
rbind(x, y)
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
# Combine the two objects with same column
devs <- rbind(tc_mem[1:2, 1:2], tc_mem[3:4, 1:2])
```
