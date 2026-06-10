# computeExpectations

This function is used to calculate expected methylation for a given
motif and sample.

## Usage

``` r
computeExpectations(binMsites, gcfreq)
```

## Arguments

- binMsites:

  Imported methylation sites with GC bin

- gcfreq:

  a `list` of GC bin frequency tables (matrices for multiple motif)

## Value

a `data.table` object with GC bin with corresponding avg methylation
