# computeExpectations

Expected methylation profile of a motif in one sample.

## Usage

``` r
computeExpectations(binMsites, gcfreq)
```

## Arguments

- binMsites:

  The matrix returned by `addGCBinToMethylome`.

- gcfreq:

  The GC bin frequency matrix of one motif.

## Value

a `data.table` with the columns `x` (distance to the site centre) and
`avg_methyl`.
