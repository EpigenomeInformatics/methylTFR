# calculate_expmeth

This function is used to calculate genome-wide expected methylation for
each motif.

## Usage

``` r
calculate_expmeth(msites, gcdist, gcfreq)
```

## Arguments

- msites:

  - imported methylation sites

- gcdist:

  - Genome wide GC distribution from (`methylTFRann`)

- gcfreq:

  - GC bin frequency table (matrix) from (`methylTFRann`)

## Value

a `data.table` object with GC bin with corresponding avg methylation
