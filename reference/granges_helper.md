# granges_helper

helper function to convert msites in EPP GRanges format

## Usage

``` r
granges_helper(grobj, chr, mscore, cov, startP, endP)
```

## Arguments

- grobj:

  GRanges object

- chr:

  chromosome

- mscore:

  methylation score

- cov:

  methylation coverage

- startP:

  start position

- endP:

  end position

## Value

a `GenomicRanges` object with EPP format
