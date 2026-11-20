# addGCBinToMethylome

Assign each methylation site to the GC bin of the genomic window it
falls in, and compute the average methylation per GC bin.

## Usage

``` r
addGCBinToMethylome(msites, gcdist, ignoreStrand = TRUE)
```

## Arguments

- msites:

  A `GRanges` object with one range per CpG (width 1) and the metadata
  columns `score` (numeric methylation fraction between 0 and 1) and
  `coverage` (read count). Usually produced by
  [`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

- gcdist:

  A `GRanges` object of 30-bp genomic windows with an integer metadata
  column `GC_bin` (1 to 5). Get it with `getGenomeGC()` from the
  annotation package.

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

## Value

A `matrix` with one row per GC bin and the columns `gcbin` and
`avg_mscore`. Pass it as `binMsites` to
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md).

## Details

The result describes how methylation depends on local GC content in this
sample.
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md)
combines it with the per-motif GC profile in `gcfreqs` to obtain the
methylation expected at the binding sites from sequence composition
alone.

## See also

[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md),
[`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md)

## Author

Irem Gunduz

## Examples

``` r
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))

bin_meth <- addGCBinToMethylome(msites, gcdist)
#> Loading required namespace: GenomeInfoDb
bin_meth
#>      gcbin avg_mscore
#> [1,]     1  0.5562500
#> [2,]     2  0.7409375
#> [3,]     3  0.8235556
#> [4,]     4  0.7962617
#> [5,]     5  0.7528302
```
