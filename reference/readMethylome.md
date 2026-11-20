# readMethylome

Import single-sample methylation calls into a `GRanges` object.

## Usage

``` r
readMethylome(filename, type, cov_threshold = 1)
```

## Arguments

- filename:

  Path to the methylation call file. Gzipped files are read directly.

- type:

  Format of the file, one of `"epp"`, `"bissnp"`, `"bismarkCytosine"`,
  `"bismarkCov"`, `"allc"` or `"encode"` (ENCODE bedMethyl).
  Case-insensitive.

- cov_threshold:

  Numeric. Sites with fewer reads than this are dropped. Default 1.

## Value

A `GRanges` object with one range per CpG and the metadata columns
`score` (methylation fraction between 0 and 1) and `coverage` (read
count). This is the `msites` input of
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md),
[`addGCBinToMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/addGCBinToMethylome.md)
and the footprint plots.

## See also

[`runMethylTFR`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFR.md),
which calls this function for every sample.

## Author

Irem Gunduz

## Examples

``` r
# Read bissnp file
bissnp_path <- system.file(
    "extdata",
    "bissnp.tsv.gz",
    package = "methylTFR"
)
bissnp <- readMethylome(bissnp_path, "bissnp")
```
