# read_methylome

read_methylome is a function to import methylation data into a GRanges
object. Bed file should be in EPP, ALLC or BisSNP format.

## Usage

``` r
read_methylome(filename, type, cov_threshold = 1)
```

## Arguments

- filename:

  filename which contains methylation data

- type:

  Type of file format. Currently supported epp, bissnp, bismarkCytosine,
  bismarkcov, allc and encode

- cov_threshold:

  numeric, coverage threshold to filter out low coverage sites, default
  is 1

## Value

a `GenomicRange` object with methylation, coverage information

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
bissnp <- read_methylome(bissnp_path, "bissnp")
```
