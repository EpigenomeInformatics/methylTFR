# rnb_sample_msites

Extract the methylation calls of a single sample from an RnBeads object
as a `GRanges` object in the layout expected by
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md).

## Usage

``` r
rnb_sample_msites(rnb_set, sites_gr, index, cov_threshold = 1, has_covg = TRUE)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- sites_gr:

  A `GRanges` object of site positions, as returned by
  `rnb_sites_to_granges`.

- index:

  Integer index of the sample to extract.

- cov_threshold:

  numeric coverage threshold.

- has_covg:

  logical, whether the object carries coverage information.

## Value

A `GRanges` object with `score` and `coverage` metadata columns,
restricted to sites with a non-missing methylation call.
