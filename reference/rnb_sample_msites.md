# rnb_sample_msites

Extract the methylation calls of a single sample from an RnBeads object
as a `GRanges` object.

## Usage

``` r
rnb_sample_msites(
  rnb_set,
  sites_gr,
  index,
  cov_threshold = 1,
  has_covg = TRUE,
  dpval_threshold = 0.05,
  has_dpval = FALSE
)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- sites_gr:

  A `GRanges` object of site positions.

- index:

  Integer index of the sample to extract.

- cov_threshold:

  numeric coverage threshold.

- has_covg:

  logical, whether coverage filtering applies.

- dpval_threshold:

  numeric detection p-value threshold.

- has_dpval:

  logical, whether detection p-value filtering applies.

## Value

A `GRanges` object restricted to valid methylation calls.
