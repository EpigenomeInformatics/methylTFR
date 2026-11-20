# calibrateDeviations

Standardise bias-corrected deviation scores against a within-sample
null, so that the resulting scores are comparable across motifs and
approximately standard normal in the absence of TF-specific signal.

## Usage

``` r
calibrateDeviations(devs, method = c("robust", "gaussian"))
```

## Arguments

- devs:

  A numeric matrix of bias-corrected deviation scores, with motifs in
  rows and samples in columns.

- method:

  Either `"robust"` (median and median absolute deviation, the default)
  or `"gaussian"` (mean and standard deviation).

## Value

A numeric matrix of the same dimensions as `devs`.

## Details

Each column (sample) is centred and scaled independently, using the
distribution of deviation scores across all motifs in that sample as the
null. With several hundred motifs, the great majority of which carry no
strong methylation signal in any single sample, this distribution is a
usable estimate of the sample's noise level.

This is deliberately different from `computeRowZScore`, which
standardises each motif against its own spread across samples. Row-wise
Z-scores are the right transform for visualising a heatmap, but they
force the standard deviation of every row to one and therefore cannot be
used to rank motifs by variability.
