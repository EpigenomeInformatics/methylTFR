# computeZScoreVariability

Identify transcription factor motifs whose methylTFR activity varies
across samples more than expected by chance.

## Usage

``` r
computeZScoreVariability(
  object,
  method = c("robust", "gaussian"),
  bootstrap = FALSE,
  niterations = 1000L,
  conf_level = 0.95,
  padjMethod = "BH"
)
```

## Arguments

- object:

  A `methylTFRdeviations` object, or a numeric matrix or data.frame of
  bias-corrected deviation scores with motifs in rows and samples in
  columns.

- method:

  Either `"robust"` (median and median absolute deviation, the default)
  or `"gaussian"` (mean and standard deviation), used to estimate the
  within-sample null.

- bootstrap:

  logical, if TRUE bootstrap confidence bounds for the variability
  estimates are computed by resampling samples with replacement.

- niterations:

  Number of bootstrap iterations, used only when `bootstrap` is TRUE.

- conf_level:

  Confidence level for the bootstrap bounds.

- padjMethod:

  Method for p-value adjustment, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default is "BH".

## Value

A `data.frame` with one row per motif and the columns `motifs`,
`variability`, `p_value` and `p_value_adjusted`, plus
`bootstrap_lower_bound` and `bootstrap_upper_bound` when `bootstrap` is
TRUE. Rows are returned in the order of the input.

## Details

Variability is defined as the standard deviation, across samples, of
calibrated deviation scores. Calibration is what makes the number
interpretable: deviation scores are first standardised against a
within-sample null estimated across all motifs (see
`calibrateDeviations`), so that a calibrated score of 1 corresponds to
one unit of sample-level noise. A motif with variability greater than 1
therefore varies across samples by more than the background spread, and
variability near 1 is what a motif carrying no signal is expected to
show.

Under the null hypothesis that a motif's calibrated scores are
independent draws from a standard normal distribution, \\(n - 1) s^2\\
follows a chi-squared distribution with \\n - 1\\ degrees of freedom,
where \\s\\ is the observed variability and \\n\\ the number of samples
with a non-missing score for that motif. P-values are obtained from the
upper tail of that distribution and adjusted for multiple testing across
motifs.

Note that this function must not be applied to the `z` assay returned by
[`deviationZScores`](https://epigenomeinformatics.github.io/methylTFR/reference/deviationZScores.md).
Those are row-wise Z-scores computed for visualisation, and their
standard deviation across samples is 1 for every motif by construction,
which would make the test degenerate. The default input is therefore the
`deviations` assay.

The approach is analogous to `chromVAR::computeVariability`, with one
difference: chromVAR calibrates each motif against a set of GC- and
width-matched background peak sets, whereas here the null is estimated
across motifs within each sample. The latter needs no additional
deviation computations and can be applied to existing results, at the
cost of assuming that most motifs are inactive in any given sample.

## See also

[`differential_deviation_test`](https://epigenomeinformatics.github.io/methylTFR/reference/differential_deviation_test.md)
for testing activity differences between predefined groups.

## Author

Irem Gunduz

## Examples

``` r
# Load example data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
devs <- cbind(tc_mem, tc_naive)

# Rank motifs by how much their activity varies across samples
var_res <- computeZScoreVariability(devs)
#> Only 10 motifs supplied. The within-sample null is estimated across motifs, so variability estimates from small motif sets should be treated as indicative only.
head(var_res[order(-var_res$variability), ])
#>         motifs variability      p_value p_value_adjusted
#> 5     MAX::MYC   2.8900591 1.461439e-12     1.461439e-11
#> 3         IRF2   1.5177307 1.389717e-02     6.948586e-02
#> 1        FOXF2   0.7950911 7.705443e-01     9.999226e-01
#> 2        FOXD1   0.7076269 8.750259e-01     9.999226e-01
#> 10 RORA(var.2)   0.6550434 9.202746e-01     9.999226e-01
#> 7         PAX6   0.6145712 9.463439e-01     9.999226e-01
```
