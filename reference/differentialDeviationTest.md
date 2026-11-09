# differentialDeviationTest

Test which motifs have different deviation scores between groups of
samples. Adapted from `chromVAR::differentialDeviations`.

## Usage

``` r
differentialDeviationTest(
  deviations,
  groups = NULL,
  motifs = rownames(deviations),
  alternative = c("two.sided", "less", "greater"),
  parametric = TRUE,
  padjMethod = "BH"
)
```

## Arguments

- deviations:

  A
  [`methylTFRdeviations`](https://epigenomeinformatics.github.io/methylTFR/reference/methylTFRdeviations-class.md)
  object, or a matrix or `data.frame` of deviation scores with motifs in
  rows and samples in columns.

- groups:

  A vector with one group label per sample. If `NULL`, the column names
  of `deviations` are used. For two groups, the factor levels set the
  direction of `mean_difference`.

- motifs:

  Motif names, by default `rownames(deviations)`.

- alternative:

  The alternative hypothesis of the two-group tests, one of
  `"two.sided"`, `"less"` or `"greater"`.

- parametric:

  Logical. If `TRUE` (default), parametric tests are used, otherwise
  non-parametric tests.

- padjMethod:

  Method for p-value adjustment, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

## Value

A `data.frame` with one row per motif and the columns `motifs`,
`p_value`, `p_value_adjusted` and `mean_difference`. For two groups,
`mean_difference` is the mean of the first group minus the mean of the
second; for more groups, it is the range of the group means.

## Details

With two groups, a Welch t-test (`parametric = TRUE`) or a Wilcoxon rank
sum test is used; with more groups, a one-way ANOVA or a Kruskal-Wallis
test.

## See also

[`computeZScoreVariability`](https://epigenomeinformatics.github.io/methylTFR/reference/computeZScoreVariability.md)
for ranking motifs without group labels.

## Examples

``` r
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
devs <- cbind(tc_mem, tc_naive)

# Construct group labels from sample names
groups <- sub("\\.bedGraph.*$", "", vapply(
    strsplit(colnames(devs), split = "_"), `[`, character(1), 1L
))

tc_result <- differentialDeviationTest(
    deviations = devs,
    groups = groups,
    alternative = "two.sided",
    parametric = TRUE,
    padjMethod = "BH"
)
head(tc_result)
#>                  motifs      p_value p_value_adjusted mean_difference
#> FOXF2             FOXF2 1.694736e-02     0.0355317162      0.02920759
#> FOXD1             FOXD1 3.644010e-01     0.4555013119     -0.01081636
#> IRF2               IRF2 1.873667e-03     0.0093683362     -0.03414276
#> MZF1(var.2) MZF1(var.2) 1.776586e-02     0.0355317162     -0.01850831
#> MAX::MYC       MAX::MYC 1.078111e-05     0.0001078111     -0.09341429
#> PPARG             PPARG 3.996324e-03     0.0133210805     -0.03428684
```
