# diff_pvalues

Test every motif for a difference between the groups, with the test
chosen from the number of groups and `parametric`.

## Usage

``` r
diff_pvalues(deviations, groups, parametric, alternative)
```

## Arguments

- deviations:

  A matrix of deviation scores, motifs in rows.

- groups:

  The group labels as a `factor`.

- parametric:

  if TRUE, use a t-test or ANOVA, otherwise a Wilcoxon or Kruskal-Wallis
  test.

- alternative:

  The alternative hypothesis of the two-group tests.

## Value

A numeric vector of p-values, one per motif.
