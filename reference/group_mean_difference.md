# group_mean_difference

Effect size per motif. For two groups, the mean of the first level minus
the mean of the second level (signed). For more groups, the range of the
group means (unsigned).

## Usage

``` r
group_mean_difference(deviations, groups)
```

## Arguments

- deviations:

  A matrix of deviation scores, motifs in rows.

- groups:

  The group labels as a `factor`.

## Value

A numeric vector with one value per motif.
