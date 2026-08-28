# resolve_diff_groups

Derive and validate the group labels used by
`differential_deviation_test`.

## Usage

``` r
resolve_diff_groups(deviations, groups)
```

## Arguments

- deviations:

  A matrix of deviation scores, motifs in rows.

- groups:

  Group labels, or NULL to take them from the column names.

## Value

The group labels as a `factor`.
