# footprint_difference

Combine the observed and expected profiles into a single corrected
curve.

## Usage

``` r
footprint_difference(plot_data, method, sample_name)
```

## Arguments

- plot_data:

  A `data.table` with the columns `x`, `avg_methyl` and `type`.

- method:

  Either `"substraction"` or `"division"`.

- sample_name:

  Sample label used in the curve label.

## Value

A `list` with the corrected `data` and the axis `lab`.
