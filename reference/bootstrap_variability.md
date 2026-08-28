# bootstrap_variability

Add bootstrap confidence bounds to a variability table.

## Usage

``` r
bootstrap_variability(z, res, niterations, conf_level)
```

## Arguments

- z:

  The calibrated Z-score matrix.

- res:

  The variability table to add the bounds to.

- niterations:

  Number of bootstrap iterations.

- conf_level:

  Confidence level of the bounds.

## Value

`res` with the two bound columns added.
