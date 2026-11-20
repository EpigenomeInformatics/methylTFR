# check_variability_inputs

Validate the inputs of `computeZScoreVariability` and return the
deviation scores as a matrix.

## Usage

``` r
check_variability_inputs(object, bootstrap, conf_level)
```

## Arguments

- object:

  A `methylTFRdeviations` object, matrix or data.frame.

- bootstrap:

  Logical, whether bootstrap bounds were requested.

- conf_level:

  Confidence level of the bootstrap bounds.

## Value

The deviation scores as a numeric matrix.
