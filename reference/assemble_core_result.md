# assemble_core_result

Close the sinks and assemble the deviations, their row-wise Z-scores and
the expected deviations into a result object.

## Usage

``` r
assemble_core_result(dev_sink, exp_sink, samples)
```

## Arguments

- dev_sink:

  The sink holding the bias-corrected deviations.

- exp_sink:

  The sink holding the expected deviations.

- samples:

  A `data.frame` with one row per sample.

## Value

a `methylTFRdeviations` object.
