# check_deviation_inputs

Validate the inputs of `computeDeviation`.

## Usage

``` r
check_deviation_inputs(motif, msites, tf_bindsites, enhancer = NULL)
```

## Arguments

- motif:

  Motif name as a character string.

- msites:

  Methylation sites as a `GRanges` object.

- tf_bindsites:

  a `GRangesList` of TF binding site positions.

- enhancer:

  a `GRanges` restricting the analysis (optional).

## Value

Invisible `NULL`. Called for the errors it raises.
