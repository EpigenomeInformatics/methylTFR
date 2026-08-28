# expected_footprint_plot

Draw the observed and expected methylation profiles.

## Usage

``` r
expected_footprint_plot(combined_data, motif, sample_name)
```

## Arguments

- combined_data:

  A `data.table` with the columns `x`, `avg_methyl` and `type`.

- motif:

  Motif name as a character string.

- sample_name:

  Sample label used in the title.

## Value

A `ggplot` object.
