# normalise_footprint_flank

Normalise a corrected footprint against its outer flanking windows, on
the same scale as the statistic itself.

## Usage

``` r
normalise_footprint_flank(difference_data, method, flankNorm)
```

## Arguments

- difference_data:

  The corrected footprint as a `data.table`.

- method:

  Either `"substraction"` or `"division"`.

- flankNorm:

  Width of the flanking window used for normalisation.

## Value

The normalised `data.table`.

## Details

Dividing a difference by its flank mean would rescale it by an arbitrary
factor, because that mean is near zero: for a typical footprint it
inflates the curve by an order of magnitude and can invert its sign.
