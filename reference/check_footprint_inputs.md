# check_footprint_inputs

Validate the inputs shared by the two footprint plotting functions.

## Usage

``` r
check_footprint_inputs(
  motif,
  tf_bindsites,
  msites,
  gc_dist,
  gcfreqs,
  enhancer = NULL
)
```

## Arguments

- motif:

  Motif name as a character string.

- tf_bindsites:

  a `GRangesList` of TF binding site positions.

- msites:

  Methylation sites as a `GRanges` object.

- gc_dist:

  a `GRanges` of the genome-wide GC distribution.

- gcfreqs:

  a `list` of GC bin frequency tables.

- enhancer:

  a `GRanges` restricting the analysis (optional).

## Value

Invisible `NULL`. Called for the errors it raises.
