# check_annotation_inputs

Validate the annotation objects required by both methylTFR entry points.

## Usage

``` r
check_annotation_inputs(tf_bindsites, gcfreqs, gc_dist, enhancer = NULL)
```

## Arguments

- tf_bindsites:

  a `GRangesList` object of TF binding site positions.

- gcfreqs:

  a `list` of GC bin frequency tables.

- gc_dist:

  a `GRanges` object of the genome-wide GC distribution.

- enhancer:

  an optional `GRanges` object of regions to restrict to.

## Value

invisible TRUE, called for the side effect of signalling errors.
