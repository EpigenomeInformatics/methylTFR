# valid_core_motifs

Drop motifs whose binding sites are empty or whose GC bin frequency
matrix is missing.

## Usage

``` r
valid_core_motifs(tf_bindsites, gcfreqs)
```

## Arguments

- tf_bindsites:

  a `GRangesList` of TF binding site positions.

- gcfreqs:

  a `list` of GC bin frequency tables.

## Value

A character vector of the motif names that can be processed.
