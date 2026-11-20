# rnb_annotation_by_id

Subset a genome-wide annotation track to the probes of an RnBeads
object, matching on probe identifier.

## Usage

``` r
rnb_annotation_by_id(rnb_set, target, assembly)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- target:

  Character scalar naming the annotation target.

- assembly:

  Character scalar naming the genome assembly.

## Value

A `data.frame` of annotation rows, or `NULL` when the identifiers cannot
be matched.

## Details

The `sites` slot holds a three-column index matrix, not a row index into
the annotation, so the track is matched by identifier rather than by
position.
