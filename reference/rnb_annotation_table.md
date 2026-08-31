# rnb_annotation_table

Look up the site or probe annotation of an RnBeads object.

## Usage

``` r
rnb_annotation_table(rnb_set, target, assembly)
```

## Arguments

- rnb_set:

  An `RnBSet` object.

- target:

  Character scalar naming the annotation target.

- assembly:

  Character scalar naming the genome assembly.

## Value

A `data.frame` with one row per site or probe.

## Details

The annotation stored in the object is used when available, otherwise
the genome-wide track registered for `target` is matched to the probes
of the object by identifier.
