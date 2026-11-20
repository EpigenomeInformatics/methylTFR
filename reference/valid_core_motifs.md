# valid_core_motifs

Drop motifs whose binding sites are empty or whose GC bin frequency
matrix is missing.

## Usage

``` r
valid_core_motifs(tf_bindsites, gcfreqs)
```

## Arguments

- tf_bindsites:

  A `GRangesList` named by motif. Each element is a `GRanges` of the
  occurrences of that motif. The names must include every name in
  `names(gcfreqs)`. Get it with
  `getTFbindsites(motifSet = "jaspar2020")` from the
  methylTFRAnnotationHg38 (or methylTFRAnnotationMm10) package, or build
  a custom set with methylTFRAnnotationBuilder.

- gcfreqs:

  A named `list` with one element per motif, named as in `tf_bindsites`.
  Each element is a numeric matrix with the five GC bins in rows and the
  positions across the binding-site window in columns. Get it with
  `getGCfreq(motifSet = ...)` from the annotation package.

## Value

A character vector of the motif names that can be processed.
