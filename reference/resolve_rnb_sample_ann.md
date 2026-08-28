# resolve_rnb_sample_ann

Default the sample annotation to the phenotype table of the RnBeads
object and check its shape.

## Usage

``` r
resolve_rnb_sample_ann(rnb_set, sample_ann, sample_ids)
```

## Arguments

- rnb_set:

  A preprocessed `RnBSet` object.

- sample_ann:

  Sample annotation supplied by the caller, or NULL.

- sample_ids:

  Character vector of sample identifiers.

## Value

The sample annotation as a `data.frame`.
