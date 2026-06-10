# compute_fp

Compute footprint returns a data.table required to create a
transcription factor footprint with label.

## Usage

``` r
compute_fp(motif_name, tf_bindsites, msites, enhancer = NULL)
```

## Arguments

- motif_name:

  - motif name eg: GATA string

- tf_bindsites:

  - a GenomicRange object contains tf binding sites positions from
    (`methylTFRann`)

- msites:

  - methylation data processed from `RnBeads`

- enhancer:

  - Specific regions such as distal motif, proximal motif

## Value

a `data.table` object to plot tf footprint
