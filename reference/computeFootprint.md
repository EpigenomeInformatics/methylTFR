# computeFootprint

Compute footprint returns a data.table required to create a
transcription factor footprint with label.

## Usage

``` r
computeFootprint(motif_name, tf_bindsites, msites, enhancer = NULL)
```

## Arguments

- motif_name:

  chracter vector containing motif name eg: GATA string

- tf_bindsites:

  TF binding sites positions from methTFRannotation package

- msites:

  Methylation data

- enhancer:

  Specific regions such as distal motif, proximal motif

## Value

a list containing plot.data and tfbs
