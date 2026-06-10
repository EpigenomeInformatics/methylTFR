# plot_footprint

Creates a footprint plot for given motifs and methylation site. It will
create png files in the specified location.

## Usage

``` r
plot_footprint(motifs, tf_bindsites, samples, img = "TF_footprint.png")
```

## Arguments

- motifs:

  - motif names as character vector

- tf_bindsites:

  - Transcript Factor binbing sites from (`methylTFRann`)

- samples:

  - Import methylation data

- img:

  - output plot png filename

## Value

image will be generated in the specified path
