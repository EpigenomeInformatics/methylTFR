# compute_deviation

compute_deviation is a function to calculate the deviation in
transcription factor

## Usage

``` r
compute_deviation(
  motif,
  msites,
  tf_bindsites,
  gcfreqs,
  gcdist,
  enhancer = NULL
)
```

## Arguments

- motif:

  - motif name

- msites:

  - imported methylation sites

- tf_bindsites:

  - a GenomicRange object contains tf binding sites positions from
    (`methylTFRann`)

- gcfreqs:

  - GC bin frequency tables (matrices for multiple motif) from
    (`methylTFRann`)

- gcdist:

  - Genome wide GC distribution from (`methylTFRann`)

- enhancer:

  - Specific regions like distal motif

## Value

a `numeric` deviation score for a given motif
