# run_methyltfr

This function is a wrapper function to calculate the deviation in
transcription factor footprint base for all given motifs per raw samples

## Usage

``` r
run_methyltfr(
  sample_ann,
  sample_dir,
  tf_bindsites = NULL,
  gcfreqs = NULL,
  gc_dist = NULL,
  sampleColName = "bedFile",
  chunkSize = 20,
  full_path = FALSE,
  annfile = NULL,
  threads = 1,
  enhancer = NULL,
  filetype = NULL,
  ignoreStrand = TRUE,
  cov_threshold = 1
)
```

## Arguments

- sample_ann:

  A tab seperated file contains sample annotations

- sample_dir:

  The directory where all bed file and annotation file stored

- tf_bindsites:

  a `GRangesList` object contains tf binding sites positions

- gcfreqs:

  a `list` of GC bin frequency tables (matrices for multiple motif)

- gc_dist:

  a `GRanges` object contains Genome wide GC distribution

- sampleColName:

  column name of the sample bed file in the annotation file

- chunkSize:

  Chunk size for parallel processing of motifs (default: 20)

- full_path:

  if TRUE, the bed file path in the annotation file is full path

- annfile:

  if provided, the sample annotation file is not read from the
  sample_dir

- threads:

  Thread count for parallel processing

- enhancer:

  a `GRanges` object specifying regions such as distal motif (optional)

- filetype:

  file type of the bed file, currently supported:
  bissnp,epp,allc,bismarkcytosine,bismarkcov,encode

- ignoreStrand:

  if TRUE, it ignores strand info from annotation

- cov_threshold:

  - numeric, coverage threshold to filter out low coverage sites,
    default is 1

## Value

a `methylTFRdeviations` object with bias-corrected deviation and
Z-scores

## See also

[`run_methylTFR_RnBeads`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methylTFR_RnBeads.md)
for running methylTFR directly on a preprocessed RnBeads object.
