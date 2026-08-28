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

## Examples

``` r
# A minimal end-to-end run on the BATF example data bundled with the
# package. The annotation objects cover a single motif, so the result
# has one row.
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

# run_methyltfr() reads per-sample calls from disk, so the bundled sites
# are written out as a bismarkCov file first.
sample_dir <- tempfile("methylTFR_example")
dir.create(sample_dir)
n_meth <- round(msites$score * msites$coverage)
write.table(
    data.frame(
        chr = as.character(GenomicRanges::seqnames(msites)),
        start = GenomicRanges::start(msites),
        end = GenomicRanges::end(msites),
        percent = msites$score * 100,
        meth = n_meth,
        unmeth = msites$coverage - n_meth
    ),
    file.path(sample_dir, "sample_1.cov"),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
write.table(
    data.frame(sampleName = "sample_1", bedFile = "sample_1.cov"),
    file.path(sample_dir, "samples.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

devs <- run_methyltfr(
    sample_ann = "samples.tsv",
    sample_dir = sample_dir,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gcdist,
    filetype = "bismarkcov"
)
#> SUCCESS [2026-08-28 14:10:41] The samples are successfully located
#> INFO [2026-08-28 14:10:42] Initializing the temp sink: methylTFR_tmp/methylTFR1b16e9154d5.h5
#> INFO [2026-08-28 14:10:42] Initializing the temp sink: methylTFR_tmp/methylTFR1b163e642a4e.h5
#> INFO [2026-08-28 14:10:42] Processing sample_1.cov
#> INFO [2026-08-28 14:10:47] Finished processing sample_1.cov
#> SUCCESS [2026-08-28 14:10:47] Computed all deviations successfully
deviations(devs)
#>      sample_1.cov
#> BATF     1.743674

unlink(sample_dir, recursive = TRUE)
```
