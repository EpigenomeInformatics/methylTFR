# run_methylTFR_RnBeads

Run the methylTFR workflow directly on a preprocessed RnBeads object,
without exporting per-sample BED files first.

## Usage

``` r
run_methylTFR_RnBeads(
  rnb_set,
  tf_bindsites = NULL,
  gcfreqs = NULL,
  gc_dist = NULL,
  chunkSize = 20,
  threads = 1,
  enhancer = NULL,
  ignoreStrand = TRUE,
  cov_threshold = 1,
  dpval_threshold = 0.05,
  sample_ann = NULL
)
```

## Arguments

- rnb_set:

  A preprocessed `RnBSet` object.

- tf_bindsites:

  a `GRangesList` of TF binding site positions.

- gcfreqs:

  a `list` of GC bin frequency tables.

- gc_dist:

  a `GRanges` of the genome-wide GC distribution.

- chunkSize:

  Chunk size for parallel processing of motifs.

- threads:

  Thread count for parallel processing.

- enhancer:

  an optional `GRanges` of regions to restrict to.

- ignoreStrand:

  if TRUE, strand information is ignored.

- cov_threshold:

  numeric, minimum coverage of a retained site.

- dpval_threshold:

  numeric, maximum detection p-value of a retained probe.

- sample_ann:

  Optional `data.frame` of sample annotation.

## Value

a `methylTFRdeviations` object with bias-corrected deviations and
Z-scores.

## Details

Methylation calls are read at single-cytosine resolution. Sequencing
sets (`RnBiseqSet`) are filtered by coverage, array sets (`RnBeadSet`)
by detection p-value. The number of sites retained per sample is
reported through logger.

## See also

[`run_methyltfr`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methyltfr.md)
for running methylTFR from per-sample BED files.

## Examples

``` r
# A minimal end-to-end run on the BATF example data bundled with the
# package. The bundled calls are wrapped in an RnBiseqSet so that the
# example exercises the same code path as a preprocessed RnBeads object.
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

if (requireNamespace("RnBeads", quietly = TRUE) &&
    requireNamespace("RnBeads.hg38", quietly = TRUE)) {
    sites <- data.frame(
        chr = as.character(GenomicRanges::seqnames(msites)),
        start = GenomicRanges::start(msites),
        strand = "*",
        stringsAsFactors = FALSE
    )
    rnb_set <- RnBeads::RnBiseqSet(
        pheno = data.frame(sampleName = "sample_1"),
        sites = sites,
        meth = matrix(msites$score, ncol = 1),
        covg = matrix(msites$coverage, ncol = 1),
        assembly = "hg38",
        summarize.regions = FALSE
    )

    devs <- run_methylTFR_RnBeads(
        rnb_set = rnb_set,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist
    )
    deviations(devs)
}
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> INFO [2026-08-31 11:37:41] Annotation target: sites | assembly: hg38
#> INFO [2026-08-31 11:37:44] Found 534 sites across 1 samples
#> INFO [2026-08-31 11:37:44] Initializing the temp sink: methylTFR_tmp/methylTFR1b555968f376.h5
#> INFO [2026-08-31 11:37:44] Initializing the temp sink: methylTFR_tmp/methylTFR1b55139fc5bc.h5
#> INFO [2026-08-31 11:37:44] Sample 1: 534 of 534 sites retained (100%)
#> INFO [2026-08-31 11:37:44] Processing 1
#> INFO [2026-08-31 11:37:48] Finished processing 1
#> SUCCESS [2026-08-31 11:37:48] Computed all deviations successfully
#>             1
#> BATF 2.009857
```
