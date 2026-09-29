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
    # identifiers.column is a session option, not a property of the object,
    # so it has to be set before the set is built or samples() falls back
    # to row numbers and methylTFR cannot name the columns it returns
    old_ids <- RnBeads::rnb.getOption("identifiers.column")
    RnBeads::rnb.options(identifiers.column = "sampleName")

    sites <- data.frame(
        chr = as.character(GenomicRanges::seqnames(msites)),
        start = GenomicRanges::start(msites),
        strand = "*",
        stringsAsFactors = FALSE
    )
    rnb_set <- RnBeads::RnBiseqSet(
        pheno = data.frame(sampleName = "sample_1"),
        sites = sites,
        meth = matrix(msites$score, ncol = 1,
                      dimnames = list(NULL, "sample_1")),
        covg = matrix(msites$coverage, ncol = 1,
                      dimnames = list(NULL, "sample_1")),
        assembly = "hg38",
        summarize.regions = FALSE
    )

    devs <- run_methylTFR_RnBeads(
        rnb_set = rnb_set,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist
    )
    RnBeads::rnb.options(identifiers.column = old_ids)
    deviations(devs)
}
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> INFO [2026-09-29 15:09:19] Annotation target: sites | assembly: hg38
#> INFO [2026-09-29 15:09:23] Found 534 sites across 1 samples
#> INFO [2026-09-29 15:09:23] Initializing the temp sink: methylTFR_tmp/methylTFR538850521060.h5
#> INFO [2026-09-29 15:09:23] Initializing the temp sink: methylTFR_tmp/methylTFR53883d749e32.h5
#> INFO [2026-09-29 15:09:23] Sample 1: 534 of 534 sites retained (100%)
#> INFO [2026-09-29 15:09:23] Processing sample_1
#> INFO [2026-09-29 15:09:26] Finished processing sample_1
#> SUCCESS [2026-09-29 15:09:26] Computed all deviations successfully
#>      sample_1
#> BATF 2.013891
```
