# run_methylTFR_RnBeads

Run the methylTFR workflow directly on a preprocessed RnBeads object,
without exporting per-sample BED files first.

This is the RnBeads-based counterpart to
[`run_methyltfr`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methyltfr.md).
Both functions share the same engine and produce numerically identical
results for the same underlying methylation calls; they differ only in
where the per-sample methylation levels come from.

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
  sample_ann = NULL
)
```

## Arguments

- rnb_set:

  A preprocessed `RnBSet` object, for example the output of
  `rnb.run.preprocessing` or a set loaded with
  [`RnBeads::load.rnb.set`](https://rdrr.io/pkg/RnBeads/man/load.rnb.set.html).

- tf_bindsites:

  a `GRangesList` object contains tf binding sites positions

- gcfreqs:

  a `list` of GC bin frequency tables (matrices for multiple motif)

- gc_dist:

  a `GRanges` object contains Genome wide GC distribution

- chunkSize:

  Chunk size for parallel processing of motifs (default: 20)

- threads:

  Thread count for parallel processing

- enhancer:

  a `GRanges` object specifying regions such as distal regulatory
  elements (optional)

- ignoreStrand:

  if TRUE, it ignores strand info from annotation

- cov_threshold:

  numeric, coverage threshold used to filter out low coverage sites,
  default is 1. Ignored for objects without coverage information.

- sample_ann:

  Optional `data.frame` of sample annotation with one row per sample,
  used as `colData`. Defaults to `RnBeads::pheno(rnb_set)`.

## Value

a `methylTFRdeviations` object with bias-corrected deviation and
Z-scores

## Details

Methylation calls are always read at single-cytosine resolution
(`type = "sites"`). Region-level summaries such as `tiling1kb` or
`distal` cannot be used, because methylTFR needs base-resolution calls
to build the footprint around each motif centre. To restrict the
analysis to a set of regulatory regions, pass those regions through the
`enhancer` argument instead.

Samples are processed one at a time and methylation levels are pulled
from the RnBeads object column by column, so disk-backed (`ff`-managed)
RnBeads sets are never loaded into memory in full.

Coverage filtering is applied only when the object carries coverage
information, which is the case for sequencing-based sets (`RnBiseqSet`).
For array-based sets `cov_threshold` is ignored and a message is
emitted.

Note that RnBeads site annotation is 1-based while
[`read_methylome`](https://epigenomeinformatics.github.io/methylTFR/reference/read_methylome.md)
reads 0-based BED coordinates as-is. The resulting one-base offset is
not corrected here, since deviation scores aggregate methylation over
windows of tens to hundreds of bases and are insensitive to a uniform
single-base shift.

## See also

[`run_methyltfr`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methyltfr.md)
for the file-based entry point.

## Author

Irem Gunduz

## Examples

``` r
# A minimal end-to-end run on the BATF example data bundled with the
# package. RnBeads and its hg38 annotation build the input object; both
# are optional dependencies.
if (requireNamespace("RnBeads", quietly = TRUE) &&
    requireNamespace("RnBeads.hg38", quietly = TRUE)) {
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))
    load(system.file(
        "extdata", "BATF_tf_bindsites.rda",
        package = "methylTFR"
    ))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

    # RnBiseqSet() takes methylation as a fraction and coverage as counts,
    # with one column per sample.
    sites <- data.frame(
        chromosome = as.character(GenomicRanges::seqnames(msites)),
        position = GenomicRanges::start(msites),
        strand = "*",
        stringsAsFactors = FALSE
    )
    rnb_set <- RnBeads::RnBiseqSet(
        pheno = data.frame(
            sampleName = "sample_1", stringsAsFactors = FALSE
        ),
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
#> INFO [2026-08-28 14:10:36] Found 534 sites across 1 samples
#> INFO [2026-08-28 14:10:36] Initializing the temp sink: methylTFR_tmp/methylTFR1b165c966482.h5
#> INFO [2026-08-28 14:10:36] Initializing the temp sink: methylTFR_tmp/methylTFR1b163f6a7d2.h5
#> INFO [2026-08-28 14:10:36] Processing 1
#> INFO [2026-08-28 14:10:41] Finished processing 1
#> SUCCESS [2026-08-28 14:10:41] Computed all deviations successfully
#>             1
#> BATF 2.009857
```
