# runMethylTFRFromRnBeads

Run the methylTFR workflow directly on a preprocessed RnBeads object,
without exporting per-sample BED files first.

## Usage

``` r
runMethylTFRFromRnBeads(
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
  sample_ann = NULL,
  tempDir = file.path(tempdir(), "methylTFR")
)
```

## Arguments

- rnb_set:

  A preprocessed `RnBSet` object, for example loaded with
  [`RnBeads::load.rnb.set()`](https://rdrr.io/pkg/RnBeads/man/load.rnb.set.html).

- tf_bindsites:

  A `GRangesList` named by motif. Each element is a `GRanges` of the
  occurrences of that motif. The names must include every name in
  `names(gcfreqs)`. Get it with
  `getTFbindsites(motifSet = "jaspar2020")` from the
  methylTFRAnnotationHg38 (or methylTFRAnnotationMm10) package, or build
  a custom set with methylTFRAnnotationBuilder.

- gcfreqs:

  A named `list` with one element per motif, named as in `tf_bindsites`.
  Each element is a numeric matrix with the five GC bins in rows and the
  positions across the binding-site window in columns. Get it with
  `getGCfreq(motifSet = ...)` from the annotation package.

- gc_dist:

  A `GRanges` object of 30-bp genomic windows with an integer metadata
  column `GC_bin` (1 to 5). Get it with `getGenomeGC()` from the
  annotation package.

- chunkSize:

  Number of motifs processed together in one block. Default 20.

- threads:

  Number of workers used to process the motifs of a chunk in parallel.
  Default 1 (serial).

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

- cov_threshold:

  Numeric, minimum coverage of a retained site. Used only for sequencing
  sets. Default 1.

- dpval_threshold:

  Numeric, maximum detection p-value of a retained probe. Used only for
  array sets. Default 0.05.

- sample_ann:

  Optional `data.frame` with one row per sample, in the same order as
  the samples in `rnb_set`. All columns are kept as `colData` of the
  result; no column is required. Default `NULL` uses
  `RnBeads::pheno(rnb_set)`.

- tempDir:

  Directory for the temporary HDF5 files written during a run. Default
  is a folder under [`tempdir()`](https://rdrr.io/r/base/tempfile.html).
  The files are deleted when the run finishes.

## Value

A
[`methylTFRdeviations`](https://epigenomeinformatics.github.io/methylTFR/reference/methylTFRdeviations-class.md)
object with the assays `deviations`, `z` and `expected`, motifs in rows
and samples in columns.

## Details

Methylation calls are read at single-cytosine resolution. Sequencing
sets (`RnBiseqSet`) are filtered by coverage, array sets (`RnBeadSet`)
by detection p-value. The number of sites retained per sample is
reported through logger.

## See also

[`runMethylTFR`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFR.md)
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

    devs <- runMethylTFRFromRnBeads(
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
#> INFO [2026-11-09 15:29:52] Annotation target: sites | assembly: hg38
#> INFO [2026-11-09 15:29:56] Found 516 sites across 1 samples
#> INFO [2026-11-09 15:29:56] Initializing the temp sink: /tmp/Rtmp8oDvL1/methylTFR/methylTFR1971273d4af5.h5
#> INFO [2026-11-09 15:29:56] Initializing the temp sink: /tmp/Rtmp8oDvL1/methylTFR/methylTFR19714c425a19.h5
#> INFO [2026-11-09 15:29:56] Sample 1: 516 of 516 sites retained (100%)
#> INFO [2026-11-09 15:29:56] Processing sample_1
#> INFO [2026-11-09 15:29:56] Finished processing sample_1
#> SUCCESS [2026-11-09 15:29:56] Computed all deviations successfully
#>       sample_1
#> BATF -0.124932
```
