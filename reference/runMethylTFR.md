# runMethylTFR

Compute bias-corrected deviation scores for every motif and every sample
listed in a sample annotation file.

## Usage

``` r
runMethylTFR(
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
  cov_threshold = 1,
  tempDir = file.path(tempdir(), "methylTFR")
)
```

## Arguments

- sample_ann:

  Name of a tab-separated (`.tsv`) or comma-separated (`.csv`) sample
  annotation file in `sample_dir`. It must have a header and contain the
  column named by `sampleColName`, which holds the methylation file
  names. All columns are kept as `colData` of the result. Not used when
  `annfile` is given.

- sample_dir:

  Directory holding the methylation call files (and, usually, the sample
  annotation file).

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

- sampleColName:

  Name of the column in `sample_ann` that holds the methylation file
  names. Default `"bedFile"`.

- chunkSize:

  Number of motifs processed together in one block. Default 20.

- full_path:

  Logical. If `TRUE`, the file names in `sampleColName` are full paths
  and `sample_dir` is not prepended. Default `FALSE`.

- annfile:

  Optional path to the sample annotation file. If given, it is used
  instead of `sample_ann` in `sample_dir`.

- threads:

  Number of workers used to process the motifs of a chunk in parallel.
  Default 1 (serial).

- enhancer:

  An optional `GRanges` object of regions, for example distal regulatory
  elements. Binding sites outside these regions are dropped. Any region
  set imported as `GRanges` works, for example the ENCODE candidate
  cis-regulatory elements (SCREEN) read with `rtracklayer::import()`.
  Default `NULL` uses all binding sites.

- filetype:

  Format of the methylation call files, one of `"bissnp"`, `"epp"`,
  `"allc"`, `"bismarkCytosine"`, `"bismarkCov"` or `"encode"`
  (case-insensitive). See
  [`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

- ignoreStrand:

  Logical. If `TRUE` (default), strand is ignored when overlapping
  methylation sites with annotations.

- cov_threshold:

  Numeric. Sites with fewer reads than this are dropped. Default 1.

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

Each sample is read with
[`readMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md),
binned by GC content with
[`addGCBinToMethylome`](https://epigenomeinformatics.github.io/methylTFR/reference/addGCBinToMethylome.md)
and scored per motif with
[`computeDeviation`](https://epigenomeinformatics.github.io/methylTFR/reference/computeDeviation.md).

## See also

[`runMethylTFRFromRnBeads`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFRFromRnBeads.md)
for running methylTFR directly on a preprocessed RnBeads object,
[`deviations`](https://epigenomeinformatics.github.io/methylTFR/reference/deviations.md)
and
[`deviationZScores`](https://epigenomeinformatics.github.io/methylTFR/reference/deviationZScores.md)
to extract the scores.

## Examples

``` r
# A minimal end-to-end run on the BATF example data bundled with the
# package. The annotation objects cover a single motif, so the result
# has one row.
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

# runMethylTFR() reads per-sample calls from disk, so the bundled sites
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

devs <- runMethylTFR(
    sample_ann = "samples.tsv",
    sample_dir = sample_dir,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gcdist,
    filetype = "bismarkCov"
)
#> SUCCESS [2026-11-09 14:50:25] The samples are successfully located
#> INFO [2026-11-09 14:50:25] Initializing the temp sink: /tmp/RtmpsM4D9Y/methylTFR/methylTFR18e3538a8a30.h5
#> INFO [2026-11-09 14:50:25] Initializing the temp sink: /tmp/RtmpsM4D9Y/methylTFR/methylTFR18e37a1ae82b.h5
#> INFO [2026-11-09 14:50:25] Processing sample_1.cov
#> INFO [2026-11-09 14:50:25] Finished processing sample_1.cov
#> SUCCESS [2026-11-09 14:50:25] Computed all deviations successfully
deviations(devs)
#>      sample_1.cov
#> BATF   -0.1167365

unlink(sample_dir, recursive = TRUE)
```
