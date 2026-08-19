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
# Not run: requires the RnBeads package, an hg38 annotation package and a
# preprocessed RnBeads set.
# \donttest{
if (requireNamespace("RnBeads", quietly = TRUE)) {
    # rnb_set <- RnBeads::load.rnb.set("reports/rnb.set_preprocessed")
    # gcfreqs <- getGCfreq(motifSet = "jaspar2020")
    # gc_dist <- getGenomeGC("hg38")
    # tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")
    #
    # deviations <- run_methylTFR_RnBeads(
    #     rnb_set = rnb_set,
    #     tf_bindsites = tf_bindsites,
    #     gcfreqs = gcfreqs,
    #     gc_dist = gc_dist,
    #     threads = 8,
    #     chunkSize = 15
    # )
}
#> NULL
# }
```
