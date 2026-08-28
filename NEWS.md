# methylTFR 0.99.6

NEW FEATURES

* Added `run_methylTFR_RnBeads()`, which runs the methylTFR workflow directly
  on a preprocessed RnBeads object instead of per-sample BED files. Samples are
  read one column at a time, so disk-backed RnBeads sets are never loaded into
  memory in full.
* Added `computeZScoreVariability()`, which ranks TF motifs by how much their
  activity varies across samples and tests each motif against a chi-squared
  null. Deviation scores are first calibrated against a within-sample null
  estimated across motifs, so that a variability above 1 is interpretable as
  "more variable than background".
* `methylTFRdeviations` objects returned by `run_methyltfr()` and
  `run_methylTFR_RnBeads()` now carry a third assay, `expected`, holding the
  GC-derived expected deviations. These were previously computed and then
  discarded.

BUG FIXES

* `run_methyltfr()` no longer rejects `.csv` sample annotation files. The
  `.tsv` branch's `else` clause caught every `.csv` file and raised an error
  after the file had already been read.
* `differential_deviation_test()` now falls back to `colnames(deviations)`
  when `groups` is NULL, instead of `colnames(groups)`, which was always NULL.
  Group labels are also validated against the number of columns.
* The invalid-input fallback for `cov_threshold` in `run_methyltfr()` now sets
  the documented default of 1 rather than 5.

DOCUMENTATION

* The heatmap in the memory T cell vignette is now drawn with
  `ComplexHeatmap` instead of `ggplot2`, with columns split by cell type
  and rows clustered.
* The compartment-agreement scatter plot in the same vignette now colours
  each motif by where it is differential: red for both compartments,
  green for CD4 only, blue for CD8 only and grey for neither.

INTERNAL

* The per-sample deviation loop shared by both entry points was factored out
  into `methyltfr_core()`, so the file-based and RnBeads-based workflows are
  guaranteed to produce identical results for the same methylation calls.

# methylTFR 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
