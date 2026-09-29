# methylTFR 0.99.9

CHANGES

* `Authors@R` lists the funders (ERA-NET Transcan-Neu III - EPILUNAR,
  grant 01KT2409; Saarland University NanoBioMed Young Investigator Grant)
  with the `fnd` role.
* `inst/CITATION` rewritten: the citation reads "Gunduz IB, Murugan SK,
  Mueller F (year). methylTFR: ... R package version x.y.z.", takes the
  version and year from DESCRIPTION, and carries the Bioconductor DOI.
* The README shows the citation and BibTeX entry.
* Added `CITATION.cff` (excluded from the build via `.Rbuildignore`).

BUG FIXES

* `rnb_sample_ids()` no longer calls `RnBeads::samples()`, which is not
  exported by all RnBeads versions (R CMD check: "Missing or unexported
  object"). Sample identifiers are read from the `identifiers.column` of
  `RnBeads::pheno()`, with the previous fallbacks unchanged.
* The vignette, README and `inst/scripts/vignette_data.R` call
  `getGenomeGC()` without an argument, which works with every version of
  the annotation packages.

# methylTFR 0.99.8

BUG FIXES

* `read_methylome(type = "encode")` returned methylation scores that were not
  fractions. `parse_encode()` divided column 11 by column 10, but in a bedMethyl
  file column 10 is read coverage and column 11 is the percentage of reads
  methylated, so a fully methylated CpG covered by five reads scored 20. The
  score is now the percentage over 100.
* `parse_encode()` skipped the first line of every file unconditionally. Files
  downloaded from ENCODE carry no header, so the first record was silently
  discarded. The header is now detected.
* `read_methylome()` stops if the parsed scores fall outside 0 to 1 rather than
  passing them to the deviation calculation.
* The `run_methylTFR_RnBeads()` example built an `RnBiseqSet` without setting
  `identifiers.column`, so `samples()` returned row numbers and the example
  stopped on the sample-identifier check. It now sets the option before
  building the set and restores it afterwards.
* The bundled `encode.tsv.gz` example has been rebuilt. Columns 10 and 11 are
  now named `coverage` and `percentMeth` rather than `blockCount` and
  `blockSizes`, and the percentages span 0 to 100, so the tests exercise the
  case the old example passed by luck.


# methylTFR 0.99.7

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
