# Test run_methyltfr

library(methylTFR)

test_that("run_methyltfr", {
  # Load test data
  load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
  load(system.file("extdata", "FOXF2_gcfreqs.rda", package = "methylTFR"))
  load(system.file("extdata", "FOXF2_tf_bindsites.rda", package = "methylTFR"))

  # Test invalid inputs
  expect_error(run_methyltfr(
    sample_ann = "sample_ann.tsv", sample_dir = ".",
    filetype = "invalid"
  ), "Please provide a valid file type")
  expect_error(run_methyltfr(
    sample_ann = "sample_ann.tsv", sample_dir = ".",
    filetype = "epp", sampleColName = NULL
  ), "Please provide a valid sample column name")
  expect_error(run_methyltfr(
    sample_ann = "sample_ann.tsv", sample_dir = ".",
    filetype = "epp"
  ), "Please load the annotation objects for given genome.")

  # Test invalid chunk size
  # expect_warning(run_methyltfr(sample_ann = "sample_ann.tsv", sample_dir = ".",
  # tf_bindsites = tf_bindsites,gcfreqs = gcfreqs, gc_dist =gcdist,filetype ="epp", chunkSize = "invalid"),
  # "Invalid chunk size detected, using default chunk size")

  # TODO add example data for testing and test the rest of the function
})
