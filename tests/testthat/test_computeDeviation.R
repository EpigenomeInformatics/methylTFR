# Test computeDeviation

library(methylTFR)

test_that("computeDeviation", {
  # Load the data
  load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
  load(system.file("extdata", "FOXF2_gcfreqs.rda", package = "methylTFR"))
  load(system.file("extdata", "FOXF2_tf_bindsites.rda", package = "methylTFR"))
  load(system.file("extdata", "example_data.rda", package = "methylTFR"))

  # Compute the deviation
  devs <- computeDeviation("FOXF2", msites, tf_bindsites, gcfreqs)

  # Check that the result is a double value
  expect_type(devs$obs_dev, "double")
  expect_type(devs$tfbs, "integer")

  # Check that the function returns an error when given incorrect input
  expect_error(
    computeDeviation(motif = NULL, msites, tf_bindsites, gcfreqs, enhancer = NULL),
    "Please provide a valid motif name"
  )
  expect_error(
    computeDeviation("FOXF2", NULL, tf_bindsites, gcfreqs),
    "Please provide a valid methylation sites with read_methylome function"
  )
  expect_error(
    computeDeviation("FOXF2", msites, NULL, gcfreqs),
    "Please provide a valid tf binding sites as GRangesList"
  )
  expect_error(
    computeDeviation("FOXF2", msites, tf_bindsites, NULL),
    "Please provide a valid GC bin frequency table list"
  )

  # Check that the function returns a warning when given incorrect options
  expect_warning(
    computeDeviation("FOXF2", msites, tf_bindsites, gcfreqs, ignoreStrand = "invalid"),
    "Found invalid strand option, using the default"
  )
  # TODO: add the exp dev then run below as a test
  # Check that the function returns the expected result when given a specific input
  # expect_equal(computeDeviation("FOXF2", msites, tf_bindsites, gcfreqs))#, binMsites = bin_meth), devs)
})
