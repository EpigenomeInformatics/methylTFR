# Test addGCbintoMethylome

library(methylTFR)

test_that("addGCbintoMethylome", {
  # load example data
  load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
  load(system.file("extdata", "example_data.rda", package = "methylTFR"))
  
  # Add GC bin
  bin_meth <- addGCBintoMethylome(msites, gcdist)
  
  # Test the output
  expect_type(bin_meth, "double")
  expect_equal(dim(bin_meth), c(5, 2))

   # Test the values in avg_mscore column
  expect_true(all(bin_meth[, "avg_mscore"] >= 0))
  expect_true(all(bin_meth[, "avg_mscore"] <= 1))

  # Test the values in gcbin column
  expect_true(all(bin_meth[, "gcbin"] %in% 1:5))

  # Test if there are any NA values
  expect_true(!any(is.na(bin_meth)))

  # Test if gcbin values are integers
  expect_true(all(floor(bin_meth[, "gcbin"]) == bin_meth[, "gcbin"]))

  # Read ALLC format
  allc_path <- system.file("extdata", "allc.tsv.gz", package = "methylTFR")
  allc <- read_methylome(allc_path, "allc")
  expect_error(addGCBintoMethylome(allc, gcdist),
  "No methylation sites found in the GC distribution")

  # Test if the function throws an error when the input is not a GRanges object
  expect_error(addGCBintoMethylome(NULL, gcdist),
  "Please provide a valid methylation sites with read_methylome function")
    expect_error(addGCBintoMethylome(msites, NULL),
  "Please provide a valid GC distribution")
  expect_warning(addGCBintoMethylome(msites, gcdist, ignoreStrand = NULL),
    "Found invalid strand option, using the default")
})