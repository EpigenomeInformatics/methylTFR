# Test read_methylome

library(methylTFR)

test_that("read_methylome", {
  # Test read_methylome for ALLC format
  allc_path <- system.file("extdata", "allc.tsv.gz", package = "methylTFR")
  allc <- read_methylome(allc_path, "allc")

  # Check the length of the allc object
  expect_equal(length(allc), 3)

  # Check the class of the allc object
  expect_s4_class(allc, "GRanges")

  # Test read_methylome for EPP format
  epp_path <- system.file("extdata", "epp.tsv.gz", package = "methylTFR")
  epp <- read_methylome(epp_path, "EPP")

  # Test read_methylome for bismarkCytosine format
  bismarkCytosine_path <- system.file("extdata", "bismarkCytosine.tsv.gz", package = "methylTFR")
  bismarkCytosine <- read_methylome(bismarkCytosine_path, "bismarkCytosine")

  # Check the class of the bismarkCytosine object
  expect_s4_class(bismarkCytosine, "GRanges")

  # Test read_methylome for bismarkCov format
  bismarkCov_path <- system.file("extdata", "bismarkCov.tsv.gz", package = "methylTFR")
  bismarkCov <- read_methylome(bismarkCov_path, "bismarkCov")

  # Check the class of the bismarkCov object
  expect_s4_class(bismarkCov, "GRanges")

  # Test read_methylome for BisSNP format
  BisSNP_path <- system.file("extdata", "bissnp.tsv.gz", package = "methylTFR")
  BisSNP <- read_methylome(BisSNP_path, "BisSNP")

  # Check the class of the BisSNP object
  expect_s4_class(BisSNP, "GRanges")

  # Test read_methylome for encode format
  encode_path <- system.file("extdata", "encode.tsv.gz", package = "methylTFR")
  encode <- read_methylome(encode_path, "encode")

  # Check the class of the encode object
  expect_s4_class(encode, "GRanges")
})
