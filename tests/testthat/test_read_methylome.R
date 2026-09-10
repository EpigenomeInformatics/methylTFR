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

    # Every row of the example is kept at the default coverage threshold
    expect_equal(length(encode), 6)

    # Column 11 is a percentage, so the score is that percentage over 100. Row
    # four is 100 percent at coverage 5: dividing column 11 by column 10, as
    # parse_encode() used to, returned 20 for it.
    expect_equal(encode$score, c(0.06, 0.03, 0, 1, 0.55, 1))
    expect_equal(encode$coverage, c(62, 62, 31, 5, 31, 10))
    expect_true(all(encode$score >= 0 & encode$score <= 1))

    # Coverage filtering reads column 10
    encode_cov <- read_methylome(encode_path, "encode", cov_threshold = 20)
    expect_equal(length(encode_cov), 4)
})

test_that("read_methylome rejects an encode file whose column 11 is not a percentage", {
    bad <- tempfile(fileext = ".tsv")
    on.exit(unlink(bad), add = TRUE)
    write.table(
        data.frame(
            chrom = "chr1", start = 1000170, end = 1000171,
            name = "x", score = 62, strand = "+",
            thickStart = 1000170, thickEnd = 1000171, itemRgb = "255,255,0",
            coverage = 5, percentMeth = 400
        ),
        bad, sep = "\t", row.names = FALSE, quote = FALSE
    )
    expect_error(read_methylome(bad, "encode"), "not a\\s+methylation percentage")
})
