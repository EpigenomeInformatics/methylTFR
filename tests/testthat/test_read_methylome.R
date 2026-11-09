# Test readMethylome

library(methylTFR)

extdataFile <- function(x) system.file("extdata", x, package = "methylTFR")

writeTsv <- function(df, header = FALSE) {
    f <- tempfile(fileext = ".tsv")
    write.table(df, f,
        sep = "\t", row.names = FALSE, col.names = header, quote = FALSE
    )
    f
}

test_that("readMethylome parses every supported format", {
    allc <- readMethylome(extdataFile("allc.tsv.gz"), "allc")
    expect_s4_class(allc, "GRanges")
    expect_equal(length(allc), 3)
    expect_equal(allc$score, round(c(1 / 2, 2 / 3, 1 / 2), 6))
    expect_equal(allc$coverage, c(2, 3, 2))

    epp <- readMethylome(extdataFile("epp.tsv.gz"), "EPP")
    expect_s4_class(epp, "GRanges")
    expect_equal(epp$score, c(1, 0.5, 1, 0.5, 0.814, 0.5))
    expect_equal(epp$coverage, c(27, 7, 20, 20, 70, 100))

    # Row four has no reads, so it is dropped
    bc <- readMethylome(extdataFile("bismarkCytosine.tsv.gz"), "bismarkCytosine")
    expect_s4_class(bc, "GRanges")
    expect_equal(length(bc), 6)
    expect_true(all(bc$coverage >= 1))

    bcov <- readMethylome(extdataFile("bismarkCov.tsv.gz"), "bismarkCov")
    expect_s4_class(bcov, "GRanges")
    expect_equal(bcov$score, c(1, 0, 1, 0, 1, 1))
    expect_true(all(as.character(GenomicRanges::strand(bcov)) == "*"))

    bissnp <- readMethylome(extdataFile("bissnp.tsv.gz"), "BisSNP")
    expect_s4_class(bissnp, "GRanges")
    expect_equal(bissnp$score, c(0.7969, 0.9062, 0.587, 0.5))
    expect_equal(bissnp$coverage, c(64, 64, 46, 4))

    encode <- readMethylome(extdataFile("encode.tsv.gz"), "encode")
    expect_s4_class(encode, "GRanges")
    expect_equal(length(encode), 6)
    # Column 11 is a percentage, so the score is that percentage over 100.
    expect_equal(encode$score, c(0.06, 0.03, 0, 1, 0.55, 1))
    expect_equal(encode$coverage, c(62, 62, 31, 5, 31, 10))
})

test_that("covThreshold filters low-coverage sites", {
    encode <- readMethylome(
        extdataFile("encode.tsv.gz"), "encode", covThreshold = 20
    )
    expect_equal(length(encode), 4)
    bissnp <- readMethylome(
        extdataFile("bissnp.tsv.gz"), "bissnp", covThreshold = 50
    )
    expect_equal(length(bissnp), 2)
})

test_that("readMethylome rejects bad arguments", {
    expect_error(
        readMethylome("no_such_file.tsv", "epp"),
        "doesn't exist"
    )
    expect_error(
        readMethylome(extdataFile("epp.tsv.gz"), "bed"),
        "is not a valid file type"
    )
    expect_error(
        readMethylome(extdataFile("epp.tsv.gz"), "epp", covThreshold = -1),
        "is not a valid coverage threshold"
    )
    expect_error(
        readMethylome(extdataFile("epp.tsv.gz"), "epp", covThreshold = "a"),
        "is not a valid coverage threshold"
    )
})

test_that("parsers reject files with too few columns", {
    short <- writeTsv(data.frame(a = "chr1", b = 1, c = 2, d = 50))
    on.exit(unlink(short), add = TRUE)
    expect_error(
        suppressMessages(readMethylome(short, "allc")),
        "at least 6 columns"
    )
    expect_error(
        suppressMessages(readMethylome(short, "bismarkCytosine")),
        "at least 5 columns"
    )
    expect_error(
        suppressMessages(readMethylome(short, "bismarkCov")),
        "at least 6 columns"
    )
    expect_error(
        suppressMessages(readMethylome(short, "encode")),
        "must contain 11 columns"
    )
})

test_that("scores outside 0 to 1 are rejected", {
    bad <- writeTsv(data.frame(
        chr = "chr1", start = 10, end = 11, pct = 150, m = 3, u = 1
    ))
    on.exit(unlink(bad), add = TRUE)
    expect_error(readMethylome(bad, "bismarkCov"), "not a fraction")
})

test_that("an encode file whose column 11 is not a percentage is rejected", {
    bad <- writeTsv(
        data.frame(
            chrom = "chr1", start = 1000170, end = 1000171,
            name = "x", score = 62, strand = "+",
            thickStart = 1000170, thickEnd = 1000171, itemRgb = "255,255,0",
            coverage = 5, percentMeth = 400
        ),
        header = TRUE
    )
    on.exit(unlink(bad), add = TRUE)
    expect_error(readMethylome(bad, "encode"), "not a\\s+methylation percentage")
})
