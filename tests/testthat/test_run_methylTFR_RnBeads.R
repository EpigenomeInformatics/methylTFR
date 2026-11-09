# Test runMethylTFRFromRnBeads and the shared engine

library(methylTFR)

test_that("runMethylTFRFromRnBeads rejects objects that are not RnBSets", {
    skip_if_not_installed("RnBeads")
    expect_error(
        runMethylTFRFromRnBeads(rnbSet = NULL),
        "Please provide a valid RnBSet object"
    )
    expect_error(
        runMethylTFRFromRnBeads(rnbSet = data.frame(a = 1)),
        "Please provide a valid RnBSet object"
    )
})

test_that("runMethylTFRFromRnBeads explains itself when RnBeads is absent", {
    skip_if(requireNamespace("RnBeads", quietly = TRUE))
    expect_error(
        runMethylTFRFromRnBeads(rnbSet = NULL),
        "RnBeads package is required"
    )
})

test_that("runMethylTFRFromRnBeads runs on a small RnBiseqSet", {
    skip_if_not_installed("RnBeads")
    skip_if_not_installed("RnBeads.hg38")
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    oldIds <- RnBeads::rnb.getOption("identifiers.column")
    RnBeads::rnb.options(identifiers.column = "sampleName")
    on.exit(RnBeads::rnb.options(identifiers.column = oldIds), add = TRUE)

    ids <- c("sample_1", "sample_2")
    sites <- data.frame(
        chr = as.character(GenomicRanges::seqnames(msites)),
        start = GenomicRanges::start(msites),
        strand = "*",
        stringsAsFactors = FALSE
    )
    meth <- cbind(msites$score, msites$score)
    covg <- cbind(msites$coverage, msites$coverage)
    dimnames(meth) <- dimnames(covg) <- list(NULL, ids)
    rnbSet <- RnBeads::RnBiseqSet(
        pheno = data.frame(sampleName = ids, group = c("a", "b")),
        sites = sites,
        meth = meth,
        covg = covg,
        assembly = "hg38",
        summarize.regions = FALSE
    )

    tmpDir <- tempfile("methylTFR_h5")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)
    devs <- runMethylTFRFromRnBeads(
        rnbSet = rnbSet,
        tfBindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcDist = gcdist,
        tempDir = tmpDir
    )
    expect_s4_class(devs, "methylTFRdeviations")
    expect_equal(colnames(deviations(devs)), ids)
    expect_equal(
        unname(deviations(devs)[, 1]), unname(deviations(devs)[, 2])
    )
    expect_true("group" %in% colnames(SummarizedExperiment::colData(devs)))
    expect_length(list.files(tmpDir, pattern = "\\.h5$"), 0)

    expect_error(
        runMethylTFRFromRnBeads(
            rnbSet = rnbSet, tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs, gcDist = gcdist,
            sampleAnn = data.frame(x = 1)
        ),
        "one row per sample"
    )
    expect_error(
        runMethylTFRFromRnBeads(
            rnbSet = rnbSet, tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs, gcDist = gcdist,
            sampleAnn = "samples.tsv"
        ),
        "sampleAnn must be a data.frame"
    )
    expect_error(
        runMethylTFRFromRnBeads(
            rnbSet = rnbSet, tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs, gcDist = gcdist,
            covThreshold = 1e6
        ),
        "No sites passed the quality thresholds"
    )
})

test_that("the shared engine validates its arguments", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))

    expect_error(
        methylTFR:::methylTFRCore(
            sampleIds = character(0),
            msitesFun = function(i) NULL,
            samples = data.frame(),
            tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gcDist = gcdist
        ),
        "No samples to process"
    )

    expect_error(
        methylTFR:::methylTFRCore(
            sampleIds = c("a", "b"),
            msitesFun = "not a function",
            samples = data.frame(id = c("a", "b")),
            tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gcDist = gcdist
        ),
        "msitesFun must be a function"
    )

    expect_error(
        methylTFR:::methylTFRCore(
            sampleIds = c("a", "b"),
            msitesFun = function(i) NULL,
            samples = data.frame(id = "a"),
            tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gcDist = gcdist
        ),
        "one row per sample"
    )

    expect_error(
        methylTFR:::methylTFRCore(
            sampleIds = "a",
            msitesFun = function(i) NULL,
            samples = data.frame(id = "a"),
            tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gcDist = gcdist,
            tempDir = c("x", "y")
        ),
        "tempDir must be a single directory path"
    )

    expect_error(
        methylTFR:::methylTFRCore(
            sampleIds = "a",
            msitesFun = function(i) data.frame(),
            samples = data.frame(id = "a"),
            tfBindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gcDist = gcdist,
            tempDir = tempfile()
        ),
        "readMethylome"
    )
})

test_that("the shared engine runs end to end on the example data", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    sampleIds <- c("sample_1", "sample_2")
    res <- methylTFR:::methylTFRCore(
        sampleIds = sampleIds,
        msitesFun = function(i) msites,
        samples = data.frame(
            sampleName = sampleIds,
            stringsAsFactors = FALSE
        ),
        tfBindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcDist = gcdist,
        chunkSize = 1,
        threads = 1,
        enhancer = range(msites)
    )

    expect_s4_class(res, "methylTFRdeviations")
    expect_equal(ncol(deviations(res)), 2L)
    expect_true(all(
        c("deviations", "z", "expected") %in%
            SummarizedExperiment::assayNames(res)
    ))
    # The same methylation calls were handed in twice, so both columns must
    # agree. This is what guarantees that the file-based and RnBeads-based
    # entry points cannot drift apart.
    expect_equal(
        unname(deviations(res)[, 1]),
        unname(deviations(res)[, 2])
    )
})
