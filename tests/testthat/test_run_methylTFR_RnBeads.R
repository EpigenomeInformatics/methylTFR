# Test run_methylTFR_RnBeads and the shared engine

library(methylTFR)

test_that("run_methylTFR_RnBeads rejects objects that are not RnBSets", {
    skip_if_not_installed("RnBeads")
    expect_error(
        run_methylTFR_RnBeads(rnb_set = NULL),
        "Please provide a valid RnBSet object"
    )
    expect_error(
        run_methylTFR_RnBeads(rnb_set = data.frame(a = 1)),
        "Please provide a valid RnBSet object"
    )
})

test_that("run_methylTFR_RnBeads explains itself when RnBeads is absent", {
    skip_if(requireNamespace("RnBeads", quietly = TRUE))
    expect_error(
        run_methylTFR_RnBeads(rnb_set = NULL),
        "RnBeads package is required"
    )
})

test_that("the shared engine validates its arguments", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))

    expect_error(
        methylTFR:::methyltfr_core(
            sample_ids = character(0),
            msites_fun = function(i) NULL,
            samples = data.frame(),
            tf_bindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gc_dist = gcdist
        ),
        "No samples to process"
    )

    expect_error(
        methylTFR:::methyltfr_core(
            sample_ids = c("a", "b"),
            msites_fun = "not a function",
            samples = data.frame(id = c("a", "b")),
            tf_bindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gc_dist = gcdist
        ),
        "msites_fun must be a function"
    )

    expect_error(
        methylTFR:::methyltfr_core(
            sample_ids = c("a", "b"),
            msites_fun = function(i) NULL,
            samples = data.frame(id = "a"),
            tf_bindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            gc_dist = gcdist
        ),
        "one row per sample"
    )
})

test_that("the shared engine runs end to end on the example data", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    sample_ids <- c("sample_1", "sample_2")
    res <- methylTFR:::methyltfr_core(
        sample_ids = sample_ids,
        msites_fun = function(i) msites,
        samples = data.frame(
            sampleName = sample_ids,
            stringsAsFactors = FALSE
        ),
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist,
        chunkSize = 1,
        threads = 1
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
