# Test plotExpectedFootprint

library(methylTFR)

load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

test_that("plotExpectedFootprint returns a ggplot", {
    p <- plotExpectedFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        sample_name = NULL,
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        returnPlotData = FALSE,
        enhancer = NULL
    )
    expect_s3_class(p, "ggplot")
})

test_that("returnPlotData = TRUE returns the plot data", {
    result <- plotExpectedFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        returnPlotData = TRUE
    )
    expect_type(result, "list")
    expect_s3_class(result$plot, "ggplot")
    expect_s3_class(result$plotDF, "data.table")
    expect_setequal(unique(result$plotDF$type), c("Expected", "Observed"))
})

test_that("plotExpectedFootprint works with an enhancer", {
    p <- plotExpectedFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        enhancer = range(msites)
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotExpectedFootprint validates its input", {
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tf_bindsites = tf_bindsites, msites = NULL,
            gc_dist = gcdist, gcfreqs = gcfreqs
        ),
        "readMethylome"
    )
    expect_error(
        plotExpectedFootprint(
            motif = NULL, tf_bindsites = tf_bindsites, msites = msites,
            gc_dist = gcdist, gcfreqs = gcfreqs
        ),
        "Please provide a valid motif name"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tf_bindsites = NULL, msites = msites,
            gc_dist = gcdist, gcfreqs = gcfreqs
        ),
        "Please provide a valid tf binding sites"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tf_bindsites = tf_bindsites, msites = msites,
            gc_dist = NULL, gcfreqs = gcfreqs
        ),
        "Please provide a valid gc_dist as GRanges"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tf_bindsites = tf_bindsites, msites = msites,
            gc_dist = gcdist, gcfreqs = NULL
        ),
        "Please provide a valid gcfreqs as list"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tf_bindsites = tf_bindsites, msites = msites,
            gc_dist = gcdist, gcfreqs = gcfreqs, enhancer = "chr1"
        ),
        "enhancer"
    )
})
