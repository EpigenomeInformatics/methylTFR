# Test plotExpectedFootprint

library(methylTFR)

load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

test_that("plotExpectedFootprint returns a ggplot", {
    p <- plotExpectedFootprint(
        motif = "BATF",
        tfBindsites = tf_bindsites,
        msites = msites,
        sampleName = NULL,
        gcDist = gcdist,
        gcfreqs = gcfreqs,
        returnPlotData = FALSE,
        enhancer = NULL
    )
    expect_s3_class(p, "ggplot")
})

test_that("returnPlotData = TRUE returns the plot data", {
    result <- plotExpectedFootprint(
        motif = "BATF",
        tfBindsites = tf_bindsites,
        msites = msites,
        gcDist = gcdist,
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
        tfBindsites = tf_bindsites,
        msites = msites,
        gcDist = gcdist,
        gcfreqs = gcfreqs,
        enhancer = range(msites)
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotExpectedFootprint validates its input", {
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tfBindsites = tf_bindsites, msites = NULL,
            gcDist = gcdist, gcfreqs = gcfreqs
        ),
        "readMethylome"
    )
    expect_error(
        plotExpectedFootprint(
            motif = NULL, tfBindsites = tf_bindsites, msites = msites,
            gcDist = gcdist, gcfreqs = gcfreqs
        ),
        "Please provide a valid motif name"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tfBindsites = NULL, msites = msites,
            gcDist = gcdist, gcfreqs = gcfreqs
        ),
        "Please provide a valid tf binding sites"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tfBindsites = tf_bindsites, msites = msites,
            gcDist = NULL, gcfreqs = gcfreqs
        ),
        "Please provide a valid gcDist as GRanges"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tfBindsites = tf_bindsites, msites = msites,
            gcDist = gcdist, gcfreqs = NULL
        ),
        "Please provide a valid gcfreqs as list"
    )
    expect_error(
        plotExpectedFootprint(
            motif = "BATF", tfBindsites = tf_bindsites, msites = msites,
            gcDist = gcdist, gcfreqs = gcfreqs, enhancer = "chr1"
        ),
        "enhancer"
    )
})
