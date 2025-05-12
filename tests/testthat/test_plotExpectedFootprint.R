# Test plotExpectedFootprint

library(methylTFR)
library(testthat)
library(ggplot2)

test_that("plotExpectedFootprint", {
    # Load example data
    load(system.file("extdata", "FOXF2_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))
    load(system.file("extdata", "FOXF2_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

    # Test: Valid inputs
    test_that("Valid inputs return a ggplot object", {
        p <- plotExpectedFootprint(
            motif = "FOXF2",
            tf_bindsites = tf_bindsites,
            msites = msites,
            sample_name = NULL,
            gc_dist = gcdist,
            gcfreqs = gcfreqs,
            returnPlotData = FALSE,
            bin_meth = NULL,
            enhancer = NULL
        )
        expect_s3_class(p, "ggplot")
    })

    # Test: Return plot data
    test_that("Valid inputs with returnPlotData = TRUE return plot data", {
        result <- plotExpectedFootprint(
            motif = "FOXF2",
            tf_bindsites = tf_bindsites,
            msites = msites,
            gc_dist = gcdist,
            gcfreqs = gcfreqs,
            returnPlotData = TRUE
        )
        expect_type(result, "list")
        expect_s3_class(result$plot, "ggplot")
        expect_s3_class(result$plotDF, "data.table")
    })

    # Test: Missing msites
    test_that("Missing msites throws an error", {
        expect_error(
            plotExpectedFootprint(
                motif = "FOXF2",
                tf_bindsites = tf_bindsites,
                msites = NULL,
                gc_dist = gcdist,
                gcfreqs = gcfreqs
            ),
            "msites must be a data frame,please provide the methylation sites"
        )
    })

    # Test: Invalid motif
    test_that("Invalid motif throws an error", {
        expect_error(
            plotExpectedFootprint(
                motif = NULL,
                tf_bindsites = tf_bindsites,
                msites = msites,
                gc_dist = gcdist,
                gcfreqs = gcfreqs
            ),
            "Please provide a valid motif name"
        )
    })

    # Test: Missing tf_bindsites
    test_that("Missing tf_bindsites throws an error", {
        expect_error(
            plotExpectedFootprint(
                motif = "FOXF2",
                tf_bindsites = NULL,
                msites = msites,
                gc_dist = gcdist,
                gcfreqs = gcfreqs
            ),
            "Please provide a valid tf binding sites as GRangesList"
        )
    })

    # Test: Missing gc_dist
    test_that("Missing gc_dist throws an error", {
        expect_error(
            plotExpectedFootprint(
                motif = "FOXF2",
                tf_bindsites = tf_bindsites,
                msites = msites,
                gc_dist = NULL,
                gcfreqs = gcfreqs
            ),
            "Please provide a valid gc_dist as GRanges"
        )
    })

    # Test: Missing gcfreqs
    test_that("Missing gcfreqs throws an error", {
        expect_error(
            plotExpectedFootprint(
                motif = "FOXF2",
                tf_bindsites = tf_bindsites,
                msites = msites,
                gc_dist = gcdist,
                gcfreqs = NULL
            ),
            "Please provide a valid gcfreqs as list"
        )
    })
})
