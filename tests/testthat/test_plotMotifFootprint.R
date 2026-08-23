# Test plotMotifFootprint

library(methylTFR)
library(testthat)
library(ggplot2)

test_that("plotMotifFootprint", {
    # Load example data
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

    # Test with valid inputs
    p <- plotMotifFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        sample_name = "Sample1",
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        enhancer = NULL,
        method = "division"
    )

    # Check if the output is a ggplot object
    expect_s3_class(p, "ggplot")

    # Test with method = "substraction"
    p_sub <- plotMotifFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        sample_name = "Sample1",
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        enhancer = NULL,
        method = "substraction"
    )

    # Check if the output is a ggplot object
    expect_s3_class(p_sub, "ggplot")
})
