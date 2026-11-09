# Test plotMotifFootprint

library(methylTFR)

load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))

test_that("plotMotifFootprint returns a ggplot for both methods", {
    for (m in c("division", "substraction")) {
        p <- plotMotifFootprint(
            motif = "BATF",
            tf_bindsites = tf_bindsites,
            msites = msites,
            sample_name = "Sample1",
            gc_dist = gcdist,
            gcfreqs = gcfreqs,
            enhancer = NULL,
            method = m
        )
        expect_s3_class(p, "ggplot")
    }
})

test_that("plotMotifFootprint works without flank normalisation", {
    p <- plotMotifFootprint(
        motif = "BATF",
        tf_bindsites = tf_bindsites,
        msites = msites,
        gc_dist = gcdist,
        gcfreqs = gcfreqs,
        flankNorm = 0
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotMotifFootprint falls back to division for a bad method", {
    expect_warning(
        p <- plotMotifFootprint(
            motif = "BATF",
            tf_bindsites = tf_bindsites,
            msites = msites,
            gc_dist = gcdist,
            gcfreqs = gcfreqs,
            method = "other"
        ),
        "default method division"
    )
    expect_s3_class(p, "ggplot")
})

test_that("normalise_footprint_flank shifts or scales the curve", {
    d <- data.table::data.table(
        x = c(-100, -60, 0, 60, 100),
        avg_methyl = c(2, 4, 1, 4, 2)
    )
    sub <- methylTFR:::normalise_footprint_flank(
        data.table::copy(d), "substraction", 10
    )
    expect_equal(sub$avg_methyl, c(0, 2, -1, 2, 0))
    div <- methylTFR:::normalise_footprint_flank(
        data.table::copy(d), "division", 10
    )
    expect_equal(div$avg_methyl, c(1, 2, 0.5, 2, 1))
    same <- methylTFR:::normalise_footprint_flank(
        data.table::copy(d), "division", 0
    )
    expect_equal(same$avg_methyl, d$avg_methyl)
})
