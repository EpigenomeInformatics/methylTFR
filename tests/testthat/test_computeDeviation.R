# Test computeDeviation

library(methylTFR)

test_that("computeDeviation", {
    # Load the data
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    # Compute binMsites
    bin_meth <- addGCBintoMethylome(msites, gcdist, TRUE)

    # Compute the deviation
    devs <- computeDeviation("BATF", msites, tf_bindsites, gcfreqs, enhancer = NULL, ignoreStrand = TRUE, bin_meth)

    # Check that the result is a double value
    expect_type(devs$dev, "double")
    expect_type(devs$exp_dev, "double")

    # Check that the function returns an error when given incorrect input
    expect_error(
        computeDeviation(motif = NULL, msites, tf_bindsites, gcfreqs, enhancer = NULL, ignoreStrand = TRUE, bin_meth),
        "Please provide a valid motif name"
    )
    expect_error(
        computeDeviation(motif = "BATF", NULL, tf_bindsites, gcfreqs, enhancer = NULL, ignoreStrand = TRUE, bin_meth),
        "Please provide a valid methylation sites with read_methylome function"
    )
    expect_error(
        computeDeviation(motif = "BATF", msites, NULL, gcfreqs, enhancer = NULL, ignoreStrand = TRUE, bin_meth),
        "Please provide a valid tf binding sites"
    )
    expect_error(
        computeDeviation("BATF", msites, tf_bindsites, NULL, enhancer = NULL, ignoreStrand = TRUE, bin_meth),
        "Please provide a valid GC bin frequency table as a matrix"
    )

    # Check that the function returns a warning when given incorrect options
    expect_warning(
        computeDeviation("BATF", msites, tf_bindsites, gcfreqs, enhancer = NULL, ignoreStrand = "invalid", bin_meth),
        "Found invalid strand option, using the default"
    )
})
