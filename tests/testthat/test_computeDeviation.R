# Test computeDeviation

library(methylTFR)

test_that("computeDeviation", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    binMeth <- addGCBinToMethylome(msites, gcdist, TRUE)
    expect_true(is.matrix(binMeth))
    expect_true(all(c("gcbin", "avg_mscore") %in% colnames(binMeth)))

    devs <- computeDeviation(
        "BATF", msites, tf_bindsites, gcfreqs,
        enhancer = NULL, ignoreStrand = TRUE, binMeth
    )
    expect_type(devs$dev, "double")
    expect_type(devs$exp_dev, "double")
    # The example data must fill the centre and both flank bins
    expect_true(is.finite(devs$dev))
    expect_true(is.finite(devs$exp_dev))

    # Explicit defaults give the same result
    devs2 <- computeDeviation(
        "BATF", msites, tf_bindsites, gcfreqs,
        binMsites = binMeth, extendBy = 130,
        windowBreaks = c(-250, -200, -25, 25, 200, 250)
    )
    expect_equal(devs, devs2)

    # Enhancer restriction
    devsEnh <- computeDeviation(
        "BATF", msites, tf_bindsites, gcfreqs,
        enhancer = range(msites), binMsites = binMeth
    )
    expect_type(devsEnh$dev, "double")

    expect_error(
        computeDeviation(
            motif = NULL, msites, tf_bindsites, gcfreqs,
            binMsites = binMeth
        ),
        "Please provide a valid motif name"
    )
    expect_error(
        computeDeviation(
            "BATF", NULL, tf_bindsites, gcfreqs,
            binMsites = binMeth
        ),
        "readMethylome"
    )
    noCov <- msites
    noCov$coverage <- NULL
    expect_error(
        computeDeviation(
            "BATF", noCov, tf_bindsites, gcfreqs,
            binMsites = binMeth
        ),
        "coverage"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, NULL, gcfreqs,
            binMsites = binMeth
        ),
        "Please provide a valid tf binding sites"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, tf_bindsites, NULL,
            binMsites = binMeth
        ),
        "getGCfreq"
    )
    expect_error(
        computeDeviation(
            "NOTAMOTIF", msites, tf_bindsites, gcfreqs,
            binMsites = binMeth
        ),
        "not in names\\(tfBindsites\\)"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, tf_bindsites, gcfreqs,
            binMsites = matrix(1:4, 2)
        ),
        "addGCBinToMethylome"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, tf_bindsites, gcfreqs,
            enhancer = "chr1", binMsites = binMeth
        ),
        "enhancer"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, tf_bindsites, gcfreqs,
            binMsites = binMeth, extendBy = -1
        ),
        "extendBy"
    )
    expect_error(
        computeDeviation(
            "BATF", msites, tf_bindsites, gcfreqs,
            binMsites = binMeth, windowBreaks = c(1, 2, 3)
        ),
        "windowBreaks"
    )

    expect_warning(
        computeDeviation(
            "BATF", msites, tf_bindsites, gcfreqs,
            ignoreStrand = "invalid", binMsites = binMeth
        ),
        "Found invalid strand option, using the default"
    )
})

test_that("addGCBinToMethylome validates its input", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    expect_error(addGCBinToMethylome(NULL, gcdist), "readMethylome")
    expect_error(addGCBinToMethylome(msites, NULL), "getGenomeGC")
    noBin <- gcdist
    S4Vectors::mcols(noBin) <- NULL
    expect_error(addGCBinToMethylome(msites, noBin), "GC_bin")
    expect_warning(
        addGCBinToMethylome(msites, gcdist, ignoreStrand = "x"),
        "invalid strand option"
    )
    farAway <- GenomicRanges::GRanges(
        "chrUn", IRanges::IRanges(1, 30), GC_bin = 1L
    )
    expect_error(
        suppressWarnings(addGCBinToMethylome(msites, farAway)),
        "No methylation sites found in the GC distribution"
    )
})

test_that("computeExpectations rejects non-matrix input", {
    expect_error(
        methylTFR:::computeExpectations(1:5, matrix(1, 5, 5)),
        "as a matrix"
    )
    expect_error(
        methylTFR:::computeExpectations(matrix(1, 5, 2), 1:5),
        "as a matrix"
    )
})
