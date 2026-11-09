# Regression tests for fixed bugs

library(methylTFR)

test_that("run_methyltfr accepts .csv sample annotation files", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))

    tmp <- tempfile("methylTFR_csv")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    ann <- file.path(tmp, "samples.csv")
    write.csv(data.frame(bedFile = "missing_sample.bed"),
        ann,
        row.names = FALSE, quote = FALSE
    )

    # The call should get past parsing and fail on the missing BED file.
    expect_error(
        run_methyltfr(
            sample_ann = "samples.csv", sample_dir = tmp,
            tf_bindsites = tf_bindsites, gcfreqs = gcfreqs,
            gc_dist = gcdist, filetype = "epp"
        ),
        "Some of the files does not exist"
    )
})

test_that("run_methyltfr reports a missing sample column clearly", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))

    tmp <- tempfile("methylTFR_col")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    ann <- file.path(tmp, "samples.tsv")
    write.table(data.frame(wrongColumn = "a.bed"),
        ann,
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    expect_error(
        run_methyltfr(
            sample_ann = "samples.tsv", sample_dir = tmp,
            tf_bindsites = tf_bindsites, gcfreqs = gcfreqs,
            gc_dist = gcdist, filetype = "epp"
        ),
        "was not found in the sample annotation file"
    )
})

test_that("differential_deviation_test falls back to column names", {
    load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
    load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
    devs <- cbind(tc_mem, tc_naive)

    mat <- deviations(devs)
    colnames(mat) <- rep(c("mem", "naive"), each = ncol(mat) / 2)

    res <- differential_deviation_test(
        deviations = mat,
        groups = NULL,
        alternative = "two.sided",
        parametric = TRUE
    )
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), nrow(mat))
    expect_true(all(res$p_value >= 0 & res$p_value <= 1, na.rm = TRUE))
})

test_that("differential_deviation_test rejects mismatched group labels", {
    load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
    mat <- deviations(tc_mem)

    expect_error(
        differential_deviation_test(
            deviations = mat,
            groups = c("a", "b"),
            alternative = "two.sided"
        ),
        "one entry per column"
    )
    expect_error(
        differential_deviation_test(
            deviations = mat,
            groups = rep("a", ncol(mat)),
            alternative = "two.sided"
        ),
        "at least two distinct groups"
    )
})

test_that("dev_helper returns NA when a needed bin is empty", {
    full <- data.table::data.table(
        x = c(-230, -100, 0, 100, 230),
        avg_methyl = c(0.8, 0.6, 0.2, 0.6, 0.8)
    )
    expect_equal(methylTFR:::dev_helper(data.table::copy(full)), 0.25)
    # The old code returned 0.857 here
    expect_true(is.na(
        methylTFR:::dev_helper(data.table::copy(full[x != 230]))
    ))
    expect_true(is.na(
        methylTFR:::dev_helper(data.table::copy(full[x != 0]))
    ))
    expect_true(is.na(
        methylTFR:::dev_helper(data.table::copy(full[x != -230]))
    ))
})

test_that("dev_helper uses custom break points", {
    d <- data.table::data.table(
        x = c(-90, 0, 90),
        avg_methyl = c(0.5, 0.25, 0.5)
    )
    expect_equal(
        methylTFR:::dev_helper(d, breaks = c(-100, -80, -10, 10, 80, 100)),
        0.5
    )
})

test_that("mean_difference is signed for two groups", {
    mat <- rbind(
        up = c(1, 1.1, 0.9, 0, 0.1, -0.1),
        down = c(0, 0.1, -0.1, 1, 1.1, 0.9)
    )
    groups <- factor(rep(c("a", "b"), each = 3), levels = c("a", "b"))
    res <- differential_deviation_test(
        mat, groups, alternative = "two.sided"
    )
    expect_equal(res$mean_difference, c(1, -1))
})

test_that("runs leave no HDF5 files behind", {
    load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
    load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
    load(system.file("extdata", "example_data.rda", package = "methylTFR"))

    tmpDir <- tempfile("methylTFR_h5")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)
    hadWdTmp <- dir.exists(file.path(getwd(), "methylTFR_tmp"))

    res <- methylTFR:::methyltfr_core(
        sample_ids = "sample_1",
        msites_fun = function(i) msites,
        samples = data.frame(sampleName = "sample_1"),
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist,
        tempDir = tmpDir
    )

    expect_s4_class(res, "methylTFRdeviations")
    expect_length(list.files(tmpDir, pattern = "\\.h5$"), 0)
    if (!hadWdTmp) {
        expect_false(dir.exists(file.path(getwd(), "methylTFR_tmp")))
    }
})
