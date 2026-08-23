# Regression tests for fixes shipped in 0.99.1

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

    # Previously the .tsv branch's else clause rejected every .csv file. The
    # call should now get past parsing and fail on the missing BED file.
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

    # groups = NULL previously resolved to colnames(groups), i.e. NULL
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
