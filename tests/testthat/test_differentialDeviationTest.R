# Test differential_deviation_test

library(methylTFR)

set.seed(11)
mat2 <- matrix(rnorm(20 * 8), nrow = 20, ncol = 8,
    dimnames = list(paste0("motif_", seq_len(20)), NULL)
)
groups2 <- rep(c("a", "b"), each = 4)
mat3 <- matrix(rnorm(20 * 9), nrow = 20, ncol = 9,
    dimnames = list(paste0("motif_", seq_len(20)), NULL)
)
groups3 <- rep(c("a", "b", "c"), each = 3)

checkResult <- function(res, n) {
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), n)
    expect_equal(
        colnames(res),
        c("motifs", "p_value", "p_value_adjusted", "mean_difference")
    )
    expect_true(all(res$p_value >= 0 & res$p_value <= 1))
    expect_true(all(res$p_value_adjusted >= res$p_value))
}

test_that("two groups use a t-test or a Wilcoxon test", {
    tRes <- differential_deviation_test(
        mat2, groups2, alternative = "two.sided", parametric = TRUE
    )
    checkResult(tRes, 20)
    expect_equal(
        tRes$p_value[1],
        t.test(mat2[1, 1:4], mat2[1, 5:8])$p.value
    )

    wRes <- differential_deviation_test(
        mat2, groups2, alternative = "less", parametric = FALSE
    )
    checkResult(wRes, 20)
    expect_equal(
        wRes$p_value[1],
        wilcox.test(mat2[1, 1:4], mat2[1, 5:8], alternative = "less")$p.value
    )
})

test_that("more than two groups use ANOVA or Kruskal-Wallis", {
    aRes <- differential_deviation_test(
        mat3, groups3, alternative = "two.sided", parametric = TRUE
    )
    checkResult(aRes, 20)
    expect_true(all(aRes$mean_difference >= 0))

    kRes <- differential_deviation_test(
        mat3, groups3, alternative = "two.sided", parametric = FALSE
    )
    checkResult(kRes, 20)
    expect_equal(
        kRes$p_value[1],
        kruskal.test(mat3[1, ], factor(groups3))$p.value
    )
})

test_that("data.frames and methylTFRdeviations objects are accepted", {
    load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
    load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
    devs <- cbind(tc_mem, tc_naive)
    groups <- rep(c("mem", "naive"), c(ncol(tc_mem), ncol(tc_naive)))

    fromObject <- differential_deviation_test(
        devs, groups, alternative = "two.sided"
    )
    fromDf <- differential_deviation_test(
        as.data.frame(deviations(devs)), groups, alternative = "two.sided"
    )
    expect_equal(fromObject$p_value, fromDf$p_value)
    expect_equal(fromObject$motifs, rownames(devs))
})

test_that("differential_deviation_test validates its input", {
    expect_error(
        differential_deviation_test(list(1, 2), groups2),
        "deviations must be a methylTFRdeviations"
    )
    expect_error(
        differential_deviation_test(mat2, groups2),
        "one of the alternatives only"
    )
    unnamed <- mat2
    colnames(unnamed) <- NULL
    expect_error(
        differential_deviation_test(unnamed, alternative = "two.sided"),
        "No group labels found"
    )
})
