# Test the methylTFRdeviations class

library(methylTFR)

load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))

test_that("accessors return the assays", {
    expect_true(is.matrix(deviations(tc_mem)) ||
        is(deviations(tc_mem), "DelayedArray"))
    expect_equal(dim(deviations(tc_mem)), dim(tc_mem))
    expect_equal(dim(deviationZScores(tc_mem)), dim(tc_mem))
    expect_error(deviations(matrix(1)), "unable to find an inherited method")
    expect_error(
        deviationZScores(data.frame(a = 1)),
        "unable to find an inherited method"
    )
})

test_that("cbind and rbind keep the class", {
    byCol <- cbind(tc_mem, tc_naive)
    expect_s4_class(byCol, "methylTFRdeviations")
    expect_equal(ncol(byCol), ncol(tc_mem) + ncol(tc_naive))
    expect_identical(cbind(tc_mem), tc_mem)

    byRow <- rbind(tc_mem[seq_len(2), ], tc_mem[3:4, ])
    expect_s4_class(byRow, "methylTFRdeviations")
    expect_equal(nrow(byRow), 4L)
    expect_identical(rbind(tc_mem), tc_mem)
})

test_that("the validity check requires the deviations and z assays", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(deviations = matrix(1, 2, 2))
    )
    expect_error(
        methods::new("methylTFRdeviations", se),
        "must contain 'deviations' and 'z'"
    )
})
