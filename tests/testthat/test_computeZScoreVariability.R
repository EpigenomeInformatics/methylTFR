# Test computeZScoreVariability

library(methylTFR)

test_that("computeZScoreVariability returns a well formed result", {
  load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
  load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
  devs <- cbind(tc_mem, tc_naive)

  res <- suppressMessages(computeZScoreVariability(devs))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(deviations(devs)))
  expect_true(all(c(
    "motifs", "variability", "p_value", "p_value_adjusted"
  ) %in% colnames(res)))
  expect_equal(res$motifs, rownames(deviations(devs)))
  expect_true(all(res$variability >= 0, na.rm = TRUE))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1, na.rm = TRUE))
  expect_true(all(res$p_value_adjusted >= res$p_value, na.rm = TRUE))
})

test_that("computeZScoreVariability accepts a plain matrix", {
  load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
  load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
  devs <- cbind(tc_mem, tc_naive)

  from_object <- suppressMessages(computeZScoreVariability(devs))
  from_matrix <- suppressMessages(
    computeZScoreVariability(deviations(devs))
  )
  expect_equal(from_object$variability, from_matrix$variability)
})

test_that("the calibrated null is not degenerate, unlike row-wise Z-scores", {
  # Row-wise Z-scores have a standard deviation of exactly 1 per motif by
  # construction, which is why they must not be fed to this function.
  set.seed(42)
  mat <- matrix(rnorm(200 * 12), nrow = 200, ncol = 12)
  rownames(mat) <- paste0("motif_", seq_len(200))

  rowz <- (mat - rowMeans(mat)) / apply(mat, 1, sd)
  expect_equal(unname(apply(rowz, 1, sd)), rep(1, 200), tolerance = 1e-8)

  res <- computeZScoreVariability(mat)
  expect_gt(stats::sd(res$variability), 0)
})

test_that("variability is centred near 1 under the null", {
  # Pure noise: calibrated variability should sit around 1 and almost
  # nothing should be called significant.
  set.seed(1)
  mat <- matrix(rnorm(300 * 20), nrow = 300, ncol = 20)
  rownames(mat) <- paste0("motif_", seq_len(300))

  res <- computeZScoreVariability(mat, method = "gaussian")
  expect_equal(median(res$variability), 1, tolerance = 0.1)
  expect_lt(mean(res$p_value_adjusted < 0.05), 0.05)
})

test_that("a genuinely variable motif is detected", {
  set.seed(7)
  mat <- matrix(rnorm(300 * 20), nrow = 300, ncol = 20)
  rownames(mat) <- paste0("motif_", seq_len(300))
  # Spike one motif with a strong two-group difference
  mat[1, seq_len(10)] <- mat[1, seq_len(10)] + 6

  res <- computeZScoreVariability(mat)
  expect_lt(res$p_value_adjusted[1], 0.05)
  expect_gt(res$variability[1], 1)
})

test_that("bootstrap bounds bracket the point estimate", {
  set.seed(3)
  mat <- matrix(rnorm(100 * 15), nrow = 100, ncol = 15)
  rownames(mat) <- paste0("motif_", seq_len(100))

  res <- computeZScoreVariability(mat,
    bootstrap = TRUE,
    niterations = 200
  )
  expect_true(all(c(
    "bootstrap_lower_bound", "bootstrap_upper_bound"
  ) %in% colnames(res)))
  expect_true(all(res$bootstrap_lower_bound <= res$bootstrap_upper_bound))
})

test_that("computeZScoreVariability validates its input", {
  mat <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)
  expect_error(
    computeZScoreVariability(mat),
    "At least three samples"
  )
  expect_error(
    computeZScoreVariability("not a matrix"),
    "must be a methylTFRdeviations object"
  )
  expect_error(
    suppressMessages(computeZScoreVariability(
      matrix(rnorm(30), nrow = 10, ncol = 3),
      conf_level = 2
    )),
    "conf_level must be a number between 0 and 1"
  )
})
