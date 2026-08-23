#' @title  anova test for multiple groups
#' @description anova test for multiple groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group
#' names or colnames(deviations)
#' @return a numeric vector of p-values
#' @importFrom stats aov
#' @keywords internal
anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- aov(devs ~ groups, tmpdf)
  p_value <- summary(res)[[1]][["Pr(>F)"]][1]
  return(p_value)
}

#' @title  t-test for two groups
#' @description t-test for two groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group
#' names or colnames(deviations)
#' @param alternative - a character string specifying
#' the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @return a numeric vector of p-values
#' @importFrom stats t.test
#' @keywords internal
t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]], splitx[[2]],
    alternative = alternative,
    paired = FALSE,
    var.equal = FALSE,
    exact = FALSE
  )$p.value)
}

#' @title wilcoxon test for two groups
#' @description wilcoxon test for two groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group
#' names or colnames(deviations)
#' @param alternative - a character string specifying
#' the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @return a numeric vector of p-values
#' @importFrom stats wilcox.test
#' @keywords internal
wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
    alternative = alternative,
    paired = FALSE
  )$p.value)
}

#' @title Kruskal-Wallis Rank Sum Test for
#'  multiple groups
#' @description Kruskal-Wallis Rank Sum Test
#'  for multiple groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group
#'  names or colnames(deviations)
#' @return a numeric vector of p-values
#' @importFrom stats kruskal.test
#' @keywords internal
kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

#' @title helper function to compute row-wise
#'  z-score of a matrix
#' @description This function takes a matrix
#'  and computes row-wise z-scores.
#' @param mat A matrix
#' @return A matrix with row-wise z-scores
#' @importFrom matrixStats rowMeans2 rowSds
#' @keywords internal
computeRowZScore <- function(mat) {
  mat <- (mat - matrixStats::rowMeans2(mat)) /
    matrixStats::rowSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
}

#' @title Helper function to compute column-wise
#' z-score of a matrix
#' @description This function takes a matrix and
#' computes column-wise z-scores.
#' @param mat A matrix
#' @return A matrix with column-wise z-scores
#' @importFrom matrixStats colMeans2 colSds
#' @keywords internal
computeColZScore <- function(mat) {
  mat <- (mat - matrixStats::colMeans2(mat)) /
    matrixStats::colSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
}
