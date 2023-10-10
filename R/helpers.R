#' @title cleanMem
#' @description cleanMem is a function to clean the memory
#' @param iter.gc - number of times to run the garbage collector
#' @return NULL
#' @keywords internal
#'
cleanMem <- function(iter.gc = 1L) {
  for (i in 1:iter.gc) {
    gc()
  }
  invisible(NULL)
}

#' @title  anova test for multiple groups
#' @description anova test for multiple groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group names or colnames(deviations)
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
#' @param groups - a character vector of group names or colnames(deviations)
#' @param alternative - a character string specifying the alternative hypothesis,
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
#' @param groups - a character vector of group names or colnames(deviations)
#' @param alternative - a character string specifying the alternative hypothesis,
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

#' @title Kruskal-Wallis Rank Sum Test for multiple groups
#' @description Kruskal-Wallis Rank Sum Test for multiple groups
#' @param x - a numeric vector of deviations
#' @param groups - a character vector of group names or colnames(deviations)
#' @return a numeric vector of p-values
#' @importFrom stats kruskal.test
#' @keywords internal
#'
kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

#' @title helper function to compute z-score of a matrix
#' @description This function takes a matrix and computes z-score.
#' @param mat A matrix.
#' @return A matrix.
#' @importFrom matrixStats rowMeans2 rowSds
#' @keywords internal
#'
computeZScore <- function(mat) {
  mat <- (mat - matrixStats::rowMeans2(mat)) / matrixStats::rowSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
}
