# anova test for multiple groups
anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- aov(devs ~ groups, tmpdf)
  p_value <- summary(res)[[1]][["Pr(>F)"]][1]
  return(p_value)
}

# t-test
t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]], splitx[[2]],
    alternative = alternative,
    paired = FALSE,
    var.equal = FALSE,
    exact = FALSE
  )$p.value)
}

# wilcoxon test for two groups
wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
    alternative = alternative,
    paired = FALSE
  )$p.value)
}

# Kruskal-Wallis Rank Sum Test for multiple groups
kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

#' differential_deviation_test
#'
#'      Differential analysis is to test which motifs are having significant
#' deviations among different cell-types.
#' Adapted from https://github.com/GreenleafLab/chromVAR/blob/master/R/differential_tests.R
#'
#' @param deviations - motif name
#' @param groups     - a character vector of group names or colnames(deviations)
#' @param motifs     - a character vector of motif names used in analysis or rownames(deviations)
#' @return a \code{data.frame} contains motifs with corresponding p-value and adj-pvalue
#' @export
#' @importFrom stats aggregate
differential_deviation_test <- function(deviations,
                                        groups = NULL,
                                        motifs = rownames(deviations),
                                        alternative = c("two.sided", "less", "greater"),
                                        parametric = TRUE,
                                        padjMethod = "BH") {
  if (is.null(groups)) {
    groups <- colnames(groups)
  } else if (length(groups) != ncol(deviations)) {
    stop(
      "invalid groups input, must be vector of lench ncol(variantion) or column",
      " name from variations dataframe"
    )
  }

  groups <- as.factor(groups)
  if (length(alternative) > 1) {
    stop(
      "Please indicate one of the alternatives only."
    )
  }
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(deviations, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(deviations, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(deviations, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(deviations, 1, kw_helper, groups)
    }
  }
  p_adj <- p.adjust(p_val, method = padjMethod)

  return(data.frame(
    motifs = motifs,
    p_value = p_val,
    adj.pvalue = p_adj
  ))
}
