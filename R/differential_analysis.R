#' @title  differential_deviation_test
#' @description Differential analysis is to test which motifs are having significant
#' deviations among different groups. Adapted from chromVAR::differential_test
#' @param deviations a \code{data.frame} contains deviations for each motif
#' @param groups a character vector of group names or colnames(deviations)
#' @param motifs a character vector of motif names used in analysis or rownames(deviations)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @param parametric logical, if TRUE, parametric tests are used,
#' otherwise non-parametric tests are used.
#' @param padjMethod method for p-value adjustment, default is "BH"
#' @return a \code{data.frame} contains motifs with corresponding p-value and
#' adjusted p-value
#' @importFrom stats aggregate p.adjust
#' @export
differential_deviation_test <- function(deviations,
                                        groups = NULL,
                                        motifs = rownames(deviations),
                                        alternative = c("two.sided", "less", "greater"),
                                        parametric = TRUE,
                                        padjMethod = "BH") {
  if (!any(class(deviations) %in% c("data.frame", "matrix", "methylTFRDeviation"))) {
    stop("deviations must be a methylTFRDeviation object, or data.frame or matrix")
  }
  if (is(deviations, "methylTFRDeviation")) {
    deviations <- deviations(deviations)
  }
  # deviations <- t(deviations)
  if (is.null(groups)) {
    groups <- colnames(groups)
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
    p_value_adjusted = p_adj
  ))
}
