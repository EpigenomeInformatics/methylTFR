#' @title resolve_diff_groups
#' @description Derive and validate the group labels used by
#' \code{differentialDeviationTest}.
#' @param deviations A matrix of deviation scores, motifs in rows.
#' @param groups Group labels, or NULL to take them from the column names.
#' @return The group labels as a \code{factor}.
#' @keywords internal
resolve_diff_groups <- function(deviations, groups) {
    if (is.null(groups)) {
        groups <- colnames(deviations)
    }
    if (is.null(groups)) {
        stop(
            "No group labels found. Provide 'groups', or supply deviations ",
            "with column names identifying the groups."
        )
    }
    groups <- as.factor(groups)
    if (length(groups) != ncol(deviations)) {
        stop("'groups' must have one entry per column of 'deviations'")
    }
    if (nlevels(groups) < 2) {
        stop("'groups' must contain at least two distinct groups")
    }
    return(groups)
}

#' @title diff_pvalues
#' @description Test every motif for a difference between the groups, with
#' the test chosen from the number of groups and \code{parametric}.
#' @param deviations A matrix of deviation scores, motifs in rows.
#' @param groups The group labels as a \code{factor}.
#' @param parametric if TRUE, use a t-test or ANOVA, otherwise a Wilcoxon
#' or Kruskal-Wallis test.
#' @param alternative The alternative hypothesis of the two-group tests.
#' @return A numeric vector of p-values, one per motif.
#' @keywords internal
diff_pvalues <- function(deviations, groups, parametric, alternative) {
    if (parametric) {
        if (nlevels(groups) == 2) {
            return(apply(deviations, 1, t_helper, groups, alternative))
        }
        return(apply(deviations, 1, anova_helper, groups))
    }
    if (nlevels(groups) == 2) {
        return(apply(deviations, 1, wilcoxon_helper, groups, alternative))
    }
    return(apply(deviations, 1, kw_helper, groups))
}

#' @title group_mean_difference
#' @description Effect size per motif. For two groups, the mean of the
#' first level minus the mean of the second level (signed). For more
#' groups, the range of the group means (unsigned).
#' @param deviations A matrix of deviation scores, motifs in rows.
#' @param groups The group labels as a \code{factor}.
#' @return A numeric vector with one value per motif.
#' @keywords internal
group_mean_difference <- function(deviations, groups) {
    group_means <- vapply(
        levels(groups),
        function(g) rowMeans(deviations[, groups == g, drop = FALSE]),
        numeric(nrow(deviations))
    )
    if (nlevels(groups) == 2) {
        return(group_means[, 1] - group_means[, 2])
    }
    return(apply(group_means, 1, function(x) max(x) - min(x)))
}

#' @title differentialDeviationTest
#' @description Test which motifs have different deviation scores between
#' groups of samples. Adapted from \code{chromVAR::differentialDeviations}.
#' @details With two groups, a Welch t-test (\code{parametric = TRUE}) or a
#' Wilcoxon rank sum test is used; with more groups, a one-way ANOVA or a
#' Kruskal-Wallis test.
#' @param deviations A \code{\linkS4class{methylTFRdeviations}} object, or
#' a matrix or \code{data.frame} of deviation scores with motifs in rows
#' and samples in columns.
#' @param groups A vector with one group label per sample. If \code{NULL},
#' the column names of \code{deviations} are used. For two groups, the
#' factor levels set the direction of \code{mean_difference}.
#' @param motifs Motif names, by default \code{rownames(deviations)}.
#' @param alternative The alternative hypothesis of the two-group tests,
#' one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param parametric Logical. If \code{TRUE} (default), parametric tests
#' are used, otherwise non-parametric tests.
#' @param padjMethod Method for p-value adjustment, passed to
#' \code{\link[stats]{p.adjust}}. Default \code{"BH"}.
#' @return A \code{data.frame} with one row per motif and the columns
#' \code{motifs}, \code{p_value}, \code{p_value_adjusted} and
#' \code{mean_difference}. For two groups, \code{mean_difference} is the
#' mean of the first group minus the mean of the second; for more groups,
#' it is the range of the group means.
#' @seealso \code{\link{computeZScoreVariability}} for ranking motifs
#' without group labels.
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @examples
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
#' devs <- cbind(tc_mem, tc_naive)
#'
#' # Construct group labels from sample names
#' groups <- sub("\\.bedGraph.*$", "", vapply(
#'     strsplit(colnames(devs), split = "_"), `[`, character(1), 1L
#' ))
#'
#' tc_result <- differentialDeviationTest(
#'     deviations = devs,
#'     groups = groups,
#'     alternative = "two.sided",
#'     parametric = TRUE,
#'     padjMethod = "BH"
#' )
#' head(tc_result)
#' @export
differentialDeviationTest <- function(
    deviations,
    groups = NULL,
    motifs = rownames(deviations),
    alternative = c("two.sided", "less", "greater"),
    parametric = TRUE,
    padjMethod = "BH"
) {
    if (!any(class(deviations) %in%
        c("data.frame", "matrix", "methylTFRdeviations"))) {
        stop(
            "deviations must be a methylTFRdeviations ",
            "object, or data.frame or matrix"
        )
    }
    if (is(deviations, "methylTFRdeviations")) {
        deviations <- deviations(deviations)
    }
    deviations <- as.matrix(deviations)
    groups <- resolve_diff_groups(deviations, groups)
    if (length(alternative) > 1) {
        stop("Please indicate one of the alternatives only.")
    }
    p_val <- diff_pvalues(deviations, groups, parametric, alternative)
    p_adj <- p.adjust(p_val, method = padjMethod)
    mean_diff <- group_mean_difference(deviations, groups)

    return(data.frame(
        motifs = motifs,
        p_value = p_val,
        p_value_adjusted = p_adj,
        mean_difference = mean_diff
    ))
}
