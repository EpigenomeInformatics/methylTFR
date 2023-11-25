#' @title cleanMem
#' @description cleanMem is a function to clean the memory
#' @param iter.gc - number of times to run the garbage collector
#' @return invisible NULL
#' @keywords internal
cleanMem <- function(iter.gc = 1L) {
  for (i in seq_along(iter.gc)) {
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
kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

#' @title helper function to compute row-wise z-score of a matrix
#' @description This function takes a matrix and computes row-wise z-scores.
#' @param mat A matrix
#' @return A matrix with row-wise z-scores
#' @importFrom matrixStats rowMeans2 rowSds
#' @keywords internal
computeRowZScore <- function(mat) {
  mat <- (mat - matrixStats::rowMeans2(mat)) / matrixStats::rowSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
}

#' @title Helper function to compute column-wise z-score of a matrix
#' @description This function takes a matrix and computes column-wise z-scores.
#' @param mat A matrix
#' @return A matrix with column-wise z-scores
#' @importFrom matrixStats colMeans2 colSds
#' @keywords internal
computeColZScore <- function(mat) {
  mat <- (mat - matrixStats::colMeans2(mat)) / matrixStats::colSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
}


#' @title compute_fp
#' @description Compute footprint returns a data.table required to create a transcription
#' factor footprint with label.
#' @param motif_name  chracter vector containing motif name eg: GATA string   
#' @param tf_bindsites  TF binding sites positions from methTFRannotation package
#' @param msites Nethylation data  
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @return a \code{data.table} object to plot tf footprint
#' @keywords internal 
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @import data.table 
computeFootprint <- function(motif_name, tf_bindsites, msites, enhancer = NULL) {
    tfbs <- tf_bindsites[[motif_name]]
    w <- width(tfbs)[1]
    tfbs <- resize(tfbs, w + 101, fix = "center")

    if (!is.null(enhancer)){
        d_hits <- findOverlaps(tfbs, enhancer)
        tfbs <- tfbs[d_hits@from]
    }
    
    hits <- findOverlaps(msites, tfbs, type = "within")
    
    mcols(tfbs)$mid_point <- round(end(tfbs) + (start(tfbs) - end(tfbs))/2)
    x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
    plot.data <- data.table(
        x = x, 
        y1 = msites[hits@from]$score, 
        y2 = msites[hits@from]$coverage
    )

    plot.data <- plot.data[, .(
        n = .N,
        avg_methyl = mean(y1),
        avg_cov = mean(y2),
        motif = motif_name,
        label = fifelse(x == max(x), as.character(motif), NA_character_)
    ), by = x]

    return(plot.data)
}
