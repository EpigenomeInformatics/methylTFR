#' @title plotZscoreVariability
#' @description plot Z-score variability of motifs across cells or samples
#' @param variability output from \code{\link{computeZScoreVariability}}
#' @param xlab label for x-axis (default is 'Sorted TFs')
#' @param n  number of top motifs to label (default is 5)
#' @param motif_names vector of motif names to use as motif_names (default is
#' variability$name)
#' @import ggplot2
#' @return ggplot object
#' @export
plotZscoreVariability <- function(variability, xlab = "Sorted TFs", n = 5,
                                  motif_names = variability$name) {
  if (!is.data.frame(variability)) {
    stop("variability must be output from computeZScoreVariability")
  }
  if (!is.character(xlab)) {
    stop("xlab must be a character value")
  }
  if (!is.numeric(n) || n %% 1 != 0 || n < 0) {
    stop("n must be a positive integer")
  }
  if (!is.character(motif_names)) {
    stop("motif_names must be a character vector")
  }
  if (length(motif_names) > nrow(variability)) {
    stop("motif_names can't have length more than the number of rows in variability")
  }
  res_df <- cbind(variability,
    rank = rank(-1 * variability$variability,
      ties.method = "random"
    ),
    annotation = motif_names
  )

  if (!"bootstrap_lower_bound" %in% colnames(variability)) {
    plot <- ggplot(res_df, aes_string(x = "rank", y = "variability")) +
      geom_point() +
      xlab(xlab) +
      ylab(ylab) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, max(res_df$variability, na.rm = TRUE) * 1.05)
      ) +
      theme_classic()
  } else {
    plot <- ggplot(res_df, aes_string(
      x = "rank",
      y = "variability",
      min = "bootstrap_lower_bound",
      max = "bootstrap_upper_bound",
      label = "annotation"
    )) +
      geom_point() +
      geom_errorbar() +
      xlab(xlab) +
      ylab("Variability") +
      theme_classic() +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, max(res_df$bootstrap_upper_bound, na.rm = TRUE) * 1.05)
      )
  }

  if (n >= 1) {
    top_df <- res_df[res_df$rank <= n, ]
    plot <- plot +
      geom_text(data = top_df, size = 3, hjust = -0.45, col = "Black")
  }
  return(plot)
}

