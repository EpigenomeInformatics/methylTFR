#' @title plotVariability
#' @description plot variability of motifs across cells or samples
#' @param variability output from \code{\link{computeVariability}}
#' @param xlab label for x-axis (default is 'Sorted TFs')
#' @param nLab  number of top motifs to label (default is 5)
#' @param motif_names vector of motif names to use as motif_names (default is
#' variability$name)
#' @import ggplot2
#' @return ggplot object
#' @export
plotVariability <- function(variability, xlab = "Sorted TFs", nLab = 5,
                            motif_names = variability$name) {
  if (!is.data.frame(variability)) {
    stop("variability must be output from computeVariability")
  }
  if (!is.character(xlab)) {
    stop("xlab must be a character value")
  }
  if (!is.numeric(nLab) || nLab %% 1 != 0 || nLab < 0) {
    stop("n must be a positive integer")
  }
  if (!is.character(motif_names)) {
    stop("motif_names must be a character vector")
  }
  if (length(motif_names) > nrow(variability)) {
    stop("motif_names can't have length more than the number of rows in variability")
  }
  res_df <- base::cbind(variability,
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

  if (nLab >= 1) {
    top_df <- res_df[res_df$rank <= nLab, ]
    plot <- plot +
      geom_text(data = top_df, size = 3, hjust = -0.45, col = "Black")
  }
  return(plot)
}

#' @title plotMotifFootprint
#' @description Creates a footprint plot for a given motif and a single methylation site.
#' @param motif Motif name as a character string
#' @param tf_bindsites Transcript Factor binding sites from the annotation package
#' @param sample Path to a single methylation data file
#' @param sample_name Optional sample name as a character string, defaults to the file name if not provided
#' @param type Type of methylation data, default is "EPP"
#' @param seed Seed for random number generation
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param gcfreqs GC frequency of the genome used to compute expected methylation
#' @param enhancer Specific region such as distal motif, proximal motif
#' @return A ggplot object of TF footprint plot for a single sample
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_line ggtitle theme_classic
plotMotifFootprint <- function(motif, tf_bindsites, sample,
                               sample_name = NULL, type = "EPP", seed = 12,
                               gc_dist, gcfreqs, enhancer = NULL) {
  # Prepare methylation sites for the single sample
  msites <- read_methylome(sample, type)

  # If sample label is not provided, use file name as label
  if (is.null(sample_name)) {
    sample_name <- basename(sample)
  }

  # Compute footprint for the motif
  plot_data <- computeFootprint(motif, tf_bindsites, msites, enhancer)

  # Compute expected footprint for the motif
  exp_data <- computeExpectedFootprint(motif, plot_data$tfbs, gcfreqs, gc_dist, seed, enhancer)

  # Add a new column to indicate whether the data is expected or observed
  exp_data[, type := "Expected"]
  plot_data$plot.data[, type := "Observed"]
  combined_data <- rbindlist(list(exp_data[, .(x, avg_methyl, type)], plot_data$plot.data[, .(x, avg_methyl, type)]))

  # Generate the footprint plot
  p1 <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
    geom_line() +
    geom_point() +
    geom_smooth(data = subset(combined_data, type == "Observed"), se = FALSE, color = "black") +
    xlab("Distance from motif center") +
    ylab("Methylation level") +
    theme_classic() +
    ggtitle(paste("TF footprint for", motif, "in", sample_name)) +
    scale_color_manual(values = c("Expected" = "blue", "Observed" = "red")) +
    theme(legend.position = "bottom")

  return(p1)
}
