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
#' @description Creates a footprint plot for given motifs and methylation site.
#' @param motifs Motif names as character vector
#' @param tf_bindsites -Transcript Factor binding sites from the annotation package
#' @param samples Sample names as character vector of paths to methylation data
#' @param sample_names Sample names as character vector
#' @param type Type of methylation data, default is "EPP"
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @return A ggplot object of TF footprint plot
#' @export
#' @importFrom ggplot2 ggplot ggsave geom_point geom_line ggtitle
#' @importFrom ggrepel geom_label_repel
plotMotifFootprint <- function(motif, tf_bindsites, samples,
                               sample_names = NULL,
                               type = "EPP", enhancer = NULL) {
  # Prepare a list of methylation sites
  msites_list <- lapply(samples, read_methylome, type = type)

  # If sample labels are not provided, use file names as labels
  if (is.null(sample_names)) {
    sample_names <- unlist(lapply(samples, basename))
  }
  names(msites_list) <- sample_names
  plot_data <- data.table()

  # Compute footprint for the single motif
  for (i in 1:length(msites_list)) {
    msites <- msites_list[[i]]
    current_plot <- computeFootprint(motif, tf_bindsites, msites,enhancer)
    current_plot[, sample_name := sample_names[i]]
    plot_data <- current_plot
  }

  # Generate the footprint plot
  p1 <- ggplot(plot_data, aes(x = x, y = avg_methyl, group = sample_name)) +
    geom_line(aes(color = sample_name)) +
    geom_point(aes(color = sample_name)) +
    xlab("Distance from motif center") +
    ylab("Methylation level") +
    theme_classic() +
    ggtitle(paste("TF footprint for", motif))

  return(p1)
}
