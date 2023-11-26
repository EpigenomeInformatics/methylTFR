#' @title computeFootprint
#' @description Compute footprint returns a data.table required to create a transcription
#' factor footprint with label.
#' @param motif_name  chracter vector containing motif name eg: GATA string
#' @param tf_bindsites  TF binding sites positions from methTFRannotation package
#' @param msites Methylation data 
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @return a \code{data.table} object to plot tf footprint
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @import data.table
computeFootprint <- function(motif_name, tf_bindsites, msites, enhancer = NULL) {
  tfbs <- tf_bindsites[[motif_name]]
  w <- width(tfbs)[1]
  tfbs <- resize(tfbs, w + 101, fix = "center")

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(tfbs, enhancer, ignore.strand = TRUE)
    tfbs <- tfbs[d_hits@from]
  }

  hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = TRUE)
  mcols(tfbs)$mid_point <- round(end(tfbs) + (start(tfbs) - end(tfbs)) / 2)
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
