#' @title computeFootprint
#' @description Compute footprint returns a data.table required to create a transcription
#' factor footprint with label.
#' @param motif_name  chracter vector containing motif name eg: GATA string
#' @param tf_bindsites  TF binding sites positions from methTFRannotation package
#' @param msites Methylation data
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @return a list containing plot.data and tfbs
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom S4Vectors mcols
#' @import data.table
computeFootprint <- function(motif_name, tf_bindsites, msites, enhancer = NULL) {
  tfbs <- tf_bindsites[[motif_name]]
  w <- width(tfbs)[1]
  tfbs <- resize(tfbs, w + 130, fix = "center")

  if (!is.null(enhancer)) {
    subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE)
  }
  hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = TRUE)
  S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) + (start(tfbs) - end(tfbs)) / 2)
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
    label = fifelse(x == max(x), as.character(motif_name), NA_character_)
  ), by = x]

  return(plot.data)
}

#' @title computeExpectedFootprint
#' @description Compute expected footprint returns a data.table required to create a transcription
#' factor footprint with label.
#' @param motif  chracter vector containing motif name eg: GATA string
#' @param gcfreqs GC frequency of the genome used to compute expected methylation
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @param msites Methylation data
#' @param bin_meth Methylation data with GC content
#' @return a \code{data.table} object to containing TF footprint
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom S4Vectors mcols
#' @import data.table
computeExpectedFootprint <- function(
    motif, gcfreqs, gc_dist,
    enhancer = NULL, msites, bin_meth) {
  gcfreq <- gcfreqs[[motif]]

  if (!is.null(enhancer)) {
    gc_dist <- suppressWarnings(subsetByOverlaps(gc_dist, enhancer, ignore.strand = TRUE))
  }
  exp_meth <- computeExpectations(bin_meth, gcfreq)
  return(exp_meth)
}

