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

  return(list(plot.data = plot.data, tfbs = NROW(tfbs)))
}

#' @title computeExpectedFootprint
#' @description Compute expected footprint returns a data.table required to create a transcription
#' factor footprint with label.
#' @param motif  chracter vector containing motif name eg: GATA string
#' @param sample_size  Number of pseudo tfbs to be generated
#' @param gcfreqs GC frequency of the genome used to compute expected methylation
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param enhancer Specific regions such as distal motif, proximal motif
#' @param seed Seed for random number generation
#' @return a \code{data.table} object to containing TF footprint
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom S4Vectors mcols
#' @import data.table
computeExpectedFootprint <- function(
    motif, sample_size, gcfreqs, gc_dist,
    seed = 12, enhancer = NULL) {
  set.seed(seed)
  gcfreq <- gcfreqs[[motif]]

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(gc_dist, enhancer, ignore.strand = TRUE)
    gc_dist <- gc_dist[d_hits@from]
  }
  gc_dist <- as.data.table(gc_dist)
  binMsites <- addGCBintoMethylome(msites, gc_dist, TRUE, sample_size, gcfreq)
  exp_meth <- lapply(binMsites, function(i) computeExpectations(i, gcfreq))

  result_list <- lapply(exp_meth, function(dt) {
    aggregated <- dt[, .(avg_meth = mean(avg_methyl)), by = x]
    return(aggregated)
  })

  final_result <- rbindlist(result_list)
  final_aggregated <- final_result[, .(avg_methyl = mean(avg_meth)), by = x]
  return(final_aggregated)
}
