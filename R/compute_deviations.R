#' @title computeExpectations
#' @description  This function is used to calculate genome-wide expected methylation
#' for each motif.
#' @param msites Imported methylation sites using \code{read_methylome} function
#' @param gcdist a \code{GRanges} object contains Genome wide GC distribution
#' @param gcfreq a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param ignoreStrand - if TRUE, it ignores strand info from annotation
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @importFrom logger log_info log_warn log_error
#' @keywords internal
computeExpectations <- function(msites, gcdist, gcfreq, ignoreStrand = TRUE) {
  if (!is.logical(ignoreStrand)) {
    logger::log_warn("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (is.null(msites) || !is(msites, "GRanges")) {
    logger::log_error("Please provide a valid methylation sites with read_methylome function")
  }
  if (is.null(gcdist) || !is(gcdist, "GRanges")) {
    logger::log_error("Please provide a valid GC distribution")
  }
  if (!is.matrix(gcfreq)) {
    logger::log_error("Please provide a valid GC bin frequency table as a matrix")
  }
  hits <- findOverlaps(msites, gcdist, type = "within", ignore.strand = ignoreStrand)
  gcmap <- data.table(
    mscore = msites[hits@from]$score,
    gcbin = gcdist[hits@to]$GC_bin
  )
  exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]
  # exp_meth <- exp_meth[, .(mscore = mean(avg_mscore)), by = gcbin]
  exp_meth <- exp_meth[order(gcbin)]
  exp_meth <- as.matrix(exp_meth)
  exp.data <- t(gcfreq) %*% exp_meth[, 2]
  mpos <- round(seq(-floor(length(exp.data) / 2), floor(length(exp.data) / 2),
    length.out = length(exp.data)
  ))
  exp.methyl <- data.table(x = mpos, avg_methyl = exp.data)
  colnames(exp.methyl) <- c("x", "avg_methyl")
  return(exp.methyl)
}

#' @title computeDeviation
#' @description computeDeviation is a function to calculate the deviation in transcription factor
#' binding sites for a given motif
#' @param motif A character vector containing motif name
#' @param msites imported methylation sites
#' @param tf_bindsites a \code{GRangesList} object contains tf binding sites positions
#' @param gcfreqs a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param gcdist a \code{GRanges} object contains Genome wide GC distribution
#' @param enhancer  a \code{GRanges} object specifying regions such as distal motif (optional)
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param intermediate if TRUE, it returns the intermediate results as well
#' @return a \code{numeric} deviation score for a given motif
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table setDT
#' @importFrom stats na.omit
#' @importFrom S4Vectors mcols
#' @import data.table
#' @examples
#' library(methylTFR)
#'
#' # Load the data
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#' load(system.file("extdata", "FOXF2_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "FOXF2_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#'
#' # Compute the deviation
#' devs <- computeDeviation("FOXF2", msites, tf_bindsites, gcfreqs, gcdist)
#' @export
computeDeviation <- function(motif, msites, tf_bindsites, gcfreqs,
                             gcdist, enhancer = NULL, ignoreStrand = TRUE,
                             intermediate = FALSE) {
  if (!is.logical(ignoreStrand)) {
    logger::log_warn("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (!is.logical(intermediate)) {
    logger::log_warn("Found invalid intermediate option, using the default")
    intermediate <- FALSE
  }
  if (is.null(motif) || !is.character(motif)) {
    logger::log_error("Please provide a valid motif name")
  }
  if (is.null(msites) || !is(msites, "GRanges")) {
    logger::log_error("Please provide a valid methylation sites with read_methylome function")
  }
  if (is.null(tf_bindsites) || !any(c(!is(tf_bindsites, "GRangesList") || !is.list(tf_bindsites))) ) {
    logger::log_error("Please provide a valid tf binding sites as GRangesList")
  }
  if (is.null(gcfreqs) || !is.list(gcfreqs)) {
    logger::log_error("Please provide a valid GC bin frequency table list")
  }
  if (is.null(gcdist) || !is(gcdist, "GRanges")) {
    logger::log_error("Please provide a valid GC distribution")
  }
  if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
    logger::log_error("Please provide a valid enhancer regions")
  }
  tfbs <- tf_bindsites[[motif]]
  gcfreq <- gcfreqs[[motif]]
  tfbs <- resize(tfbs, width(tfbs)[1] + 101, fix = "center")

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(tfbs, enhancer, ignore.strand = ignoreStrand)
    tfbs <- tfbs[d_hits@from]
  }

  exp_meth <- computeExpectations(msites, gcdist, gcfreq, ignoreStrand)
  hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = ignoreStrand)

  S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) + ((start(tfbs) - end(tfbs)) / 2))
  x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
  sum_meth <- data.table::data.table(
    x = x,
    avg_methyl = msites[hits@from]$score # ,
    # avg_cov = msites[hits@from]$coverage
  )

  # Convert sum_meth and exp_meth to data.tables
  setDT(sum_meth)
  setDT(exp_meth)
  # sum_meth <- sum_meth[, .(n = .N, avg_methyl = mean(avg_methyl)), by = x]

  # GC bias correction
  obs_dev <- dev_helper(sum_meth)
  exp_dev <- dev_helper(exp_meth)
  dev <- obs_dev - exp_dev

  if (intermediate) {
    return(data.table::data.table(exp_dev, dev))
  } else {
    return(dev)
  }
}

#' @title dev_helper
#' @description Helper function for computeDeviation
#' @param data A data.table object with columns 'x', 'avg_methyl'
#' @return A numeric value
#' @keywords internal
dev_helper <- function(data) {
  data[, cuts := cut(x, c(-250, -200, -25, 25, 200, 250))]
  interval_mean <- data[, .(n = .N, mean = mean(avg_methyl)), by = cuts]
  interval_mean <- na.omit(interval_mean[order(cuts)])
  num_intervals <- nrow(interval_mean)
  if (num_intervals > 0) {
    avg_first_last <- ((interval_mean$mean[1] + interval_mean$mean[num_intervals]) / 2)
    var <- interval_mean$mean[(num_intervals + 1) %/% 2] / avg_first_last
  } else {
    var <- NA
  }
  return(var)
}
