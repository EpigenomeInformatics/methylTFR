#' @title computeExpectations
#' @description  This function is used to calculate expected methylation for a given
#' motif and sample.
#' @param binMsites Imported methylation sites with GC bin
#' @param gcfreq a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @importFrom logger log_info log_warn log_error
#' @keywords internal
computeExpectations <- function(binMsites, gcfreq) {
  if (!is.matrix(binMsites)) {
    stop("Please provide a valid GC bin frequency table as a matrix")
  }
  if (!is.matrix(gcfreq)) {
    stop("Please provide a valid GC bin frequency table as a matrix")
  }
  exp.data <- t(gcfreq) %*% binMsites[, 2]
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
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param enhancer  a \code{GRanges} object specifying regions such as distal motif (optional)
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param exp_dev a numeric vector of expected deviations computed by computeExpectedDeviation
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
#' # Add GC bin
#' bin_meth <- computeExpectedDeviation("FOXF2", msites, gcfreqs, gc_dist,TRUE)
#'
#' # Compute the deviation
#' devs <- computeDeviation("FOXF2", msites, tf_bindsites, gcfreqs)
#' @export
computeDeviation <- function(motif, msites, tf_bindsites, gcfreqs, gc_dist,
                             enhancer = NULL, ignoreStrand = TRUE,exp_dev) {
  if (!is.logical(ignoreStrand)) {
    warning("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (is.null(motif) || !is.character(motif)) {
    stop("Please provide a valid motif name")
  }
  if (is.null(msites) || !is(msites, "GRanges")) {
    stop("Please provide a valid methylation sites with read_methylome function")
  }
  if (is.null(tf_bindsites) || !any(c(!is(tf_bindsites, "GRangesList") ||
    !is.list(tf_bindsites)))) {
    stop("Please provide a valid tf binding sites as GRangesList")
  }
  if (is.null(gcfreqs) || !is.list(gcfreqs)) {
    stop("Please provide a valid GC bin frequency table list")
  }
  if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
    stop("Please provide a valid enhancer regions")
  }
  tfbs <- tf_bindsites[[motif]]
  gcfreq <- gcfreqs[[motif]]
  tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(tfbs, enhancer, ignore.strand = ignoreStrand)
    tfbs <- tfbs[d_hits@from]
  }
  exp_dev <- unlist(exp_dev)
  hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = ignoreStrand)
  if (length(hits@from) == 0) {
    stop(paste0("No methylation sites found in the", motif, " binding sites"))
  }
  #binMsites <- addGCBintoMethylome(msites, gc_dist,ignoreStrand, 10000,threads,gcfreq) 
  #exp_meth <- computeExpectations(binMsites, gcfreq)
  # WARNING: This is a trial implementation
  #exp_meth <- lapply(binMsites, function(x) computeExpectations(x, gcfreq))#, mc.cores = threads)

  S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) + ((start(tfbs) - end(tfbs)) / 2))
  sum_meth <- data.table::data.table(
    x = start(msites[hits@from]) - tfbs[hits@to]$mid_point,
    avg_methyl = msites[hits@from]$score
  )
  # GC bias correction
  obs_dev <- dev_helper(sum_meth)
  #exp_dev <- dev_helper(exp_meth)
  #dev <- (obs_dev - exp_dev)

  # WARNING: This is a trial implementation
  dev <- (obs_dev - mean(exp_dev))/sd(exp_dev)
  obs_dev <- (obs_dev - mean(exp_dev))
  return(data.table::data.table(obs_dev, dev)) 
  # WARNING: Changed exp_dev to obs_dev
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
