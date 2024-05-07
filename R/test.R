#' @title computeExpectations helper
#' @description Compute the expected methylation score for each GC bin
#' @param msites - GRanges object - methylation sites
#' @param gcdist - GRanges object - GC distribution
#' @param ignoreStrand - logical - ignore strand information
#' @param gcShuffleSize - integer - number of GC sites to sample
#' @param gcfreq - matrix - GC frequency matrix
#' @return - list of matrices
#' @keywords internal
compute_exp_meth <- function(msites, gcdist, gcShuffleSize, ignoreStrand, gcfreq) {
  sampled_gcdist <- sample(gcdist, size = gcShuffleSize, replace = TRUE, prob = gcfreq)
  hits <- findOverlaps(msites, sampled_gcdist, type = "within", ignore.strand = ignoreStrand)
  if (length(hits@from) == 0) {
    stop("No methylation sites found in the GC distribution")
  }
  gcmap <- data.table(
    mscore = msites[hits@from]$score,
    gcbin = sampled_gcdist[hits@to]$GC_bin
  )
  exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]
  exp_meth <- exp_meth[order(gcbin)]
  return(as.matrix(exp_meth))
}

#' @title addGCBintoMethylome
#' @description  This function is used to add GC bin to the methylation sites
#' and calculate the average methylation for each bin.
#' @param msites Imported methylation sites using \code{read_methylome} function
#' @param gcdist a \code{GRanges} object contains Genome wide GC distribution
#' @param ignoreStrand - if TRUE, it ignores strand info from annotation
#' @param sample_size - number of GC sites to sample, default is 10000
#' @param gcfreq - a list of GC bin frequency tables (matrices for multiple motif)
#' @return a list of matrices
#' @importFrom logger log_info log_warn log_error
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom parallel mclapply
#' @export
#' @author Irem Gunduz
addGCBintoMethylome <- function(msites, gcdist, ignoreStrand = TRUE,
 sample_size = 10000,gcfreq) {
  if (!is.logical(ignoreStrand)) {
    warning("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (is.null(msites) || !is(msites, "GRanges")) {
    stop("Please provide a valid methylation sites with read_methylome function")
  }
  if (is.null(gcdist) || !is(gcdist, "GRanges")) {
    stop("Please provide a valid GC distribution")
  } 
  #gcfreq_resized <- gcfreq[gcdist$GC_bin] # This will select only the first column
  exp_meth_list <- lapply(250:(NCOL(gcfreq) - 250),  # This will focus on motif center
  function(x) compute_exp_meth(msites, gcdist, sample_size,
  ignoreStrand,gcfreq[gcdist$GC_bin,x])) # Test if this makes function slower
  return(exp_meth_list)
}

#' @title computeExpectedDeviation
#' @description computeExpectedDeviation is a function to calculate the expected deviation in transcription factor
#' binding sites for a given motif
#' @param motif A character vector containing motif name
#' @param msites imported methylation sites
#' @param gcfreqs a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param enhancer  a \code{GRanges} object specifying regions such as distal motif (optional)
#' @return a \code{numeric} deviation score for a given motif
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table setDT
#' @export 
computeExpectedDeviation <- function(motif, msites, gcfreqs, gc_dist,
                             ignoreStrand = TRUE,enhancer=NULL) {
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
  if (is.null(gcfreqs) || !is.list(gcfreqs)) {
    stop("Please provide a valid GC bin frequency table list")
  }
  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(gc_dist, enhancer, ignore.strand = ignoreStrand)
    gc_dist <- gc_dist[d_hits@from]
  }
  gcfreq <- gcfreqs[[motif]]
  binMsites <- addGCBintoMethylome(msites, gc_dist,ignoreStrand,10000,gcfreq) 
  exp_meth <- lapply(binMsites, function(x) computeExpectations(x, gcfreq))
  exp_dev <- lapply(exp_meth, dev_helper)
  return(unlist(exp_dev))
}