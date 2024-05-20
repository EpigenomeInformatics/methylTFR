#' @title computeExpectations helper
#' @description Compute the expected methylation score for each GC bin
#' @param msites - GRanges object - methylation sites
#' @param gcdist - GRanges object - GC distribution
#' @param ignoreStrand - logical - ignore strand information
#' @param gcShuffleSize - integer - number of GC sites to sample
#' @param gcfreq - matrix - GC frequency matrix
#' @return - list of matrices
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table rbindlist
#' @importFrom dplyr sample_n
#' @keywords internal

compute_exp_meth <- function(msites, gcdist, gcShuffleSize, ignoreStrand, gcfreq) {
  n_samples <- round(gcfreq * gcShuffleSize)
  sampled_gcdist <- lapply(seq_along(n_samples), function(i) {
    sample_n(gcdist[gcdist$GC_bin == i], n_samples[i], replace = FALSE, prob = NULL)
  })
  sampled_gcdist <- rbindlist(sampled_gcdist)
  sampled_gcdist <- GRanges(
    seqnames = sampled_gcdist$seqnames,
    ranges = IRanges(
      start = sampled_gcdist$start,
      end = sampled_gcdist$end
    ),
    strand = sampled_gcdist$strand,
    GC_bias = sampled_gcdist$GC_bias,
    GC_bin = sampled_gcdist$GC_bin
  )
  hits <- findOverlaps(msites, sampled_gcdist, type = "within", ignore.strand = ignoreStrand)

  if (length(hits@from) == 0) {
    return(as.matrix(data.frame(gcbin = 1:5, avg_mscore = rep(0, 5))))
  }
  gcmap <- data.table(
    mscore = msites[hits@from]$score,
    gcbin = sampled_gcdist[hits@to]$GC_bin
  )
  exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]

  # Set gcbin as a key column in exp_meth
  setkey(exp_meth, gcbin)

  # Join exp_meth with a data table of all gcbin values and replace NA with 0
  exp_meth <- exp_meth[.(1:5), on = .(gcbin)]
  exp_meth[is.na(exp_meth)] <- 0
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
#' @keywords internal
#' @author Irem Gunduz
addGCBintoMethylome <- function(
    msites, gcdist, ignoreStrand = TRUE,
    sample_size = 25000, gcfreq) {
  if (!is.logical(ignoreStrand)) {
    warning("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (is.null(msites) || !is(msites, "GRanges")) {
    stop("Please provide a valid methylation sites with read_methylome function")
  }
  # Find overlaps between msites and gcdist
  # hits <- findOverlaps(msites, gcdist, type = "within", ignore.strand = ignoreStrand)
  # gcdist <- gcdist[hits@to] # To assure every GC_bin is present in the final result
  exp_meth_list <- lapply(
    251:(NCOL(gcfreq) - 250), # This will focus on motif center
    function(x) {
      compute_exp_meth(
        msites, gcdist, sample_size,
        ignoreStrand, gcfreq[, x]
      )
    }
  ) # This will calculate the expected methylation for each GC bin
  return(exp_meth_list)
}

#' @title computeExpectations
#' @description  This function is used to calculate expected methylation for a given
#' motif and sample.
#' @param binMsites_sub Imported methylation sites with GC bin
#' @param gcfreq a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @importFrom logger log_info log_warn log_error
#' @keywords internal
computeExpectations <- function(binMsites_sub, gcfreq) {
  if (!is.matrix(binMsites_sub)) {
    stop("Please provide a valid GC bin frequency table as a matrix")
  }
  if (!is.matrix(gcfreq)) {
    stop("Please provide a valid GC bin frequency table as a matrix")
  }
  exp.data <- t(gcfreq) %*% binMsites_sub[, 2]
  mpos <- round(seq(-floor(length(exp.data) / 2), floor(length(exp.data) / 2),
    length.out = length(exp.data)
  ))
  exp.methyl <- data.table(x = mpos, avg_methyl = exp.data)
  colnames(exp.methyl) <- c("x", "avg_methyl")
  return(exp.methyl)
}

#' @title computeExpectedDeviation
#' @description computeExpectedDeviation is a function to calculate the expected deviation in transcription factor
#' binding sites for a given motif
#' @param motif A character vector containing motif name
#' @param msites imported methylation sites
#' @param gcfreqs a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param sample_size number of GC sites to sample, default is 25000
#' @return a \code{numeric} deviation score for a given motif
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table setDT
#' @export
computeExpectedDeviation <- function(motif, msites, gcfreqs, gc_dist,
                                     ignoreStrand = TRUE, sample_size = 25000) {
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
  gcfreq <- gcfreqs[[motif]]
  binMsites <- addGCBintoMethylome(msites, gc_dist, ignoreStrand, sample_size, gcfreq)
  exp_meth <- lapply(binMsites, function(i) computeExpectations(i, gcfreq))
  exp_dev <- lapply(exp_meth, dev_helper)
  return(unlist(exp_dev))
}
