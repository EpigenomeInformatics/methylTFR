#' @title addGCBinToMethylome
#' @description Assign each methylation site to the GC bin of the genomic
#' window it falls in, and compute the average methylation per GC bin.
#' @details The result describes how methylation depends on local GC
#' content in this sample. \code{\link{computeDeviation}} combines it with
#' the per-motif GC profile in \code{gcfreqs} to obtain the methylation
#' expected at the binding sites from sequence composition alone.
#' @inheritParams computeDeviation
#' @param gcdist A \code{GRanges} object of 30-bp genomic windows with an
#' integer metadata column \code{GC_bin} (1 to 5). Get it with
#' \code{getGenomeGC()} from the annotation package.
#' @return A \code{matrix} with one row per GC bin and the columns
#' \code{gcbin} and \code{avg_mscore}. Pass it as \code{binMsites} to
#' \code{\link{computeDeviation}}.
#' @seealso \code{\link{computeDeviation}}, \code{\link{readMethylome}}
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @export
#' @examples
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#'
#' bin_meth <- addGCBinToMethylome(msites, gcdist)
#' bin_meth
#' @author Irem Gunduz
addGCBinToMethylome <- function(msites, gcdist, ignoreStrand = TRUE) {
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    check_msites(msites)
    check_gc_dist(gcdist)
    hits <- findOverlaps(msites, gcdist, ignore.strand = ignoreStrand)
    if (length(hits@from) == 0) {
        stop("No methylation sites found in the GC distribution")
    }
    gcmap <- data.table(
        mscore = msites[hits@from]$score,
        gcbin = gcdist[hits@to]$GC_bin
    )
    exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]
    exp_meth <- exp_meth[order(gcbin)]
    return(as.matrix(exp_meth))
}


#' @title computeExpectations
#' @description Expected methylation profile of a motif in one sample.
#' @param binMsites The matrix returned by \code{addGCBinToMethylome}.
#' @param gcfreq The GC bin frequency matrix of one motif.
#' @return a \code{data.table} with the columns \code{x} (distance to the
#' site centre) and \code{avg_methyl}.
#' @importFrom data.table data.table
#' @keywords internal
computeExpectations <- function(binMsites, gcfreq) {
    if (!is.matrix(binMsites)) {
        stop("Please provide a valid GC bin frequency table as a matrix")
    }
    if (!is.matrix(gcfreq)) {
        stop("Please provide a valid GC bin frequency table as a matrix")
    }
    exp_data <- t(gcfreq) %*% binMsites[, 2]
    mpos <- round(seq(-floor(length(exp_data) / 2),
        floor(length(exp_data) / 2),
        length.out = length(exp_data)
    ))
    exp.methyl <- data.table(x = mpos, avg_methyl = exp_data)
    colnames(exp.methyl) <- c("x", "avg_methyl")
    return(exp.methyl)
}
