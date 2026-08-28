#' @title addGCBintoMethylome
#' @description  This function is used to add GC bin
#'  to the methylation sites
#' and calculate the average methylation for each bin.
#' @param msites Imported methylation sites using
#' \code{read_methylome} function
#' @param gcdist a \code{GRanges} object contains
#'  Genome wide GC distribution
#' @param ignoreStrand - if TRUE, it ignores strand
#'  info from annotation
#' @return a \code{matrix} object with GC bin with
#' corresponding avg methylation
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @importFrom logger log_info log_warn log_error
#' @export
#' @examples
#' library(methylTFR)
#'
#' # Load the data
#' load(system.file("extdata",
#'     "gcdist_subset.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "example_data.rda",
#'     package = "methylTFR"
#' ))
#'
#' # Add GC bin
#' bin_meth <- addGCBintoMethylome(msites, gcdist)
#' @author Irem Gunduz
addGCBintoMethylome <- function(
    msites,
    gcdist,
    ignoreStrand = TRUE
) {
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    if (is.null(msites) || !is(msites, "GRanges")) {
        stop(
            "Please provide a valid methylation sites with ",
            "read_methylome function"
        )
    }
    if (is.null(gcdist) || !is(gcdist, "GRanges")) {
        stop("Please provide a valid GC distribution")
    }
    hits <- findOverlaps(msites, gcdist,
        ignore.strand = ignoreStrand
    )
    if (length(hits@from) == 0) {
        stop("No methylation sites found in the GC distribution")
    }
    gcmap <- data.table(
        mscore = msites[hits@from]$score,
        gcbin = gcdist[hits@to]$GC_bin
    )
    exp_meth <- gcmap[, .(avg_mscore = mean(mscore)),
        by = gcbin
    ]
    exp_meth <- exp_meth[order(gcbin)]
    return(as.matrix(exp_meth))
}


#' @title computeExpectations
#' @description  This function is used to calculate expected
#'  methylation for a given
#' motif and sample.
#' @param binMsites Imported methylation sites with GC bin
#' @param gcfreq a \code{list} of GC bin frequency tables
#' (matrices for multiple motif)
#' @return a \code{data.table} object with GC bin with
#' corresponding avg methylation
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
    mpos <- round(seq(-floor(length(exp.data) / 2),
        floor(length(exp.data) / 2),
        length.out = length(exp.data)
    ))
    exp.methyl <- data.table(
        x = mpos,
        avg_methyl = exp.data
    )
    colnames(exp.methyl) <- c("x", "avg_methyl")
    return(exp.methyl)
}
