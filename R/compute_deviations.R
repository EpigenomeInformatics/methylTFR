#' @title computeDeviation
#' @description computeDeviation is a function to calculate
#'  the deviation in transcription factor
#' binding sites for a given motif
#' @param motif A character vector containing motif name
#' @param msites imported methylation sites
#' @param tf_bindsites a \code{GRangesList} object contains
#'  tf binding sites positions
#' @param enhancer  a \code{GRanges} object specifying regions
#' such as distal motif (optional)
#' @param ignoreStrand if TRUE, it ignores strand info from
#' annotation
#' @param binMsites A matrix object with GC bin with corresponding
#'  avg methylation
#' @param gcfreqs a \code{list} of GC bin frequency tables
#'  (matrices for multiple motif)
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
#' load(system.file("extdata",
#'     "BATF_tf_bindsites.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "example_data.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "BATF_gcfreqs.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "gcdist_subset.rda",
#'     package = "methylTFR"
#' ))
#'
#' # Compute binMsites
#' bin_meth <- addGCBintoMethylome(msites, gcdist, TRUE)
#'
#' # Compute the deviation
#' devs <- computeDeviation("BATF",
#'     msites,
#'     tf_bindsites,
#'     gcfreqs,
#'     enhancer = NULL,
#'     ignoreStrand = TRUE,
#'     bin_meth
#' )
#' @export
computeDeviation <- function(
        motif, msites, tf_bindsites, gcfreqs,
        enhancer = NULL, ignoreStrand = TRUE,
        binMsites) {
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }
    if (is.null(msites) || !is(msites, "GRanges")) {
        stop("Please provide a valid methylation
        sites with read_methylome function")
    }
    if (is.null(tf_bindsites) ||
        !any(c(!is(tf_bindsites, "GRangesList") ||
            !is.list(tf_bindsites)))) {
        stop("Please provide a valid tf binding sites as GRangesList")
    }
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("Please provide a valid enhancer regions")
    }
    tfbs <- tf_bindsites[[motif]]
    tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
    gcfreq <- gcfreqs[[motif]]
    if (!is.null(enhancer)) {
        tfbs <- subsetByOverlaps(tfbs, enhancer,
            ignore.strand = ignoreStrand
        )
    }
    hits <- findOverlaps(msites, tfbs,
        type = "within",
        ignore.strand = ignoreStrand
    )
    if (length(hits@from) == 0) {
        stop(paste0(
            "No methylation sites found in the",
            motif, " binding sites"
        ))
    }
    exp_meth <- computeExpectations(binMsites, gcfreq)

    S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) +
        ((start(tfbs) - end(tfbs)) / 2))
    sum_meth <- data.table(
        x = start(msites[hits@from]) - tfbs[hits@to]$mid_point,
        avg_methyl = msites[hits@from]$score
    )

    # GC bias correction
    obs_dev <- dev_helper(sum_meth)
    exp_dev <- dev_helper(exp_meth)
    dev <- obs_dev - exp_dev
    return(data.table(dev, exp_dev))
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
