#' @title computeFootprint
#' @description Observed methylation profile around the binding sites of
#' one motif.
#' @inheritParams computeDeviation
#' @param motifName Name of the motif.
#' @return a \code{data.table} with one row per distance to the site
#' centre.
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols
#' @import data.table
computeFootprint <- function(
    motifName, tfBindsites, msites, enhancer = NULL, extendBy = 130
) {
    tfbs <- tfBindsites[[motifName]]
    tfbs <- resize(tfbs, width(tfbs)[1] + extendBy, fix = "center")

    if (!is.null(enhancer)) {
        tfbs <- subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE)
    }

    hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = TRUE)
    S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) +
        (start(tfbs) - end(tfbs)) / 2)
    x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
    plotData <- data.table(
        x = x,
        y1 = msites[hits@from]$score,
        y2 = msites[hits@from]$coverage
    )

    plotData <- plotData[, .(
        n = .N,
        avg_methyl = mean(y1),
        avg_cov = mean(y2),
        motif = motifName,
        label = fifelse(
            x == max(x),
            as.character(motifName), NA_character_
        )
    ), by = x]

    return(plotData)
}

#' @title computeExpectedFootprint
#' @description Methylation profile expected from GC content around the
#' binding sites of one motif.
#' @inheritParams computeDeviation
#' @return a \code{data.table} with the columns \code{x} and
#' \code{avg_methyl}.
#' @keywords internal
#' @importFrom IRanges subsetByOverlaps
#' @import data.table
computeExpectedFootprint <- function(
    motif, gcfreqs, gcDist, enhancer = NULL, msites
) {
    gcfreq <- gcfreqs[[motif]]

    if (!is.null(enhancer)) {
        gcDist <- subsetByOverlaps(gcDist, enhancer, ignore.strand = TRUE)
    }
    binMeth <- addGCBinToMethylome(msites, gcDist)
    return(computeExpectations(binMeth, gcfreq))
}
