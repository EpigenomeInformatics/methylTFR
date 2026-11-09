#' @title computeFootprint
#' @description Observed methylation profile around the binding sites of
#' one motif.
#' @inheritParams computeDeviation
#' @param motif_name Name of the motif.
#' @return a \code{data.table} with one row per distance to the site
#' centre.
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols
#' @import data.table
computeFootprint <- function(
    motif_name, tf_bindsites, msites, enhancer = NULL, extendBy = 130
) {
    tfbs <- tf_bindsites[[motif_name]]
    tfbs <- resize(tfbs, width(tfbs)[1] + extendBy, fix = "center")

    if (!is.null(enhancer)) {
        tfbs <- subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE)
    }

    hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = TRUE)
    S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) +
        (start(tfbs) - end(tfbs)) / 2)
    x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
    plot_data <- data.table(
        x = x,
        y1 = msites[hits@from]$score,
        y2 = msites[hits@from]$coverage
    )

    plot_data <- plot_data[, .(
        n = .N,
        avg_methyl = mean(y1),
        avg_cov = mean(y2),
        motif = motif_name,
        label = fifelse(
            x == max(x),
            as.character(motif_name), NA_character_
        )
    ), by = x]

    return(plot_data)
}

#' @title computeExpectedFootprint
#' @description Methylation profile expected from GC content around the
#' binding sites of one motif.
#' @inheritParams computeDeviation
#' @inheritParams run_methyltfr
#' @return a \code{data.table} with the columns \code{x} and
#' \code{avg_methyl}.
#' @keywords internal
#' @importFrom IRanges subsetByOverlaps
#' @import data.table
computeExpectedFootprint <- function(
    motif, gcfreqs, gc_dist, enhancer = NULL, msites
) {
    gcfreq <- gcfreqs[[motif]]

    if (!is.null(enhancer)) {
        gc_dist <- subsetByOverlaps(gc_dist, enhancer, ignore.strand = TRUE)
    }
    bin_meth <- addGCBintoMethylome(msites, gc_dist)
    return(computeExpectations(bin_meth, gcfreq))
}
