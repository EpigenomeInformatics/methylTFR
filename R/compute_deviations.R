#' @title check_deviation_inputs
#' @description Validate the inputs of \code{computeDeviation}.
#' @inheritParams computeDeviation
#' @param motif Motif name as a character string.
#' @param extendBy Number of base pairs added to each binding site.
#' @param windowBreaks Break points of the distance bins.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
check_deviation_inputs <- function(
    motif, msites, tf_bindsites, gcfreqs, binMsites, enhancer = NULL,
    extendBy = 130, windowBreaks = c(-250, -200, -25, 25, 200, 250)
) {
    if (is.null(motif) || !is.character(motif) || length(motif) != 1) {
        stop("Please provide a valid motif name")
    }
    check_msites(msites)
    check_tf_bindsites(tf_bindsites)
    check_gcfreqs(gcfreqs)
    check_enhancer(enhancer)
    if (!motif %in% names(tf_bindsites)) {
        stop("Motif ", motif, " is not in names(tf_bindsites)")
    }
    if (!motif %in% names(gcfreqs)) {
        stop("Motif ", motif, " is not in names(gcfreqs)")
    }
    if (!is.matrix(binMsites) ||
        !all(c("gcbin", "avg_mscore") %in% colnames(binMsites))) {
        stop(
            "binMsites must be the matrix returned by ",
            "addGCBintoMethylome()"
        )
    }
    if (!is.numeric(extendBy) || length(extendBy) != 1 || extendBy < 0) {
        stop("extendBy must be a single non-negative number")
    }
    if (!is.numeric(windowBreaks) || length(windowBreaks) != 6 ||
        is.unsorted(windowBreaks, strictly = TRUE)) {
        stop("windowBreaks must be six strictly increasing numbers")
    }
    invisible(NULL)
}

#' @title computeDeviation
#' @description Compute the bias-corrected deviation score of one motif in
#' one sample.
#' @details
#' Each binding site is widened by \code{extendBy} base pairs and every CpG
#' inside it is placed by its distance to the site centre. The distances
#' are cut into five bins at \code{windowBreaks}: an outer left flank, an
#' inner left region, the centre, an inner right region and an outer right
#' flank. The observed score is the mean methylation in the centre bin
#' divided by the mean of the two outer flank bins. The same ratio is
#' computed for the methylation expected from GC content (see
#' \code{\link{addGCBintoMethylome}}), and the deviation is the observed
#' minus the expected ratio. If the centre or either outer flank bin holds
#' no CpG, the score is \code{NA}.
#'
#' Every CpG counts once, whatever its read coverage: the bin means are
#' unweighted means of the per-CpG methylation fractions. Coverage only
#' enters through the \code{cov_threshold} filter applied when the data are
#' read.
#' @inheritParams run_methyltfr
#' @param msites A \code{GRanges} object with one range per CpG (width 1)
#' and the metadata columns \code{score} (numeric methylation fraction
#' between 0 and 1) and \code{coverage} (read count). Usually produced by
#' \code{\link{read_methylome}}.
#' @param binMsites The matrix returned by \code{\link{addGCBintoMethylome}},
#' with the columns \code{gcbin} and \code{avg_mscore}.
#' @param motif Name of the motif, one of \code{names(tf_bindsites)}.
#' @param extendBy Number of base pairs added to the width of each binding
#' site before overlapping it with the methylation sites, split evenly
#' between both sides. Default 130.
#' @param windowBreaks Six strictly increasing distances from the site
#' centre, in base pairs, that define the five distance bins. Bins 1 and 5
#' are the outer flanks and bin 3 is the centre. Default
#' \code{c(-250, -200, -25, 25, 200, 250)}.
#' @return A \code{data.table} with one row and the columns \code{dev}
#' (bias-corrected deviation) and \code{exp_dev} (expected deviation).
#' @seealso \code{\link{run_methyltfr}} to score all motifs and samples at
#' once, \code{\link{plotMotifFootprint}} to see the profile behind the
#' score.
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom data.table data.table setDT
#' @importFrom S4Vectors mcols
#' @import data.table
#' @examples
#' load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#' bin_meth <- addGCBintoMethylome(msites, gcdist, TRUE)
#'
#' computeDeviation(
#'     motif = "BATF",
#'     msites = msites,
#'     tf_bindsites = tf_bindsites,
#'     gcfreqs = gcfreqs,
#'     enhancer = NULL,
#'     ignoreStrand = TRUE,
#'     binMsites = bin_meth
#' )
#' @export
computeDeviation <- function(
    motif, msites, tf_bindsites, gcfreqs, enhancer = NULL,
    ignoreStrand = TRUE, binMsites, extendBy = 130,
    windowBreaks = c(-250, -200, -25, 25, 200, 250)
) {
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    check_deviation_inputs(
        motif, msites, tf_bindsites, gcfreqs, binMsites, enhancer,
        extendBy, windowBreaks
    )
    tfbs <- tf_bindsites[[motif]]
    tfbs <- resize(tfbs, width(tfbs)[1] + extendBy, fix = "center")
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
        stop(
            "No methylation sites found in the ",
            motif, " binding sites"
        )
    }
    exp_meth <- computeExpectations(binMsites, gcfreq)

    S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) +
        ((start(tfbs) - end(tfbs)) / 2))
    sum_meth <- data.table(
        x = start(msites[hits@from]) - tfbs[hits@to]$mid_point,
        avg_methyl = msites[hits@from]$score
    )

    # GC bias correction
    obs_dev <- dev_helper(sum_meth, windowBreaks)
    exp_dev <- dev_helper(exp_meth, windowBreaks)
    dev <- obs_dev - exp_dev
    return(data.table(dev, exp_dev))
}

#' @title dev_helper
#' @description Ratio of the mean methylation in the centre bin to the mean
#' of the two outer flank bins.
#' @param data A \code{data.table} with the columns \code{x} (distance to
#' the site centre) and \code{avg_methyl}.
#' @param breaks Six strictly increasing break points defining five bins.
#' @return A numeric value, or \code{NA} if the centre bin or either outer
#' flank bin is empty.
#' @keywords internal
dev_helper <- function(data, breaks = c(-250, -200, -25, 25, 200, 250)) {
    data[, cuts := cut(x, breaks, labels = FALSE)]
    bin_means <- data[!is.na(cuts), .(avg = mean(avg_methyl)), by = cuts]
    centre <- bin_means[cuts == 3L, avg]
    left <- bin_means[cuts == 1L, avg]
    right <- bin_means[cuts == 5L, avg]
    if (!length(centre) || !length(left) || !length(right)) {
        return(NA_real_)
    }
    centre / ((left + right) / 2)
}
