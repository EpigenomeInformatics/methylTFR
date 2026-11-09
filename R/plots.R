#' @title checkFootprintInputs
#' @description Validate the inputs shared by the two footprint plotting
#' functions.
#' @inheritParams computeDeviation
#' @param motif Motif name as a character string.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
checkFootprintInputs <- function(
    motif, tfBindsites, msites, gcDist, gcfreqs, enhancer = NULL
) {
    checkMsites(msites)
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }
    checkTfBindsites(tfBindsites)
    checkEnhancer(enhancer)
    checkGcDist(gcDist)
    checkGcfreqs(gcfreqs)
    invisible(NULL)
}

#' @title expectedFootprintPlot
#' @description Draw the observed and expected methylation profiles.
#' @param combinedData A \code{data.table} with the columns \code{x},
#' \code{avg_methyl} and \code{type}.
#' @param motif Motif name as a character string.
#' @param sampleName Sample label used in the title.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab ggtitle
#' @importFrom ggplot2 theme_classic scale_color_manual theme
#' @keywords internal
expectedFootprintPlot <- function(combinedData, motif, sampleName) {
    ggplot(combinedData, aes(
        x = x,
        y = avg_methyl,
        color = type
    )) +
        geom_line() +
        geom_point() +
        xlab("Distance from motif center") +
        ylab("Methylation level") +
        theme_classic() +
        ggtitle(paste("TF footprint for", motif, "in", sampleName)) +
        scale_color_manual(values = c(
            "Expected" = "blue",
            "Observed" = "red"
        )) +
        theme(legend.position = "bottom")
}

#' @title plotExpectedFootprint
#' @description Plot the observed methylation profile around the binding
#' sites of one motif together with the profile expected from GC content.
#' @inheritParams computeDeviation
#' @param motif Name of the motif, one of \code{names(tfBindsites)}.
#' @param sampleName Optional sample label used in the plot title.
#' Default \code{"sample"}.
#' @param returnPlotData Logical. If \code{TRUE}, return the plot data as
#' well as the plot.
#' @return A \code{ggplot} object, or, if \code{returnPlotData} is
#' \code{TRUE}, a list with the elements \code{plot} and \code{plotDF}.
#' @seealso \code{\link{plotMotifFootprint}} for the bias-corrected
#' profile, \code{\link{computeDeviation}} for the score it summarises.
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_line ggtitle theme_classic
#' @examples
#' load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#' p <- plotExpectedFootprint(
#'     motif = "BATF",
#'     tfBindsites = tf_bindsites,
#'     msites = msites,
#'     sampleName = "Sample1",
#'     gcDist = gcdist,
#'     gcfreqs = gcfreqs,
#'     enhancer = NULL,
#'     returnPlotData = FALSE
#' )
plotExpectedFootprint <- function(
    motif, tfBindsites, msites,
    sampleName = NULL, gcDist, gcfreqs,
    enhancer = NULL, returnPlotData = FALSE
) {
    if (is.null(sampleName)) {
        sampleName <- "sample"
    }
    checkFootprintInputs(
        motif, tfBindsites, msites, gcDist, gcfreqs, enhancer
    )

    plotData <- computeFootprint(motif, tfBindsites, msites, enhancer)
    expData <- computeExpectedFootprint(
        motif, gcfreqs, gcDist, enhancer, msites
    )

    expData[, type := "Expected"]
    plotData[, type := "Observed"]
    combinedData <- rbindlist(list(
        expData[, .(x, avg_methyl, type)],
        plotData[, .(x, avg_methyl, type)]
    ))

    p1 <- expectedFootprintPlot(combinedData, motif, sampleName)
    if (returnPlotData) {
        return(list(plot = p1, plotDF = combinedData))
    }
    return(p1)
}

#' @title footprintDifference
#' @description Combine the observed and expected profiles into a single
#' corrected curve.
#' @param plotData A \code{data.table} with the columns \code{x},
#' \code{avg_methyl} and \code{type}.
#' @param method Either \code{"substraction"} or \code{"division"}.
#' @param sampleName Sample label used in the curve label.
#' @return A \code{list} with the corrected \code{data} and the axis
#' \code{lab}.
#' @keywords internal
footprintDifference <- function(plotData, method, sampleName) {
    if (method == "substraction") {
        differenceData <- plotData[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
        ), by = x]
        differenceData[, type := paste("Obs. sub. Exp.", sampleName)]
        lab <- "(Observed - Expected)"
    } else {
        differenceData <- plotData[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]
        ), by = x]
        differenceData[, type := paste("Obs. div. Exp.", sampleName)]
        lab <- "(Observed / Expected)"
    }
    return(list(data = differenceData, lab = lab))
}

#' @title normaliseFootprintFlank
#' @description Normalise a corrected footprint against its outer flanking
#' windows, on the same scale as the statistic itself.
#' @details Dividing a difference by its flank mean would rescale it by an
#' arbitrary factor, because that mean is near zero: for a typical footprint
#' it inflates the curve by an order of magnitude and can invert its sign.
#' @param differenceData The corrected footprint as a \code{data.table}.
#' @param method Either \code{"substraction"} or \code{"division"}.
#' @param flankNorm Width of the flanking window used for normalisation.
#' @return The normalised \code{data.table}.
#' @keywords internal
normaliseFootprintFlank <- function(differenceData, method, flankNorm) {
    if (is.null(flankNorm) || flankNorm <= 0) {
        return(differenceData)
    }
    flank <- max(abs(differenceData$x), na.rm = TRUE)
    idx <- abs(differenceData$x) >= flank - flankNorm
    normFactor <- mean(differenceData$avg_methyl[idx], na.rm = TRUE)
    if (method == "substraction") {
        differenceData[, avg_methyl := avg_methyl - normFactor]
    } else {
        differenceData[, avg_methyl := avg_methyl / normFactor]
    }
    return(differenceData)
}

#' @title plotMotifFootprint
#' @description Plot the bias-corrected methylation profile around the
#' binding sites of one motif in one sample.
#' @inheritParams plotExpectedFootprint
#' @param method How observed and expected profiles are combined, either
#' \code{"substraction"} (observed minus expected) or \code{"division"}
#' (observed over expected, the default).
#' @param flankNorm Width in base pairs of the outer flanking window used
#' to normalise the curve. Default 50; \code{0} turns normalisation off.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotExpectedFootprint}} for the two profiles
#' before correction.
#' @importFrom ggplot2 aes xlab ylab theme xlim
#' @importFrom ggplot2 ggplot geom_line ggtitle theme_classic
#' @export
#' @examples
#' load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#' p <- plotMotifFootprint(
#'     motif = "BATF",
#'     tfBindsites = tf_bindsites,
#'     msites = msites,
#'     sampleName = "Sample1",
#'     gcDist = gcdist,
#'     gcfreqs = gcfreqs,
#'     enhancer = NULL,
#'     method = "division"
#' )
plotMotifFootprint <- function(
    motif, tfBindsites, msites,
    sampleName = NULL, gcDist, gcfreqs,
    enhancer = NULL, method = "division",
    flankNorm = 50
) {
    if (is.null(method) || !method %in% c("substraction", "division")) {
        method <- "division"
        warning("Invalid method, using the default method division.")
    }
    if (is.null(sampleName)) {
        sampleName <- "sample"
    }
    plotData <- plotExpectedFootprint(
        motif, tfBindsites, msites,
        sampleName = sampleName,
        gcDist = gcDist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE
    )$plotDF

    diff <- footprintDifference(plotData, method, sampleName)
    differenceData <- normaliseFootprintFlank(diff$data, method, flankNorm)

    ggplot(differenceData, aes(
        x = x,
        y = avg_methyl,
        color = type
    )) +
        geom_line() +
        xlab("Distance from motif center") +
        ylab(paste0("Methylation difference ", diff$lab)) +
        theme_classic() +
        ggtitle(paste(
            "TF footprint difference for", motif, "in", sampleName
        )) +
        theme(legend.position = "bottom") +
        xlim(-200, 200)
}
