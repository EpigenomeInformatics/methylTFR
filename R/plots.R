#' @title check_footprint_inputs
#' @description Validate the inputs shared by the two footprint plotting
#' functions.
#' @inheritParams computeDeviation
#' @inheritParams runMethylTFR
#' @param motif Motif name as a character string.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
check_footprint_inputs <- function(
    motif, tf_bindsites, msites, gc_dist, gcfreqs, enhancer = NULL
) {
    check_msites(msites)
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }
    check_tf_bindsites(tf_bindsites)
    check_enhancer(enhancer)
    check_gc_dist(gc_dist)
    check_gcfreqs(gcfreqs)
    invisible(NULL)
}

#' @title expected_footprint_plot
#' @description Draw the observed and expected methylation profiles.
#' @param combined_data A \code{data.table} with the columns \code{x},
#' \code{avg_methyl} and \code{type}.
#' @param motif Motif name as a character string.
#' @param sample_name Sample label used in the title.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab ggtitle
#' @importFrom ggplot2 theme_classic scale_color_manual theme
#' @keywords internal
expected_footprint_plot <- function(combined_data, motif, sample_name) {
    ggplot(combined_data, aes(
        x = x,
        y = avg_methyl,
        color = type
    )) +
        geom_line() +
        geom_point() +
        xlab("Distance from motif center") +
        ylab("Methylation level") +
        theme_classic() +
        ggtitle(paste("TF footprint for", motif, "in", sample_name)) +
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
#' @inheritParams runMethylTFR
#' @param motif Name of the motif, one of \code{names(tf_bindsites)}.
#' @param sample_name Optional sample label used in the plot title.
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
#'     tf_bindsites = tf_bindsites,
#'     msites = msites,
#'     sample_name = "Sample1",
#'     gc_dist = gcdist,
#'     gcfreqs = gcfreqs,
#'     enhancer = NULL,
#'     returnPlotData = FALSE
#' )
plotExpectedFootprint <- function(
    motif, tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, returnPlotData = FALSE
) {
    if (is.null(sample_name)) {
        sample_name <- "sample"
    }
    check_footprint_inputs(
        motif, tf_bindsites, msites, gc_dist, gcfreqs, enhancer
    )

    plot_data <- computeFootprint(motif, tf_bindsites, msites, enhancer)
    exp_data <- computeExpectedFootprint(
        motif, gcfreqs, gc_dist, enhancer, msites
    )

    exp_data[, type := "Expected"]
    plot_data[, type := "Observed"]
    combined_data <- rbindlist(list(
        exp_data[, .(x, avg_methyl, type)],
        plot_data[, .(x, avg_methyl, type)]
    ))

    p1 <- expected_footprint_plot(combined_data, motif, sample_name)
    if (returnPlotData) {
        return(list(plot = p1, plotDF = combined_data))
    }
    return(p1)
}

#' @title footprint_difference
#' @description Combine the observed and expected profiles into a single
#' corrected curve.
#' @param plot_data A \code{data.table} with the columns \code{x},
#' \code{avg_methyl} and \code{type}.
#' @param method Either \code{"substraction"} or \code{"division"}.
#' @param sample_name Sample label used in the curve label.
#' @return A \code{list} with the corrected \code{data} and the axis
#' \code{lab}.
#' @keywords internal
footprint_difference <- function(plot_data, method, sample_name) {
    if (method == "substraction") {
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
        ), by = x]
        difference_data[, type := paste("Obs. sub. Exp.", sample_name)]
        lab <- "(Observed - Expected)"
    } else {
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]
        ), by = x]
        difference_data[, type := paste("Obs. div. Exp.", sample_name)]
        lab <- "(Observed / Expected)"
    }
    return(list(data = difference_data, lab = lab))
}

#' @title normalise_footprint_flank
#' @description Normalise a corrected footprint against its outer flanking
#' windows, on the same scale as the statistic itself.
#' @details Dividing a difference by its flank mean would rescale it by an
#' arbitrary factor, because that mean is near zero: for a typical footprint
#' it inflates the curve by an order of magnitude and can invert its sign.
#' @param difference_data The corrected footprint as a \code{data.table}.
#' @param method Either \code{"substraction"} or \code{"division"}.
#' @param flankNorm Width of the flanking window used for normalisation.
#' @return The normalised \code{data.table}.
#' @keywords internal
normalise_footprint_flank <- function(difference_data, method, flankNorm) {
    if (is.null(flankNorm) || flankNorm <= 0) {
        return(difference_data)
    }
    flank <- max(abs(difference_data$x), na.rm = TRUE)
    idx <- abs(difference_data$x) >= flank - flankNorm
    norm_factor <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
    if (method == "substraction") {
        difference_data[, avg_methyl := avg_methyl - norm_factor]
    } else {
        difference_data[, avg_methyl := avg_methyl / norm_factor]
    }
    return(difference_data)
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
#'     tf_bindsites = tf_bindsites,
#'     msites = msites,
#'     sample_name = "Sample1",
#'     gc_dist = gcdist,
#'     gcfreqs = gcfreqs,
#'     enhancer = NULL,
#'     method = "division"
#' )
plotMotifFootprint <- function(
    motif, tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, method = "division",
    flankNorm = 50
) {
    if (is.null(method) || !method %in% c("substraction", "division")) {
        method <- "division"
        warning("Invalid method, using the default method division.")
    }
    if (is.null(sample_name)) {
        sample_name <- "sample"
    }
    plot_data <- plotExpectedFootprint(
        motif, tf_bindsites, msites,
        sample_name = sample_name,
        gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE
    )$plotDF

    diff <- footprint_difference(plot_data, method, sample_name)
    difference_data <- normalise_footprint_flank(diff$data, method, flankNorm)

    ggplot(difference_data, aes(
        x = x,
        y = avg_methyl,
        color = type
    )) +
        geom_line() +
        xlab("Distance from motif center") +
        ylab(paste0("Methylation difference ", diff$lab)) +
        theme_classic() +
        ggtitle(paste(
            "TF footprint difference for", motif, "in", sample_name
        )) +
        theme(legend.position = "bottom") +
        xlim(-200, 200)
}
