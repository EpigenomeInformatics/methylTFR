#' @title check_footprint_inputs
#' @description Validate the inputs shared by the two footprint plotting
#' functions.
#' @param motif Motif name as a character string.
#' @param tf_bindsites a \code{GRangesList} of TF binding site positions.
#' @param msites Methylation sites as a \code{GRanges} object.
#' @param gc_dist a \code{GRanges} of the genome-wide GC distribution.
#' @param gcfreqs a \code{list} of GC bin frequency tables.
#' @param enhancer a \code{GRanges} restricting the analysis (optional).
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
check_footprint_inputs <- function(
    motif, tf_bindsites, msites, gc_dist, gcfreqs, enhancer = NULL
) {
    if (is.null(msites)) {
        stop(
            "msites must be a data frame, ",
            "please provide the methylation sites"
        )
    }
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }
    if (is.null(tf_bindsites) ||
        !any(c(!is(tf_bindsites, "GRangesList") ||
            !is.list(tf_bindsites)))) {
        stop("Please provide a valid tf binding sites as GRangesList")
    }
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("Please provide a valid enhancer regions")
    }
    if (is.null(gc_dist) || !is(gc_dist, "GRanges")) {
        stop("Please provide a valid gc_dist as GRanges")
    }
    if (is.null(gcfreqs) || !is(gcfreqs, "list")) {
        stop("Please provide a valid gcfreqs as list")
    }
    if (!is(msites, "GRanges")) {
        stop(
            "Please provide a valid methylation ",
            "sites with read_methylome function"
        )
    }
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
        ggtitle(paste(
            "TF footprint for", motif, "in",
            sample_name
        )) +
        scale_color_manual(values = c(
            "Expected" = "blue",
            "Observed" = "red"
        )) +
        theme(legend.position = "bottom")
}

#' @title plotExpectedFootprint
#' @description Creates a footprint plot of expected
#' vs observed methylation
#' for a given motif and a single methylation site.
#' @param motif Motif name as a character string
#' @param tf_bindsites Transcript Factor binding sites
#' from the annotation package
#' @param msites Methylation sites as a data frame
#' @param sample_name Optional sample name as a character
#'  string, defaults
#' to the file name if not provided
#' @param gc_dist a \code{GRanges} object contains Genome
#'  wide GC distribution
#' @param gcfreqs GC frequency of the genome used to
#' compute expected methylation
#' @param enhancer Specific region such as distal motif,
#'  proximal motif
#' @param returnPlotData Logical indicating whether to
#'  return the plot data
#' @return A ggplot object of TF footprint plot for a
#'  single sample
#' @export
#' @importFrom methods is
#' @importFrom ggplot2 aes xlab ylab scale_color_manual theme xlim
#' @examples
#' library(ggplot2)
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
#' # Plot the expected footprint
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
#'
#' @importFrom ggplot2 ggplot geom_point geom_line ggtitle theme_classic
plotExpectedFootprint <- function(
    motif, tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, returnPlotData = FALSE
) {
    # If sample label is not provided, use file name as label
    if (is.null(sample_name)) {
        sample_name <- "sample"
    }
    check_footprint_inputs(
        motif, tf_bindsites, msites, gc_dist, gcfreqs, enhancer
    )

    # Compute footprint for the motif
    plot_data <- computeFootprint(
        motif,
        tf_bindsites, msites, enhancer
    )

    # Compute expected footprint for the motif
    exp_data <- computeExpectedFootprint(
        motif,
        gcfreqs, gc_dist, enhancer, msites
    )

    # Add a new column to indicate whether the data is expected or observed
    exp_data[, type := "Expected"]
    plot_data[, type := "Observed"]
    combined_data <- rbindlist(list(
        exp_data[, .(x, avg_methyl, type)],
        plot_data[, .(x, avg_methyl, type)]
    ))

    p1 <- expected_footprint_plot(combined_data, motif, sample_name)
    if (returnPlotData) {
        return(list(plot = p1, plotDF = combined_data))
    } else {
        return(p1)
    }
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
        # Calculate observed - expected methylation
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
        ),
        by = x
        ]
        difference_data[, type := paste(
            "Obs. sub. Exp.",
            sample_name
        )]
        lab <- "(Observed - Expected)"
    } else {
        # Calculate observed / expected methylation
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]
        ),
        by = x
        ]
        difference_data[, type := paste(
            "Obs. div. Exp.",
            sample_name
        )]
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
#' @description Creates a footprint plot of bias
#' corrected methylation for a given motif and sample
#' @param motif Motif name as a character string
#' @param tf_bindsites Transcript Factor binding
#' sites from the annotation package
#' @param msites Methylation sites as a data frame
#' @param sample_name Optional sample name as a
#' character string, defaults to the file name if not provided
#' @param gc_dist a \code{GRanges} object contains
#'  Genome wide GC distribution
#' @param gcfreqs GC frequency of the genome used
#' to compute expected methylation
#' @param enhancer Specific region such as distal motif,
#'  proximal motif
#' @param method Method to calculate the difference,
#' either "substraction" or "division"
#' @param flankNorm Numeric value indicating the flank normalization
#' distance, default is 50
#' @return A ggplot object of TF footprint difference plot
#'  for a single sample
#' @importFrom methods is
#' @importFrom ggplot2 aes xlab ylab scale_color_manual theme xlim
#' @export
#' @examples
#' library(ggplot2)
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
#' # Plot the expected footprint
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
#'
#' @importFrom ggplot2 ggplot geom_point geom_line ggtitle theme_classic
plotMotifFootprint <- function(
    motif,
    tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, method = "division",
    flankNorm = 50
) {
    if (is.null(method) ||
        !method %in% c("substraction", "division")) {
        method <- "division"
        warning("method is not provided,using default method substraction.")
    }
    # Get the expected methylation for the motif
    plot_data <- plotExpectedFootprint(
        motif, tf_bindsites, msites,
        sample_name = sample_name,
        gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE
    )$plotDF

    diff <- footprint_difference(plot_data, method, sample_name)
    difference_data <- diff$data
    lab <- diff$lab
    difference_data <- normalise_footprint_flank(
        difference_data, method, flankNorm
    )

    p_combined <- ggplot(difference_data, aes(
        x = x,
        y = avg_methyl,
        color = type
    )) +
        geom_line() + # Add lines for each type
        xlab("Distance from motif center") +
        ylab(paste0("Methylation difference ", lab)) +
        theme_classic() +
        ggtitle(paste(
            "TF footprint difference for",
            motif, "in ",
            sample_name
        )) +
        theme(legend.position = "bottom") +
        xlim(-200, 200) # Set the x-axis limits

    return(p_combined)
}
