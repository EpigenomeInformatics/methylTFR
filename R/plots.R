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
#' @param bin_meth a \code{matrix} object containing GC
#' bin with corresponding
#' avg methylation
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
#'     "FOXF2_tf_bindsites.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "example_data.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "FOXF2_gcfreqs.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "gcdist_subset.rda",
#'     package = "methylTFR"
#' ))
#'
#' # Plot the expected footprint
#' p <- plotExpectedFootprint(
#'     motif = "FOXF2",
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
plotExpectedFootprint <- function(motif, tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, returnPlotData = FALSE,
    bin_meth = NULL) {
    if (is.null(msites)) {
        stop("msites must be a data frame,
        please provide the methylation sites")
    }
    # If sample label is not provided, use file name as label
    if (is.null(sample_name)) {
        sample_name <- "sample"
    }
    # Check if motif is valid
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }

    # Check if tf_bindsites is a GRangesList
    if (is.null(tf_bindsites) ||
        !any(c(!is(tf_bindsites, "GRangesList") ||
            !is.list(tf_bindsites)))) {
        stop("Please provide a valid tf
        binding sites as GRangesList")
    }
    # Check if enhancer is a GRanges
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("Please provide a valid enhancer regions")
    }

    # Check if gc_dist is a GRanges
    if (is.null(gc_dist) || !is(gc_dist, "GRanges")) {
        stop("Please provide a valid gc_dist as GRanges")
    }

    # Check if gcfreqs is a list
    if (is.null(gcfreqs) || !is(gcfreqs, "list")) {
        stop("Please provide a valid gcfreqs as list")
    }
    # Check if msites is a GRanges
    if (is.null(msites) || !is(msites, "GRanges")) {
        stop("Please provide a valid methylatio
        n sites with read_methylome function")
    }
    # Check if bin_meth is a matrix
    if (is.null(bin_meth) || !is.matrix(bin_meth)) {
        bin_meth <- addGCBintoMethylome(msites, gc_dist, TRUE)
    }
    # Compute footprint for the motif
    plot_data <- computeFootprint(
        motif,
        tf_bindsites, msites, enhancer
    )

    # Compute expected footprint for the motif
    exp_data <- computeExpectedFootprint(
        motif,
        gcfreqs, gc_dist, enhancer, msites, bin_meth
    )

    # Add a new column to indicate whether the data is expected or observed
    exp_data[, type := "Expected"]
    plot_data[, type := "Observed"]
    combined_data <- rbindlist(list(
        exp_data[, .(x, avg_methyl, type)],
        plot_data[, .(x, avg_methyl, type)]
    ))

    # Generate the footprint plot
    p1 <- ggplot(combined_data, aes(
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
    if (returnPlotData) {
        return(list(plot = p1, plotDF = combined_data))
    } else {
        return(p1)
    }
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
#' @param bin_meth a \code{matrix} object containing GC bin with
#'  corresponding avg methylation
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
#'     "FOXF2_tf_bindsites.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "example_data.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "FOXF2_gcfreqs.rda",
#'     package = "methylTFR"
#' ))
#' load(system.file("extdata",
#'     "gcdist_subset.rda",
#'     package = "methylTFR"
#' ))
#'
#' # Plot the expected footprint
#' p <- plotMotifFootprint(
#'     motif = "FOXF2",
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
plotMotifFootprint <- function(motif,
    tf_bindsites, msites,
    sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, method = "division",
    bin_meth = NULL) {
    if (is.null(method) ||
        !method %in% c("substraction", "division")) {
        method <- "division"
        warning("method is not provided,
        using default method substraction.")
    }
    # Get the expected methylation for the motif
    plot_data <- plotExpectedFootprint(
        motif, tf_bindsites, msites,
        sample_name = sample_name,
        gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE,
        bin_meth = bin_meth
    )$plotDF

    if (method == "substraction") {
        # Calculate observed - expected methylation
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
        ),
        by = x
        ]
        difference_data[, type := paste(
            "Observed substraction Expected",
            sample_name
        )]
        lab <- "(Observed - Expected)"
    } else if (method == "division") {
        # Calculate observed / expected methylation
        difference_data <- plot_data[, .(
            avg_methyl =
                avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]
        ),
        by = x
        ]
        difference_data[, type := paste(
            "Observed divided Expected",
            sample_name
        )]
        lab <- "(Observed / Expected)"
    }

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
