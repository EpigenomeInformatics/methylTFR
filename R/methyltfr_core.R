#' @title check_annotation_inputs
#' @description Validate the annotation objects required by both methylTFR
#' entry points.
#' @param tf_bindsites a \code{GRangesList} object of TF binding site
#' positions.
#' @param gcfreqs a \code{list} of GC bin frequency tables.
#' @param gc_dist a \code{GRanges} object of the genome-wide GC distribution.
#' @param enhancer an optional \code{GRanges} object of regions to restrict to.
#' @return invisible TRUE, called for the side effect of signalling errors.
#' @importFrom methods is
#' @keywords internal
check_annotation_inputs <- function(
    tf_bindsites, gcfreqs, gc_dist,
    enhancer = NULL
) {
    if (any(vapply(
        list(tf_bindsites, gcfreqs, gc_dist), is.null, logical(1)
    ))) {
        stop("Please load the annotation objects for given genome.")
    }
    if (!is(tf_bindsites, "GRangesList") && !is.list(tf_bindsites)) {
        stop("tf_bindsites must be a GRangesList object")
    }
    if (!is(gcfreqs, "list")) {
        stop("gcfreqs must be a list object")
    }
    if (!is(gc_dist, "GRanges")) {
        stop("gc_dist must be a GRanges object")
    }
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("enhancer must be a GRanges object")
    }
    invisible(TRUE)
}


#' @title check_run_options
#' @description Validate the run-time options shared by both methylTFR entry
#' points, substituting documented defaults where a value is unusable.
#' @param chunkSize Chunk size for parallel processing of motifs.
#' @param threads Thread count for parallel processing.
#' @param ignoreStrand if TRUE, strand information is ignored.
#' @param cov_threshold numeric coverage threshold.
#' @return A named list with the validated values.
#' @keywords internal
check_run_options <- function(
    chunkSize = 20, threads = 1,
    ignoreStrand = TRUE, cov_threshold = 1
) {
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    if (!is.numeric(chunkSize) || chunkSize < 1) {
        warning("Invalid chunk size detected, using default chunk size")
        chunkSize <- 20
    }
    if (!is.numeric(threads) || threads < 1) {
        warning("Invalid thread count detected, using default thread count")
        threads <- 1
    }
    if (!is.numeric(cov_threshold) || cov_threshold < 0) {
        warning("Invalid cov_threshold detected, using default cov_threshold")
        cov_threshold <- 1
    }
    list(
        chunkSize = chunkSize, threads = threads,
        ignoreStrand = ignoreStrand, cov_threshold = cov_threshold
    )
}


#' @title read_sample_annotation
#' @description Read a methylTFR sample annotation table.
#' @param annfile Path to a \code{.csv} or \code{.tsv} annotation file.
#' @param sampleColName Name of the column holding the per-sample file names.
#' @return A \code{data.frame} of sample annotation.
#' @importFrom utils read.table
#' @keywords internal
read_sample_annotation <- function(annfile, sampleColName) {
    if (!file.exists(annfile)) {
        stop(sprintf(
            "%s does not exist, please check the file path !!", annfile
        ))
    }
    if (endsWith(annfile, ".csv")) {
        samples <- read.table(annfile,
            header = TRUE, sep = ",", stringsAsFactors = FALSE
        )
    } else if (endsWith(annfile, ".tsv")) {
        samples <- read.table(annfile,
            header = TRUE, sep = "\t", stringsAsFactors = FALSE
        )
    } else {
        stop(
            "Please provide a valid annotation file with a ",
            ".csv or .tsv extension"
        )
    }
    if (!sampleColName %in% colnames(samples)) {
        stop(sprintf(
            "Column '%s' was not found in the sample annotation file",
            sampleColName
        ))
    }
    samples
}


#' @title methyltfr_core
#' @description Internal engine shared by \code{\link{run_methyltfr}} and
#' \code{\link{run_methylTFR_RnBeads}}. It validates the motif set, allocates
#' the on-disk sinks, iterates over samples, computes per-motif deviations in
#' chunks and assembles the resulting \code{methylTFRdeviations} object.
#'
#' The only difference between the two public entry points is where the
#' per-sample methylation calls come from. That difference is isolated in the
#' \code{msites_fun} argument, so both entry points share identical numerical
#' behaviour.
#' @param sample_ids A character vector of sample identifiers. Used for the
#' column names of the resulting object and to size the sinks.
#' @param msites_fun A function of a single integer \code{i} returning a
#' \code{GRanges} object of methylation calls for sample \code{i}, with a
#' numeric \code{score} metadata column holding methylation levels in
#' \code{[0, 1]}.
#' @param samples A \code{data.frame} of sample annotation with one row per
#' entry of \code{sample_ids}, used as \code{colData}.
#' @param tf_bindsites a \code{GRangesList} object containing TF binding site
#' positions.
#' @param gcfreqs a \code{list} of GC bin frequency tables.
#' @param gc_dist a \code{GRanges} object containing the genome-wide GC
#' distribution.
#' @param chunkSize Chunk size for parallel processing of motifs.
#' @param threads Thread count for parallel processing.
#' @param enhancer a \code{GRanges} object restricting the analysis to a set of
#' regions such as distal regulatory elements (optional).
#' @param ignoreStrand if TRUE, strand information is ignored.
#' @return a \code{methylTFRdeviations} object with bias-corrected deviations,
#' row-wise Z-scores and expected deviations.
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom parallel mclapply
#' @importFrom logger log_info
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray DelayedArray close
#' @importFrom methods as new is
#' @keywords internal
methyltfr_core <- function(
    sample_ids,
    msites_fun,
    samples,
    tf_bindsites,
    gcfreqs,
    gc_dist,
    chunkSize = 20,
    threads = 1,
    enhancer = NULL,
    ignoreStrand = TRUE
) {
    if (!is.character(sample_ids) || length(sample_ids) == 0) {
        stop("No samples to process.")
    }
    if (!is.function(msites_fun)) {
        stop("msites_fun must be a function of a single sample index.")
    }
    if (nrow(samples) != length(sample_ids)) {
        stop("Sample annotation must have one row per sample.")
    }

    motifs <- names(gcfreqs)

    # Validate motifs: discard if TFBS is empty or matrix is missing
    valid_motifs <- vapply(motifs, function(m) {
        has_tfbs <- !is.null(tf_bindsites[[m]]) && length(tf_bindsites[[m]]) > 0
        has_matrix <- !is.null(gcfreqs[[m]])
        return(has_tfbs && has_matrix)
    }, logical(1))

    if (any(!valid_motifs)) {
        num_discarded <- sum(!valid_motifs)
        log_info(
            "Discarding ", num_discarded,
            " motifs due to empty TFBS or missing matrix."
        )
        motifs <- motifs[valid_motifs]
    }

    if (length(motifs) == 0) {
        stop("No valid motifs remaining after validation.")
    }

    # Split the motifs into chunks
    numChunks <- ceiling(length(motifs) / chunkSize)
    motif_chunks <- split(motifs, rep(seq_len(numChunks),
        each = chunkSize,
        length.out = length(motifs)
    ))

    # Create the temp sinks
    dev_sink <- create_sink(sample_ids, motifs)
    exp_sink <- create_sink(sample_ids, motifs)

    # Set the grids
    dev_grid <- set_grid(sample_ids, motif_chunks)
    exp_grid <- set_grid(sample_ids, motif_chunks)

    if (!is.null(enhancer)) {
        gc_dist <- subsetByOverlaps(gc_dist,
            enhancer,
            ignore.strand = ignoreStrand
        )
    }

    for (i in seq_along(sample_ids)) {
        sample_name <- sample_ids[i]
        msites <- msites_fun(i)
        if (!is(msites, "GRanges")) {
            stop(
                "msites_fun did not return a GRanges object for sample ",
                sample_name
            )
        }
        log_info("Processing ", sample_name)
        bin_meth <- addGCBintoMethylome(
            msites,
            gc_dist, ignoreStrand
        )

        # Process motifs in chunks
        for (j in seq_along(motif_chunks)) {
            chunk_motifs <- motif_chunks[[j]]

            sample_deviations <- mclapply(chunk_motifs,
                computeDeviation,
                msites = msites,
                tf_bindsites = tf_bindsites,
                gcfreqs = gcfreqs,
                binMsites = bin_meth,
                enhancer = enhancer,
                mc.cores = threads,
                ignoreStrand = ignoreStrand
            )
            names(sample_deviations) <- chunk_motifs

            # Write the block to the sink
            write_block_to_sink(
                lapply(sample_deviations, function(x) x$dev),
                dev_grid, i, j, dev_sink
            )
            write_block_to_sink(
                lapply(sample_deviations, function(x) x$exp_dev),
                exp_grid, i, j, exp_sink
            )
            rm(sample_deviations)
        }
        rm(msites)
        cleanMem()
        log_info("Finished processing ", sample_name)
    }
    log_success("Computed all deviations successfully")

    # Close the sinks
    DelayedArray::close(dev_sink)
    DelayedArray::close(exp_sink)
    deviation <- as.matrix(t(as(dev_sink, "DelayedArray")))
    exp_dev <- as.matrix(t(as(exp_sink, "DelayedArray")))

    se <- SummarizedExperiment(
        assays = list(
            deviations = deviation,
            z = computeRowZScore(deviation),
            expected = exp_dev
        ),
        colData = samples,
        rowData = DataFrame(motifs = row.names(deviation))
    )
    return(new("methylTFRdeviations", se))
}
