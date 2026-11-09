#' @title check_msites
#' @description Check that methylation sites have the structure returned by
#' \code{readMethylome}.
#' @param msites The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @keywords internal
check_msites <- function(msites) {
    if (is.null(msites) || !is(msites, "GRanges")) {
        stop(
            "msites must be a GRanges object of methylation sites; ",
            "create it with readMethylome()"
        )
    }
    missing_cols <- setdiff(
        c("score", "coverage"), colnames(S4Vectors::mcols(msites))
    )
    if (length(missing_cols) > 0) {
        stop(
            "msites is missing the metadata column(s) ",
            toString(missing_cols), "; create it with readMethylome()"
        )
    }
    invisible(NULL)
}

#' @title check_tf_bindsites
#' @description Check that TF binding sites are a named list of ranges.
#' @param tf_bindsites The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
check_tf_bindsites <- function(tf_bindsites) {
    if (is.null(tf_bindsites) ||
        !(is(tf_bindsites, "GRangesList") || is.list(tf_bindsites))) {
        stop(
            "Please provide a valid tf binding sites as GRangesList; ",
            "get it with getTFbindsites() from the annotation package"
        )
    }
    if (is.null(names(tf_bindsites))) {
        stop("tf_bindsites must be named by motif")
    }
    invisible(NULL)
}

#' @title check_gcfreqs
#' @description Check that GC frequency tables are a named list.
#' @param gcfreqs The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
check_gcfreqs <- function(gcfreqs) {
    if (is.null(gcfreqs) || !is(gcfreqs, "list")) {
        stop(
            "Please provide a valid gcfreqs as list; ",
            "get it with getGCfreq() from the annotation package"
        )
    }
    if (is.null(names(gcfreqs))) {
        stop("gcfreqs must be named by motif")
    }
    invisible(NULL)
}

#' @title check_gc_dist
#' @description Check that the genome-wide GC distribution is a
#' \code{GRanges} object with a \code{GC_bin} column.
#' @param gc_dist The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @keywords internal
check_gc_dist <- function(gc_dist) {
    if (is.null(gc_dist) || !is(gc_dist, "GRanges")) {
        stop(
            "Please provide a valid gc_dist as GRanges; ",
            "get it with getGenomeGC() from the annotation package"
        )
    }
    if (!"GC_bin" %in% colnames(S4Vectors::mcols(gc_dist))) {
        stop(
            "gc_dist must have a GC_bin metadata column; ",
            "get it with getGenomeGC() from the annotation package"
        )
    }
    invisible(NULL)
}

#' @title check_enhancer
#' @description Check the optional region set.
#' @param enhancer The object to check, or \code{NULL}.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
check_enhancer <- function(enhancer) {
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("Please provide a valid enhancer regions as GRanges")
    }
    invisible(NULL)
}

#' @title check_annotation_inputs
#' @description Validate the annotation objects required by both methylTFR
#' entry points.
#' @inheritParams runMethylTFR
#' @return invisible TRUE, called for the side effect of signalling errors.
#' @keywords internal
check_annotation_inputs <- function(
    tf_bindsites, gcfreqs, gc_dist, enhancer = NULL
) {
    if (any(vapply(
        list(tf_bindsites, gcfreqs, gc_dist), is.null, logical(1)
    ))) {
        stop("Please load the annotation objects for given genome.")
    }
    check_tf_bindsites(tf_bindsites)
    check_gcfreqs(gcfreqs)
    check_gc_dist(gc_dist)
    check_enhancer(enhancer)
    no_sites <- setdiff(names(gcfreqs), names(tf_bindsites))
    if (length(no_sites) > 0) {
        stop(
            "names(gcfreqs) and names(tf_bindsites) do not match. ",
            "No binding sites for: ", toString(utils::head(no_sites, 5)),
            if (length(no_sites) > 5) ", ..." else "",
            ". Use getTFbindsites() and getGCfreq() with matching motif sets."
        )
    }
    invisible(TRUE)
}


#' @title check_run_options
#' @description Validate the run-time options shared by both methylTFR entry
#' points, substituting documented defaults where a value is unusable.
#' @inheritParams runMethylTFR
#' @return A named list with the validated values.
#' @keywords internal
check_run_options <- function(
    chunkSize = 20, threads = 1, ignoreStrand = TRUE, cov_threshold = 1
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


#' @title check_core_inputs
#' @description Validate the arguments shared by both methylTFR entry points.
#' @param sample_ids A character vector of sample identifiers.
#' @param msites_fun A function of a single integer sample index.
#' @param samples A \code{data.frame} with one row per sample.
#' @param tempDir Directory for the temporary HDF5 files.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
check_core_inputs <- function(sample_ids, msites_fun, samples, tempDir) {
    if (!is.character(sample_ids) || length(sample_ids) == 0) {
        stop("No samples to process.")
    }
    if (!is.function(msites_fun)) {
        stop("msites_fun must be a function of a single sample index.")
    }
    if (nrow(samples) != length(sample_ids)) {
        stop("Sample annotation must have one row per sample.")
    }
    if (!is.character(tempDir) || length(tempDir) != 1 || !nzchar(tempDir)) {
        stop("tempDir must be a single directory path")
    }
    invisible(NULL)
}

#' @title valid_core_motifs
#' @description Drop motifs whose binding sites are empty or whose GC bin
#' frequency matrix is missing.
#' @inheritParams runMethylTFR
#' @return A character vector of the motif names that can be processed.
#' @importFrom logger log_info
#' @keywords internal
valid_core_motifs <- function(tf_bindsites, gcfreqs) {
    motifs <- names(gcfreqs)
    valid_motifs <- vapply(motifs, function(m) {
        has_tfbs <- !is.null(tf_bindsites[[m]]) &&
            length(tf_bindsites[[m]]) > 0
        has_matrix <- !is.null(gcfreqs[[m]])
        return(has_tfbs && has_matrix)
    }, logical(1))

    if (any(!valid_motifs)) {
        log_info(
            "Discarding ", sum(!valid_motifs),
            " motifs due to empty TFBS or missing matrix."
        )
        motifs <- motifs[valid_motifs]
    }

    if (length(motifs) == 0) {
        stop("No valid motifs remaining after validation.")
    }
    return(motifs)
}

#' @title bpparam_from_threads
#' @description Build the \pkg{BiocParallel} back-end used to spread the
#' motifs of one chunk over workers.
#' @details A forking back-end is used where the platform supports it and a
#' socket back-end on Windows, so that \code{threads} has the same meaning
#' on every platform. \code{threads = 1} runs serially in the current
#' process.
#' @param threads Thread count for parallel processing.
#' @return A \code{BiocParallelParam} object.
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#' @keywords internal
bpparam_from_threads <- function(threads) {
    if (is.null(threads) || !is.numeric(threads) || threads <= 1) {
        return(SerialParam())
    }
    if (.Platform$OS.type == "windows") {
        return(SnowParam(workers = threads))
    }
    return(MulticoreParam(workers = threads))
}

#' @title process_core_sample
#' @description Compute and write the deviations of one sample, one motif
#' chunk at a time.
#' @inheritParams runMethylTFR
#' @param index Integer index of the sample within \code{sample_ids}.
#' @param sample_ids A character vector of sample identifiers.
#' @param msites_fun A function of a single integer sample index.
#' @param motif_chunks A \code{list} of character vectors of motif names.
#' @param dev_grid,exp_grid The grids the blocks are written on.
#' @param dev_sink,exp_sink The sinks the blocks are written to.
#' @param BPPARAM A \code{BiocParallelParam} object.
#' @return Invisible \code{NULL}. Called for its effect on the sinks.
#' @importFrom BiocParallel bplapply
#' @importFrom logger log_info
#' @keywords internal
process_core_sample <- function(
    index, sample_ids, msites_fun, motif_chunks, tf_bindsites, gcfreqs,
    gc_dist, dev_grid, exp_grid, dev_sink, exp_sink, BPPARAM, enhancer,
    ignoreStrand
) {
    sample_name <- sample_ids[index]
    msites <- msites_fun(index)
    check_msites(msites)
    log_info("Processing ", sample_name)
    bin_meth <- addGCBinToMethylome(msites, gc_dist, ignoreStrand)

    # Process motifs in chunks
    for (j in seq_along(motif_chunks)) {
        chunk_motifs <- motif_chunks[[j]]

        sample_deviations <- bplapply(chunk_motifs,
            computeDeviation,
            msites = msites,
            tf_bindsites = tf_bindsites,
            gcfreqs = gcfreqs,
            binMsites = bin_meth,
            enhancer = enhancer,
            ignoreStrand = ignoreStrand,
            BPPARAM = BPPARAM
        )
        names(sample_deviations) <- chunk_motifs

        write_block_to_sink(
            lapply(sample_deviations, function(x) x$dev),
            dev_grid, index, j, dev_sink
        )
        write_block_to_sink(
            lapply(sample_deviations, function(x) x$exp_dev),
            exp_grid, index, j, exp_sink
        )
    }
    log_info("Finished processing ", sample_name)
    invisible(NULL)
}

#' @title assemble_core_result
#' @description Close the sinks and assemble the deviations, their row-wise
#' Z-scores and the expected deviations into a result object.
#' @param dev_sink The sink holding the bias-corrected deviations.
#' @param exp_sink The sink holding the expected deviations.
#' @param samples A \code{data.frame} with one row per sample.
#' @return a \code{methylTFRdeviations} object.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray DelayedArray close
#' @importFrom methods as new
#' @keywords internal
assemble_core_result <- function(dev_sink, exp_sink, samples) {
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

#' @title methyltfr_core
#' @description Internal engine shared by \code{\link{runMethylTFR}} and
#' \code{\link{runMethylTFRFromRnBeads}}. It validates the motif set,
#' allocates the on-disk sinks, iterates over samples, computes per-motif
#' deviations in chunks and assembles the resulting
#' \code{methylTFRdeviations} object.
#'
#' The only difference between the two public entry points is where the
#' per-sample methylation calls come from. That difference is isolated in
#' the \code{msites_fun} argument, so both entry points share identical
#' numerical behaviour.
#' @inheritParams runMethylTFR
#' @param sample_ids A character vector of sample identifiers. Used for the
#' column names of the resulting object and to size the sinks.
#' @param msites_fun A function of a single integer \code{i} returning the
#' methylation sites of sample \code{i} as a \code{GRanges} with the
#' metadata columns \code{score} and \code{coverage}.
#' @param samples A \code{data.frame} of sample annotation with one row per
#' entry of \code{sample_ids}, used as \code{colData}.
#' @return a \code{methylTFRdeviations} object with bias-corrected
#' deviations, row-wise Z-scores and expected deviations.
#' @importFrom IRanges subsetByOverlaps
#' @importFrom logger log_success
#' @keywords internal
methyltfr_core <- function(
    sample_ids, msites_fun, samples, tf_bindsites, gcfreqs, gc_dist,
    chunkSize = 20, threads = 1, enhancer = NULL, ignoreStrand = TRUE,
    tempDir = file.path(tempdir(), "methylTFR")
) {
    check_core_inputs(sample_ids, msites_fun, samples, tempDir)
    motifs <- valid_core_motifs(tf_bindsites, gcfreqs)
    BPPARAM <- bpparam_from_threads(threads)

    # Split the motifs into chunks
    numChunks <- ceiling(length(motifs) / chunkSize)
    motif_chunks <- split(motifs, rep(seq_len(numChunks),
        each = chunkSize,
        length.out = length(motifs)
    ))

    # Temporary HDF5 sinks, removed when the run ends
    dev_sink <- create_sink(sample_ids, motifs, tempDir = tempDir)
    exp_sink <- create_sink(sample_ids, motifs, tempDir = tempDir)
    on.exit(unlink(c(dev_sink@filepath, exp_sink@filepath)), add = TRUE)

    dev_grid <- set_grid(sample_ids, motif_chunks)
    exp_grid <- set_grid(sample_ids, motif_chunks)

    if (!is.null(enhancer)) {
        gc_dist <- subsetByOverlaps(gc_dist,
            enhancer,
            ignore.strand = ignoreStrand
        )
    }

    for (i in seq_along(sample_ids)) {
        process_core_sample(
            index = i, sample_ids = sample_ids, msites_fun = msites_fun,
            motif_chunks = motif_chunks, tf_bindsites = tf_bindsites,
            gcfreqs = gcfreqs, gc_dist = gc_dist, dev_grid = dev_grid,
            exp_grid = exp_grid, dev_sink = dev_sink,
            exp_sink = exp_sink, BPPARAM = BPPARAM,
            enhancer = enhancer, ignoreStrand = ignoreStrand
        )
    }
    log_success("Computed all deviations successfully")

    return(assemble_core_result(dev_sink, exp_sink, samples))
}
