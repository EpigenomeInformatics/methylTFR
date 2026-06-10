#' @title run_methyltfr
#' @description This function is a wrapper function to
#' calculate the deviation
#' in transcription factor
#' footprint base for all given motifs per raw samples
#' @param sample_ann A tab seperated file contains
#' sample annotations
#' @param sample_dir The directory where all bed file
#'  and annotation file stored
#' @param threads Thread count for parallel processing
#' @param enhancer a \code{GRanges} object specifying
#' regions such as distal motif (optional)
#' @param tf_bindsites a \code{GRangesList} object contains
#'  tf binding sites positions
#' @param gcfreqs a \code{list} of GC bin frequency tables
#'  (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains
#' Genome wide GC distribution
#' @param sampleColName column name of the sample bed
#'  file in the annotation file
#' @param chunkSize Chunk size for parallel processing
#'  of motifs (default: 20)
#' @param full_path if TRUE, the bed file path in the
#' annotation file is full path
#' @param annfile if provided, the sample annotation file
#'  is not read from the sample_dir
#' @param filetype file type of the bed file,
#' currently supported:
#'  bissnp,epp,allc,bismarkcytosine,bismarkcov,encode
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param cov_threshold - numeric, coverage threshold to
#'  filter out low coverage sites,
#' default is 1
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom data.table data.table
#' @importFrom parallel mclapply
#' @importFrom logger log_info log_error
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray DelayedArray close
#' @importFrom utils read.table
#' @importFrom methods as new
#' @importFrom stats sd
#' @import logger
#' @return a \code{methylTFRdeviations} object with
#' bias-corrected deviation and Z-scores
#' @export
run_methyltfr <- function(
        sample_ann, sample_dir, tf_bindsites = NULL,
        gcfreqs = NULL, gc_dist = NULL,
        sampleColName = "bedFile", chunkSize = 20,
        full_path = FALSE, annfile = NULL, threads = 1,
        enhancer = NULL, filetype = NULL,
        ignoreStrand = TRUE, cov_threshold = 1) {
    if (!tolower(filetype) %in% c(
        "bissnp", "epp", "allc", "bismarkcytosine",
        "bismarkcov", "encode"
    )) {
        stop("Please provide a valid file type")
    }
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    if (is.null(sampleColName) || !is.character(sampleColName)) {
        stop("Please provide a valid sample column name")
    }
    if (!is.numeric(chunkSize) || chunkSize < 1) {
        warning("Invalid chunk size detected, using default chunk size")
        chunkSize <- 20
    }
    if (any(sapply(list(tf_bindsites, gcfreqs, gc_dist), is.null))) {
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
    if (is.null(annfile) || !is.character(annfile)) {
        if (is.null(sample_ann) || !is.character(sample_ann)) {
            stop("Please provide a valid sample annotation file")
        }
        if (is.null(sample_dir) || !is.character(sample_dir)) {
            stop("Please provide a valid sample directory")
        }
        if (!dir.exists(sample_dir)) {
            stop("Sample directory does not exist,
            please check the directory path")
        }
    }
    if (!is.numeric(threads) || threads < 1) {
        warning("Invalid thread count detected,
        using default thread count")
        threads <- 1
    }
    if (!is.numeric(cov_threshold)) {
        warning("Invalid cov_threshold detected,
         using default cov_threshold")
        cov_threshold <- 5
    }
    if (!is.logical(full_path)) {
        stop("Invalid full path flag detected,
        please provide a valid logical value")
    }
    if (is.null(annfile) || !is.character(annfile)) {
        annfile <- file.path(sample_dir, sample_ann)
    }
    if (!file.exists(annfile)) {
        stop(" %s does not exist, please check
        the file path !!", annfile)
    }
    if (endsWith(annfile, ".csv")) {
        samples <- read.table(annfile,
            header = TRUE, sep = ",", stringsAsFactors = FALSE
        )
    }
    if (endsWith(annfile, ".tsv")) {
        samples <- read.table(annfile,
            header = TRUE, sep = "\t", stringsAsFactors = FALSE
        )
    } else {
        stop("Please provide a valid annotation
         file with .csv or .tsv extension")
    }
    if (full_path) {
        files_list <- samples[, sampleColName]
    } else {
        files_list <- file.path(
            sample_dir,
            samples[, sampleColName]
        )
    }
    if (!all(file.exists(files_list))) {
        stop("Some of the files does not exist,
        please check the file path!")
    }
    log_success("The samples are successfully located")

    motifs <- names(gcfreqs)

    # Validate motifs: discard if TFBS is empty or matrix is missing
    valid_motifs <- vapply(motifs, function(m) {
        has_tfbs <- !is.null(tf_bindsites[[m]]) && length(tf_bindsites[[m]]) > 0
        has_matrix <- !is.null(gcfreqs[[m]])
        return(has_tfbs && has_matrix)
    }, logical(1))

    if (any(!valid_motifs)) {
        num_discarded <- sum(!valid_motifs)
        log_info(paste("Discarding", num_discarded, 
                       "motifs due to empty TFBS or missing matrix."))
        motifs <- motifs[valid_motifs]
    }

    if (length(motifs) == 0) {
        stop("No valid motifs remaining after validation.")
    }

    # Split the motifs into chunks
    numChunks <- ceiling(length(motifs) / chunkSize)
    motif_chunks <- split(motifs, rep(1:numChunks,
        each = chunkSize,
        length.out = length(motifs)
    ))

    # Create a temp sinks
    dev_sink <- create_sink(files_list, motifs)
    z_sink <- create_sink(files_list, motifs)

    # set the grid
    dev_grid <- set_grid(files_list, motif_chunks)
    z_grid <- set_grid(files_list, motif_chunks)

    if (!is.null(enhancer)) {
        gc_dist <- subsetByOverlaps(gc_dist,
            enhancer,
            ignore.strand = ignoreStrand
        )
    }

    for (i in seq_along(files_list)) {
        bedfile <- files_list[i]
        sample_name <- basename(bedfile)
        msites <- read_methylome(bedfile,
            type = filetype, cov_threshold
        )
        log_info("Processing ", sample_name)
        bin_meth <- addGCBintoMethylome(
            msites,
            gc_dist, ignoreStrand
        )

        # Process motifs in chunks
        for (j in seq_along(motif_chunks)) {
            # Get the current chunk
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
                z_grid, i, j, z_sink
            )
            rm(sample_deviations)
        }
        rm(msites)
        cleanMem()
        log_info("Finished processing ", sample_name)
    }
    log_success("Computed all deviations successfully")

    # Close the sink
    DelayedArray::close(dev_sink)
    DelayedArray::close(z_sink)
    deviation <- as.matrix(t(as(dev_sink, "DelayedArray")))
    exp_dev <- as.matrix(t(as(z_sink, "DelayedArray")))

    # Compute the sd and normalize the deviation
    se <- SummarizedExperiment(
        assays = list(
            deviations = deviation,
            z = computeRowZScore(deviation) # ,
            # exp_dev = exp_dev
        ),
        colData = samples,
        rowData = DataFrame(motifs = row.names(deviation))
    )
    return(new("methylTFRdeviations", se))
}