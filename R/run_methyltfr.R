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
#' @seealso \code{\link{run_methylTFR_RnBeads}} for running methylTFR
#' directly on a preprocessed RnBeads object.
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
    if (any(vapply(
        list(tf_bindsites, gcfreqs, gc_dist), is.null,
        logical(1)
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
        cov_threshold <- 1
    }
    if (!is.logical(full_path)) {
        stop("Invalid full path flag detected,
        please provide a valid logical value")
    }
    if (is.null(annfile) || !is.character(annfile)) {
        annfile <- file.path(sample_dir, sample_ann)
    }
    if (!file.exists(annfile)) {
        stop(sprintf(
            "%s does not exist, please check the file path !!",
            annfile
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
        stop("Please provide a valid annotation
         file with .csv or .tsv extension")
    }
    if (!sampleColName %in% colnames(samples)) {
        stop(sprintf(
            "Column '%s' was not found in the sample annotation file",
            sampleColName
        ))
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

    # Per-sample reader handed to the shared engine
    msites_fun <- function(i) {
        read_methylome(files_list[i],
            type = filetype,
            cov_threshold = cov_threshold
        )
    }

    methyltfr_core(
        sample_ids = basename(files_list),
        msites_fun = msites_fun,
        samples = samples,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        chunkSize = chunkSize,
        threads = threads,
        enhancer = enhancer,
        ignoreStrand = ignoreStrand
    )
}
