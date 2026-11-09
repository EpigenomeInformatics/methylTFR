#' @title check_file_run_inputs
#' @description Validate the file-specific arguments of
#' \code{run_methyltfr}.
#' @param filetype File type of the methylation call files.
#' @param sampleColName Column name holding the file names.
#' @param full_path if TRUE, the annotation file holds full paths.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
check_file_run_inputs <- function(filetype, sampleColName, full_path) {
    if (is.null(filetype) || !tolower(filetype) %in% c(
        "bissnp", "epp", "allc", "bismarkcytosine",
        "bismarkcov", "encode"
    )) {
        stop("Please provide a valid file type")
    }
    if (is.null(sampleColName) || !is.character(sampleColName)) {
        stop("Please provide a valid sample column name")
    }
    if (!is.logical(full_path)) {
        stop(
            "Invalid full path flag detected, ",
            "please provide a valid logical value"
        )
    }
    invisible(NULL)
}

#' @title resolve_annotation_file
#' @description Work out the path of the sample annotation file.
#' @param annfile Explicit path to the annotation file, or NULL.
#' @param sample_ann Name of the annotation file inside \code{sample_dir}.
#' @param sample_dir Directory holding the methylation call files.
#' @return The path of the annotation file as a character string.
#' @keywords internal
resolve_annotation_file <- function(annfile, sample_ann, sample_dir) {
    if (!is.null(annfile) && is.character(annfile)) {
        return(annfile)
    }
    if (is.null(sample_ann) || !is.character(sample_ann)) {
        stop("Please provide a valid sample annotation file")
    }
    if (is.null(sample_dir) || !is.character(sample_dir)) {
        stop("Please provide a valid sample directory")
    }
    if (!dir.exists(sample_dir)) {
        stop(
            "Sample directory does not exist, ",
            "please check the directory path"
        )
    }
    return(file.path(sample_dir, sample_ann))
}

#' @title locate_sample_files
#' @description Build and check the list of per-sample methylation files.
#' @param samples The sample annotation as a \code{data.frame}.
#' @param sample_dir Directory holding the methylation call files.
#' @param sampleColName Column name holding the file names.
#' @param full_path if TRUE, the annotation file holds full paths.
#' @return A character vector of existing file paths.
#' @importFrom logger log_success
#' @keywords internal
locate_sample_files <- function(samples, sample_dir, sampleColName, full_path) {
    if (full_path) {
        files_list <- samples[, sampleColName]
    } else {
        files_list <- file.path(sample_dir, samples[, sampleColName])
    }
    if (!all(file.exists(files_list))) {
        stop(
            "Some of the files does not exist, ",
            "please check the file path!"
        )
    }
    log_success("The samples are successfully located")
    return(files_list)
}

#' @title run_methyltfr
#' @description Compute bias-corrected deviation scores for every motif and
#' every sample listed in a sample annotation file.
#' @details Each sample is read with \code{\link{read_methylome}}, binned by
#' GC content with \code{\link{addGCBintoMethylome}} and scored per motif
#' with \code{\link{computeDeviation}}.
#' @param tf_bindsites A \code{GRangesList} named by motif. Each element is a
#' \code{GRanges} of the occurrences of that motif. The names must include
#' every name in \code{names(gcfreqs)}. Get it with
#' \code{getTFbindsites(motifSet = "jaspar2020")} from the
#' \pkg{methylTFRAnnotationHg38} (or \pkg{methylTFRAnnotationMm10}) package,
#' or build a custom set with \pkg{methylTFRAnnotationBuilder}.
#' @param gcfreqs A named \code{list} with one element per motif, named as
#' in \code{tf_bindsites}. Each element is a numeric matrix with the five GC
#' bins in rows and the positions across the binding-site window in
#' columns. Get it with \code{getGCfreq(motifSet = ...)} from the
#' annotation package.
#' @param gc_dist A \code{GRanges} object of 30-bp genomic windows with an
#' integer metadata column \code{GC_bin} (1 to 5). Get it with
#' \code{getGenomeGC()} from the annotation package.
#' @param enhancer An optional \code{GRanges} object of regions, for example
#' distal regulatory elements. Binding sites outside these regions are
#' dropped. Any region set imported as \code{GRanges} works, for example
#' the ENCODE candidate cis-regulatory elements (SCREEN) read with
#' \code{rtracklayer::import()}. Default \code{NULL} uses all binding sites.
#' @param ignoreStrand Logical. If \code{TRUE} (default), strand is ignored
#' when overlapping methylation sites with annotations.
#' @param chunkSize Number of motifs processed together in one block.
#' Default 20.
#' @param threads Number of workers used to process the motifs of a chunk in
#' parallel. Default 1 (serial).
#' @param cov_threshold Numeric. Sites with fewer reads than this are
#' dropped. Default 1.
#' @param tempDir Directory for the temporary HDF5 files written during a
#' run. Default is a folder under \code{tempdir()}. The files are deleted
#' when the run finishes.
#' @param sample_ann Name of a tab-separated (\code{.tsv}) or comma-separated
#' (\code{.csv}) sample annotation file in \code{sample_dir}. It must have
#' a header and contain the column named by \code{sampleColName}, which
#' holds the methylation file names. All columns are kept as
#' \code{colData} of the result. Not used when \code{annfile} is given.
#' @param annfile Optional path to the sample annotation file. If given,
#' it is used instead of \code{sample_ann} in \code{sample_dir}.
#' @param sample_dir Directory holding the methylation call files (and,
#' usually, the sample annotation file).
#' @param sampleColName Name of the column in \code{sample_ann} that holds
#' the methylation file names. Default \code{"bedFile"}.
#' @param full_path Logical. If \code{TRUE}, the file names in
#' \code{sampleColName} are full paths and \code{sample_dir} is not
#' prepended. Default \code{FALSE}.
#' @param filetype Format of the methylation call files, one of
#' \code{"bissnp"}, \code{"epp"}, \code{"allc"}, \code{"bismarkCytosine"},
#' \code{"bismarkCov"} or \code{"encode"} (case-insensitive). See
#' \code{\link{read_methylome}}.
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom data.table data.table
#' @importFrom logger log_info log_error
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray DelayedArray close
#' @importFrom utils read.table
#' @importFrom methods as new
#' @importFrom stats sd
#' @import logger
#' @return A \code{\linkS4class{methylTFRdeviations}} object with the assays
#' \code{deviations}, \code{z} and \code{expected}, motifs in rows and
#' samples in columns.
#' @seealso \code{\link{run_methylTFR_RnBeads}} for running methylTFR
#' directly on a preprocessed RnBeads object, \code{\link{deviations}} and
#' \code{\link{deviationZScores}} to extract the scores.
#' @examples
#' # A minimal end-to-end run on the BATF example data bundled with the
#' # package. The annotation objects cover a single motif, so the result
#' # has one row.
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#' # run_methyltfr() reads per-sample calls from disk, so the bundled sites
#' # are written out as a bismarkCov file first.
#' sample_dir <- tempfile("methylTFR_example")
#' dir.create(sample_dir)
#' n_meth <- round(msites$score * msites$coverage)
#' write.table(
#'     data.frame(
#'         chr = as.character(GenomicRanges::seqnames(msites)),
#'         start = GenomicRanges::start(msites),
#'         end = GenomicRanges::end(msites),
#'         percent = msites$score * 100,
#'         meth = n_meth,
#'         unmeth = msites$coverage - n_meth
#'     ),
#'     file.path(sample_dir, "sample_1.cov"),
#'     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
#' )
#' write.table(
#'     data.frame(sampleName = "sample_1", bedFile = "sample_1.cov"),
#'     file.path(sample_dir, "samples.tsv"),
#'     sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' devs <- run_methyltfr(
#'     sample_ann = "samples.tsv",
#'     sample_dir = sample_dir,
#'     tf_bindsites = tf_bindsites,
#'     gcfreqs = gcfreqs,
#'     gc_dist = gcdist,
#'     filetype = "bismarkCov"
#' )
#' deviations(devs)
#'
#' unlink(sample_dir, recursive = TRUE)
#' @export
run_methyltfr <- function(
    sample_ann, sample_dir, tf_bindsites = NULL,
    gcfreqs = NULL, gc_dist = NULL,
    sampleColName = "bedFile", chunkSize = 20,
    full_path = FALSE, annfile = NULL, threads = 1,
    enhancer = NULL, filetype = NULL,
    ignoreStrand = TRUE, cov_threshold = 1,
    tempDir = file.path(tempdir(), "methylTFR")
) {
    check_file_run_inputs(filetype, sampleColName, full_path)
    check_annotation_inputs(tf_bindsites, gcfreqs, gc_dist, enhancer)
    opts <- check_run_options(chunkSize, threads, ignoreStrand, cov_threshold)

    annfile <- resolve_annotation_file(annfile, sample_ann, sample_dir)
    samples <- read_sample_annotation(annfile, sampleColName)
    files_list <- locate_sample_files(
        samples, sample_dir, sampleColName, full_path
    )

    # Per-sample reader handed to the shared engine
    msites_fun <- function(i) {
        read_methylome(files_list[i],
            type = filetype,
            cov_threshold = opts$cov_threshold
        )
    }

    methyltfr_core(
        sample_ids = basename(files_list),
        msites_fun = msites_fun,
        samples = samples,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        chunkSize = opts$chunkSize,
        threads = opts$threads,
        enhancer = enhancer,
        ignoreStrand = opts$ignoreStrand,
        tempDir = tempDir
    )
}
