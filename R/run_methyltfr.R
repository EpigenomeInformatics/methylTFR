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
#'     filetype = "bismarkcov"
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
    ignoreStrand = TRUE, cov_threshold = 1
) {
    if (!tolower(filetype) %in% c(
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
    check_annotation_inputs(tf_bindsites, gcfreqs, gc_dist, enhancer)
    opts <- check_run_options(
        chunkSize, threads, ignoreStrand, cov_threshold
    )

    if (is.null(annfile) || !is.character(annfile)) {
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
        annfile <- file.path(sample_dir, sample_ann)
    }

    samples <- read_sample_annotation(annfile, sampleColName)

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
        ignoreStrand = opts$ignoreStrand
    )
}
