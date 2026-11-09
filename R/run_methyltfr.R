#' @title checkFileRunInputs
#' @description Validate the file-specific arguments of
#' \code{runMethylTFR}.
#' @param filetype File type of the methylation call files.
#' @param sampleColName Column name holding the file names.
#' @param fullPath if TRUE, the annotation file holds full paths.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
checkFileRunInputs <- function(filetype, sampleColName, fullPath) {
    if (is.null(filetype) || !tolower(filetype) %in% c(
        "bissnp", "epp", "allc", "bismarkcytosine",
        "bismarkcov", "encode"
    )) {
        stop("Please provide a valid file type")
    }
    if (is.null(sampleColName) || !is.character(sampleColName)) {
        stop("Please provide a valid sample column name")
    }
    if (!is.logical(fullPath)) {
        stop(
            "Invalid full path flag detected, ",
            "please provide a valid logical value"
        )
    }
    invisible(NULL)
}

#' @title resolveAnnotationFile
#' @description Work out the path of the sample annotation file.
#' @param sampleAnn File name inside \code{sampleDir}, or a path to the
#' file.
#' @param sampleDir Directory holding the methylation call files.
#' @return The path of the annotation file as a character string.
#' @keywords internal
resolveAnnotationFile <- function(sampleAnn, sampleDir) {
    if (is.null(sampleAnn) || !is.character(sampleAnn) ||
        length(sampleAnn) != 1) {
        stop("Please provide a valid sample annotation file")
    }
    if (is.null(sampleDir) || !is.character(sampleDir)) {
        stop("Please provide a valid sample directory")
    }
    if (!dir.exists(sampleDir)) {
        stop(
            "Sample directory does not exist, ",
            "please check the directory path"
        )
    }
    inDir <- file.path(sampleDir, sampleAnn)
    if (file.exists(inDir)) {
        return(inDir)
    }
    return(sampleAnn)
}

#' @title locateSampleFiles
#' @description Build and check the list of per-sample methylation files.
#' @param samples The sample annotation as a \code{data.frame}.
#' @param sampleDir Directory holding the methylation call files.
#' @param sampleColName Column name holding the file names.
#' @param fullPath if TRUE, the annotation file holds full paths.
#' @return A character vector of existing file paths.
#' @importFrom logger log_success
#' @keywords internal
locateSampleFiles <- function(samples, sampleDir, sampleColName, fullPath) {
    if (fullPath) {
        filesList <- samples[, sampleColName]
    } else {
        filesList <- file.path(sampleDir, samples[, sampleColName])
    }
    if (!all(file.exists(filesList))) {
        stop(
            "Some of the files does not exist, ",
            "please check the file path!"
        )
    }
    log_success("The samples are successfully located")
    return(filesList)
}

#' @title runMethylTFR
#' @description Compute bias-corrected deviation scores for every motif and
#' every sample listed in a sample annotation file.
#' @details Each sample is read with \code{\link{readMethylome}}, binned by
#' GC content with \code{\link{addGCBinToMethylome}} and scored per motif
#' with \code{\link{computeDeviation}}.
#' @param tfBindsites A \code{GRangesList} named by motif. Each element is a
#' \code{GRanges} of the occurrences of that motif. The names must include
#' every name in \code{names(gcfreqs)}. Get it with
#' \code{getTFbindsites(motifSet = "jaspar2020")} from the
#' \pkg{methylTFRAnnotationHg38} (or \pkg{methylTFRAnnotationMm10}) package,
#' or build a custom set with \pkg{methylTFRAnnotationBuilder}.
#' @param gcfreqs A named \code{list} with one element per motif, named as
#' in \code{tfBindsites}. Each element is a numeric matrix with the five GC
#' bins in rows and the positions across the binding-site window in
#' columns. Get it with \code{getGCfreq(motifSet = ...)} from the
#' annotation package.
#' @param gcDist A \code{GRanges} object of 30-bp genomic windows with an
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
#' @param covThreshold Numeric. Sites with fewer reads than this are
#' dropped. Default 1.
#' @param tempDir Directory for the temporary HDF5 files written during a
#' run. Default is a folder under \code{tempdir()}. The files are deleted
#' when the run finishes.
#' @param sampleAnn Name of a tab-separated (\code{.tsv}) or comma-separated
#' (\code{.csv}) sample annotation file in \code{sampleDir}, or a path to
#' it. It must have a header and contain the column named by
#' \code{sampleColName}, which holds the methylation file names. All
#' columns are kept as \code{colData} of the result.
#' @param sampleDir Directory holding the methylation call files (and,
#' usually, the sample annotation file).
#' @param sampleColName Name of the column in \code{sampleAnn} that holds
#' the methylation file names. Default \code{"bedFile"}.
#' @param fullPath Logical. If \code{TRUE}, the file names in
#' \code{sampleColName} are full paths and \code{sampleDir} is not
#' prepended. Default \code{FALSE}.
#' @param filetype Format of the methylation call files, one of
#' \code{"bissnp"}, \code{"epp"}, \code{"allc"}, \code{"bismarkCytosine"},
#' \code{"bismarkCov"} or \code{"encode"} (case-insensitive). See
#' \code{\link{readMethylome}}.
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
#' @seealso \code{\link{runMethylTFRFromRnBeads}} for running methylTFR
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
#' # runMethylTFR() reads per-sample calls from disk, so the bundled sites
#' # are written out as a bismarkCov file first.
#' sampleDir <- tempfile("methylTFR_example")
#' dir.create(sampleDir)
#' nMeth <- round(msites$score * msites$coverage)
#' write.table(
#'     data.frame(
#'         chr = as.character(GenomicRanges::seqnames(msites)),
#'         start = GenomicRanges::start(msites),
#'         end = GenomicRanges::end(msites),
#'         percent = msites$score * 100,
#'         meth = nMeth,
#'         unmeth = msites$coverage - nMeth
#'     ),
#'     file.path(sampleDir, "sample_1.cov"),
#'     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
#' )
#' write.table(
#'     data.frame(sampleName = "sample_1", bedFile = "sample_1.cov"),
#'     file.path(sampleDir, "samples.tsv"),
#'     sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' devs <- runMethylTFR(
#'     sampleAnn = "samples.tsv",
#'     sampleDir = sampleDir,
#'     tfBindsites = tf_bindsites,
#'     gcfreqs = gcfreqs,
#'     gcDist = gcdist,
#'     filetype = "bismarkCov"
#' )
#' deviations(devs)
#'
#' unlink(sampleDir, recursive = TRUE)
#' @export
runMethylTFR <- function(
    sampleAnn, sampleDir, tfBindsites = NULL,
    gcfreqs = NULL, gcDist = NULL,
    sampleColName = "bedFile", chunkSize = 20,
    fullPath = FALSE, threads = 1,
    enhancer = NULL, filetype = NULL,
    ignoreStrand = TRUE, covThreshold = 1,
    tempDir = file.path(tempdir(), "methylTFR")
) {
    checkFileRunInputs(filetype, sampleColName, fullPath)
    checkAnnotationInputs(tfBindsites, gcfreqs, gcDist, enhancer)
    opts <- checkRunOptions(chunkSize, threads, ignoreStrand, covThreshold)

    annFile <- resolveAnnotationFile(sampleAnn, sampleDir)
    samples <- readSampleAnnotation(annFile, sampleColName)
    filesList <- locateSampleFiles(
        samples, sampleDir, sampleColName, fullPath
    )

    # Per-sample reader handed to the shared engine
    msitesFun <- function(i) {
        readMethylome(filesList[i],
            type = filetype,
            covThreshold = opts$covThreshold
        )
    }

    methylTFRCore(
        sampleIds = basename(filesList),
        msitesFun = msitesFun,
        samples = samples,
        tfBindsites = tfBindsites,
        gcfreqs = gcfreqs,
        gcDist = gcDist,
        chunkSize = opts$chunkSize,
        threads = opts$threads,
        enhancer = enhancer,
        ignoreStrand = opts$ignoreStrand,
        tempDir = tempDir
    )
}
