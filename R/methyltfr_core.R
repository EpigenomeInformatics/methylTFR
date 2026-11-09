#' @title checkMsites
#' @description Check that methylation sites have the structure returned by
#' \code{readMethylome}.
#' @param msites The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @keywords internal
checkMsites <- function(msites) {
    if (is.null(msites) || !is(msites, "GRanges")) {
        stop(
            "msites must be a GRanges object of methylation sites; ",
            "create it with readMethylome()"
        )
    }
    missingCols <- setdiff(
        c("score", "coverage"), colnames(S4Vectors::mcols(msites))
    )
    if (length(missingCols) > 0) {
        stop(
            "msites is missing the metadata column(s) ",
            toString(missingCols), "; create it with readMethylome()"
        )
    }
    invisible(NULL)
}

#' @title checkTfBindsites
#' @description Check that TF binding sites are a named list of ranges.
#' @param tfBindsites The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
checkTfBindsites <- function(tfBindsites) {
    if (is.null(tfBindsites) ||
        !(is(tfBindsites, "GRangesList") || is.list(tfBindsites))) {
        stop(
            "Please provide a valid tf binding sites as GRangesList; ",
            "get it with getTFbindsites() from the annotation package"
        )
    }
    if (is.null(names(tfBindsites))) {
        stop("tfBindsites must be named by motif")
    }
    invisible(NULL)
}

#' @title checkGcfreqs
#' @description Check that GC frequency tables are a named list.
#' @param gcfreqs The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
checkGcfreqs <- function(gcfreqs) {
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

#' @title checkGcDist
#' @description Check that the genome-wide GC distribution is a
#' \code{GRanges} object with a \code{GC_bin} column.
#' @param gcDist The object to check.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @keywords internal
checkGcDist <- function(gcDist) {
    if (is.null(gcDist) || !is(gcDist, "GRanges")) {
        stop(
            "Please provide a valid gcDist as GRanges; ",
            "get it with getGenomeGC() from the annotation package"
        )
    }
    if (!"GC_bin" %in% colnames(S4Vectors::mcols(gcDist))) {
        stop(
            "gcDist must have a GC_bin metadata column; ",
            "get it with getGenomeGC() from the annotation package"
        )
    }
    invisible(NULL)
}

#' @title checkEnhancer
#' @description Check the optional region set.
#' @param enhancer The object to check, or \code{NULL}.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
checkEnhancer <- function(enhancer) {
    if (!is.null(enhancer) && !is(enhancer, "GRanges")) {
        stop("Please provide a valid enhancer regions as GRanges")
    }
    invisible(NULL)
}

#' @title checkAnnotationInputs
#' @description Validate the annotation objects required by both methylTFR
#' entry points.
#' @inheritParams runMethylTFR
#' @return invisible TRUE, called for the side effect of signalling errors.
#' @keywords internal
checkAnnotationInputs <- function(
    tfBindsites, gcfreqs, gcDist, enhancer = NULL
) {
    if (any(vapply(
        list(tfBindsites, gcfreqs, gcDist), is.null, logical(1)
    ))) {
        stop("Please load the annotation objects for given genome.")
    }
    checkTfBindsites(tfBindsites)
    checkGcfreqs(gcfreqs)
    checkGcDist(gcDist)
    checkEnhancer(enhancer)
    noSites <- setdiff(names(gcfreqs), names(tfBindsites))
    if (length(noSites) > 0) {
        stop(
            "names(gcfreqs) and names(tfBindsites) do not match. ",
            "No binding sites for: ", toString(utils::head(noSites, 5)),
            if (length(noSites) > 5) ", ..." else "",
            ". Use getTFbindsites() and getGCfreq() with matching motif sets."
        )
    }
    invisible(TRUE)
}


#' @title checkRunOptions
#' @description Validate the run-time options shared by both methylTFR entry
#' points, substituting documented defaults where a value is unusable.
#' @inheritParams runMethylTFR
#' @return A named list with the validated values.
#' @keywords internal
checkRunOptions <- function(
    chunkSize = 20, threads = 1, ignoreStrand = TRUE, covThreshold = 1
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
    if (!is.numeric(covThreshold) || covThreshold < 0) {
        warning("Invalid covThreshold detected, using default covThreshold")
        covThreshold <- 1
    }
    list(
        chunkSize = chunkSize, threads = threads,
        ignoreStrand = ignoreStrand, covThreshold = covThreshold
    )
}


#' @title readSampleAnnotation
#' @description Read a methylTFR sample annotation table.
#' @param annFile Path to a \code{.csv} or \code{.tsv} annotation file.
#' @param sampleColName Name of the column holding the per-sample file names.
#' @return A \code{data.frame} of sample annotation.
#' @importFrom utils read.table
#' @keywords internal
readSampleAnnotation <- function(annFile, sampleColName) {
    if (!file.exists(annFile)) {
        stop(sprintf(
            "%s does not exist, please check the file path !!", annFile
        ))
    }
    if (endsWith(annFile, ".csv")) {
        samples <- read.table(annFile,
            header = TRUE, sep = ",", stringsAsFactors = FALSE
        )
    } else if (endsWith(annFile, ".tsv")) {
        samples <- read.table(annFile,
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


#' @title checkCoreInputs
#' @description Validate the arguments shared by both methylTFR entry points.
#' @param sampleIds A character vector of sample identifiers.
#' @param msitesFun A function of a single integer sample index.
#' @param samples A \code{data.frame} with one row per sample.
#' @param tempDir Directory for the temporary HDF5 files.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @keywords internal
checkCoreInputs <- function(sampleIds, msitesFun, samples, tempDir) {
    if (!is.character(sampleIds) || length(sampleIds) == 0) {
        stop("No samples to process.")
    }
    if (!is.function(msitesFun)) {
        stop("msitesFun must be a function of a single sample index.")
    }
    if (nrow(samples) != length(sampleIds)) {
        stop("Sample annotation must have one row per sample.")
    }
    if (!is.character(tempDir) || length(tempDir) != 1 || !nzchar(tempDir)) {
        stop("tempDir must be a single directory path")
    }
    invisible(NULL)
}

#' @title validCoreMotifs
#' @description Drop motifs whose binding sites are empty or whose GC bin
#' frequency matrix is missing.
#' @inheritParams runMethylTFR
#' @return A character vector of the motif names that can be processed.
#' @importFrom logger log_info
#' @keywords internal
validCoreMotifs <- function(tfBindsites, gcfreqs) {
    motifs <- names(gcfreqs)
    validMotifs <- vapply(motifs, function(m) {
        hasTfbs <- !is.null(tfBindsites[[m]]) &&
            length(tfBindsites[[m]]) > 0
        hasMatrix <- !is.null(gcfreqs[[m]])
        return(hasTfbs && hasMatrix)
    }, logical(1))

    if (any(!validMotifs)) {
        log_info(
            "Discarding ", sum(!validMotifs),
            " motifs due to empty TFBS or missing matrix."
        )
        motifs <- motifs[validMotifs]
    }

    if (length(motifs) == 0) {
        stop("No valid motifs remaining after validation.")
    }
    return(motifs)
}

#' @title bpparamFromThreads
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
bpparamFromThreads <- function(threads) {
    if (is.null(threads) || !is.numeric(threads) || threads <= 1) {
        return(SerialParam())
    }
    if (.Platform$OS.type == "windows") {
        return(SnowParam(workers = threads))
    }
    return(MulticoreParam(workers = threads))
}

#' @title processCoreSample
#' @description Compute and write the deviations of one sample, one motif
#' chunk at a time.
#' @inheritParams runMethylTFR
#' @param index Integer index of the sample within \code{sampleIds}.
#' @param sampleIds A character vector of sample identifiers.
#' @param msitesFun A function of a single integer sample index.
#' @param motifChunks A \code{list} of character vectors of motif names.
#' @param devGrid,expGrid The grids the blocks are written on.
#' @param devSink,expSink The sinks the blocks are written to.
#' @param BPPARAM A \code{BiocParallelParam} object.
#' @return Invisible \code{NULL}. Called for its effect on the sinks.
#' @importFrom BiocParallel bplapply
#' @importFrom logger log_info
#' @keywords internal
processCoreSample <- function(
    index, sampleIds, msitesFun, motifChunks, tfBindsites, gcfreqs,
    gcDist, devGrid, expGrid, devSink, expSink, BPPARAM, enhancer,
    ignoreStrand
) {
    sampleName <- sampleIds[index]
    msites <- msitesFun(index)
    checkMsites(msites)
    log_info("Processing ", sampleName)
    binMeth <- addGCBinToMethylome(msites, gcDist, ignoreStrand)

    # Process motifs in chunks
    for (j in seq_along(motifChunks)) {
        chunkMotifs <- motifChunks[[j]]

        sampleDeviations <- bplapply(chunkMotifs,
            computeDeviation,
            msites = msites,
            tfBindsites = tfBindsites,
            gcfreqs = gcfreqs,
            binMsites = binMeth,
            enhancer = enhancer,
            ignoreStrand = ignoreStrand,
            BPPARAM = BPPARAM
        )
        names(sampleDeviations) <- chunkMotifs

        writeBlockToSink(
            lapply(sampleDeviations, function(x) x$dev),
            devGrid, index, j, devSink
        )
        writeBlockToSink(
            lapply(sampleDeviations, function(x) x$exp_dev),
            expGrid, index, j, expSink
        )
        rm(sampleDeviations)
    }
    rm(msites)
    cleanMem()
    log_info("Finished processing ", sampleName)
    invisible(NULL)
}

#' @title assembleCoreResult
#' @description Close the sinks and assemble the deviations, their row-wise
#' Z-scores and the expected deviations into a result object.
#' @param devSink The sink holding the bias-corrected deviations.
#' @param expSink The sink holding the expected deviations.
#' @param samples A \code{data.frame} with one row per sample.
#' @return a \code{methylTFRdeviations} object.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray DelayedArray close
#' @importFrom methods as new
#' @keywords internal
assembleCoreResult <- function(devSink, expSink, samples) {
    DelayedArray::close(devSink)
    DelayedArray::close(expSink)
    deviation <- as.matrix(t(as(devSink, "DelayedArray")))
    expDev <- as.matrix(t(as(expSink, "DelayedArray")))

    se <- SummarizedExperiment(
        assays = list(
            deviations = deviation,
            z = computeRowZScore(deviation),
            expected = expDev
        ),
        colData = samples,
        rowData = DataFrame(motifs = row.names(deviation))
    )
    return(new("methylTFRdeviations", se))
}

#' @title methylTFRCore
#' @description Internal engine shared by \code{\link{runMethylTFR}} and
#' \code{\link{runMethylTFRFromRnBeads}}. It validates the motif set,
#' allocates the on-disk sinks, iterates over samples, computes per-motif
#' deviations in chunks and assembles the resulting
#' \code{methylTFRdeviations} object.
#'
#' The only difference between the two public entry points is where the
#' per-sample methylation calls come from. That difference is isolated in
#' the \code{msitesFun} argument, so both entry points share identical
#' numerical behaviour.
#' @inheritParams runMethylTFR
#' @param sampleIds A character vector of sample identifiers. Used for the
#' column names of the resulting object and to size the sinks.
#' @param msitesFun A function of a single integer \code{i} returning the
#' methylation sites of sample \code{i} as a \code{GRanges} with the
#' metadata columns \code{score} and \code{coverage}.
#' @param samples A \code{data.frame} of sample annotation with one row per
#' entry of \code{sampleIds}, used as \code{colData}.
#' @return a \code{methylTFRdeviations} object with bias-corrected
#' deviations, row-wise Z-scores and expected deviations.
#' @importFrom IRanges subsetByOverlaps
#' @importFrom logger log_success
#' @keywords internal
methylTFRCore <- function(
    sampleIds, msitesFun, samples, tfBindsites, gcfreqs, gcDist,
    chunkSize = 20, threads = 1, enhancer = NULL, ignoreStrand = TRUE,
    tempDir = file.path(tempdir(), "methylTFR")
) {
    checkCoreInputs(sampleIds, msitesFun, samples, tempDir)
    motifs <- validCoreMotifs(tfBindsites, gcfreqs)
    BPPARAM <- bpparamFromThreads(threads)

    # Split the motifs into chunks
    numChunks <- ceiling(length(motifs) / chunkSize)
    motifChunks <- split(motifs, rep(seq_len(numChunks),
        each = chunkSize,
        length.out = length(motifs)
    ))

    # Temporary HDF5 sinks, removed when the run ends
    devSink <- createSink(sampleIds, motifs, tempDir = tempDir)
    expSink <- createSink(sampleIds, motifs, tempDir = tempDir)
    on.exit(unlink(c(devSink$path, expSink$path)), add = TRUE)

    devGrid <- setGrid(sampleIds, motifChunks)
    expGrid <- setGrid(sampleIds, motifChunks)

    if (!is.null(enhancer)) {
        gcDist <- subsetByOverlaps(gcDist,
            enhancer,
            ignore.strand = ignoreStrand
        )
    }

    for (i in seq_along(sampleIds)) {
        processCoreSample(
            index = i, sampleIds = sampleIds, msitesFun = msitesFun,
            motifChunks = motifChunks, tfBindsites = tfBindsites,
            gcfreqs = gcfreqs, gcDist = gcDist, devGrid = devGrid,
            expGrid = expGrid, devSink = devSink$sink,
            expSink = expSink$sink, BPPARAM = BPPARAM,
            enhancer = enhancer, ignoreStrand = ignoreStrand
        )
    }
    log_success("Computed all deviations successfully")

    return(assembleCoreResult(devSink$sink, expSink$sink, samples))
}
