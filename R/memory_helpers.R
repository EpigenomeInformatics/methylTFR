#' @title createSink
#' @description Create a temporary HDF5 sink for storing methylTFR results.
#' @param filesList A character vector of sample names.
#' @param motifs A character vector of motifs.
#' @param tempDir Directory the HDF5 file is written to. Created if it does
#' not exist.
#' @param pattern Prefix of the temporary file name.
#' @param fileext File extension of the temporary file.
#' @param verbose A logical indicating whether to print messages.
#' @return A list with the elements \code{sink}, the
#' \code{HDF5RealizationSink}, and \code{path}, the file it writes to.
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom logger log_info
#' @keywords internal
createSink <- function(
    filesList, motifs, tempDir = file.path(tempdir(), "methylTFR"),
    pattern = "methylTFR", fileext = ".h5", verbose = TRUE
) {
    if (!dir.exists(tempDir)) {
        dir.create(tempDir, recursive = TRUE, showWarnings = FALSE)
    }
    path <- tempfile(pattern = pattern, tmpdir = tempDir, fileext = fileext)

    sink <- HDF5Array::HDF5RealizationSink(
        dim = c(length(filesList), length(motifs)),
        dimnames = list(basename(filesList), motifs),
        type = "double",
        filepath = path,
        name = "methylTFRmat", level = 6
    )

    if (verbose) {
        logger::log_info(paste0("Initializing the temp sink: ", path))
    }

    return(list(sink = sink, path = path))
}

#' @title setGrid
#' @description Set the grid for the methylTFR sink
#' @param filesList A character vector of file names
#' @param motifChunks A list of motif chunks
#' @return A methylTFR grid
#' @keywords internal
#' @importFrom DelayedArray ArbitraryArrayGrid
setGrid <- function(filesList, motifChunks) {
    DelayedArray::ArbitraryArrayGrid(list(
        cumsum(lengths(filesList)),
        cumsum(lengths(motifChunks))
    ))
}

#' @title writeBlockToSink
#' @description Write the block to the methylTFR sink
#' @param devValues A numeric vector of deviation values
#' @param grid A methylTFR grid
#' @param i An integer specifying the row index
#' @param j An integer specifying the column index
#' @param sink A methylTFR sink
#' @importFrom DelayedArray write_block
#' @return Invisible \code{NULL}. Called for its effect on the sink.
#' @keywords internal
writeBlockToSink <- function(devValues, grid, i, j, sink) {
    DelayedArray::write_block(
        block = as.matrix(t(unlist(devValues))),
        viewport = grid[[as.integer(i), as.integer(j)]],
        sink = sink
    )
    rm(devValues)
    cleanMem()
}

#' @title cleanMem
#' @description cleanMem is a function to clean the memory
#' @param iter.gc - number of times to run the garbage collector
#' @return invisible NULL
#' @keywords internal
cleanMem <- function(iter.gc = 1L) {
    for (i in seq_len(iter.gc)) {
        gc()
    }
    invisible(NULL)
}
