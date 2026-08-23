#' @title create_sink
#' @description Create a temp sink for storing methylTFR results
#' @param files_list A character vector of file names
#' @param motifs A character vector of motifs
#' @param temp_dir A character vector specifying the temp directory
#' @param pattern A character vector specifying the pattern for the temp file
#' @param fileext A character vector specifying the file extension for the temp file
#' @param verbose A logical indicating whether to print messages
#' @return A methylTFR sink
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom logger log_info
#' @keywords internal
create_sink <- function(
  files_list, motifs, temp_dir = "methylTFR_tmp", pattern = "methylTFR",
  fileext = ".h5", verbose = TRUE
) {
  # Create a temp sink
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }
  tempfile <- tempfile(pattern = pattern, tmpdir = temp_dir, fileext = fileext)

  # Create a sink for each region type
  sink <- HDF5Array::HDF5RealizationSink(
    dim = c(length(files_list), length(motifs)),
    dimnames = list(basename(files_list), motifs),
    type = "double",
    filepath = tempfile,
    name = paste0("methylTFRmat"), level = 6
  )

  if (verbose) {
    logger::log_info(paste0("Initializing the temp sink: ", tempfile))
  }

  return(sink)
}

#' @title set_grid
#' @description Set the grid for the methylTFR sink
#' @param files_list A character vector of file names
#' @param motif_chunks A list of motif chunks
#' @return A methylTFR grid
#' @keywords internal
#' @importFrom DelayedArray ArbitraryArrayGrid
set_grid <- function(files_list, motif_chunks) {
  # set the grid
  grid <- DelayedArray::ArbitraryArrayGrid(list(
    cumsum(lengths(files_list)),
    cumsum(lengths(motif_chunks))
  ))

  return(grid)
}

#' @title write_block_to_sink
#' @description Write the block to the methylTFR sink
#' @param dev_values A numeric vector of deviation values
#' @param grid A methylTFR grid
#' @param i An integer specifying the row index
#' @param j An integer specifying the column index
#' @param sink A methylTFR sink
#' @importFrom DelayedArray write_block
#' @return sink
#' @keywords internal
write_block_to_sink <- function(dev_values, grid, i, j, sink) {
  # Write the block to the sink
  sink <- DelayedArray::write_block(
    block = as.matrix(t(unlist(dev_values))),
    viewport = grid[[as.integer(i), as.integer(j)]],
    sink = sink
  )
  rm(dev_values)
  cleanMem()
}

#' @title cleanMem
#' @description cleanMem is a function to clean the memory
#' @param iter.gc - number of times to run the garbage collector
#' @return invisible NULL
#' @keywords internal
cleanMem <- function(iter.gc = 1L) {
  for (i in seq_along(iter.gc)) {
    gc()
  }
  invisible(NULL)
}
