#' @title run_methyltfr
#' @description This function is a wrapper function to calculate the deviation in transcription factor
#'  footprint base for all given motifs using parallel package
#' @param sample_ann A tab seperated file contains sample annotations
#' @param sample_dir The directory where all bed file and annotation file stored
#' @param threads Thread count for parallel processing
#' @param enhancer a \code{GRanges} object specifying regions such as distal motif (optional)
#' @param tf_bindsites a \code{GRangesList} object contains tf binding sites positions
#' @param gcfreqs a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param sampleColName column name of the sample bed file in the annotation file
#' @param chunkSize Chunk size for parallel processing of motifs (default: 20)
#' @param full_path if TRUE, the bed file path in the annotation file is full path
#' @param annfile if provided, the sample annotation file is not read from the sample_dir
#' @param filetype file type of the bed file, currently supported: bissnp,epp,allc,bismarkcytosine,bismarkcov,encode
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom parallel mclapply
#' @importFrom DelayedArray write_block close ArbitraryArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom logger log_info log_error
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
#' @importFrom methods as new
#' @importFrom stats sd
#' @return a \code{methylTFRdeviations} object with deviation and zscore
#' @export
run_methyltfr <- function(
    sample_ann, sample_dir, tf_bindsites = NULL,
    gcfreqs = NULL, gc_dist = NULL, sampleColName = "bedFile", chunkSize = 20,
    full_path = FALSE, annfile = NULL, threads = 1, enhancer = NULL,
    filetype = NULL, ignoreStrand = TRUE) {

  if (!tolower(filetype) %in% c("bissnp", "epp", "allc", "bismarkcytosine", "bismarkcov", "encode")) {
    logger::log_error("Please provide a valid file type")
  }
  if (!is.logical(ignoreStrand)) {
    logger::log_warn("Found invalid strand option, using the default")
    ignoreStrand <- TRUE
  }
  if (is.null(sampleColName) || !is.character(sampleColName)) {
    logger::log_error("Please provide a valid sample column name")
  }
  if (!is.numeric(chunkSize) || chunkSize < 1) {
    logger::log_warn("Invalid chunk size detected, using default chunk size")
    chunkSize <- 20
  }
  if (any(sapply(list(tf_bindsites, gcfreqs, gc_dist), is.null))) {
    logger::log_error("Please load the annotation objects for given genome.")
  }
  if (is.null(annfile) || !is.character(annfile)) {
    if (is.null(sample_ann) || !is.character(sample_ann)) {
      logger::log_error("Please provide a valid sample annotation file")
    }
    if (is.null(sample_dir) || !is.character(sample_dir)) {
      logger::log_error("Please provide a valid sample directory")
    }
    if (!dir.exists(sample_dir)) {
      logger::log_error("Sample directory does not exist, please check the directory path")
    }
  }
  if (!is.numeric(threads) || threads < 1) {
    logger::log_warn("Invalid thread count detected, using default thread count")
    threads <- 1
  }
  if (!is.logical(full_path)) {
    logger::log_error("Invalid full path flag detected, please provide a valid logical value")
  }
  if (is.null(annfile) || !is.character(annfile)) {
    annfile <- file.path(sample_dir, sample_ann)
  }
  if (!file.exists(annfile)) {
    logger::log_error(" %s does not exist, please check the file path !!", annfile)
  }
  if (endsWith(annfile, ".csv")) {
    samples <- read.table(annfile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  }
  if (endsWith(annfile, ".tsv")) {
    samples <- read.table(annfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    logger::log_error("Please provide a valid annotation file with .csv or .tsv extension")
  }
  if (full_path) {
    files_list <- samples[, sampleColName]
  } else {
    files_list <- file.path(sample_dir, samples[, sampleColName])
  }
  if (!all(file.exists(files_list))) {
    logger::log_error("Some of the files does not exist, please check the file path!")
  }
  logger::log_success("The samples are successfully located")

  motifs <- names(gcfreqs)
  # Split the motifs into chunks
  numChunks <- ceiling(length(motifs) / chunkSize)
  motif_chunks <- split(motifs, rep(1:numChunks, each = chunkSize, length.out = length(motifs)))

  # Create a temp sinks
  dev_sink <- create_sink(files_list, motifs)
  z_sink <- create_sink(files_list, motifs)

  # set the grid
  dev_grid <- set_grid(files_list, motif_chunks)
  z_grid <- set_grid(files_list, motif_chunks)

  for (i in seq_along(files_list)) {
    bedfile <- files_list[i]
    sample_name <- basename(bedfile)
    msites <- read_methylome(bedfile, type = filetype)
    logger::log_info("Processing ", sample_name)

    # Process motifs in chunks
    for (j in seq_along(motif_chunks)) {
      # Get the current chunk
      chunk_motifs <- motif_chunks[[j]]

      sample_deviations <- parallel::mclapply(chunk_motifs,
        computeDeviation,
        msites = msites,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcdist = gc_dist,
        enhancer = enhancer,
        mc.cores = threads,
        ignoreStrand = ignoreStrand,
        intermediate = TRUE
      )

      # Write the block to the sink
      write_block_to_sink(lapply(sample_deviations, function(x) x$obs_dev),dev_grid, i, j, dev_sink)
      write_block_to_sink(lapply(sample_deviations, function(x) x$exp_dev),z_grid, i, j, z_sink)
      rm(sample_deviations)
    }
    rm(msites)
    cleanMem()
    logger::log_info("Finished processing ", sample_name)
  }
  logger::log_success("Computed all deviations successfully")

  # Close the sink
  DelayedArray::close(dev_sink)
  DelayedArray::close(z_sink)
  deviation <- as.matrix(t(as(dev_sink, "DelayedArray")))
  exp_dev <- as.matrix(t(as(z_sink, "DelayedArray")))

  # Compute the sd and normalize the deviation
  #norm_dev <- deviation - matrixStats::colMeans2(exp_dev, na.rm = TRUE)
  sd <- apply(exp_dev, 2, sd, na.rm = TRUE)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = deviation, z = as.matrix(deviation / sd)),
    colData = samples, rowData = DataFrame(motifs = row.names(deviation)))
  return(new("methylTFRdeviations", se)) 
  #return(as.matrix(deviation))
}