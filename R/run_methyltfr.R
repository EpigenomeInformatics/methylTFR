#' @title run_methyltfr
#' @description This function is a wrapper function to calculate the deviation in transcription factor
#' footprint base for all given motifs per raw samples
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
#' @importFrom logger log_info log_error
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
#' @importFrom methods as new
#' @importFrom stats sd
#' @return a \code{methylTFRdeviations} object with bias-corrected deviation and Z-scores
#' @export
run_methyltfr <- function(
    sample_ann, sample_dir, tf_bindsites = NULL,
    gcfreqs = NULL, gc_dist = NULL, sampleColName = "bedFile", chunkSize = 20,
    full_path = FALSE, annfile = NULL, threads = 1, enhancer = NULL,
    filetype = NULL, ignoreStrand = TRUE) {
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
      stop("Sample directory does not exist, please check the directory path")
    }
  }
  if (!is.numeric(threads) || threads < 1) {
    warning("Invalid thread count detected, using default thread count")
    threads <- 1
  }
  if (!is.logical(full_path)) {
    stop("Invalid full path flag detected, please provide a valid logical value")
  }
  if (is.null(annfile) || !is.character(annfile)) {
    annfile <- file.path(sample_dir, sample_ann)
  }
  if (!file.exists(annfile)) {
    stop(" %s does not exist, please check the file path !!", annfile)
  }
  if (endsWith(annfile, ".csv")) {
    samples <- read.table(annfile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  }
  if (endsWith(annfile, ".tsv")) {
    samples <- read.table(annfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    stop("Please provide a valid annotation file with .csv or .tsv extension")
  }
  if (full_path) {
    files_list <- samples[, sampleColName]
  } else {
    files_list <- file.path(sample_dir, samples[, sampleColName])
  }
  if (!all(file.exists(files_list))) {
    stop("Some of the files does not exist, please check the file path!")
  }
  logger::log_success("The samples are successfully located")

  motifs <- names(gcfreqs)
  # Split the motifs into chunks
  numChunks <- ceiling(length(motifs) / chunkSize)
  motif_chunks <- split(motifs, rep(1:numChunks,
    each = chunkSize,
    length.out = length(motifs)
  ))

  # Create a temp sinks
  dev_sink <- methylTFR:::create_sink(files_list, motifs)
  z_sink <- methylTFR:::create_sink(files_list, motifs)

  # set the grid
  dev_grid <- methylTFR:::set_grid(files_list, motif_chunks)
  z_grid <- methylTFR:::set_grid(files_list, motif_chunks)

  for (i in seq_along(files_list)) {
    bedfile <- files_list[i]
    sample_name <- basename(bedfile)
    msites <- read_methylome(bedfile, type = filetype)
    logger::log_info("Processing ", sample_name)

    # Process motifs in chunks
    for (j in seq_along(motif_chunks)) {
      # Get the current chunk
      chunk_motifs <- motif_chunks[[j]]

      # Compute expected deviations
      exp_dev <- parallel::mclapply(chunk_motifs, computeExpectedDeviation,
        msites = msites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        enhancer = enhancer,
        ignoreStrand = ignoreStrand,
        mc.cores = threads
      )
      # Compute bias corrected deviations
      sample_deviations <- parallel::mclapply(chunk_motifs,
        computeDeviation,
        msites = msites,
        gc_dist = gc_dist,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        enhancer = enhancer,
        mc.cores = threads,
        ignoreStrand = ignoreStrand,
        exp_dev = exp_dev
      )

      # Write the block to the sink
      methylTFR:::write_block_to_sink(
        lapply(sample_deviations, function(x) x$dev),
        dev_grid, i, j, dev_sink
      )
      methylTFR:::write_block_to_sink(
        lapply(sample_deviations, function(x) x$obs_dev),
        z_grid, i, j, z_sink
      )
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
  obs_dev <- as.matrix(t(as(z_sink, "DelayedArray")))

  # create summarized experiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      deviations = obs_dev, z = deviation
    ),
    colData = samples, rowData = DataFrame(motifs = row.names(deviation))
  )
  return(new("methylTFRdeviations", se))
}
