#' @title run_methylTFR_RnBeads
#' @description Run methylTFR on a \code{RnBeads} object
#' @param rnb_set a \code{RnBiseqSet}, \code{RnBeadRawSet} or \code{RnBeadSet} object
#' @param tf_bindsites Transcript Factor binding sites from the annotation package
#' @param threads Thread count for parallel processing
#' @param enhancer a \code{GRanges} object specifying enhancer such as distal regions (optional)
#' @param tf_bindsites a \code{GRangesList} object contains tf binding sites positions
#' @param gcfreqs a \code{list} of GC bin frequency tables (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains Genome wide GC distribution
#' @param chunkSize Chunk size for parallel processing of motifs (default: 20)
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @importFrom logger log_info log_success
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray close
#' @importFrom methods is
#' @importFrom RnBeads annotation
#' @author Irem B. Gunduz
#' @export
#' @return a \code{methylTFRdeviations} object with bias-corrected deviation and Z-scores
#'
run_methylTFR_RnBeads <- function(rnb_set, tf_bindsites, threads = 1,
                                  enhancer = NULL, gcfreqs = NULL, gc_dist = NULL,
                                  chunkSize = 20, ignoreStrand = TRUE) {
  if (any(sapply(list(tf_bindsites, gcfreqs, gc_dist), is.null))) {
    stop("Please load the annotation objects for given genome.")
  }
  if (!is(rnb_set, "RnBiseqSet") &&
    !is(rnb_set, "RnBeadRawSet") &&
    !is(rnb_set, "RnBeadSet")) {
    stop("rnb_set must be a RnBeads object")
  }
  if (!is(tf_bindsites, "GRangesList")) {
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
  if (!is.numeric(threads) || threads < 1) {
    RnBeads::logger.warning("Invalid thread count detected, using default thread as 1")
    threads <- 1
  }
  if (!is.numeric(chunkSize) || chunkSize < 1) {
    RnBeads::logger.warning("Invalid chunk size detected, using default chunk size")
    chunkSize <- 20
  }
  if (!is.logical(ignoreStrand)) {
    stop("ignoreStrand must be a logical value")
  }

  # Get the sample names
  sample_names <- samples(rnb_set) # @meth.sites)
  if (is.null(sample_names) || length(sample_names) == 0) {
    stop("No sample names found in the RnBeads object")
  }
  RnBeads::logger.info("Retriveing methylation sites")
  sites <- RnBeads::annotation(rnb_set, type = "sites", add.names = FALSE, include.regions = FALSE)[1:4]

  # Split the motifs into chunks
  motifs <- names(gcfreqs)
  numChunks <- ceiling(length(motifs) / chunkSize)
  motif_chunks <- split(motifs, rep(1:numChunks,
    each = chunkSize,
    length.out = length(motifs)
  ))

  # Create a temp sinks
  dev_sink <- methylTFR:::create_sink(sample_names, motifs, verbose = FALSE)
  z_sink <- methylTFR:::create_sink(sample_names, motifs, verbose = FALSE)

  # set the grid
  dev_grid <- methylTFR:::set_grid(sample_names, motif_chunks)
  z_grid <- methylTFR:::set_grid(sample_names, motif_chunks)


  if (!is.null(enhancer)) {
    d_hits <- suppressWarnings(findOverlaps(gc_dist, enhancer, ignore.strand = ignoreStrand))
    gc_dist <- gc_dist[d_hits@from]
  }
  gc_dist <- as.data.table(gc_dist)

  RnBeads::logger.start("Computing deviations")
  for (i in seq_along(sample_names)) {
    msites <- get_rnbeads_sample(i, rnb_set, sample_names, sites)
    RnBeads::logger.info(paste0("Processing sample ", sample_names[i]))

    # Process motifs in chunks
    for (j in seq_along(motif_chunks)) {
      # Get the current chunk
      chunk_motifs <- motif_chunks[[j]]

      # Compute bias corrected deviations
      sample_deviations <- parallel::mclapply(chunk_motifs,
        computeDeviation,
        msites = msites,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        enhancer = enhancer,
        mc.cores = threads,
        ignoreStrand = ignoreStrand
      )
      names(sample_deviations) <- chunk_motifs


      # Compute expected deviations
      exp_dev <- parallel::mclapply(chunk_motifs, computeExpectedDeviation,
        msites = msites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        sample_deviations = sample_deviations,
        ignoreStrand = ignoreStrand,
        mc.cores = threads
      )

      # Write the block to the sink
      methylTFR:::write_block_to_sink(
        lapply(exp_dev, function(x) x$obs_dev),
        dev_grid, i, j, dev_sink
      )
      methylTFR:::write_block_to_sink(
        lapply(exp_dev, function(x) x$z_score),
        z_grid, i, j, z_sink
      )
      rm(exp_dev, sample_deviations)
    }
    rm(msites)
    methylTFR:::cleanMem()
    RnBeads::logger.info(paste0("Completed processing ", sample_names[i]))
  }

  RnBeads::logger.completed()

  # Close the sink
  DelayedArray::close(dev_sink)
  DelayedArray::close(z_sink)
  obs_dev <- as.matrix(t(as(dev_sink, "DelayedArray")))
  z_dev <- as.matrix(t(as(z_sink, "DelayedArray")))

  # Compute the sd and normalize the deviation
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      deviations = obs_dev, z = z_dev
    ),
    colData = rnb_set@pheno,
    rowData = DataFrame(motifs = row.names(obs_dev))
  )
  return(new("methylTFRdeviations", se))
}

#' @title get_rnbeads_sample
#' @description Get the methylation sites for a sample from a \code{RnBiseqSet} object
#' @param i index of the sample
#' @param rnb_set a \code{RnBiseqSet} object
#' @param sample_names a \code{character} vector of sample names
#' @param sites a \code{GRanges} object of methylation sites
#' @importFrom RnBeads covg meth
#' @importFrom data.table data.table
#' @importFrom GenomicRanges granges
#' @importFrom logger log_info
#' @return a \code{GRanges} object of methylation sites
#' @author Irem Gunduz
#' @keywords internal
get_rnbeads_sample <- function(i, rnb_set, sample_names, sites) {
  sample_name <- basename(sample_names[i])
  mm <- RnBeads::meth(rnb_set, type = "sites", row.names = FALSE, i = NULL, j = i)
  msites <- data.table::data.table(sites, mm)
  cov <- RnBeads::covg(rnb_set, type = "sites", row.names = FALSE, i = NULL, j = i)
  msites <- data.table::data.table(msites, cov)
  colnames(msites) <- c("chr", "start", "end", "strand", "score", "cov")
  msites <- methylTFR:::granges_helper(msites, msites$chr, msites$score, msites$cov, msites$start, msites$end)
  return(msites)
}
