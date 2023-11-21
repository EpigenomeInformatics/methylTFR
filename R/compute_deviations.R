#' @title calculate_expmeth
#' @description  This function is used to calculate genome-wide expected methylation for each motif.
#' @param msites - imported methylation sites
#' @param gcdist - Genome wide GC distribution
#' @param gcfreq - GC bin frequency table (matrix)
#' @param ignoreStrand - if TRUE, it ignores strand info from annotation
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @keywords internal
calculate_expmeth <- function(msites, gcdist, gcfreq, ignoreStrand) {
  hits <- findOverlaps(msites, gcdist, type = "within", ignore.strand = ignoreStrand)
  gcmap <- data.table(
    mscore = msites[hits@from]$score,
    gcbin = gcdist[hits@to]$GC_bin
  )
  exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]
  # exp_meth <- exp_meth[, .(mscore = mean(avg_mscore)), by = gcbin]
  exp_meth <- exp_meth[order(gcbin)]
  exp_meth <- as.matrix(exp_meth)
  exp.data <- t(gcfreq) %*% exp_meth[, 2]
  mpos <- round(seq(-floor(length(exp.data) / 2), floor(length(exp.data) / 2), length.out = length(exp.data)))
  exp.methyl <- data.table(x = mpos, avg_methyl = exp.data)
  colnames(exp.methyl) <- c("x", "exp_avg_methyl")
  return(exp.methyl)
}


#' @title compute_deviation
#' @description compute_deviation is a function to calculate the deviation in transcription factor
#' @param motif        - motif name
#' @param msites       - imported methylation sites
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions
#' @param gcfreqs      - GC bin frequency tables (matrices for multiple motif)
#' @param gcdist       - Genome wide GC distribution
#' @param enhancer     - Specific regions like distal motif
#' @param ignoreStrand - if TRUE, it ignores strand info from annotation
#' @return a \code{numeric} deviation score for a given motif
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table setDT
#' @importFrom stats na.omit
#' @importFrom S4Vectors mcols
#' @import data.table
#' @keywords internal
compute_deviation <- function(motif, msites, tf_bindsites, gcfreqs,
                              gcdist, enhancer, ignoreStrand) {
  tfbs <- tf_bindsites[[motif]]
  gcfreq <- gcfreqs[[motif]]
  w <- width(tfbs)[1] # TODO remove these lines and make sure annot is 500bp
  tfbs <- resize(tfbs, w + 101, fix = "center")

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(tfbs, enhancer, ignoreStrand)
    tfbs <- tfbs[d_hits@from]
  }

  exp_meth <- calculate_expmeth(msites, gcdist, gcfreq, ignoreStrand)
  hits <- findOverlaps(msites, tfbs, type = "within", ignore.strand = ignoreStrand)

  S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) + ((start(tfbs) - end(tfbs)) / 2))
  x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
  sum_meth <- data.table::data.table(
    x = x,
    avg_methyl = msites[hits@from]$score,
    avg_cov = msites[hits@from]$coverage
  )

  # Convert sum_meth and exp_meth to data.tables
  setDT(sum_meth)
  setDT(exp_meth)

  # Group by 'x' and summarize in sum_meth
  sum_meth <- sum_meth[, .(n = .N, avg_methyl = mean(avg_methyl), avg_cov = mean(avg_cov)), by = x]

  # GC bias correction
  dt_join <- merge(sum_meth, exp_meth, by = "x")

  # Calculate 'diff'
  dt_join[, diff := abs(avg_methyl - exp_avg_methyl)]

  # Create intervals and calculate mean
  dt_join[, cuts := cut(x, c(-200, -100, -10, 10, 100, 200))]
  interval_mean <- dt_join[, .(n = .N, mean = mean(diff)), by = cuts]
  # interval_mean$mean <- round(interval_mean$mean,6)
  interval_mean <- na.omit(interval_mean[order(cuts)])

  num_intervals <- nrow(interval_mean)
  if (num_intervals > 0) {
    avg_first_last <- ((interval_mean$mean[1] + interval_mean$mean[num_intervals]) / 2)
    var <- interval_mean$mean[(num_intervals + 1) %/% 2] / avg_first_last
  } else {
    var <- NA
  }
  return(var)
}



#' @title run_methyltfr
#' @description This function is a wrapper function to calculate the deviation in transcription factor
#'  footprint base for all given motifs using parallel package
#' @param sample_ann   - a tab seperated file contains sample annotations
#' @param sample_dir   - directory where all bed file and annotation file stored
#' @param threads      - thread count for parallel processing
#' @param enhancer     - Specific regions like distal motif
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions
#' @param gcfreqs      - GC bin frequency tables (matrices for multiple motif)
#' @param gc_dist       - Genome wide GC distribution
#' @param sampleColName - column name of the sample bed file in the annotation file
#' @param chunkSize     - chunk size for parallel processing of motifs
#' @param full_path     - if TRUE, the bed file path in the annotation file is full path
#' @param annfile       - if provided, the sample annotation file is not read from the sample_dir
#' @param filetype      - file type of the bed file, currently supported: bissnp,epp,allc,bismark
#' @param ignoreStrand - if TRUE, it ignores strand info from annotation
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom parallel mclapply
#' @importFrom DelayedArray write_block close ArbitraryArrayGrid
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom logger log_info log_error
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
#' @importFrom methods as
#' @return deviation score matrix for all samples and motifs
#' @export
run_methyltfr <- function(
    sample_ann, sample_dir, tf_bindsites = NULL,
    gcfreqs = NULL, gc_dist = NULL, sampleColName = "bedFile", chunkSize = 20,
    full_path = FALSE, annfile = NULL, threads = 1, enhancer = NULL,
    filetype = NULL, ignoreStrand = TRUE) {
  if (!tolower(filetype) %in% c("bissnp", "epp", "allc")) {
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
  logger::log_success("The samples and annotation package are loaded successfully")

  motifs <- names(gcfreqs)
  # Split the motifs into chunks
  numChunks <- ceiling(length(motifs) / chunkSize)
  motif_chunks <- split(motifs, rep(1:numChunks, each = chunkSize, length.out = length(motifs)))

  # Create a temp sink
  if (!dir.exists("methylTFR_tmp")) {
    dir.create("methylTFR_tmp")
  }
  tempfile <- tempfile(pattern = "methylTFR", tmpdir = "methylTFR_tmp", fileext = ".h5")
  # Create a sink for each region type
  sink <- HDF5Array::HDF5RealizationSink(
    dim = c(length(files_list), length(motifs)),
    dimnames = list(basename(files_list), motifs),
    type = "double",
    filepath = tempfile,
    name = paste0("methylTFRmat"), level = 6
  )
  logger::log_info(paste0("Initializing the temp sink: ", tempfile))
  # set the grid
  grid <- DelayedArray::ArbitraryArrayGrid(list(
    cumsum(lengths(files_list)),
    cumsum(lengths(motif_chunks))
  ))

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
        compute_deviation,
        msites = msites,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcdist = gc_dist,
        enhancer = enhancer,
        mc.cores = threads,
        ignoreStrand = ignoreStrand
      )

      # Write the block to the sink
      DelayedArray::write_block(
        block = as.matrix(t(unlist(sample_deviations))),
        viewport = grid[[as.integer(i), as.integer(j)]], sink = sink
      )
      rm(sample_deviations)
      cleanMem()
    }
    rm(msites)
    cleanMem()
    logger::log_info("Finished processing ", sample_name)
  }
  logger::log_success("Computed all deviations successfully")
  # Close the sink
  DelayedArray::close(sink)
  deviation <- t(as(sink, "DelayedArray"))
  # file.remove(tempfile) # TODO find a better solution for this

  # se <- SummarizedExperiment::SummarizedExperiment(
  # assays = list(deviations = as.matrix(deviation), z = computeZScore(deviation)),
  # colData = samples, rowData = DataFrame(motifs = row.names(deviation))
  # )
  # return(se) # new("methylTFRDeviations", se)) #TODO
  return(as.matrix(deviation))
}
