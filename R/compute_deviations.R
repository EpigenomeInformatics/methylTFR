#' @title read_methylome
#' @description read_methylome is a function to import methylation data into R Genomic
#'  Range object. Bed file should be processed in the pipeline developed by
#'  Fabian Müller and Christoph Bock. (EPP format)
#' @param  - bedfilename which contains methylation data EPP format
#' @return \code{GenomicRange} object with methylation, coverage information
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom stringr str_split str_replace
read_methylome <- function(filename, type = c("EPP", "BisSNP")) {
  if (!file.exists(filename)) {
    stop(paste(filename, " doesn't exist or path is incorrect !!"))
  }
  type <- match.arg(type)
  if (type == "EPP") {
    # Parse EPP file format
    msites <- fread(filename, header = FALSE, showProgress = FALSE)
    mcov <- unlist(stringr::str_split(msites$V4, "/"))
    mcov <- as.numeric(stringr::str_replace(mcov, "'", ""))
    cov <- mcov[seq_along(mcov) %% 2 == 0]
    mscore <- round(msites$V5 / 1000, 3)
  }
  if (type == "BisSNP") {
    # Parse BisSNP Tab-Separated file format
    msites <- fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    mscore <- round(msites$V4 / 100, 3)
    cov <- msites$V5
  }
  msites <- GenomicRanges::GRanges(
    seqnames = msites$V1,
    ranges = IRanges::IRanges(start = msites$V2, end = msites$V2 + 2),
    strand = msites$V6,
    score = mscore,
    methylation = msites$V4,
    coverage = cov
  )
  return(msites)
}


#' calculate_expmeth
#'
#'      This function is used to calculate genome-wide expected methylation
#'  for each motif.
#'
#' @param msites - imported methylation sites
#' @param gcdist - Genome wide GC distribution from (\code{methylTFRann})
#' @param gcfreq - GC bin frequency table (matrix) from (\code{methylTFRann})
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
calculate_expmeth <- function(msites, gcdist, gcfreq) {
  hits <- findOverlaps(msites, gcdist, type = "within")
  gcmap <- data.table(
    mscore = msites[hits@from]$score,
    gcbin = gcdist[hits@to]$GC_bin
  )
  exp_meth <- gcmap[, .(avg_mscore = mean(mscore)), by = gcbin]
  exp_meth <- exp_meth[, .(mscore = mean(avg_mscore)), by = gcbin]
  exp_meth <- exp_meth[order(gcbin)]
  exp_meth <- as.matrix(exp_meth)
  exp.data <- t(gcfreq) %*% exp_meth[, 2]
  mpos <- seq(-floor(length(exp.data) / 2), floor(length(exp.data) / 2))
  exp.methyl <- data.table(x = mpos, avg_methyl = exp.data)
  colnames(exp.methyl) <- c("x", "exp_avg_methyl")
  return(exp.methyl)
}


#' @title compute_deviation
#' @description compute_deviation is a function to calculate the deviation in transcription factor
#' @param motifs        - list of motifs
#' @param msites       - imported methylation sites
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions from (\code{methylTFRann})
#' @param gcfreqs      - GC bin frequency tables (matrices for multiple motif) from (\code{methylTFRann})
#' @param gcdist       - Genome wide GC distribution from (\code{methylTFRann})
#' @param enhancer     - Specific regions like distal motif
#' @return deviation score for a given motif
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end mcols
#' @importFrom data.table data.table setDT
#' @importFrom dplyr n left_join
compute_deviation <- function(motif, msites, tf_bindsites, gcfreqs,
                              gcdist, enhancer = NULL) {
  tfbs <- tf_bindsites[[motif]]
  gcfreq <- gcfreqs[[motif]]
  w <- width(tfbs)[1]
  tfbs <- resize(tfbs, w + 101, fix = "center")

  if (!is.null(enhancer)) {
    d_hits <- findOverlaps(tfbs, enhancer)
    tfbs <- tfbs[d_hits@from]
  }

  exp_meth <- calculate_expmeth(msites, gcdist, gcfreq)
  hits <- findOverlaps(msites, tfbs, type = "within")

  mcols(tfbs)$mid_point <- round(end(tfbs) + ((start(tfbs) - end(tfbs)) / 2))
  x <- start(msites[hits@from]) - tfbs[hits@to]$mid_point
  plot.data <- data.table::data.table(
    x = x,
    y1 = msites[hits@from]$score,
    y2 = msites[hits@from]$coverage
  )

  # Convert plot.data and exp_meth to data.tables
  setDT(plot.data)
  setDT(exp_meth)

  # Group by 'x' and summarize in plot.data
  plot.data <- plot.data[, .(n = .N, avg_methyl = mean(y1), avg_cov = mean(y2)), by = x]

  # GC bias correction
  dt_join <- merge(plot.data, exp_meth, by = "x")

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
#' @param genome       - human genome version default: hg38
#' @param threads      - thread count for parallel processing
#' @param enhancer     - Specific regions like distal motif
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions from (\code{methylTFRann})
#' @param gcfreqs      - GC bin frequency tables (matrices for multiple motif) from (\code{methylTFRann})
#' @param gc_dist       - Genome wide GC distribution from (\code{methylTFRann})
#' @param sampleColName - column name of the sample bed file in the annotation file
#' @param chunkSize     - chunk size for parallel processing of motifs
#' @param full_path     - if TRUE, the bed file path in the annotation file is full path
#' @param annfile       - if provided, the sample annotation file is not read from the sample_dir
#' @return deviation score matrix for all samples
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
#' @importFrom parallel mclapply
#' @importFrom logger log_info log_error log_appender appender_file
run_methyltfr <- function(
    sample_ann, sample_dir, genome = "hg38", tf_bindsites = NULL,
    gcfreqs = NULL, gc_dist = NULL, sampleColName = "bedFile", chunkSize = 20,
    full_path = FALSE, annfile = NULL, threads = 1, enhancer = NULL, filetype = c("EPP", "BisSNP")) {
  if (is.null(genome)) {
    logger::log_error("Please provide the genome version !!")
  }
  if (genome != "hg38" && any(sapply(list(tf_bindsites, gcfreqs, gc_dist), is.null))) {
    logger::log_error("Please provide the tf_bindsites, gcfreqs and gc_dist for the provided genome!")
  }
  if (any(sapply(list(tf_bindsites, gcfreqs, gc_dist), is.null))) {
    if (genome == "hg38") {
      logger::log_info("Loading the hg38 package")
      if (!require("methylTFRann")) {
        logger::log_error("methylTFRann package is not installed in your environment !!")
      }
      tf_bindsites <- getTFbindsites()
      gc_dist <- getGenomeGC()
      gcfreqs <- getGCfreq()
    }
  }
  # TODO add type checker for sample_ann
  if (!is.null(annfile) & is.character(annfile)) {
    sample_ann <- annfile
  } else {
    annfile <- file.path(sample_dir, sample_ann)
  }
  if (!file.exists(annfile)) {
    logger::log_error(" %s does not exist, please check the file path !!", annfile)
  }
  if (endsWith(annfile, ".csv")) {
    samples <- read.table(annfile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  } else {
    samples <- read.table(annfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  if (full_path) {
    files_list <- samples[, sampleColName]
  } else {
    files_list <- file.path(sample_dir, samples[, sampleColName])
  }
  if (!all(file.exists(files_list))) {
    logger::log_error("Some of the bed files are not exist, please check the file path !!")
  }
  logger::log_info("The samples and annotation package are loaded successfully !!")

  motifs <- names(gcfreqs)
  # deviation <- matrix(NA, nrow = length(motifs), ncol = length(files_list))

  # Split the motifs into chunks
  numChunks <- ceiling(length(motifs) / chunkSize)
  motif_chunks <- split(motifs, rep(1:numChunks, each = chunkSize, length.out = length(motifs)))

  tempfile <- tempfile(pattern = "methylTFR", tmpdir = tempdir(), fileext = ".h5")
  # Create a sink for each region type
  sink <- HDF5Array::HDF5RealizationSink(
    dim = c(length(files_list), length(motifs)),
    dimnames = list(basename(files_list), motifs),
    type = "double",
    filepath = tempfile,
    name = paste0("methylTFRmat"), level = 6
  )
  logger::log_info(paste0("Initializing the sink file in temp: ", tempfile))
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

      sample_deviations <- mclapply(chunk_motifs,
        compute_deviation,
        msites = msites,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcdist = gc_dist,
        enhancer = enhancer,
        mc.cores = threads
      )
      # row_indices <- match(chunk_motifs, motifs)

      # Write the block to the sink
      DelayedArray::write_block(
        block = as.matrix(t(unlist(sample_deviations))),
        viewport = grid[[as.integer(i), as.integer(j)]], sink = sink
      )
      # deviation[row_indices, i]  <- unlist(sample_deviations)
      rm(sample_deviations)
      methylTFR:::cleanMem()
    }
    rm(msites)
    methylTFR:::cleanMem()
    logger::log_info("Finished processing ", sample_name)
  }

  # Close the sink
  DelayedArray::close(sink)
  deviation <- t(as(sink, "DelayedArray"))
  deviation <- as.matrix(deviation)
  file.remove(tempfile) # TODO find a better solution for this
  # rownames(deviation) <- motifs
  # colnames(deviation) <- basename(files_list)
  return(deviation)
}

#' @title cleanMem
#' @description cleanMem is a function to clean the memory
#' @param iter.gc - number of times to run the garbage collector
#' @return NULL
#' @keywords internal
#'
cleanMem <- function(iter.gc = 1L) {
  for (i in 1:iter.gc) {
    gc()
  }
  invisible(NULL)
}
