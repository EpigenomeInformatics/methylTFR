#' @title compute_gc
#' @description Process motifs to matrix containing GC frequencies per motif
#' @param motif Motif name
#' @importFrom Biostrings DNAString letterFrequencyInSlidingView
#' @param X DNA sequence
#' @return matrix contatining gc frequency table
#' @author Sarath Kumar
#' @keywords internal
compute_gc <- function(X) {
  x <- DNAString(as.character(X))
  # center around
  nucfreqs <- letterFrequencyInSlidingView(x,
    view.width = 30,
    letters = c("A", "C", "G", "T")
  )
  gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  return(gc)
}

#' @title convert_to_bins
#' @description Convert GC frequencies to bins
#' @param x GC frequencies
#' @param gc_bin GC bin
#' @return table containing GC bins
#' @keywords internal
#' @author Sarath Kumar
convert_to_bins <- function(x, gc_bin) {
  bin <- cbind(1:length(x), findInterval(x, gc_bin, rightmost.closed = TRUE))
  return(bin)
}

#' @title convert_to_matrix
#' @description Convert GC bins to matrix
#' @param bins GC bins
#' @return matrix contatining GC bins
#' @keywords internal
#' @author Sarath Kumar
convert_to_matrix <- function(bins) {
  mat <- matrix(0, 5, dim(bins)[1])
  bins <- na.omit(bins)
  mat[cbind(bins[, 2], bins[, 1])] <- 1
  return(mat)
}


#' @title get_gcdist
#' @description helper function to calculate GC distribution
#' @param chr Chromosome name
#' @param chr_len Chromosome length
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @return matrix contatining GC distribution
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom IRanges IRanges
#' @keywords internal
get_gcdist <- function(chr, chr_len, genome) {
  qgr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = seq(1, chr_len[chr], 30), width = 30)
  )
  # last interval might be out of chromosome length
  qgr <- qgr[-length(qgr)]
  dna_seq <- getSeq(genome, qgr)
  # calculate GC content
  nucfreqs <- letterFrequency(dna_seq, c("A", "C", "G", "T"))
  gc_tmp <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  gc_tmp <- na.omit(gc_tmp)

  return(gc_tmp)
}


#' @title compute_gc_genome
#' @description helper function to calculate GC distribution of a chromosome
#' @param chr Chromosome name
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @keywords internal
#' @return data.frame contatining GC distribution
compute_gc_genome <- function(chr) {
  qgr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = seq(1, chr_len[chr], 30), width = 30)
  )
  qgr <- qgr[-length(qgr)]
  dna_seq <- getSeq(genome, qgr)
  nucfreqs <- letterFrequency(dna_seq, c("A", "C", "G", "T"))
  gc_tmp <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  gc_bin <- seq(0, 1, length.out = 6)
  gcbin <- findInterval(gc_tmp, gc_bin, rightmost.closed = TRUE)
  values(qgr) <- DataFrame(GC_bias = gc_tmp, GC_bin = gcbin)
  return(qgr)
}

#' @title calculate_GCdist
#' @description calculating GC distribution and GC frequency table
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @param threads Number of threads
#' @param onlyMain Only main chromosomes
#' @param includeSexChr keep the sex chromosomes (chrX, chrY, chrM)
#' @importFrom parallel mclapply
#' @importFrom GenomeInfoDb seqlengths
#' @return matrix containing GC distribution
#' @export
#' @examples
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' gc_dist <- calculate_GCdist(BSgenome.Hsapiens.UCSC.hg38, threads = 1)
#' }
calculate_GCdist <- function(genome, threads = 1, onlyMain = TRUE, includeSexChr = TRUE) {
  if (!is.numeric(threads)) {
    stop("Number of threads must be numeric!")
  }
  if (!is.logical(onlyMain)) {
    stop("onlyMain must be logical!")
  }
  if (!is.logical(includeSexChr)) {
    stop("includeSexChr must be logical!")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be a BSgenome object")
  }
  chr_len <- seqlengths(genome)
  if (onlyMain) {
    filter <- ifelse(includeSexChr, "^chr[0-9MXY]+$", "^chr[0-9]+$")
    chr_len <- chr_len[grepl(filter, names(chr_len))]
  }
  chr_names <- names(chr_len)

  gc_dist <- parallel::mclapply(chr_names, get_gcdist,
    chr_len = chr_len,
    genome = genome,
    mc.cores = threads
  )
  return(unlist(gc_dist))
}


#' @title processMotifs2GCMatrix
#' @description Process motifs to matrix containing GC frequencies per motif
#' @param motif Motif name
#' @param gc_bin GC bin
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @param tf_bindsites TF bindsites
#' @return Matrix containing GC frequencies per motif
#' @importFrom logger log_info
#' @importFrom Biostrings getSeq
#' @importFrom methods is
#' @author Irem Gunduz
#' @export
#' @examples
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(JASPAR2020)
#' library(methylTFR)
#'
#' motifPFMatrixList <- getMatrixSet(
#'   x = JASPAR2020,
#'   opts = list(species = 9606, all_versions = FALSE, collection = "CORE")
#' )
#' tf_bindsites <- motifBSFromPFMatrixList(motifPFMatrixList[1], BSgenome.Hsapiens.UCSC.hg38, 1)
#' gc_dist <- calculate_gcdist(genome = BSgenome.Hsapiens.UCSC.hg38, threads = 1)
#' gc_bin <- quantile(gc_dist, probs = seq(0, 1, 1 / 5))
#' gc_matrix <- processMotifs2GCMatrix(
#'   tf_bindsites[[1]], names(tf_bindsites)[1],
#'   gc_bin, BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
processMotifs2GCMatrix <- function(tf_bindsites, motif, gc_bin, genome) {
  if (!class(tf_bindsites) %in% c("list", "GRangesList")) {
    stop("tf_bindsites must be a list object")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be a BSgenome object")
  }
  if (!is.character(motif)) {
    stop("motif must be a character")
  }
  if (!is.numeric(gc_bin)) {
    stop("gc_bin must be numeric")
  }
  dna_seq <- getSeq(genome, tf_bindsites[[motif]])
  logger::log_info(paste("Processing compute gc .. ", motif))
  motif_gc <- lapply(dna_seq, compute_gc)
  logger::log_info(paste("Processing convert to bins .. ", motif))
  gcbins <- lapply(motif_gc, convert_to_bins, gc_bin)
  logger::log_info(paste("Processing convert to matrix ... ", motif))
  gcmat <- lapply(gcbins, convert_to_matrix)
  m_gcfreq <- Reduce("+", gcmat, accumulate = FALSE)
  logger::log_info(paste("Normalizing matrix...", motif))
  normalized_matrix <- sweep(as.matrix(m_gcfreq), 2, colSums(as.matrix(m_gcfreq)), FUN = "/")
  return(normalized_matrix)
}


#' @title computeGenomeWideGC
#' @description Compute genome-wide GC distribution for given genome
#' @param genome BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38)
#' @param onlyMain Only main chromosomes
#' @param includeSexChr keep the sex chromosomes (chrX, chrY, chrM), only valid when onlyMain is TRUE
#' @param num_cores Number of cores
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom parallel mclapply
#' @importFrom stats complete.cases
#' @importFrom methods is
#' @return GC distribution GRanges object
#' @author Irem Gunduz
#' @export
#' @examples
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' gc_dist <- calculate_GCdist(BSgenome.Hsapiens.UCSC.hg38, threads = 1)
#' }
computeGenomeWideGC <- function(genome, onlyMain = TRUE, includeSexChr = TRUE, num_cores = 1) {
  if (!is(genome, "BSgenome")) {
    stop("genome must be a BSgenome object")
  }
  chr_len <- seqlengths(genome)
  chr_names <- names(chr_len)
  if (!is.logical(onlyMain)) {
    stop("onlyMain must be logical!")
  }
  if (!is.logical(includeSexChr)) {
    stop("includeSexChr must be logical!")
  }
  if (!is.numeric(num_cores)) {
    stop("Number of cores must be numeric!")
  }
  if (onlyMain) {
    filter <- ifelse(includeSexChr, "^chr[0-9MXY]+$", "^chr[0-9]+$")
    chr_names <- chr_names[grepl(filter, chr_names)]
  }
  t_qgr <- parallel::mclapply(chr_names, compute_gc_genome, mc.cores = num_cores)
  t_qgr <- do.call(c, t_qgr)
  t_qgr <- t_qgr[complete.cases(t_qgr$GC_bias, t_qgr$GC_bin), ]
  return(t_qgr)
}

#' @title motifBSFromPFMatrixList
#' @description This function extracts motif binding sites from a PFMatrixList using a specified genome.
#' @param motifPFMatrixList PFMatrixList object containing motif profiles
#' @param genome BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38)
#' @param threads Number of parallel threads for processing
#' @param onlyMainChr Logical, indicating whether to include only main chromosomes
#' @param includeSexChr Logical, indicating whether to include sex chromosomes (valid when onlyMainChr is TRUE)
#' @return A GRangesList object containing motif binding sites
#' @importFrom GenomicRanges GRanges resize GRangesList
#' @importFrom parallel mclapply
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom BSgenome getSeq
#' @importFrom motifmatchr matchMotifs
#' @importFrom methods is
#' @author Irem Gunduz
#' @export
#' @examples 
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(JASPAR2020)
#' motifPFMatrixList <- TFBSTools::getMatrixSet(
#'   x = JASPAR2020,
#'   opts = list(species = 9606, all_versions = FALSE, collection = "CORE")
#' )
#' result <- processMotifData(motifPFMatrixList[1], BSgenome.Hsapiens.UCSC.hg38, 1)
#' }
motifBSFromPFMatrixList <- function(
    motifPFMatrixList, genome, threads = 2,
    onlyMainChr = TRUE, includeSexChr = TRUE) {
  seqNames <- seqnames(genome)
  if (!is(genome, "BSgenome")) {
    stop("genome must be a BSgenome object")
  }
  if (!is(motifPFMatrixList, "PFMatrixList")) {
    stop("motifPFMatrixList must be a PFMatrixList object")
  }
  if (!is.numeric(threads)) {
    stop("Number of threads must be numeric!")
  }
  if (!is.logical(onlyMainChr)) {
    stop("onlyMainChr must be logical!")
  }
  if (!is.logical(includeSexChr)) {
    stop("includeSexChr must be logical!")
  }
  filter <- ifelse(onlyMainChr, ifelse(includeSexChr, "chr[0-9MXY]+$", "chr[0-9]+$"), seqNames)
  seqNames <- seqNames[grepl(filter, seqNames)]
  motif_names <- names(motifPFMatrixList)
  if (!is.numeric(threads)) {
    stop("Number of threads must be numeric!")
  }
  result_list <- lapply(motif_names, function(mo) {
    logger::log_info(paste0("Started processing ", mo))
    motif_data <- motifPFMatrixList[[mo]]
    tf_binding <- parallel::mclapply(seqNames, function(chr) {
      chr_seq <- getSeq(genome, chr)
      motif_ix <- matchMotifs(motif_data, chr_seq, out = "positions")
      motif_ix <- unlist(motif_ix[[1]])

      curr_motif <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = start(motif_ix), width = unique(width(motif_ix))),
        strand = unlist(motif_ix@elementMetadata["strand"])
      )
      values(curr_motif) <- motif_ix@elementMetadata["score"]
      curr_motif <- resize(curr_motif, motif_width + 400, fix = "center")
      return(curr_motif)
    }, mc.cores = threads)

    tf_binding <- do.call(c, tf_binding)
    logger::log_success(paste0("Finished processing ", mo))
    return(tf_binding)
  })

  logger::log_info("Finished processing all motifs, converting results into GRangesList")
  tf_bindsites <- GRangesList(result_list)
  names(tf_bindsites) <- motif_names
  return(tf_bindsites)
}
