#' @title Process motifs to matrix containing GC frequencies per motif
#' @description Process motifs to matrix containing GC frequencies per motif
#' @param motif Motif name
#' @importFrom Biostrings DNAString letterFrequencyInSlidingView
#' @param X DNA sequence
#' @return gc frequency
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

#' @title Convert GC frequencies to bins
#' @description Convert GC frequencies to bins
#' @param x GC frequencies
#' @param gc_bin GC bin
#' @return GC bins
#' @keywords internal
#' @author Sarath Kumar
convert_to_bins <- function(x, gc_bin) {
  bin <- cbind(1:length(x), findInterval(x, gc_bin, rightmost.closed = TRUE))
  return(bin)
}

#' @title Convert GC bins to matrix
#' @description Convert GC bins to matrix
#' @param bins GC bins
#' @return GC matrix
#' @keywords internal
#' @author Sarath Kumar
convert_to_matrix <- function(bins) {
  mat <- matrix(0, 5, dim(bins)[1])
  bins <- na.omit(bins)
  mat[cbind(bins[, 2], bins[, 1])] <- 1
  return(mat)
}


#' @title helper function to calculate GC distribution
#' @description helper function to calculate GC distribution
#' @param chr Chromosome name
#' @param chr_len Chromosome length
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @return GC distribution
#' @keywords internal
get_gcdist <- function(chr, chr_len, genome) {
  require(GenomicRanges)
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


#' @title helper function to calculate GC distribution of a chromosome
#' @description helper function to calculate GC distribution of a chromosome
#' @param chr Chromosome name
#' @keywords internal
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

#' @title calculating GC distribution and GC frequency table
#' @description calculating GC distribution and GC frequency table
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @param threads Number of threads
#' @param onlyMain Only main chromosomes
#' @param includeSexChr keep the sex chromosomes (chrX, chrY, chrM)
#' @importFrom parallel mclapply
#' @return GC distribution
#' @export
calculate_gcdist <- function(genome, threads = 4, onlyMain = TRUE, includeSexChr = TRUE) {
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


#' @title Process motifs to matrix containing GC frequencies per motif
#' @description Process motifs to matrix containing GC frequencies per motif
#' @param motif Motif name
#' @param gc_bin GC bin
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @param tf_bindsites TF bindsites
#' @return Matrix containing GC frequencies per motif
#' @importFrom logger log_info
#' @export
#' @author Irem Gunduz
processMotifs2GCMatrix <- function(tf_bindsites, motif, gc_bin, genome) {
  tfbs <- tf_bindsites[[motif]]
  dna_seq <- getSeq(genome, tfbs)
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


#' @title Compute genome-wide GC distribution for given genome
#' @description Compute genome-wide GC distribution for given genome
#' @param genome Genome (e.g. BSgenome.Hsapiens.UCSC.hg38)
#' @param onlyMain Only main chromosomes
#' @param includeSexChr keep the sex chromosomes (chrX, chrY, chrM), only valid when onlyMain is TRUE
#' @param num_cores Number of cores
#' @importFrom GenomicRanges GRanges
#' @return GC distribution GRanges object
#' @export
#' @author Irem Gunduz
computeGenomeWideGC <- function(genome, onlyMain = TRUE, includeSexChr = TRUE, num_cores = 1) {
  chr_len <- seqlengths(genome)
  chr_names <- names(chr_len)

  if (onlyMain) {
    filter <- ifelse(includeSexChr, "^chr[0-9MXY]+$", "^chr[0-9]+$")
    chr_names <- chr_names[grepl(filter, chr_names)]
  }
  t_qgr <- parallel::mclapply(chr_names, compute_gc_genome, mc.cores = num_cores)
  t_qgr <- do.call(c, t_qgr)
  t_qgr <- t_qgr[complete.cases(t_qgr$GC_bias, t_qgr$GC_bin), ]
  return(t_qgr)
}

