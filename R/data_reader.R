#' @title read_methylome
#' @description read_methylome is a function to import methylation data into a GRanges
#' object. Bed file should be in EPP,ALLC or BisSNP format.
#' @param filename - filename which contains methylation data
#' @param type Type of file format. Currently supported epp, bissnp, bismarkCytosine,
#' bismarkcov,allc and encode
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import data.table
#' @import R.utils
#' @importFrom data.table := .N fread
#' @importFrom stringr str_split str_replace
#' @return a \code{GenomicRange} object with methylation, coverage information
#' @author Irem Gunduz
#' @examples
#' # Read bissnp file
#' bissnp_path <- system.file("extdata", "bissnp.tsv.gz", package = "methylTFR")
#' bissnp <- read_methylome(bissnp_path, "bissnp")
#' @export
read_methylome <- function(filename, type) {
  type <- tolower(type)
  if (!file.exists(filename)) {
    stop(paste(filename, " doesn't exist or path is incorrect!"))
  }
  if (!type %in% c("bissnp", "epp", "allc", "bismarkcytosine", "bismarkcov", "encode")) {
    stop(paste(type, " is not a valid file type!"))
  }
  if (type == "epp") {
    # Parse EPP file format
    msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    mcov <- unlist(stringr::str_split(msites$V4, "/"))
    mcov <- as.numeric(stringr::str_replace(mcov, "'", ""))
    cov <- mcov[seq_along(mcov) %% 2 == 0]
    mscore <- round(msites$V5 / 1000, 6)
    msites <- msites[, .(V1, V2, V3, V6), ]
  }
  if (type == "bissnp") {
    # Parse BisSNP Tab-Separated file format
    msites <- data.table::fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    mscore <- round(msites$V4 / 100, 6)
    cov <- msites$V5
    msites <- msites[, .(V1, V2, V3, V6), ]
  }
  if (type == "allc") {
    # Parse ALLC Tab-Separated file format
    allc <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(allc) < 6) {
      logger::log_warn(paste0(filename, " is an invalid allc file!"))
      stop("allc file must contain at least 6 columns!")
    }
    cov <- allc$V6
    mscore <- round(allc$V5 / cov, 6)
    msites <- allc[, .(V1, V2, V2, V3), ]
  }
  if (type == "bismarkcytosine") {
    # Parse bismarkCytosine file format
    msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(msites) < 5) {
      logger::log_warn(paste0(filename, " is an invalid bismarkCytosine file!"))
      stop("bismarkCytosine file must contain at least 5 columns!")
    }
    cov <- (msites$V5 + msites$V4) # methylated + unmethylated
    mscore <- round(msites$V4 / cov, 6)
    msites <- msites[, .(V1, V2, V2, V3), ]
  }
  if (type == "bismarkcov") {
    # Parse bismarkCov file format
    msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(msites) < 6) {
      logger::log_warn(paste0(filename, " is an invalid bismarkCov file!"))
      stop("bismarkCov file must contain at least 6 columns!")
    }
    mscore <- round(msites$V4 / 100, 6)
    cov <- msites$V5 + msites$V6
    msites <- msites[, .(V1, V2, V3), ]
    msites$strand <- rep("*", nrow(msites))
  }
  if (type == "encode") {
    # Parse encode file format
    msites <- data.table::fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    if (ncol(msites) < 11) {
      logger::log_warn(paste0(filename, " is an invalid encode file!"))
      stop("encode file must contain 11 columns!")
    }
    mscore <- round(msites$V11 / msites$V10, 6)
    cov <- msites$V10
    msites <- msites[, .(V1, V2, V3, V6), ]
  }
  colnames(msites) <- c("chr", "start", "end", "strand")
  # convert to GRanges object
  gr_obj <- granges_helper(
    grobj = msites,
    chr = msites$chr,
    mscore = mscore,
    cov = cov,
    startP = msites$start,
    endP = msites$end
  )
  return(gr_obj)
}


#' @title granges_helper
#' @description  helper function to convert msites in EPP GRanges format
#' @param  grobj - GRanges object
#' @param  chr - chromosome
#' @param  mscore - methylation score
#' @param  cov - methylation coverage
#' @param  startP - start position
#' @param  endP - end position
#' @return a \code{GenomicRanges} object with EPP format
#' @keywords internal
granges_helper <- function(
    grobj, chr, mscore, cov, # meth,
    startP, endP) {
  gr_obj <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = startP, end = endP),
    strand = grobj$strand,
    score = mscore,
    coverage = cov
  )
  if (any(is.na(gr_obj$score))) {
    # Remove for NaN values
    gr_obj <- gr_obj[!is.nan(gr_obj$score)]
  }
  return(gr_obj)
}
