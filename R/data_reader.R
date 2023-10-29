#' @title read_methylome
#' @description read_methylome is a function to import methylation data into R Genomic
#'  Range object. Bed file should be processed in the pipeline developed by
#'  Fabian Müller and Christoph Bock. (EPP format)
#' @param  filename - filename which contains methylation data
#' @param  type - type of file format (epp, bissnp)
#' @return \code{GenomicRange} object with methylation, coverage information
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom stringr str_split str_replace
read_methylome <- function(filename, type) {
  type <- tolower(type)
  if (!file.exists(filename)) {
    stop(paste(filename, " doesn't exist or path is incorrect!"))
  }
  if (!type %in% c("bissnp", "epp")) {
    stop(paste(type, " is not a valid file type!"))
  }
  if (type == "epp") {
    # Parse EPP file format
    msites <- fread(filename, header = FALSE, showProgress = FALSE)
    mcov <- unlist(stringr::str_split(msites$V4, "/"))
    mcov <- as.numeric(stringr::str_replace(mcov, "'", ""))
    cov <- mcov[seq_along(mcov) %% 2 == 0]
    mscore <- round(msites$V5 / 1000, 3)
  }
  if (type == "bissnp") {
    # Parse BisSNP Tab-Separated file format
    msites <- fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    mscore <- round(msites$V4 / 100, 3)
    cov <- msites$V5
  }
  msites <- GenomicRanges::GRanges(
    seqnames = msites$V1,
    ranges = IRanges::IRanges(start = msites$V2, end = msites$V3),
    strand = msites$V6,
    score = mscore,
    methylation = msites$V4,
    coverage = cov
  )
  return(msites)
}
