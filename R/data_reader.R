#' @title read_methylome
#' @description read_methylome is a function to import methylation data into
#' a GRanges object. Bed file should be in EPP, ALLC or BisSNP format.
#' @param filename filename which contains methylation data
#' @param type Type of file format. Currently supported epp, bissnp,
#' bismarkCytosine, bismarkcov, allc and encode
#' @param cov_threshold numeric, coverage threshold to filter out low coverage
#' sites, default is 1
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
#' bissnp_path <- system.file(
#'   "extdata",
#'   "bissnp.tsv.gz",
#'   package = "methylTFR"
#' )
#' bissnp <- read_methylome(bissnp_path, "bissnp")
#' @export
read_methylome <- function(filename, type, cov_threshold = 1) {
  type <- tolower(type)
  if (!file.exists(filename)) {
    stop(filename, " doesn't exist or path is incorrect!")
  }
  valid_types <- c(
    "bissnp", "epp", "allc", "bismarkcytosine", "bismarkcov", "encode"
  )
  if (!type %in% valid_types) {
    stop(type, " is not a valid file type!")
  }
  if (!is.numeric(cov_threshold) || cov_threshold < 0) {
    stop(cov_threshold, " is not a valid coverage threshold!")
  }

  parsed_data <- switch(type,
    "epp" = parse_epp(filename),
    "bissnp" = parse_bissnp(filename),
    "allc" = parse_allc(filename),
    "bismarkcytosine" = parse_bismarkcytosine(filename),
    "bismarkcov" = parse_bismarkcov(filename),
    "encode" = parse_encode(filename)
  )

  msites <- parsed_data$msites
  colnames(msites) <- c("chr", "start", "end", "strand")

  gr_obj <- granges_helper(
    grobj = msites,
    chr = msites$chr,
    mscore = parsed_data$mscore,
    cov = parsed_data$cov,
    startP = msites$start,
    endP = msites$end
  )

  # Filter low coverage sites
  gr_obj <- gr_obj[gr_obj$coverage >= cov_threshold]
  return(gr_obj)
}

#' @keywords internal
parse_epp <- function(filename) {
  msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
  mcov <- unlist(stringr::str_split(msites$V4, "/"))
  mcov <- as.numeric(stringr::str_replace(mcov, "'", ""))
  cov <- mcov[seq_along(mcov) %% 2 == 0]
  mscore <- round(msites$V5 / 1000, 6)
  msites <- msites[, .(V1, V2, V3, V6)]
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @keywords internal
parse_bissnp <- function(filename) {
  msites <- data.table::fread(
    filename,
    header = FALSE, skip = 1, showProgress = FALSE
  )
  mscore <- round(msites$V4 / 100, 6)
  cov <- msites$V5
  msites <- msites[, .(V1, V2, V3, V6)]
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @keywords internal
parse_allc <- function(filename) {
  allc <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
  if (ncol(allc) < 6) {
    logger::log_warn(sprintf("%s is an invalid allc file!", filename))
    stop("allc file must contain at least 6 columns!")
  }
  cov <- allc$V6
  mscore <- round(allc$V5 / cov, 6)
  msites <- allc[, .(V1, V2, V2, V3)]
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @keywords internal
parse_bismarkcytosine <- function(filename) {
  msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
  if (ncol(msites) < 5) {
    logger::log_warn(
      sprintf("%s is an invalid bismarkCytosine file!", filename)
    )
    stop("bismarkCytosine file must contain at least 5 columns!")
  }
  cov <- (msites$V5 + msites$V4)
  mscore <- round(msites$V4 / cov, 6)
  msites <- msites[, .(V1, V2, V2, V3)]
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @keywords internal
parse_bismarkcov <- function(filename) {
  msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
  if (ncol(msites) < 6) {
    logger::log_warn(sprintf("%s is an invalid bismarkCov file!", filename))
    stop("bismarkCov file must contain at least 6 columns!")
  }
  mscore <- round(msites$V4 / 100, 6)
  cov <- msites$V5 + msites$V6
  msites <- msites[, .(V1, V2, V3)]
  msites$strand <- rep("*", nrow(msites))
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @keywords internal
parse_encode <- function(filename) {
  msites <- data.table::fread(
    filename,
    header = FALSE, skip = 1, showProgress = FALSE
  )
  if (ncol(msites) < 11) {
    logger::log_warn(sprintf("%s is an invalid encode file!", filename))
    stop("encode file must contain 11 columns!")
  }
  mscore <- round(msites$V11 / msites$V10, 6)
  cov <- msites$V10
  msites <- msites[, .(V1, V2, V3, V6)]
  list(msites = msites, mscore = mscore, cov = cov)
}

#' @title granges_helper
#' @description helper function to convert msites in EPP GRanges format
#' @param grobj GRanges object
#' @param chr chromosome
#' @param mscore methylation score
#' @param cov methylation coverage
#' @param startP start position
#' @param endP end position
#' @return a \code{GenomicRanges} object with EPP format
#' @keywords internal
granges_helper <- function(grobj, chr, mscore, cov, startP, endP) {
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
