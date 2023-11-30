#' @title read_methylome
#' @description read_methylome is a function to import methylation data into R Genomic
#'  Range object. Bed file should be in EPP,ALLC or BisSNP format.
#' @param  filename - filename which contains methylation data
#' @param  type Type of file format. Currently supported epp, bissnp,bismarkCytosine,bismarkcov,allc and encode
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom stringr str_split str_replace
#' @return a \code{GenomicRange} object with methylation, coverage information
#' @author Irem Gunduz
#' @examples
#' # Read EPP file
#' epp_path <- system.file("extdata", "epp.tsv.gz", package = "methylTFR")
#' epp <- read_methylome(epp_path, "EPP")
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
    mscore <- round(msites$V5 / 1000, 3)
    msites <- msites[, .(V1, V2, V3, V6)]
    # colnames(msites) <- c("chr", "start", "end", "strand")
  }
  if (type == "bissnp") {
    # Parse BisSNP Tab-Separated file format
    msites <- data.table::fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    mscore <- round(msites$V4 / 100, 3)
    cov <- msites$V5
    msites <- msites[, .(V1, V2, V3, V6)]
    # colnames(msites) <- c("chr", "start", "end", "strand")
  }
  if (type == "allc") {
    # Parse ALLC Tab-Separated file format
    allc <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(allc) < 6) {
      logger::log_warn(paste0(filename, " is an invalid allc file!"))
      stop("allc file must contain at least 6 columns!")
    }
    cov <- allc$V6
    mscore <- round(allc$V5 / cov, 3)
    msites <- allc[, .(V1, V2, V2, V3)]
    # colnames(msites) <- c("chr", "start", "end", "strand") # , "meth")
  }
  if (type == "bismarkcytosine") {
    # Parse bismarkCytosine file format
    msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(msites) < 8) {
      logger::log_warn(paste0(filename, " is an invalid bismarkCytosine file!"))
      stop("bismarkCytosine file must contain at least 8 columns!")
    }
    mscore <- round(msites$V5 / msites$V6, 3)
    cov <- msites$V6
    msites <- msites[, .(V1, V3, V3, V4)]
  }
  if (type == "bismarkcov") {
    # Parse bismarkCov file format
    msites <- data.table::fread(filename, header = FALSE, showProgress = FALSE)
    if (ncol(msites) < 6) {
      logger::log_warn(paste0(filename, " is an invalid bismarkCov file!"))
      stop("bismarkCov file must contain at least 6 columns!")
    }
    mscore <- round(msites$V4 / 100, 3)
    cov <- msites$V5 + msites$V6
    msites <- msites[, .(V1, V2, V3)]
    msites$strand <- rep("*", nrow(msites))
  }
  if (type == "encode") {
    # Parse encode file format
    msites <- data.table::fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
    if (ncol(msites) < 11) {
      logger::log_warn(paste0(filename, " is an invalid encode file!"))
      stop("encode file must contain 11 columns!")
    }
    mscore <- round(msites$V11 / msites$V10, 3)
    cov <- msites$V10
    msites <- msites[, .(V1, V2, V3, V6)]
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

#' @title EPP_helper
#' @description  helper function to convert msites in EPP GRanges format
#' @param  granges - GRanges object
#' @return a \code{GenomicRanges} object with EPP format
#' @keywords internal
EPP_helper <- function(granges) {
  tiles <- data.table::as.data.table(granges)
  tiles[, score2 := round((round(tiles$score * tiles$coverage, 0) / tiles$coverage) * 1000, 0)]
  tiles[, score := paste0(round(tiles$score * tiles$coverage, 0), "/", tiles$coverage)]
  tiles <- tiles[, .(seqnames, start, end, score, score2, strand)]
  return(tiles)
}

#' @title convertToEPP
#' @description  convertToEPP is a function to convert different objects to EPP format
#' and save as bed file
#' @param obj - object to be converted. Currently support GRanges, GRangesList,
#' CompressedGRangesList, data.frame, methylRawList, methylRawListDB, methylRaw objects
#' @param save - logical, whether to save the converted object as bed file
#' @param filePath - character, path to save the converted object
#' @param threads - integer, number of threads to use
#' @param verbose - logical, whether to print out the progress
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom data.table fwrite as.data.table
#' @author Irem Gunduz
#' @export
#' @return a list of \code{GenomicRanges} object with EPP format or if save TRUE, invisible NULL
convertToEPP <- function(obj, save = FALSE, filePath = NULL, threads = 1, verbose = TRUE) {
  if (!class(obj) %in% c(
    "data.frame", "GRanges",
    "GRangesList", "CompressedGRangesList", "methylRawList", "methylRawListDB"
  )) {
    stop(paste("Object class", class(obj), "is not supported!"))
  }
  if (!is.logical(save)) {
    stop("save must be logical, please set save=TRUE or save=FALSE")
  }
  if (save & c(is.null(filePath) || !is.character(filePath))) {
    stop("filePath must be specified chracter if save=TRUE")
  }
  if (!is.numeric(threads)) {
    stop("threads must be numeric")
  }
  if (!is.logical(verbose)) {
    verbose <- FALSE
  }
  if (is.data.frame(obj)) {
    if (any(colnames(obj) %in% c("chr", "start", "end", "strand"))) {
      obj <- GenomicRanges::makeGRangesFromDataFrame(obj, keep.extra.columns = TRUE)
    } else {
      stop("data.frame must contain columns chr,start,end,strand,score,coverage,methylation")
    }
  }
  if (is(obj, "GRanges")) {
    tiles <- EPP_helper(obj)
  }
  if (class(obj) %in% c("GRangesList", "CompressedGRangesList")) {
    tiles <- parallel::mclapply(obj, function(x) {
      EPP_helper(x)
    }, mc.cores = threads)
  }
  if (class(obj) %in% c("methylRawList", "methylRawListDB")) {
    obj <- methylKitToEPP(obj)
    tiles <- parallel::mclapply(obj, function(x) {
      EPP_helper(x)
    }, mc.cores = threads)
  }
  if (!save) {
    return(tiles)
  } else {
    for (sample in names(tiles)) {
      if (!file.exists(paste0(filePath, "/", sample, ".bed"))) {
        if (verbose) {
          logger::log_info(paste0("Writing EPP file for ", sample))
        }
        data.table::fwrite(data.table::as.data.table(tiles[[sample]]),
          paste0(filePath, "/", sample, ".bed"), TRUE, FALSE,
          sep = "\t", row.names = FALSE, col.names = FALSE, showProgress = FALSE
        )
        if (verbose) {
          logger::log_success(paste0("EPP file for ", sample, " is written!"))
        }
      } else {
        logger::log_warn(paste0("File already exists in filePath for sample ", sample, "!"))
      }
    }
    return(invisible(NULL))
  }
}

#' @title convert methylKit object to EPP format
#' @description convert methylKit object to EPP format
#' @param  mKit - methylKit object
#' @return a \code{GenomicRanges} or  object with EPP format
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom methylKit getSampleID
#' @keywords internal
methylKitToEPP <- function(mKit) {
  if (is(mKit, "methylRaw")) {
    mscore <- round(x$numCs / x$coverage, 2)
    meth <- paste0(x$numCs, "/", x$coverage)
    gr_obj <- granges_helper(
      grobj = mKit,
      chr = mKit$chr,
      mscore = mscore,
      cov = mKit$coverage,
      startP = mKit$start,
      endP = mKit$end
    )
  }
  if (is(mKit, "methylRawList") || is(mKit, "methylRawListDB")) {
    gr_obj <- lapply(mKit, function(x) {
      mscore <- round(x$numCs / x$coverage, 2)
      meth <- paste0(x$numCs, "/", x$coverage)
      return(granges_helper(
        grobj = x,
        chr = x$chr,
        mscore = mscore,
        cov = x$coverage,
        startP = x$start,
        endP = x$end
      ))
    })
    gr_obj <- GRangesList(gr_obj)
    names(gr_obj) <- getSampleID(mKit)
  }
  return(gr_obj)
}
