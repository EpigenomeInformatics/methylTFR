#' @title readMethylome
#' @description Import single-sample methylation calls into a
#' \code{GRanges} object.
#' @param filename Path to the methylation call file. Gzipped files are
#' read directly.
#' @param type Format of the file, one of \code{"epp"}, \code{"bissnp"},
#' \code{"bismarkCytosine"}, \code{"bismarkCov"}, \code{"allc"} or
#' \code{"encode"} (ENCODE bedMethyl). Case-insensitive.
#' @param cov_threshold Numeric. Sites with fewer reads than this are
#' dropped. Default 1.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import data.table
#' @import R.utils
#' @importFrom data.table := .N fread
#' @importFrom stringr str_split str_replace
#' @return A \code{GRanges} object with one range per CpG and the metadata
#' columns \code{score} (methylation fraction between 0 and 1) and
#' \code{coverage} (read count). This is the \code{msites} input of
#' \code{\link{computeDeviation}}, \code{\link{addGCBinToMethylome}} and
#' the footprint plots.
#' @seealso \code{\link{runMethylTFR}}, which calls this function for every
#' sample.
#' @author Irem Gunduz
#' @examples
#' # Read bissnp file
#' bissnp_path <- system.file(
#'     "extdata",
#'     "bissnp.tsv.gz",
#'     package = "methylTFR"
#' )
#' bissnp <- readMethylome(bissnp_path, "bissnp")
#' @export
readMethylome <- function(filename, type, cov_threshold = 1) {
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
    gr_obj <- gr_obj[!is.na(gr_obj$coverage) &
        gr_obj$coverage >= cov_threshold]

    # Scores must be fractions
    if (length(gr_obj) > 0) {
        rng <- range(gr_obj$score, na.rm = TRUE)
        if (is.finite(rng[1]) && (rng[1] < 0 || rng[2] > 1)) {
            stop(
                "Methylation scores parsed from ", filename, " span ",
                signif(rng[1], 3), " to ", signif(rng[2], 3),
                ", which is not a fraction. Check that '", type,
                "' is the right format for this file."
            )
        }
    }
    return(gr_obj)
}

#' @title Parse an EPP file
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
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

#' @title Parse a BisSNP file
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
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

#' @title Parse an ALLC file
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
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

#' @title Parse a Bismark cytosine report
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
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

#' @title Parse a Bismark coverage file
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
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

#' @title Parse an ENCODE bedMethyl file
#' @param filename Path to the file.
#' @return A list with the elements msites, mscore and cov.
#' @keywords internal
parse_encode <- function(filename) {
    # bedMethyl: column 10 is coverage, column 11 is percent methylated
    msites <- data.table::fread(
        filename,
        header = "auto", showProgress = FALSE
    )
    if (ncol(msites) < 11) {
        logger::log_warn(sprintf("%s is an invalid encode file!", filename))
        stop("encode file must contain 11 columns!")
    }
    cov <- as.numeric(msites[[10]])
    pct <- as.numeric(msites[[11]])
    if (any(pct > 100, na.rm = TRUE) || any(pct < 0, na.rm = TRUE)) {
        stop(
            "Column 11 of ", filename, " is outside 0-100, so it is not a ",
            "methylation percentage. Check that this is a bedMethyl file."
        )
    }
    mscore <- round(pct / 100, 6)
    msites <- data.table::data.table(
        V1 = msites[[1]], V2 = msites[[2]], V3 = msites[[3]], V6 = msites[[6]]
    )
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
