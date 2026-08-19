#' @title run_methylTFR_RnBeads
#' @description Run the methylTFR workflow directly on a preprocessed
#' \pkg{RnBeads} object, without exporting per-sample BED files first.
#'
#' This is the RnBeads-based counterpart to \code{\link{run_methyltfr}}. Both
#' functions share the same engine and produce numerically identical results
#' for the same underlying methylation calls; they differ only in where the
#' per-sample methylation levels come from.
#'
#' @details
#' Methylation calls are always read at single-cytosine resolution
#' (\code{type = "sites"}). Region-level summaries such as \code{tiling1kb} or
#' \code{distal} cannot be used, because methylTFR needs base-resolution calls
#' to build the footprint around each motif centre. To restrict the analysis to
#' a set of regulatory regions, pass those regions through the \code{enhancer}
#' argument instead.
#'
#' Samples are processed one at a time and methylation levels are pulled from
#' the RnBeads object column by column, so disk-backed (\code{ff}-managed)
#' RnBeads sets are never loaded into memory in full.
#'
#' Coverage filtering is applied only when the object carries coverage
#' information, which is the case for sequencing-based sets
#' (\code{RnBiseqSet}). For array-based sets \code{cov_threshold} is ignored
#' and a message is emitted.
#'
#' Note that RnBeads site annotation is 1-based while
#' \code{\link{read_methylome}} reads 0-based BED coordinates as-is. The
#' resulting one-base offset is not corrected here, since deviation scores
#' aggregate methylation over windows of tens to hundreds of bases and are
#' insensitive to a uniform single-base shift.
#'
#' @param rnb_set A preprocessed \code{RnBSet} object, for example the output
#' of \code{rnb.run.preprocessing} or a set loaded with
#' \code{RnBeads::load.rnb.set}.
#' @param tf_bindsites a \code{GRangesList} object contains
#'  tf binding sites positions
#' @param gcfreqs a \code{list} of GC bin frequency tables
#'  (matrices for multiple motif)
#' @param gc_dist a \code{GRanges} object contains
#' Genome wide GC distribution
#' @param chunkSize Chunk size for parallel processing
#'  of motifs (default: 20)
#' @param threads Thread count for parallel processing
#' @param enhancer a \code{GRanges} object specifying
#' regions such as distal regulatory elements (optional)
#' @param ignoreStrand if TRUE, it ignores strand info from annotation
#' @param cov_threshold numeric, coverage threshold used to filter out low
#' coverage sites, default is 1. Ignored for objects without coverage
#' information.
#' @param sample_ann Optional \code{data.frame} of sample annotation with one
#' row per sample, used as \code{colData}. Defaults to
#' \code{RnBeads::pheno(rnb_set)}.
#' @return a \code{methylTFRdeviations} object with
#' bias-corrected deviation and Z-scores
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom logger log_info log_warn
#' @importFrom methods is
#' @seealso \code{\link{run_methyltfr}} for the file-based entry point.
#' @author Irem Gunduz
#' @examples
#' # Not run: requires the RnBeads package, an hg38 annotation package and a
#' # preprocessed RnBeads set.
#' \donttest{
#' if (requireNamespace("RnBeads", quietly = TRUE)) {
#'     # rnb_set <- RnBeads::load.rnb.set("reports/rnb.set_preprocessed")
#'     # gcfreqs <- getGCfreq(motifSet = "jaspar2020")
#'     # gc_dist <- getGenomeGC("hg38")
#'     # tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")
#'     #
#'     # deviations <- run_methylTFR_RnBeads(
#'     #     rnb_set = rnb_set,
#'     #     tf_bindsites = tf_bindsites,
#'     #     gcfreqs = gcfreqs,
#'     #     gc_dist = gc_dist,
#'     #     threads = 8,
#'     #     chunkSize = 15
#'     # )
#' }
#' }
#' @export
run_methylTFR_RnBeads <- function(
        rnb_set, tf_bindsites = NULL,
        gcfreqs = NULL, gc_dist = NULL,
        chunkSize = 20, threads = 1,
        enhancer = NULL, ignoreStrand = TRUE,
        cov_threshold = 1, sample_ann = NULL) {
    if (!requireNamespace("RnBeads", quietly = TRUE)) {
        stop(
            "The RnBeads package is required for run_methylTFR_RnBeads(). ",
            "Install it with ",
            "BiocManager::install('RnBeads'), or export your samples and use ",
            "run_methyltfr() instead."
        )
    }
    if (is.null(rnb_set) || !is(rnb_set, "RnBSet")) {
        stop("Please provide a valid RnBSet object")
    }
    if (!is.logical(ignoreStrand)) {
        warning("Found invalid strand option, using the default")
        ignoreStrand <- TRUE
    }
    if (!is.numeric(chunkSize) || chunkSize < 1) {
        warning("Invalid chunk size detected, using default chunk size")
        chunkSize <- 20
    }
    if (!is.numeric(threads) || threads < 1) {
        warning("Invalid thread count detected, using default thread count")
        threads <- 1
    }
    if (!is.numeric(cov_threshold) || cov_threshold < 0) {
        warning("Invalid cov_threshold detected, using default cov_threshold")
        cov_threshold <- 1
    }
    if (any(vapply(
        list(tf_bindsites, gcfreqs, gc_dist), is.null,
        logical(1)
    ))) {
        stop("Please load the annotation objects for given genome.")
    }
    if (!is(tf_bindsites, "GRangesList") && !is.list(tf_bindsites)) {
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

    sample_ids <- as.character(RnBeads::samples(rnb_set))
    if (length(sample_ids) == 0) {
        stop("The RnBSet object does not contain any samples")
    }

    if (is.null(sample_ann)) {
        sample_ann <- as.data.frame(RnBeads::pheno(rnb_set),
            stringsAsFactors = FALSE
        )
    }
    if (!is.data.frame(sample_ann)) {
        stop("sample_ann must be a data.frame")
    }
    if (nrow(sample_ann) != length(sample_ids)) {
        stop("sample_ann must have one row per sample in the RnBSet object")
    }

    sites_gr <- rnb_sites_to_granges(rnb_set, ignoreStrand)
    has_covg <- rnb_has_coverage(rnb_set)
    if (!has_covg) {
        log_warn(
            "The RnBSet object does not carry coverage information; ",
            "cov_threshold is ignored."
        )
    }
    log_info(
        "Found ", length(sites_gr), " sites across ",
        length(sample_ids), " samples"
    )

    msites_fun <- function(i) {
        rnb_sample_msites(
            rnb_set = rnb_set,
            sites_gr = sites_gr,
            index = i,
            cov_threshold = cov_threshold,
            has_covg = has_covg
        )
    }

    methyltfr_core(
        sample_ids = sample_ids,
        msites_fun = msites_fun,
        samples = sample_ann,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        chunkSize = chunkSize,
        threads = threads,
        enhancer = enhancer,
        ignoreStrand = ignoreStrand
    )
}


#' @title rnb_sites_to_granges
#' @description Build a \code{GRanges} object of the site annotation of an
#' RnBeads object. The order of the ranges matches the row order of the
#' methylation matrix returned by \code{RnBeads::meth}.
#' @param rnb_set An \code{RnBSet} object.
#' @param ignoreStrand if TRUE, all ranges are returned with strand \code{"*"}.
#' @return A \code{GRanges} object with one range per site.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
rnb_sites_to_granges <- function(rnb_set, ignoreStrand = TRUE) {
    ann <- RnBeads::annotation(rnb_set, type = "sites")
    if (is.null(ann) || nrow(ann) == 0) {
        stop("The RnBSet object does not contain any site annotation")
    }
    required <- c("Chromosome", "Start")
    missing_cols <- setdiff(required, colnames(ann))
    if (length(missing_cols) > 0) {
        stop(
            "Unexpected RnBeads site annotation, missing column(s): ",
            paste(missing_cols, collapse = ", ")
        )
    }
    ends <- if ("End" %in% colnames(ann)) ann$End else ann$Start
    strands <- "*"
    if (!ignoreStrand && "Strand" %in% colnames(ann)) {
        strands <- as.character(ann$Strand)
        strands[is.na(strands) | !strands %in% c("+", "-")] <- "*"
    }
    GenomicRanges::GRanges(
        seqnames = as.character(ann$Chromosome),
        ranges = IRanges::IRanges(
            start = as.integer(ann$Start),
            end = as.integer(ends)
        ),
        strand = strands
    )
}


#' @title rnb_has_coverage
#' @description Test whether an RnBeads object carries coverage information.
#' @param rnb_set An \code{RnBSet} object.
#' @return A logical scalar.
#' @keywords internal
rnb_has_coverage <- function(rnb_set) {
    # Try the subsetting form first so that large disk-backed sets are not
    # materialised, then fall back for RnBeads versions without the j argument.
    res <- tryCatch(
        !is.null(RnBeads::covg(rnb_set, type = "sites", j = 1L)),
        error = function(e) NULL
    )
    if (is.null(res)) {
        res <- tryCatch(
            !is.null(RnBeads::covg(rnb_set, type = "sites")),
            error = function(e) FALSE
        )
    }
    return(isTRUE(res))
}


#' @title rnb_sample_msites
#' @description Extract the methylation calls of a single sample from an
#' RnBeads object as a \code{GRanges} object in the layout expected by
#' \code{\link{computeDeviation}}.
#' @param rnb_set An \code{RnBSet} object.
#' @param sites_gr A \code{GRanges} object of site positions, as returned by
#' \code{rnb_sites_to_granges}.
#' @param index Integer index of the sample to extract.
#' @param cov_threshold numeric coverage threshold.
#' @param has_covg logical, whether the object carries coverage information.
#' @return A \code{GRanges} object with \code{score} and \code{coverage}
#' metadata columns, restricted to sites with a non-missing methylation call.
#' @importFrom logger log_warn
#' @keywords internal
rnb_sample_msites <- function(rnb_set, sites_gr, index,
    cov_threshold = 1, has_covg = TRUE) {
    index <- as.integer(index)
    mvals <- rnb_column(RnBeads::meth, rnb_set, index)
    if (length(mvals) != length(sites_gr)) {
        stop(
            "The number of methylation values does not match the number of ",
            "annotated sites; the RnBSet object appears to be inconsistent."
        )
    }
    keep <- !is.na(mvals)
    if (has_covg) {
        cvals <- rnb_column(RnBeads::covg, rnb_set, index)
        if (length(cvals) != length(sites_gr)) {
            stop(
                "The number of coverage values does not match the number of ",
                "annotated sites."
            )
        }
        keep <- keep & !is.na(cvals) & cvals >= cov_threshold
    } else {
        cvals <- rep(NA_real_, length(sites_gr))
    }
    if (!any(keep)) {
        stop(
            "No sites passed the coverage threshold for sample index ", index
        )
    }
    gr <- sites_gr[keep]
    gr$score <- as.numeric(mvals[keep])
    gr$coverage <- as.numeric(cvals[keep])
    return(gr)
}


#' @title rnb_column
#' @description Extract a single sample column from an RnBeads accessor,
#' falling back to full extraction on RnBeads versions that do not support
#' column subsetting.
#' @param accessor An RnBeads accessor function, either \code{RnBeads::meth} or
#' \code{RnBeads::covg}.
#' @param rnb_set An \code{RnBSet} object.
#' @param index Integer index of the sample to extract.
#' @return A numeric vector with one value per site.
#' @keywords internal
rnb_column <- function(accessor, rnb_set, index) {
    vals <- tryCatch(
        accessor(rnb_set, type = "sites", j = index),
        error = function(e) NULL
    )
    if (is.null(vals)) {
        vals <- accessor(rnb_set, type = "sites")
        if (is.matrix(vals) || is.data.frame(vals)) {
            vals <- vals[, index]
        }
    }
    return(as.numeric(as.vector(vals)))
}
