#' @title check_rnb_inputs
#' @description Check that \pkg{RnBeads} is available and that the object
#' handed in is an \code{RnBSet}.
#' @param rnb_set The object to validate.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
check_rnb_inputs <- function(rnb_set) {
    if (!requireNamespace("RnBeads", quietly = TRUE)) {
        stop(
            "The RnBeads package is required for run_methylTFR_RnBeads(). ",
            "Install it with BiocManager::install('RnBeads'), or export ",
            "your samples and use run_methyltfr() instead."
        )
    }
    if (is.null(rnb_set) || !is(rnb_set, "RnBSet")) {
        stop("Please provide a valid RnBSet object")
    }
    invisible(NULL)
}

#' @title resolve_rnb_sample_ann
#' @description Default the sample annotation to the phenotype table of the
#' RnBeads object and check its shape.
#' @param rnb_set A preprocessed \code{RnBSet} object.
#' @param sample_ann Sample annotation supplied by the caller, or NULL.
#' @param sample_ids Character vector of sample identifiers.
#' @return The sample annotation as a \code{data.frame}.
#' @keywords internal
resolve_rnb_sample_ann <- function(rnb_set, sample_ann, sample_ids) {
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
    return(sample_ann)
}

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
#' # A minimal end-to-end run on the BATF example data bundled with the
#' # package. RnBeads and its hg38 annotation build the input object; both
#' # are optional dependencies.
#' if (requireNamespace("RnBeads", quietly = TRUE) &&
#'     requireNamespace("RnBeads.hg38", quietly = TRUE)) {
#'     load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#'     load(system.file(
#'         "extdata", "BATF_tf_bindsites.rda",
#'         package = "methylTFR"
#'     ))
#'     load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#'     load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#'     # RnBiseqSet() takes methylation as a fraction and coverage as counts,
#'     # with one column per sample.
#'     sites <- data.frame(
#'         chromosome = as.character(GenomicRanges::seqnames(msites)),
#'         position = GenomicRanges::start(msites),
#'         strand = "*",
#'         stringsAsFactors = FALSE
#'     )
#'     rnb_set <- RnBeads::RnBiseqSet(
#'         pheno = data.frame(
#'             sampleName = "sample_1", stringsAsFactors = FALSE
#'         ),
#'         sites = sites,
#'         meth = matrix(msites$score, ncol = 1),
#'         covg = matrix(msites$coverage, ncol = 1),
#'         assembly = "hg38",
#'         summarize.regions = FALSE
#'     )
#'
#'     devs <- run_methylTFR_RnBeads(
#'         rnb_set = rnb_set,
#'         tf_bindsites = tf_bindsites,
#'         gcfreqs = gcfreqs,
#'         gc_dist = gcdist
#'     )
#'     deviations(devs)
#' }
#' @export
run_methylTFR_RnBeads <- function(
    rnb_set, tf_bindsites = NULL, gcfreqs = NULL, gc_dist = NULL,
    chunkSize = 20, threads = 1, enhancer = NULL, ignoreStrand = TRUE,
    cov_threshold = 1, sample_ann = NULL
) {
    check_rnb_inputs(rnb_set)
    check_annotation_inputs(tf_bindsites, gcfreqs, gc_dist, enhancer)
    opts <- check_run_options(chunkSize, threads, ignoreStrand, cov_threshold)

    sample_ids <- rnb_sample_ids(rnb_set)

    sample_ann <- resolve_rnb_sample_ann(rnb_set, sample_ann, sample_ids)

    sites_gr <- rnb_sites_to_granges(rnb_set, opts$ignoreStrand)
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
            cov_threshold = opts$cov_threshold,
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
        chunkSize = opts$chunkSize,
        threads = opts$threads,
        enhancer = enhancer,
        ignoreStrand = opts$ignoreStrand
    )
}


#' @title rnb_sample_ids
#' @description Determine the sample identifiers of an RnBeads object.
#' @details RnBeads exports \code{samples} with \code{exportMethods} rather
#' than \code{export}, so the generic is not reachable as
#' \code{RnBeads::samples} and referring to it that way fails
#' \code{R CMD check}. The identifiers are therefore taken from the column
#' names of the methylation matrix, falling back to the row names of the
#' sample annotation.
#' @param rnb_set An \code{RnBSet} object.
#' @return A character vector of sample identifiers.
#' @keywords internal
rnb_sample_ids <- function(rnb_set) {
    ids <- tryCatch(
        colnames(RnBeads::meth(rnb_set, type = "sites", i = 1L)),
        error = function(e) NULL
    )
    if (is.null(ids) || length(ids) == 0) {
        ids <- tryCatch(rownames(RnBeads::pheno(rnb_set)),
            error = function(e) NULL
        )
    }
    if (is.null(ids) || length(ids) == 0) {
        nsamples <- tryCatch(nrow(RnBeads::pheno(rnb_set)),
            error = function(e) 0L
        )
        if (!is.null(nsamples) && nsamples > 0) {
            ids <- paste0("sample_", seq_len(nsamples))
        }
    }
    if (is.null(ids) || length(ids) == 0) {
        stop(
            "Could not determine sample identifiers from the RnBSet object; ",
            "pass them explicitly via sample_ann."
        )
    }
    return(as.character(ids))
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
rnb_sample_msites <- function(
    rnb_set, sites_gr, index,
    cov_threshold = 1, has_covg = TRUE
) {
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
