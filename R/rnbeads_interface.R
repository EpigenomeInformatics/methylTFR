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
#' @details
#' Methylation calls are read at single-cytosine resolution. Sequencing
#' sets (\code{RnBiseqSet}) are filtered by coverage, array sets
#' (\code{RnBeadSet}) by detection p-value. The number of sites retained
#' per sample is reported through \pkg{logger}.
#'
#' @param rnb_set A preprocessed \code{RnBSet} object.
#' @param tf_bindsites a \code{GRangesList} of TF binding site positions.
#' @param gcfreqs a \code{list} of GC bin frequency tables.
#' @param gc_dist a \code{GRanges} of the genome-wide GC distribution.
#' @param chunkSize Chunk size for parallel processing of motifs.
#' @param threads Thread count for parallel processing.
#' @param enhancer an optional \code{GRanges} of regions to restrict to.
#' @param ignoreStrand if TRUE, strand information is ignored.
#' @param cov_threshold numeric, minimum coverage of a retained site.
#' @param dpval_threshold numeric, maximum detection p-value of a
#' retained probe.
#' @param sample_ann Optional \code{data.frame} of sample annotation.
#' @return a \code{methylTFRdeviations} object with bias-corrected
#' deviations and Z-scores.
#' @importFrom logger log_info
#' @importFrom methods is
#' @export
run_methylTFR_RnBeads <- function(
    rnb_set, tf_bindsites = NULL, gcfreqs = NULL, gc_dist = NULL,
    chunkSize = 20, threads = 1, enhancer = NULL, ignoreStrand = TRUE,
    cov_threshold = 1, dpval_threshold = 0.05, sample_ann = NULL
) {
    check_rnb_inputs(rnb_set)
    check_annotation_inputs(tf_bindsites, gcfreqs, gc_dist, enhancer)
    opts <- check_run_options(chunkSize, threads, ignoreStrand, cov_threshold)

    sample_ids <- rnb_sample_ids(rnb_set)
    sample_ann <- resolve_rnb_sample_ann(rnb_set, sample_ann, sample_ids)
    sites_gr <- rnb_sites_to_granges(rnb_set, opts$ignoreStrand)
    mode <- rnb_quality_mode(rnb_set, dpval_threshold)
    log_info(
        "Found ", length(sites_gr), " sites across ",
        length(sample_ids), " samples"
    )

    msites_fun <- function(i) {
        rnb_sample_msites(
            rnb_set = rnb_set, sites_gr = sites_gr, index = i,
            cov_threshold = opts$cov_threshold, has_covg = mode$has_covg,
            dpval_threshold = dpval_threshold, has_dpval = mode$has_dpval
        )
    }

    methyltfr_core(
        sample_ids = sample_ids, msites_fun = msites_fun,
        samples = sample_ann, tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs, gc_dist = gc_dist, chunkSize = opts$chunkSize,
        threads = opts$threads, enhancer = enhancer,
        ignoreStrand = opts$ignoreStrand
    )
}


#' @title rnb_sample_ids
#' @description Determine the sample identifiers of an RnBeads object.
#' @param rnb_set An \code{RnBSet} object.
#' @return A character vector of sample identifiers.
#' @keywords internal
rnb_sample_ids <- function(rnb_set) {
    ids <- tryCatch(
        colnames(RnBeads::meth(rnb_set, type = "sites", i = 1L)),
        error = function(e) NULL
    )
    if (length(ids) == 0) {
        ids <- tryCatch(
            rownames(RnBeads::pheno(rnb_set)),
            error = function(e) NULL
        )
    }
    if (length(ids) == 0) {
        nsamples <- tryCatch(
            nrow(RnBeads::pheno(rnb_set)),
            error = function(e) 0L
        )
        if (length(nsamples) == 1 && nsamples > 0) {
            ids <- paste0("sample_", seq_len(nsamples))
        }
    }
    if (length(ids) == 0) {
        stop("Could not determine sample identifiers from the RnBSet object.")
    }
    return(as.character(ids))
}


#' @title rnb_annotation_target
#' @description Resolve the annotation target of an RnBeads object.
#' @param rnb_set An \code{RnBSet} object.
#' @return A character scalar naming the annotation target, \code{"sites"}
#' for sequencing sets and the array platform for array sets.
#' @importFrom methods is
#' @keywords internal
rnb_annotation_target <- function(rnb_set) {
    if (is(rnb_set, "RnBeadSet")) rnb_set@target else "sites"
}


#' @title rnb_annotation_table
#' @description Look up the site or probe annotation of an RnBeads object.
#' @details The annotation stored in the object is used when available,
#' otherwise the genome-wide track registered for \code{target} is
#' subset to the sites of the object.
#' @param rnb_set An \code{RnBSet} object.
#' @param target Character scalar naming the annotation target.
#' @param assembly Character scalar naming the genome assembly.
#' @return A \code{data.frame} with one row per site or probe.
#' @keywords internal
rnb_annotation_table <- function(rnb_set, target, assembly) {
    pull <- function(expr) tryCatch(expr, error = function(e) NULL)

    ann <- pull(RnBeads::annotation(rnb_set, type = target))
    if (is.null(ann)) {
        ann <- pull(RnBeads::annotation(rnb_set, type = "sites"))
    }
    if (is.null(ann) && target != "sites") {
        ann <- pull({
            track <- RnBeads::rnb.get.annotation(target, assembly)
            RnBeads::rnb.annotation2data.frame(track)[rnb_set@sites, ]
        })
    }
    if (is.null(ann) || nrow(ann) == 0) {
        stop(
            "Could not extract coordinates for target '", target,
            "' and assembly '", assembly, "'. Load the matching RnBeads ",
            "annotation, or a custom annotation for this array, before ",
            "calling methylTFR."
        )
    }
    return(ann)
}


#' @title rnb_sites_to_granges
#' @description Build a \code{GRanges} object of the site or probe
#' annotation of an RnBeads object.
#' @param rnb_set An \code{RnBSet} object.
#' @param ignoreStrand if TRUE, all ranges are returned with strand
#' \code{"*"}.
#' @return A \code{GRanges} object with one range per site or probe.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom logger log_info
#' @importFrom methods is
#' @keywords internal
rnb_sites_to_granges <- function(rnb_set, ignoreStrand = TRUE) {
    assembly <- rnb_set@assembly
    target <- rnb_annotation_target(rnb_set)
    log_info("Annotation target: ", target, " | assembly: ", assembly)
    invisible(requireNamespace(paste0("RnBeads.", assembly), quietly = TRUE))

    ann <- rnb_annotation_table(rnb_set, target, assembly)
    missing_cols <- setdiff(c("Chromosome", "Start"), colnames(ann))
    if (length(missing_cols) > 0) {
        stop("RnBeads annotation misses column(s): ", toString(missing_cols))
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
            start = as.integer(ann$Start), end = as.integer(ends)
        ),
        strand = strands
    )
}


#' @title rnb_has_coverage
#' @description Test whether an RnBeads object carries coverage
#' information.
#' @param rnb_set An \code{RnBSet} object.
#' @return A logical scalar.
#' @keywords internal
rnb_has_coverage <- function(rnb_set) {
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


#' @title rnb_has_dpval
#' @description Test whether an RnBeads object carries detection
#' p-values.
#' @param rnb_set An \code{RnBSet} object.
#' @return A logical scalar.
#' @keywords internal
rnb_has_dpval <- function(rnb_set) {
    res <- tryCatch(
        !is.null(RnBeads::dpval(rnb_set, type = "sites", j = 1L)),
        error = function(e) NULL
    )
    if (is.null(res)) {
        res <- tryCatch(
            !is.null(RnBeads::dpval(rnb_set, type = "sites")),
            error = function(e) FALSE
        )
    }
    return(isTRUE(res))
}


#' @title rnb_quality_mode
#' @description Decide which per-site quality filter applies to an
#' RnBeads object and report it.
#' @param rnb_set An \code{RnBSet} object.
#' @param dpval_threshold numeric detection p-value threshold.
#' @return A list with the logical flags \code{has_covg} and
#' \code{has_dpval}.
#' @importFrom logger log_info log_warn
#' @keywords internal
rnb_quality_mode <- function(rnb_set, dpval_threshold) {
    has_covg <- rnb_has_coverage(rnb_set)
    has_dpval <- rnb_has_dpval(rnb_set)
    if (has_dpval) {
        log_info("Filtering probes at detection p-value <= ", dpval_threshold)
    } else if (!has_covg) {
        log_warn(
            "RnBSet carries neither coverage nor detection p-values, ",
            "no quality filtering is applied"
        )
    }
    list(has_covg = has_covg, has_dpval = has_dpval)
}


#' @title rnb_column
#' @description Extract a single sample column from an RnBeads accessor.
#' @param accessor An RnBeads accessor function.
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


#' @title rnb_sample_msites
#' @description Extract the methylation calls of a single sample from an
#' RnBeads object as a \code{GRanges} object.
#' @param rnb_set An \code{RnBSet} object.
#' @param sites_gr A \code{GRanges} object of site positions.
#' @param index Integer index of the sample to extract.
#' @param cov_threshold numeric coverage threshold.
#' @param has_covg logical, whether coverage filtering applies.
#' @param dpval_threshold numeric detection p-value threshold.
#' @param has_dpval logical, whether detection p-value filtering applies.
#' @return A \code{GRanges} object restricted to valid methylation calls.
#' @importFrom logger log_info
#' @keywords internal
rnb_sample_msites <- function(
    rnb_set, sites_gr, index, cov_threshold = 1, has_covg = TRUE,
    dpval_threshold = 0.05, has_dpval = FALSE
) {
    index <- as.integer(index)
    mvals <- rnb_column(RnBeads::meth, rnb_set, index)
    if (length(mvals) != length(sites_gr)) {
        stop("Methylation values do not match annotated sites.")
    }

    keep <- !is.na(mvals)
    if (has_covg) {
        cvals <- rnb_column(RnBeads::covg, rnb_set, index)
        keep <- keep & !is.na(cvals) & cvals >= cov_threshold
    } else {
        cvals <- rep(NA_real_, length(sites_gr))
    }
    if (has_dpval) {
        dvals <- rnb_column(RnBeads::dpval, rnb_set, index)
        keep <- keep & !is.na(dvals) & dvals <= dpval_threshold
    }

    log_info(
        "Sample ", index, ": ", sum(keep), " of ", length(keep),
        " sites retained (", round(100 * mean(keep), 2), "%)"
    )
    if (!any(keep)) {
        stop("No sites passed the quality thresholds for sample index ", index)
    }

    gr <- sites_gr[keep]
    gr$score <- as.numeric(mvals[keep])
    gr$coverage <- as.numeric(cvals[keep])
    return(gr)
}
