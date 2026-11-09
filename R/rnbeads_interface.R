#' @title checkRnbInputs
#' @description Check that \pkg{RnBeads} is available and that the object
#' handed in is an \code{RnBSet}.
#' @param rnbSet The object to validate.
#' @return Invisible \code{NULL}. Called for the errors it raises.
#' @importFrom methods is
#' @keywords internal
checkRnbInputs <- function(rnbSet) {
    if (!requireNamespace("RnBeads", quietly = TRUE)) {
        stop(
            "The RnBeads package is required for runMethylTFRFromRnBeads(). ",
            "Install it with BiocManager::install('RnBeads'), or export ",
            "your samples and use runMethylTFR() instead."
        )
    }
    if (is.null(rnbSet) || !is(rnbSet, "RnBSet")) {
        stop("Please provide a valid RnBSet object")
    }
    invisible(NULL)
}


#' @title resolveRnbSampleAnn
#' @description Default the sample annotation to the phenotype table of the
#' RnBeads object and check its shape.
#' @param rnbSet A preprocessed \code{RnBSet} object.
#' @param sampleAnn Sample annotation supplied by the caller, or NULL.
#' @param sampleIds Character vector of sample identifiers.
#' @return The sample annotation as a \code{data.frame}.
#' @keywords internal
resolveRnbSampleAnn <- function(rnbSet, sampleAnn, sampleIds) {
    if (is.null(sampleAnn)) {
        sampleAnn <- as.data.frame(RnBeads::pheno(rnbSet),
            stringsAsFactors = FALSE
        )
    }
    if (!is.data.frame(sampleAnn)) {
        stop("sampleAnn must be a data.frame")
    }
    if (nrow(sampleAnn) != length(sampleIds)) {
        stop("sampleAnn must have one row per sample in the RnBSet object")
    }
    return(sampleAnn)
}


#' @title runMethylTFRFromRnBeads
#' @description Run the methylTFR workflow directly on a preprocessed
#' \pkg{RnBeads} object, without exporting per-sample BED files first.
#'
#' @details
#' Methylation calls are read at single-cytosine resolution. Sequencing
#' sets (\code{RnBiseqSet}) are filtered by coverage, array sets
#' (\code{RnBeadSet}) by detection p-value. The number of sites retained
#' per sample is reported through \pkg{logger}.
#'
#' @inheritParams runMethylTFR
#' @param rnbSet A preprocessed \code{RnBSet} object, for example loaded
#' with \code{RnBeads::load.rnb.set()}.
#' @param covThreshold Numeric, minimum coverage of a retained site. Used
#' only for sequencing sets. Default 1.
#' @param dpvalThreshold Numeric, maximum detection p-value of a retained
#' probe. Used only for array sets. Default 0.05.
#' @param sampleAnn Optional \code{data.frame} with one row per sample, in
#' the same order as the samples in \code{rnbSet}. All columns are kept as
#' \code{colData} of the result; no column is required. Default
#' \code{NULL} uses \code{RnBeads::pheno(rnbSet)}.
#' @return A \code{\linkS4class{methylTFRdeviations}} object with the assays
#' \code{deviations}, \code{z} and \code{expected}, motifs in rows and
#' samples in columns.
#' @examples
#' # A minimal end-to-end run on the BATF example data bundled with the
#' # package. The bundled calls are wrapped in an RnBiseqSet so that the
#' # example exercises the same code path as a preprocessed RnBeads object.
#' load(system.file("extdata", "example_data.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
#' load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
#' load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
#'
#' if (requireNamespace("RnBeads", quietly = TRUE) &&
#'     requireNamespace("RnBeads.hg38", quietly = TRUE)) {
#'     # identifiers.column is a session option, not a property of the object,
#'     # so it has to be set before the set is built or samples() falls back
#'     # to row numbers and methylTFR cannot name the columns it returns
#'     old_ids <- RnBeads::rnb.getOption("identifiers.column")
#'     RnBeads::rnb.options(identifiers.column = "sampleName")
#'
#'     sites <- data.frame(
#'         chr = as.character(GenomicRanges::seqnames(msites)),
#'         start = GenomicRanges::start(msites),
#'         strand = "*",
#'         stringsAsFactors = FALSE
#'     )
#'     rnbSet <- RnBeads::RnBiseqSet(
#'         pheno = data.frame(sampleName = "sample_1"),
#'         sites = sites,
#'         meth = matrix(msites$score, ncol = 1,
#'                       dimnames = list(NULL, "sample_1")),
#'         covg = matrix(msites$coverage, ncol = 1,
#'                       dimnames = list(NULL, "sample_1")),
#'         assembly = "hg38",
#'         summarize.regions = FALSE
#'     )
#'
#'     devs <- runMethylTFRFromRnBeads(
#'         rnbSet = rnbSet,
#'         tfBindsites = tf_bindsites,
#'         gcfreqs = gcfreqs,
#'         gcDist = gcdist
#'     )
#'     RnBeads::rnb.options(identifiers.column = old_ids)
#'     deviations(devs)
#' }
#' @seealso \code{\link{runMethylTFR}} for running methylTFR from
#' per-sample BED files.
#' @importFrom logger log_info
#' @importFrom methods is
#' @export
runMethylTFRFromRnBeads <- function(
    rnbSet, tfBindsites = NULL, gcfreqs = NULL, gcDist = NULL,
    chunkSize = 20, threads = 1, enhancer = NULL, ignoreStrand = TRUE,
    covThreshold = 1, dpvalThreshold = 0.05, sampleAnn = NULL,
    tempDir = file.path(tempdir(), "methylTFR")
) {
    checkRnbInputs(rnbSet)
    checkAnnotationInputs(tfBindsites, gcfreqs, gcDist, enhancer)
    opts <- checkRunOptions(chunkSize, threads, ignoreStrand, covThreshold)

    sampleIds <- rnbSampleIds(rnbSet)
    sampleAnn <- resolveRnbSampleAnn(rnbSet, sampleAnn, sampleIds)
    sitesGr <- rnbSitesToGRanges(rnbSet, opts$ignoreStrand)
    mode <- rnbQualityMode(rnbSet, dpvalThreshold)
    log_info(
        "Found ", length(sitesGr), " sites across ",
        length(sampleIds), " samples"
    )

    msitesFun <- function(i) {
        rnbSampleMsites(
            rnbSet = rnbSet, sitesGr = sitesGr, index = i,
            covThreshold = opts$covThreshold, hasCovg = mode$hasCovg,
            dpvalThreshold = dpvalThreshold, hasDpval = mode$hasDpval
        )
    }

    methylTFRCore(
        sampleIds = sampleIds, msitesFun = msitesFun,
        samples = sampleAnn, tfBindsites = tfBindsites,
        gcfreqs = gcfreqs, gcDist = gcDist, chunkSize = opts$chunkSize,
        threads = opts$threads, enhancer = enhancer,
        ignoreStrand = opts$ignoreStrand, tempDir = tempDir
    )
}


#' @title rnbSampleIds
#' @description Determine the sample identifiers of an RnBeads object.
#' @param rnbSet An \code{RnBSet} object.
#' @return A character vector of sample identifiers.
#' @keywords internal
rnbSampleIds <- function(rnbSet) {
    # Same as RnBeads' samples(): the column named by the
    # "identifiers.column" option. samples() is not exported by all
    # RnBeads versions, so it is not called directly.
    ids <- tryCatch({
        id_col <- RnBeads::rnb.getOption("identifiers.column")
        ph <- RnBeads::pheno(rnbSet)
        if (!is.null(id_col) && id_col %in% colnames(ph)) ph[[id_col]]
    }, error = function(e) NULL)
    if (length(ids) == 0) {
        ids <- tryCatch(
            colnames(RnBeads::meth(rnbSet, type = "sites", i = 1L)),
            error = function(e) NULL
        )
    }
    if (length(ids) == 0) {
        ids <- tryCatch(
            rownames(RnBeads::pheno(rnbSet)),
            error = function(e) NULL
        )
    }
    if (length(ids) == 0) {
        stop("Could not determine sample identifiers from the RnBSet object.")
    }
    ids <- as.character(ids)
    if (any(grepl("\\.idat$|_Grn|_Red", ids, ignore.case = TRUE))) {
        stop(
            "Sample identifiers look like IDAT channel files, not samples: '",
            ids[1], "'. The RnBSet was imported with the raw channel files ",
            "as samples, so each array appears twice and the methylation ",
            "values are not comparable. Re-import with a sample sheet and ",
            "an identifiers.column before running methylTFR."
        )
    }
    if (identical(ids, as.character(seq_along(ids)))) {
        stop(
            "Sample identifiers are the row numbers 1..", length(ids),
            ", so the RnBSet carries no usable sample names in this ",
            "session. identifiers.column is a session option rather than ",
            "a property of the saved object, so set it before loading, ",
            "for example rnb.options(identifiers.column = 'Sample.ID'), ",
            "or assign colnames on the methylation matrix yourself."
        )
    }
    if (anyDuplicated(ids)) {
        stop(
            "The RnBSet has duplicated sample identifiers, for example '",
            ids[anyDuplicated(ids)], "'. Resolve or average the duplicates ",
            "before running methylTFR."
        )
    }
    return(ids)
}


#' @title rnbAnnotationTarget
#' @description Resolve the annotation target of an RnBeads object.
#' @param rnbSet An \code{RnBSet} object.
#' @return A character scalar naming the annotation target, \code{"sites"}
#' for sequencing sets and the array platform for array sets.
#' @importFrom methods is
#' @keywords internal
rnbAnnotationTarget <- function(rnbSet) {
    if (is(rnbSet, "RnBeadSet")) rnbSet@target else "sites"
}


#' @title rnbProbeIds
#' @description Read the probe identifiers of an RnBeads object.
#' @details The \code{sites} slot of an \code{RnBSet} is an index matrix
#' whose row names carry the probe identifiers.
#' @param rnbSet An \code{RnBSet} object.
#' @return A character vector of probe identifiers, or \code{NULL}.
#' @keywords internal
rnbProbeIds <- function(rnbSet) {
    ids <- tryCatch(rownames(rnbSet@sites), error = function(e) NULL)
    if (length(ids) == 0) {
        ids <- tryCatch(
            rownames(RnBeads::meth(rnbSet, type = "sites")),
            error = function(e) NULL
        )
    }
    if (length(ids) == 0) NULL else as.character(ids)
}


#' @title rnbAnnotationById
#' @description Subset a genome-wide annotation track to the probes of an
#' RnBeads object, matching on probe identifier.
#' @details The \code{sites} slot holds a three-column index matrix, not a
#' row index into the annotation, so the track is matched by identifier
#' rather than by position.
#' @param rnbSet An \code{RnBSet} object.
#' @param target Character scalar naming the annotation target.
#' @param assembly Character scalar naming the genome assembly.
#' @return A \code{data.frame} of annotation rows, or \code{NULL} when the
#' identifiers cannot be matched.
#' @keywords internal
rnbAnnotationById <- function(rnbSet, target, assembly) {
    ids <- rnbProbeIds(rnbSet)
    if (is.null(ids)) {
        return(NULL)
    }
    track <- RnBeads::rnb.get.annotation(target, assembly)
    full <- RnBeads::rnb.annotation2data.frame(track)
    idx <- match(ids, rownames(full))
    if (anyNA(idx)) {
        return(NULL)
    }
    full[idx, , drop = FALSE]
}


#' @title rnbAnnotationTable
#' @description Look up the site or probe annotation of an RnBeads object.
#' @details The annotation stored in the object is used when available,
#' otherwise the genome-wide track registered for \code{target} is
#' matched to the probes of the object by identifier.
#' @param rnbSet An \code{RnBSet} object.
#' @param target Character scalar naming the annotation target.
#' @param assembly Character scalar naming the genome assembly.
#' @return A \code{data.frame} with one row per site or probe.
#' @keywords internal
rnbAnnotationTable <- function(rnbSet, target, assembly) {
    pull <- function(expr) tryCatch(expr, error = function(e) NULL)

    ann <- pull(RnBeads::annotation(rnbSet, type = target))
    if (is.null(ann)) {
        ann <- pull(RnBeads::annotation(rnbSet, type = "sites"))
    }
    if (is.null(ann) && target != "sites") {
        ann <- pull(rnbAnnotationById(rnbSet, target, assembly))
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


#' @title rnbSitesToGRanges
#' @description Build a \code{GRanges} object of the site or probe
#' annotation of an RnBeads object.
#' @param rnbSet An \code{RnBSet} object.
#' @param ignoreStrand if TRUE, all ranges are returned with strand
#' \code{"*"}.
#' @return A \code{GRanges} object with one range per site or probe.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom logger log_info
#' @importFrom methods is
#' @keywords internal
rnbSitesToGRanges <- function(rnbSet, ignoreStrand = TRUE) {
    assembly <- rnbSet@assembly
    target <- rnbAnnotationTarget(rnbSet)
    log_info("Annotation target: ", target, " | assembly: ", assembly)
    invisible(requireNamespace(paste0("RnBeads.", assembly), quietly = TRUE))

    ann <- rnbAnnotationTable(rnbSet, target, assembly)
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


#' @title rnbHasCoverage
#' @description Test whether an RnBeads object carries coverage
#' information.
#' @param rnbSet An \code{RnBSet} object.
#' @return A logical scalar.
#' @keywords internal
rnbHasCoverage <- function(rnbSet) {
    res <- tryCatch(
        !is.null(RnBeads::covg(rnbSet, type = "sites", j = 1L)),
        error = function(e) NULL
    )
    if (is.null(res)) {
        res <- tryCatch(
            !is.null(RnBeads::covg(rnbSet, type = "sites")),
            error = function(e) FALSE
        )
    }
    return(isTRUE(res))
}


#' @title rnbHasDpval
#' @description Test whether an RnBeads object carries detection
#' p-values.
#' @param rnbSet An \code{RnBSet} object.
#' @return A logical scalar.
#' @keywords internal
rnbHasDpval <- function(rnbSet) {
    res <- tryCatch(
        !is.null(RnBeads::dpval(rnbSet, type = "sites", j = 1L)),
        error = function(e) NULL
    )
    if (is.null(res)) {
        res <- tryCatch(
            !is.null(RnBeads::dpval(rnbSet, type = "sites")),
            error = function(e) FALSE
        )
    }
    return(isTRUE(res))
}


#' @title rnbQualityMode
#' @description Decide which per-site quality filter applies to an
#' RnBeads object and report it.
#' @param rnbSet An \code{RnBSet} object.
#' @param dpvalThreshold numeric detection p-value threshold.
#' @return A list with the logical flags \code{hasCovg} and
#' \code{hasDpval}.
#' @importFrom logger log_info log_warn
#' @keywords internal
rnbQualityMode <- function(rnbSet, dpvalThreshold) {
    hasCovg <- rnbHasCoverage(rnbSet)
    hasDpval <- rnbHasDpval(rnbSet)
    if (hasDpval) {
        log_info("Filtering probes at detection p-value <= ", dpvalThreshold)
    } else if (!hasCovg) {
        log_warn(
            "RnBSet carries neither coverage nor detection p-values, ",
            "no quality filtering is applied"
        )
    }
    list(hasCovg = hasCovg, hasDpval = hasDpval)
}


#' @title rnbColumn
#' @description Extract a single sample column from an RnBeads accessor.
#' @param accessor An RnBeads accessor function.
#' @param rnbSet An \code{RnBSet} object.
#' @param index Integer index of the sample to extract.
#' @return A numeric vector with one value per site.
#' @keywords internal
rnbColumn <- function(accessor, rnbSet, index) {
    vals <- tryCatch(
        accessor(rnbSet, type = "sites", j = index),
        error = function(e) NULL
    )
    if (is.null(vals)) {
        vals <- accessor(rnbSet, type = "sites")
        if (is.matrix(vals) || is.data.frame(vals)) {
            vals <- vals[, index]
        }
    }
    return(as.numeric(as.vector(vals)))
}


#' @title rnbSampleMsites
#' @description Extract the methylation calls of a single sample from an
#' RnBeads object as a \code{GRanges} object.
#' @param rnbSet An \code{RnBSet} object.
#' @param sitesGr A \code{GRanges} object of site positions.
#' @param index Integer index of the sample to extract.
#' @param covThreshold numeric coverage threshold.
#' @param hasCovg logical, whether coverage filtering applies.
#' @param dpvalThreshold numeric detection p-value threshold.
#' @param hasDpval logical, whether detection p-value filtering applies.
#' @return A \code{GRanges} object restricted to valid methylation calls.
#' @importFrom logger log_info
#' @keywords internal
rnbSampleMsites <- function(
    rnbSet, sitesGr, index, covThreshold = 1, hasCovg = TRUE,
    dpvalThreshold = 0.05, hasDpval = FALSE
) {
    index <- as.integer(index)
    mvals <- rnbColumn(RnBeads::meth, rnbSet, index)
    if (length(mvals) != length(sitesGr)) {
        stop("Methylation values do not match annotated sites.")
    }

    keep <- !is.na(mvals)
    if (hasCovg) {
        cvals <- rnbColumn(RnBeads::covg, rnbSet, index)
        keep <- keep & !is.na(cvals) & cvals >= covThreshold
    } else {
        cvals <- rep(NA_real_, length(sitesGr))
    }
    if (hasDpval) {
        dvals <- rnbColumn(RnBeads::dpval, rnbSet, index)
        keep <- keep & !is.na(dvals) & dvals <= dpvalThreshold
    }

    log_info(
        "Sample ", index, ": ", sum(keep), " of ", length(keep),
        " sites retained (", round(100 * mean(keep), 2), "%)"
    )
    if (!any(keep)) {
        stop("No sites passed the quality thresholds for sample index ", index)
    }

    gr <- sitesGr[keep]
    gr$score <- as.numeric(mvals[keep])
    gr$coverage <- as.numeric(cvals[keep])
    return(gr)
}
