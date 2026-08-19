#!/usr/bin/env Rscript

#####################################################################
# tcell_data.R
#
# Regenerates inst/extdata/immuneDeviations.rda, the example dataset
# used by the T cell memory vignette.
#
# Source: pseudobulk methylTFR deviation scores from the single-cell
# immune methylome atlas, scored against the JASPAR2020 motif set
# restricted to distal regulatory regions ("jaspar2020_distal"). One
# pseudobulk per donor sample per cell type.
#
# Data citation:
#   Gunduz IB, Wei B, Chen DC, Wang W, Hariharan M, Norell T, et al.
#   Dissecting epigenome dynamics in human immune cells upon viral and
#   chemical exposure by multimodal single-cell profiling.
#   bioRxiv 2025.09.09.675101. doi:10.1101/2025.09.09.675101
#
# Conventions follow src/integration/10_zdiff.R of the atlas
# manuscript, so that anything computed from this dataset lines up with
# the published figures:
#
#   * cell type is taken from the sample name, as the token before the
#     first underscore, with a trailing ".bedGraph" removed
#   * differential testing is two-sided and parametric, BH-adjusted
#   * differences are oriented NAIVE MINUS MEMORY, so a positive value
#     means higher deviation in the naive state
#
# The full set is seven cell types x 632 motifs x 38 samples, too large
# to ship in a Bioconductor package. It is reduced on the SAMPLE axis
# only: a fixed number of donor samples per cell type. Every motif is
# retained.
#
# Subsetting motifs would be the more obvious way to shrink the object,
# but selecting motifs by variability and then demonstrating
# computeZScoreVariability() on the result is circular, and pre-filtering
# to variable motifs inflates the within-sample null that the function
# calibrates against. Cutting samples leaves the motif distribution
# untouched.
#
# Run from the package root:
#     Rscript inst/scripts/tcell_data.R
# or, to point it somewhere else:
#     MTFR_PSEUDOBULK_DIR=/path/to/pseudobulks \
#       Rscript inst/scripts/tcell_data.R
#
# Environment variables:
#     MTFR_PSEUDOBULK_DIR  directory holding *_deviations.RDS
#                          (defaults to the cluster path below)
#     MTFR_N_PER_TYPE      samples kept per cell type (default 15)
#     MTFR_EXTDATA_DIR     output directory (default inst/extdata)
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
    library(methylTFR)
    library(SummarizedExperiment)
})

## ------------------------------------------------------------------
## Configuration
## ------------------------------------------------------------------

# Defaults to the cluster location of the per-sample pseudobulks used in
# 10_zdiff.R. Override with MTFR_PSEUDOBULK_DIR to run it anywhere else.
default_pseudobulk_dir <- file.path(
    "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/",
    "data", "sample_pseudobulks"
)
pseudobulk_dir <- Sys.getenv("MTFR_PSEUDOBULK_DIR", default_pseudobulk_dir)
if (!nzchar(pseudobulk_dir) || !dir.exists(pseudobulk_dir)) {
    stop(
        "Pseudobulk directory not found: ", pseudobulk_dir,
        ". Set MTFR_PSEUDOBULK_DIR to the directory holding the ",
        "*_deviations.RDS files."
    )
}
n_per_type <- as.integer(Sys.getenv("MTFR_N_PER_TYPE", "15"))
extdata_dir <- Sys.getenv("MTFR_EXTDATA_DIR", file.path("inst", "extdata"))

# The four T cell subsets carry the naive-to-memory contrast; the other
# three provide the lineage contrast that makes an unsupervised
# variability screen informative.
keep_types <- c(
    "Tc-Naive", "Tc-Mem", "Th-Naive", "Th-Mem",
    "B-cell", "NK-cell", "Monocyte"
)

# Naive first, memory second. The order matters: differences are taken
# as naive minus memory.
subset_pairs <- list(
    Tc = c(naive = "Tc-Naive", memory = "Tc-Mem"),
    Th = c(naive = "Th-Naive", memory = "Th-Mem")
)

dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## Helpers, matching 10_zdiff.R
## ------------------------------------------------------------------

# Cell type is the token before the first underscore of the sample name,
# e.g. "Tc-Mem_OP_S4_Long_D1.bedGraph.bed" -> "Tc-Mem".
get_groupname <- function(x) {
    sub("\\.bedGraph.*$", "", vapply(
        strsplit(x, split = "_"), `[`, character(1), 1L
    ))
}

## ------------------------------------------------------------------
## 1. Load and join the pseudobulk objects
## ------------------------------------------------------------------

files <- file.path(pseudobulk_dir, paste0(keep_types, "_deviations.RDS"))
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
    stop(
        "Missing pseudobulk file(s): ",
        paste(basename(missing), collapse = ", ")
    )
}

objs <- lapply(files, readRDS)
names(objs) <- keep_types

motifs <- rownames(objs[[1]])
if (!all(vapply(objs, function(x) identical(rownames(x), motifs),
    logical(1)
))) {
    stop("The pseudobulk objects do not share a common motif set.")
}

dev_mat <- do.call(base::cbind, lapply(objs, function(x) {
    as.matrix(methylTFR::deviations(x))
}))

# Derive the cell type the same way 10_zdiff.R does, then check it
# against the colData that the pseudobulk objects already carry.
cell_type <- get_groupname(colnames(dev_mat))

from_coldata <- unlist(lapply(objs, function(x) {
    as.character(as.data.frame(colData(x))$cell_type)
}), use.names = FALSE)
if (!identical(cell_type, from_coldata)) {
    warning(
        "Cell types parsed from the sample names disagree with the ",
        "cell_type column of the colData; using the parsed names."
    )
}

sample_annot <- do.call(base::rbind, unname(lapply(objs, function(x) {
    as.data.frame(colData(x))[, c("CommonMinID", "condition")]
})))
rownames(sample_annot) <- colnames(dev_mat)
sample_annot$cell_type <- factor(cell_type, levels = keep_types)

message(sprintf(
    "Loaded %d motifs x %d samples across %d cell types",
    nrow(dev_mat), ncol(dev_mat), length(keep_types)
))

## ------------------------------------------------------------------
## 2. Drop motifs that are undefined in any sample
## ------------------------------------------------------------------
## A motif has no deviation score in a sample where too few of its
## binding sites are covered, which is expected for sparse pseudobulks.
## This is the only motif-level filter applied.

finite_motif <- apply(dev_mat, 1, function(r) all(is.finite(r)))
if (any(!finite_motif)) {
    message(sprintf(
        "Dropping %d motif(s) undefined in at least one sample: %s",
        sum(!finite_motif),
        paste(rownames(dev_mat)[!finite_motif], collapse = ", ")
    ))
    dev_mat <- dev_mat[finite_motif, , drop = FALSE]
}

## ------------------------------------------------------------------
## 3. Reduce the sample axis
## ------------------------------------------------------------------

selected <- unlist(lapply(keep_types, function(ct) {
    idx <- which(sample_annot$cell_type == ct)
    if (length(idx) <= n_per_type) idx else sort(sample(idx, n_per_type))
}))

dev_mat <- dev_mat[, selected, drop = FALSE]
sample_annot <- sample_annot[selected, , drop = FALSE]

message(sprintf(
    "Keeping %d motifs x %d samples (%d per cell type)",
    nrow(dev_mat), ncol(dev_mat), n_per_type
))

## ------------------------------------------------------------------
## 4. Row-wise Z-scores
## ------------------------------------------------------------------
## Recomputed on the subset rather than carried over, because row-wise
## Z-scores depend on which samples are present. This is the same
## transform 10_zdiff.R applies via computeRowZScore, so
## deviationZScores() on the shipped object can be used directly for a
## Z-score difference without reaching for an internal function.

row_sd <- apply(dev_mat, 1, stats::sd)
z_mat <- (dev_mat - rowMeans(dev_mat)) / row_sd
z_mat[!is.finite(z_mat)] <- 0

## ------------------------------------------------------------------
## 5. Sanity checks
## ------------------------------------------------------------------
## Confirms that the shipped object still supports both analyses the
## vignette performs: an unsupervised variability screen across all cell
## types, and the naive-to-memory contrast within each T compartment.

variability <- computeZScoreVariability(dev_mat, method = "robust")
variability <- variability[order(-variability$variability), ]
message(sprintf(
    "  variability: %d motifs adj. p < 0.05; top: %s (%.2f)",
    sum(variability$p_value_adjusted < 0.05, na.rm = TRUE),
    variability$motifs[1], variability$variability[1]
))

for (nm in names(subset_pairs)) {
    pair <- subset_pairs[[nm]]
    keep <- sample_annot$cell_type %in% pair
    grp <- factor(as.character(sample_annot$cell_type[keep]),
        levels = c(pair[["naive"]], pair[["memory"]])
    )
    res <- differential_deviation_test(
        deviations = dev_mat[, keep, drop = FALSE],
        groups = grp,
        alternative = "two.sided",
        parametric = TRUE,
        padjMethod = "BH"
    )
    is_naive <- grp == pair[["naive"]]
    sub_z <- z_mat[, keep, drop = FALSE]
    # naive minus memory, as in 10_zdiff.R
    res$zdiff <- rowMeans(sub_z[, is_naive, drop = FALSE]) -
        rowMeans(sub_z[, !is_naive, drop = FALSE])
    res <- res[order(res$p_value_adjusted), ]
    message(sprintf(
        "  %s: %d/%d motifs adj. p < 0.05; strongest %s (zdiff %+.2f)",
        paste(rev(pair), collapse = " vs "),
        sum(res$p_value_adjusted < 0.05, na.rm = TRUE), nrow(res),
        res$motifs[1], res$zdiff[1]
    ))
}

## ------------------------------------------------------------------
## 6. Assemble a methylTFRdeviations object
## ------------------------------------------------------------------

se <- SummarizedExperiment(
    assays = list(deviations = dev_mat, z = z_mat),
    colData = S4Vectors::DataFrame(sample_annot),
    rowData = S4Vectors::DataFrame(motifs = rownames(dev_mat))
)
immuneDeviations <- methods::new("methylTFRdeviations", se)

metadata(immuneDeviations) <- list(
    motifSet = "jaspar2020_distal",
    genome = "hg38",
    source = paste(
        "Pseudobulk single-cell methylomes of human immune cells;",
        "see inst/scripts/tcell_data.R"
    ),
    citation = paste(
        "Gunduz IB, Wei B, Chen DC, Wang W, Hariharan M, Norell T,",
        "et al. Dissecting epigenome dynamics in human immune cells",
        "upon viral and chemical exposure by multimodal single-cell",
        "profiling. bioRxiv 2025.09.09.675101.",
        "doi:10.1101/2025.09.09.675101"
    ),
    contrastOrientation = "naive minus memory"
)

## ------------------------------------------------------------------
## 7. Save
## ------------------------------------------------------------------

out_file <- file.path(extdata_dir, "immuneDeviations.rda")
save(immuneDeviations, file = out_file, compress = "xz")

message(sprintf(
    "\nWrote %s\n  %d motifs x %d samples\n  %.1f KB",
    out_file, nrow(immuneDeviations), ncol(immuneDeviations),
    file.size(out_file) / 1024
))
print(table(colData(immuneDeviations)$cell_type))
