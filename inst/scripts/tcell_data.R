#!/usr/bin/env Rscript

#####################################################################
# tcell_data.R
#
# Regenerates inst/extdata/tcellMemory.rda, the example dataset used by
# the T cell memory vignette.
#
# Source: pseudobulk methylTFR deviation scores from the single-cell
# immune methylome atlas, scored against the JASPAR2020 motif set
# restricted to distal regulatory regions ("jaspar2020_distal"). One
# pseudobulk per donor sample per cell type.
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
#   * a motif counts as changed when abs(difference) > 0.5 and the
#     adjusted p-value is below 0.05
#
# The full set is seven cell types x 632 motifs x 38 samples, too large
# to ship in a Bioconductor package. This script cuts both axes: the
# four T cell subsets only, and the most variable N motifs among them.
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
#     MTFR_N_MOTIFS        number of motifs to retain (default 200)
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
    "/icbb/projects/igunduz/mTFR_bias_fix_v3",
    "all_pseudobulks_310824", "jaspar2020_distal"
)
pseudobulk_dir <- Sys.getenv("MTFR_PSEUDOBULK_DIR", default_pseudobulk_dir)
if (!nzchar(pseudobulk_dir) || !dir.exists(pseudobulk_dir)) {
    stop(
        "Pseudobulk directory not found: ", pseudobulk_dir,
        ". Set MTFR_PSEUDOBULK_DIR to the directory holding the ",
        "*_deviations.RDS files."
    )
}
n_motifs <- as.integer(Sys.getenv("MTFR_N_MOTIFS", "200"))
extdata_dir <- Sys.getenv("MTFR_EXTDATA_DIR", file.path("inst", "extdata"))

# Naive first, memory second, in both compartments. The order matters:
# differences are taken as naive minus memory.
subset_pairs <- list(
    Tc = c(naive = "Tc-Naive", memory = "Tc-Mem"),
    Th = c(naive = "Th-Naive", memory = "Th-Mem")
)
keep_types <- unname(unlist(subset_pairs))

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
    "Loaded %d motifs x %d samples across %d subsets",
    nrow(dev_mat), ncol(dev_mat), length(keep_types)
))

## ------------------------------------------------------------------
## 2. Drop motifs that are undefined in any sample
## ------------------------------------------------------------------
## A motif has no deviation score in a sample where too few of its
## binding sites are covered, which is expected for sparse pseudobulks.

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
## 3. Keep the most variable motifs
## ------------------------------------------------------------------
## Selection is on variability across the four subsets, so the retained
## motifs are the ones that carry a naive-to-memory or CD4-to-CD8
## signal. Motifs that never move teach the reader nothing.

variability <- computeZScoreVariability(dev_mat, method = "robust")
variability <- variability[order(-variability$variability), ]
n_keep <- min(n_motifs, nrow(variability))
selected <- variability$motifs[seq_len(n_keep)]

message(sprintf(
    "Retaining the %d most variable motifs (variability %.2f to %.2f)",
    n_keep, variability$variability[n_keep], variability$variability[1]
))

dev_mat <- dev_mat[selected, , drop = FALSE]

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
## 5. Sanity check: the naive-to-memory contrast must survive the cut
## ------------------------------------------------------------------
## Reproduces the 10_zdiff.R contrast on the reduced matrix, oriented
## naive minus memory, purely as a check that the shipped dataset still
## shows what the vignette claims it shows.
##
## Note the two scales. P-values come from the raw bias-corrected
## deviations, whose differences here span roughly +/- 0.15. The 0.5
## effect-size cutoff used in 10_zdiff.R applies to the Z-score
## difference, not to that raw difference, so it is applied to the z
## assay below. Mixing the two would silently return nothing.

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
    sub_dev <- dev_mat[, keep, drop = FALSE]
    sub_z <- z_mat[, keep, drop = FALSE]
    # naive minus memory, as in 10_zdiff.R
    res$diff <- rowMeans(sub_dev[, is_naive, drop = FALSE]) -
        rowMeans(sub_dev[, !is_naive, drop = FALSE])
    res$zdiff <- rowMeans(sub_z[, is_naive, drop = FALSE]) -
        rowMeans(sub_z[, !is_naive, drop = FALSE])
    res <- res[order(res$p_value_adjusted), ]
    n_changed <- sum(
        abs(res$zdiff) > 0.5 & res$p_value_adjusted < 0.05,
        na.rm = TRUE
    )
    message(sprintf(
        "  %s: %d/%d motifs with adj. p < 0.05; %d also |zdiff| > 0.5",
        paste(rev(pair), collapse = " vs "),
        sum(res$p_value_adjusted < 0.05, na.rm = TRUE), nrow(res), n_changed
    ))
    message(sprintf(
        "    strongest: %s (diff = %+.3f, zdiff = %+.2f, adj. p = %.2e)",
        res$motifs[1], res$diff[1], res$zdiff[1], res$p_value_adjusted[1]
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
tcellMemory <- methods::new("methylTFRdeviations", se)

metadata(tcellMemory) <- list(
    motifSet = "jaspar2020_distal",
    genome = "hg38",
    source = paste(
        "Pseudobulk single-cell methylomes of human immune cells;",
        "see inst/scripts/tcell_data.R"
    ),
    contrastOrientation = "naive minus memory"
)

## ------------------------------------------------------------------
## 7. Save
## ------------------------------------------------------------------

out_file <- file.path(extdata_dir, "tcellMemory.rda")
save(tcellMemory, file = out_file, compress = "xz")

message(sprintf(
    "\nWrote %s\n  %d motifs x %d samples\n  %.1f KB",
    out_file, nrow(tcellMemory), ncol(tcellMemory),
    file.size(out_file) / 1024
))
print(table(colData(tcellMemory)$cell_type))

