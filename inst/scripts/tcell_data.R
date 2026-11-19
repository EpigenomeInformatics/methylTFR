#!/usr/bin/env Rscript

#####################################################################
# tcell_data.R
#
# Regenerates inst/extdata/immuneDeviations.rda, the example dataset
# used by the T cell memory vignette.
#
# Source: pseudobulk deviation scores (jaspar2020_distal) from
#   Gunduz IB et al. bioRxiv 2025.09.09.675101.
#   doi:10.1101/2025.09.09.675101
#
# The source data are under controlled access, so the shipped samples
# are simulated from per-cell-type motif means and SDs; all motifs are
# kept. Differences are naive minus memory, as in 10_zdiff.R.
#
# Run from the package root:
#     MTFR_PSEUDOBULK_DIR=/path/to/pseudobulks \
#       Rscript inst/scripts/tcell_data.R
#
# Environment variables:
#     MTFR_PSEUDOBULK_DIR  directory holding *_deviations.RDS (required)
#     MTFR_N_PER_TYPE      simulated samples per cell type (default 10)
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

# Directory holding the per-sample pseudobulks used in 10_zdiff.R.
pseudobulk_dir <- Sys.getenv("MTFR_PSEUDOBULK_DIR")
if (!nzchar(pseudobulk_dir) || !dir.exists(pseudobulk_dir)) {
    stop(
        "Pseudobulk directory not found: ", pseudobulk_dir,
        ". Set MTFR_PSEUDOBULK_DIR to the directory holding the ",
        "*_deviations.RDS files."
    )
}
n_per_type <- as.integer(Sys.getenv("MTFR_N_PER_TYPE", "10"))
extdata_dir <- Sys.getenv("MTFR_EXTDATA_DIR", file.path("inst", "extdata"))

# T cell subsets plus three other lineages
keep_types <- c(
    "Tc-Naive", "Tc-Mem", "Th-Naive", "Th-Mem",
    "B-cell", "NK-cell", "Monocyte"
)

# Naive first: differences are naive minus memory
subset_pairs <- list(
    Tc = c(naive = "Tc-Naive", memory = "Tc-Mem"),
    Th = c(naive = "Th-Naive", memory = "Th-Mem")
)

dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## Helpers, matching 10_zdiff.R
## ------------------------------------------------------------------

# Cell type is the token before the first underscore of the sample name,
# e.g. "Tc-Mem_01.bedGraph.bed" -> "Tc-Mem".
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

# Check parsed cell types against colData
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

# Only the cell type is kept; donor IDs and conditions are not shipped
sample_annot <- data.frame(
    cell_type = factor(cell_type, levels = keep_types),
    row.names = colnames(dev_mat)
)

message(sprintf(
    "Loaded %d motifs x %d samples across %d cell types",
    nrow(dev_mat), ncol(dev_mat), length(keep_types)
))

## ------------------------------------------------------------------
## 2. Drop motifs that are undefined in any sample
## ------------------------------------------------------------------

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
## 3. Simulate samples from per-cell-type summary statistics
## ------------------------------------------------------------------
## The source data are under controlled access (dbGaP phs003204.v1.p1),
## so no per-sample values are shipped. For each cell type, the mean and
## SD of every motif are taken across all source samples, and
## n_per_type samples are drawn from a normal distribution with those
## parameters.

min_source <- 10
n_source <- table(sample_annot$cell_type)
stopifnot(all(n_source >= min_source))

type_mean <- vapply(keep_types, function(ct) {
    rowMeans(dev_mat[, sample_annot$cell_type == ct, drop = FALSE])
}, numeric(nrow(dev_mat)))
type_sd <- vapply(keep_types, function(ct) {
    apply(dev_mat[, sample_annot$cell_type == ct, drop = FALSE], 1, stats::sd)
}, numeric(nrow(dev_mat)))
rownames(type_mean) <- rownames(type_sd) <- rownames(dev_mat)

sim <- lapply(keep_types, function(ct) {
    matrix(
        stats::rnorm(nrow(type_mean) * n_per_type,
            mean = type_mean[, ct], sd = type_sd[, ct]
        ),
        nrow = nrow(type_mean),
        dimnames = list(
            rownames(type_mean),
            sprintf("%s_%02d", ct, seq_len(n_per_type))
        )
    )
})
# Three significant digits keep the .rda file small
sim_mat <- signif(do.call(base::cbind, sim), 3)

sample_annot <- data.frame(
    cell_type = factor(rep(keep_types, each = n_per_type), levels = keep_types),
    row.names = colnames(sim_mat)
)
stopifnot(
    identical(colnames(sample_annot), "cell_type"),
    !any(colnames(sim_mat) %in% colnames(dev_mat))
)

message(sprintf(
    "Simulated %d motifs x %d samples (%d per cell type) from %d source samples",
    nrow(sim_mat), ncol(sim_mat), n_per_type, sum(n_source)
))

# The source matrix is not used past this point
dev_mat <- sim_mat
rm(sim_mat, sim, objs)

## ------------------------------------------------------------------
## 4. Row-wise Z-scores
## ------------------------------------------------------------------
## Computed on the simulated samples

row_sd <- apply(dev_mat, 1, stats::sd)
z_mat <- (dev_mat - rowMeans(dev_mat)) / row_sd
z_mat[!is.finite(z_mat)] <- 0
z_mat <- round(z_mat, 3)

## ------------------------------------------------------------------
## 5. Sanity checks
## ------------------------------------------------------------------

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
    res <- differentialDeviationTest(
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
        "Simulated from per-cell-type motif means and SDs of pseudobulk",
        "single-cell methylomes of human immune cells;",
        "see inst/scripts/tcell_data.R"
    ),
    simulated = TRUE,
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
