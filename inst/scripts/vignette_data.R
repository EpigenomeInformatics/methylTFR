# Example WGBS data from ENCODE (ENCSR663MXB, GRCh38), open access.
# https://www.encodeproject.org/experiments/ENCSR663MXB/
# Run from the package root.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
    library(methylTFRAnnotationHg38)
    library(jsonlite)
})

extdata <- file.path("inst", "extdata")

saveCompressed <- function(..., file) {
    path <- file.path(extdata, file)
    save(..., file = path)
    tools::resaveRdaFiles(path, compress = "auto")
}

# Download the CpG bedMethyl file
encodeUrl <- "https://www.encodeproject.org"
meta <- fromJSON(
    paste0(encodeUrl, "/experiments/ENCSR663MXB/?format=json"),
    simplifyVector = FALSE
)
bedFile <- Filter(function(f) {
    identical(f$file_format, "bed") &&
        identical(f$output_type, "methylation state at CpG") &&
        identical(f$assembly, "GRCh38") &&
        identical(f$status, "released")
}, meta$files)[[1]]
localBed <- file.path(tempdir(), paste0(bedFile$accession, ".bed.gz"))
if (!file.exists(localBed)) {
    download.file(paste0(encodeUrl, bedFile$href), localBed, mode = "wb")
}

# Load the annotation
gcfreqs <- getGCfreq(motifSet = "jaspar2020")
gcdist <- getGenomeGC()
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

batf <- tf_bindsites[["BATF"]]
batfWindows <- resize(batf, width(batf)[1] + 130, fix = "center")
msitesAll <- readMethylome(localBed, "encode")

############################################################################
# Example data for computeDeviation and runMethylTFR
############################################################################
# First 1000 chr1 CpGs in BATF windows
msites <- msitesAll[seqnames(msitesAll) == "chr1"]
msites <- subsetByOverlaps(msites, batfWindows, ignore.strand = TRUE)
msites <- msites[seq_len(min(1000, length(msites)))]
gcdist_sub <- subsetByOverlaps(gcdist[seqnames(gcdist) == "chr1"], msites)

# All example CpGs are on chr1, so only chr1 sites are kept
tf_bindsites <- list(BATF = batf[seqnames(batf) == "chr1"])
gcfreqs <- list(BATF = gcfreqs[["BATF"]])

saveCompressed(msites, file = "example_data.rda")
saveCompressed(tf_bindsites, file = "BATF_tf_bindsites.rda")
saveCompressed(gcfreqs, file = "BATF_gcfreqs.rda")
gcdist <- gcdist_sub
saveCompressed(gcdist, file = "gcdist_subset.rda")

############################################################################
# Example data for plotMotifFootprint
############################################################################
# chr1 CpGs in BATF windows, and only the GC windows they fall in
msites_sub <- msitesAll[seqnames(msitesAll) == "chr1"]
msites_sub <- subsetByOverlaps(msites_sub, batfWindows, ignore.strand = TRUE)
gcdist <- subsetByOverlaps(getGenomeGC(), msites_sub, ignore.strand = TRUE)

saveCompressed(msites_sub, file = "msites_sub.rda")
saveCompressed(gcdist, file = "gcdist_BATF.rda")

############################################################################
# Example deviations simulated from Gunduz 2025 summary statistics
############################################################################
# The source data are under controlled access (dbGaP phs003204.v1.p1), so
# samples are simulated from the per-motif mean and SD of each cell type.
# Set MTFR_PSEUDOBULK_DIR to the folder with the *_deviations.RDS files
pseudobulkDir <- Sys.getenv("MTFR_PSEUDOBULK_DIR")
if (nzchar(pseudobulkDir)) {
    simulateType <- function(ct, nMotifs = 10, nSamples = 5) {
        src <- readRDS(file.path(pseudobulkDir, paste0(ct, "_deviations.RDS")))
        dev <- as.matrix(deviations(src))[seq_len(nMotifs), , drop = FALSE]
        stopifnot(ncol(dev) >= 10, all(is.finite(dev)))
        ids <- sprintf("%s_%02d", ct, seq_len(nSamples))
        sim <- matrix(
            rnorm(nMotifs * nSamples,
                mean = rowMeans(dev), sd = apply(dev, 1, sd)
            ),
            nrow = nMotifs, dimnames = list(rownames(dev), ids)
        )
        z <- (sim - rowMeans(sim)) / apply(sim, 1, sd)
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(deviations = sim, z = z),
            colData = S4Vectors::DataFrame(
                cell_type = rep(ct, nSamples), row.names = ids
            )
        )
        S4Vectors::metadata(se) <- list(simulated = TRUE)
        methods::new("methylTFRdeviations", se)
    }
    set.seed(42)
    tc_mem <- simulateType("Tc-Mem")
    tc_naive <- simulateType("Tc-Naive")
    saveCompressed(tc_mem, file = "tc_mem.rda")
    saveCompressed(tc_naive, file = "tc_naive.rda")
}
