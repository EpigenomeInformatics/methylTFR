# Creates the example data in inst/extdata from open-access data.
#
# Methylation calls: ENCODE whole-genome bisulfite sequencing, experiment
# ENCSR663MXB (GRCh38), CpG bedMethyl file. ENCODE data are released without
# restrictions: https://www.encodeproject.org/experiments/ENCSR663MXB/
#
# Annotations: methylTFRAnnotationHg38
# (https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38).
#
suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
    library(methylTFRAnnotationHg38)
    library(jsonlite)
})

extdataDir <- Sys.getenv("MTFR_EXTDATA_DIR", file.path("inst", "extdata"))
dir.create(extdataDir, showWarnings = FALSE, recursive = TRUE)

saveCompressed <- function(..., file) {
    path <- file.path(extdataDir, file)
    save(..., file = path)
    tools::resaveRdaFiles(path, compress = "auto")
}

############################################################################
# Download the ENCODE CpG bedMethyl file
############################################################################

encodeUrl <- "https://www.encodeproject.org"
experiment <- "ENCSR663MXB"
meta <- fromJSON(sprintf(
    "%s/experiments/%s/?format=json", encodeUrl, experiment
), simplifyVector = FALSE)

isCpgBed <- function(f) {
    identical(f$file_format, "bed") &&
        identical(f$output_type, "methylation state at CpG") &&
        identical(f$assembly, "GRCh38") &&
        identical(f$status, "released")
}
cpgFiles <- Filter(isCpgBed, meta$files)
if (length(cpgFiles) == 0) {
    stop("No released GRCh38 CpG bedMethyl file found for ", experiment)
}
bedFile <- cpgFiles[[1]]
message("Using ", bedFile$accession, " from ", experiment)

localBed <- file.path(tempdir(), paste0(bedFile$accession, ".bed.gz"))
if (!file.exists(localBed)) {
    utils::download.file(
        paste0(encodeUrl, bedFile$href), localBed, mode = "wb"
    )
}

############################################################################
# Example data for computeDeviation, run_methyltfr and related functions
############################################################################

gcfreqs <- getGCfreq(motifSet = "jaspar2020")
gcdist <- getGenomeGC()
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

msitesAll <- read_methylome(localBed, "encode")
msitesAll <- msitesAll[seqnames(msitesAll) == "chr1"]

# Keep chr1 CpGs within the BATF binding-site windows used by
# computeDeviation(), and the first 1000 of them
batf <- tf_bindsites[["BATF"]]
batfWindows <- resize(batf, width(batf)[1] + 130, fix = "center")
msites <- subsetByOverlaps(msitesAll, batfWindows, ignore.strand = TRUE)
msites <- msites[seq_len(min(1000, length(msites)))]

gcdist <- gcdist[seqnames(gcdist) == "chr1"]
gcdist <- subsetByOverlaps(gcdist, msites)

tf_bindsites <- list(BATF = batf)
gcfreqs <- list(BATF = gcfreqs[["BATF"]])

saveCompressed(msites, file = "example_data.rda")
saveCompressed(tf_bindsites, file = "BATF_tf_bindsites.rda")
saveCompressed(gcfreqs, file = "BATF_gcfreqs.rda")
saveCompressed(gcdist, file = "gcdist_subset.rda")

############################################################################
# Example data for plotMotifFootprint and related functions
############################################################################

msitesAll <- read_methylome(localBed, "encode")
gcdist <- getGenomeGC()

msites_sub <- subsetByOverlaps(msitesAll, batfWindows, ignore.strand = TRUE)
msites_sub <- msites_sub[!seqnames(msites_sub) %in% c("chrX", "chrY")]
gcdist <- subsetByOverlaps(gcdist, batfWindows, ignore.strand = TRUE)

saveCompressed(msites_sub, file = "msites_sub.rda")
saveCompressed(gcdist, file = "gcdist_BATF.rda")

############################################################################
# Example deviations from the immune methylome atlas (Gunduz et al. 2025)
############################################################################
# These are per-motif pseudobulk deviation scores, not CpG calls. Set
# MTFR_PSEUDOBULK_DIR to the directory holding the *_deviations.RDS files.

pseudobulkDir <- Sys.getenv("MTFR_PSEUDOBULK_DIR")
if (nzchar(pseudobulkDir)) {
    tc_mem <- readRDS(file.path(pseudobulkDir, "Tc-Mem_deviations.RDS"))
    tc_naive <- readRDS(file.path(pseudobulkDir, "Tc-Naive_deviations.RDS"))
    tc_mem <- tc_mem[1:10, 1:5]
    tc_naive <- tc_naive[1:10, 1:5]
    saveCompressed(tc_mem, file = "tc_mem.rda")
    saveCompressed(tc_naive, file = "tc_naive.rda")
} else {
    message("MTFR_PSEUDOBULK_DIR not set; tc_mem.rda and tc_naive.rda kept")
}
