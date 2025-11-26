# Example WGBS data downloaded from European Genome-Phenome Archive. 
# It is converted to GRanges format, so it can be directly used as example data.
# source from https://ega-archive.org/studies/EGAS00001001624

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(GenomicRanges)
  library(methylTFRAnnotationHg38)
  library(logger)
  library(rtracklayer)
})

# Load the annotation
gcfreqs <-  getGCfreq(motifSet = "jaspar2020_distal")
gcdist <- getGenomeGC("hg38")
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

# Read the data as GRanges object
filename <- "/icbb/projects/skumar/memoryTcells/bed/51_Hf03_BlEM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"
msites <- read_methylome(filename,"bissnp")

#Subset data and gcdist to chr1
msites <- msites[seqnames(msites) == "chr1"]
gcdist <- gcdist[seqnames(gcdist) == "chr1"]

#Subset to first 1000 rows
msites <- msites[1:1000]
gcdist <- subsetByOverlaps(gcdist, msites)

# Extract BATF annotations
tf_bindsites <- tf_bindsites[["BATF"]]
tf_bindsites <- list(tf_bindsites)
names(tf_bindsites) <- "BATF"
gcfreqs <- gcfreqs[["BATF"]]
gcfreqs <- list(gcfreqs)
names(gcfreqs) <- "BATF"

# Save the data as RDA files
save(msites, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/example_data.rda")
save(tf_bindsites, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/BATF_tf_bindsites.rda")
save(gcfreqs, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/BATF_gcfreqs.rda")
save(gcdist, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/gcdist_subset.rda")

# Compresse the data
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/example_data.rda", "auto")
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/BATF_tf_bindsites.rda", "auto")
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/BATF_gcfreqs.rda", "auto")
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/gcdist_subset.rda", "auto")

############################################################################
# Example data for plotMotifFootprint and related functions
############################################################################
# Load the full dataset
msites <- read_methylome(filename,"bissnp",1)
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
gcdist <- getGenomeGC("hg38")
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

# Check the overlaps with gc_freqs
motif <- "BATF"
tfbs <- tf_bindsites[[motif]]
tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
tfbs <- subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE)
gcfreq <- gcfreqs[[motif]]

# Save the subsetted data
msites_sub <- subsetByOverlaps(msites, tfbs, ignore.strand = TRUE)
gcdist <- subsetByOverlaps(gcdist, tfbs, ignore.strand = TRUE)

# Remove x chromosome
msites_sub <- msites_sub[seqnames(msites_sub) != "chrX"]
save(msites_sub, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/msites_sub.rda")
save(gcdist, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/gcdist_BATF.rda")

# Compress the data
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/msites_sub.rda", "auto")
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/gcdist_BATF.rda", "auto")

############################################################################
# Example deviations from Gunduz 2025 paper
############################################################################

# Load methylTFRdeviations objects
tc_mem <- readRDS("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks/Tc-Mem_deviations.RDS")
tc_naive <- readRDS("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks/Tc-Naive_deviations.RDS")

# Merge the two objects
#devs <- cbind(tc_mem[1:10,1:5], tc_naive[1:10,1:5])
tc_mem <- tc_mem[1:10,1:5]
tc_naive <- tc_naive[1:10,1:5]

# Save the subsetted objects
save(tc_mem, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/tc_mem.rda")
save(tc_naive, file = "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/tc_naive.rda")

# Compress the data
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/tc_mem.rda", "auto")
tools::resaveRdaFiles("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/tc_naive.rda", "auto")

############################################################################

