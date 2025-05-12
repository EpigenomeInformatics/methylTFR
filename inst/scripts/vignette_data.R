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
gcdist <- getGenomeGC()
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

# Extract FOXF2 annotations
tf_bindsites <- tf_bindsites[1]
tf_bindsites <- list(unlist(tf_bindsites)[1:1000])
names(tf_bindsites) <- "FOXF2"
gcfreqs <- gcfreqs[1]

# Save the data as RDA files
save(msites, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/example_data.rda")
save(tf_bindsites, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/FOXF2_tf_bindsites.rda")
save(gcfreqs, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/FOXF2_gcfreqs.rda")
save(gcdist, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/gcdist_subset.rda")

# Compresse the data
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/example_data.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/FOXF2_tf_bindsites.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/FOXF2_gcfreqs.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/gcdist_subset.rda", "auto")

############################################################################
# Example data for plotMotifFootprint and related functions
############################################################################
# Load the full dataset
msites <- read_methylome(filename,"bissnp")

# Check the overlaps with gc_freqs
motif <- "FOXF2"
tfbs <- tf_bindsites[[motif]]
tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
gcfreq <- gcfreqs[[motif]]
hits <- suppressWarnings(findOverlaps(msites, tfbs, type = "within", ignore.strand = TRUE))

# Save bin_meth to make it faster
bin_meth <- addGCBintoMethylome(msites, gcdist,TRUE) 
save(bin_meth, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/bin_meth.rda")

# Save the subsetted data
msites_sub <- msites[hits@from]
save(msites_sub, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/msites_sub.rda")

# Compress the data
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/msites_sub.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/bin_meth.rda", "auto")

############################################################################
# Example deviations from Gunduz 2025 paper
############################################################################

# Load methylTFRdeviations objects
tc_mem <- readRDS("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks/Tc-Mem_deviations.RDS")
tc_naive <- readRDS("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks/Tc-Naive_deviations.RDS")

# Merge the two objects
devs <- cbind(tc_mem, tc_naive)

# Save the merged object
save(devs, file = "/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/tcells_deviations.rda")

# Compress the data
tools::resaveRdaFiles("/icbb/projects/igunduz/irem_github/methylTFR/inst/extdata/tcells_deviations.rda", "auto")

############################################################################
