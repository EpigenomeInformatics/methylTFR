# Example WGBS data downloaded from European Genome-Phenome Archive. 
#It is converted to GRanges format, so it can be directly used as example data.
# source from https://ega-archive.org/studies/EGAS00001001624

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(GenomicRanges)
  library(methylTFRAnnotationHg38)
  library(logger)
})

# Load the annotation
gcfreqs <-  getGCfreq(motifSet = "jaspar2020")
gcdist <- getGenomeGC()
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")

# Read the data as GRanges object
filename <- "/icbb/projects/skumar/memoryTcells/bed/51_Hf03_BlEM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"
msites <- read_methylome(filename,"bissnp")

#Subset data and gcdist to chr1
msites <- msites[seqnames(msites) == "chr1"]
gcdist <- gcdist[seqnames(gcdist) == "chr1"]

# Extract FOXF2 annotations
tf_bindsites <- tf_bindsites[1]
gcfreqs <- gcfreqs[1]

# Save the data as RDA files
save(msites, file = "/icbb/projects/igunduz/methylTFR/inst/extdata/example_data.rda")
save(tf_bindsites, file = "/icbb/projects/igunduz/methylTFR/inst/extdata/FOXF2_tf_bindsites.rda")
save(gcfreqs, file = "/icbb/projects/igunduz/methylTFR/inst/extdata/FOXF2_gcfreqs.rda")
save(gcdist, file = "/icbb/projects/igunduz/methylTFR/inst/extdata/gcdist_subset.rda")

# Compresse the data
tools::resaveRdaFiles("/icbb/projects/igunduz/methylTFR/inst/extdata/example_data.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/methylTFR/inst/extdata/FOXF2_tf_bindsites.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/methylTFR/inst/extdata/FOXF2_gcfreqs.rda", "auto")
tools::resaveRdaFiles("/icbb/projects/igunduz/methylTFR/inst/extdata/gcdist_subset.rda", "auto")
