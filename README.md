# methylTFR

`methylTFR` is an R-package to analyze the methylation data from Infinium 450K microarray and bisulfite sequencing protocols. This tool aims to indentify the methylation patterns or variation on transcription factor binding sites in each individual cells or samples. This tool is still under development and testing stage.

# Requirements

- R (>= 3.5.0)
- GenomicRanges
- data.table
- dplyr
- stringr

```R

# To install dependencies 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages(c("data.table", "dplyr", "stringr"))

```

# Installation

```sh

git clone git@github.com:EpigenomeInformatics/methylTFR.git
cd methylTFR
Rscript -e "library(devtools); devtools::install('.')"

```

# QuickStart

```R

library(methylTFRann)
library(methylTFR)

library(parallel)
library(data.table)

# load annotation files from methylTFRann
tf_bindsites <- getTFbindsites()
gc_dist <- getGenomeGC()
motif_gcfreq <- getGCfreq()

motif_list <- names(motif_gcfreq)

# read files to be processed
files_list <- list.files(path="/data/blueprint/bed", pattern='*.bed', full.names=T)
deviation_distal <- data.frame(motifs = names(motif_gcfreq))
# process each file and calculate variability score for each samples or cells
for (fname in files_list) {
    basename = unlist(str_split(fname, "/"))[10]
    prefix = str_replace(basename, ".bed", "")
    if (! prefix %in% colnames(deviation_distal)){
        print(paste0(fname, " Processing ..."))
        msites <- read_methylome(fname)
        assign(prefix, mclapply(motif_list, compute_variability, msites = msites, tf_bindsites = tf_bindsites, gcfreqs = motif_gcfreq, gcdist = gc_dist, mc.cores = 16))
        deviation_distal[prefix] = unlist(get(prefix))
        print("Done")
        save(deviation_distal, file = "/data/gc_corrected_distal_deviation_all.Rds")
    }
}

```
