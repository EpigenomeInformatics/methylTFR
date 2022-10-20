
library(methylTFRann)
library(methylTFR)

library(parallel)
library(data.table)

# from methylTFRann
tf_bindsites <- getTFbindsites()
gc_dist <- getGenomeGC()
motif_gcfreq <- getGCfreq()

motif_list <- names(motif_gcfreq)

# files to be processed
files_list <- list.files(path="/data/blueprint/bed", pattern='*.bed', full.names=T)
deviation_distal <- data.frame(motifs = names(motif_gcfreq))
for (fname in files_list) {
    basename = unlist(str_split(fname, "/"))[10]
    prefix = str_replace(basename, ".bed", "")
    if (! prefix %in% colnames(deviation_distal)){
        print(paste0(fname, " Processing ..."))
        msites <- read_methylome(fname)
        assign(prefix, mclapply(motif_list, compute_deviation, msites = msites, tf_bindsites = tf_bindsites, gcfreqs = motif_gcfreq, gcdist = gc_dist, mc.cores = 16))
        deviation_distal[prefix] = unlist(get(prefix))
        print("Done")
        save(deviation_distal, file = "/data/gc_corrected_distal_deviation_all.Rds")
    }
}