#' read_methylome
#' 
#'  read_methylome is a function to import methylation data into R Genomic 
#'  Range object. Bed file should be processed in the pipeline developed by 
#'  Fabian Müller and Christoph Bock. (EPP format)
#' @param  - bedfilename which contains methylation data EPP format 
#' @return \code{GenomicRange} object with methylation, coverage information
#' @export 
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom stringr str_split str_replace
read_methylome <- function(filename){
    if (!file.exists(filename)) {
        stop(paste(filename, " doesn't exist or path is incorrect !!"))
    }

    msites <- fread(filename, header = FALSE)
    mcov <- unlist(str_split(msites$V4, "/"))
    mcov <- as.numeric(str_replace(mcov, "'", ""))

    row_odd <- seq_len(length(mcov)) %% 2 

    msites <- GRanges(seqnames = msites$V1, 
                    ranges = IRanges(start = msites$V2, width = 3), 
                    strand = msites$V6,
                    T = mcov[row_odd == 0],
                    M = mcov[row_odd == 1],
                    score = msites$V5,
                    methylation = msites$V4,
                    coverage = mcov[row_odd == 0])
    return(msites)
}


#' calculate_expmeth
#' 
#'      This function is used to calculate genome-wide expected methylation 
#'  for each motif. 
#'
#' @param msites - imported methylation sites 
#' @param gcdist - Genome wide GC distribution from (\code{methylTFRann})
#' @param gcfreq - GC bin frequency table (matrix) from (\code{methylTFRann})
#' @return a \code{data.table} object with GC bin with corresponding avg methylation
#' @export 
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
calculate_expmeth <- function(msites, gcdist, gcfreq){
    hits <- findOverlaps(msites, gcdist, type = "within")
    gcmap <- data.table(mscore = msites[hits@from]$score, 
                        gcbin = gcdist[hits@to]$GC_bin)
    exp_meth <- gcmap %>% 
                group_by(gcbin) %>% 
                summarise(avg_mscore = mean(mscore))
    exp_meth <- as.matrix(exp_meth %>% 
                          group_by(gcbin) %>% 
                          summarise(mscore = mean(avg_mscore)/1000))
    exp.data <- t(gcfreq) %*% exp_meth[, 2]
    mpos <- seq(-floor(length(exp.data)/2), floor(length(exp.data)/2))
    exp.methyl <-  data.table(x = mpos, avg_methyl = exp.data)
    colnames(exp.methyl) <- c("x", "exp_avg_methyl")
    return(exp.methyl)
}


#' compute_deviation
#' 
#'      This function is function to calculate the deviation in transcription factor
#'  footprint base for a given motif
#'
#' @param motifs        - list of motifs 
#' @param msites       - imported methylation sites
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions from (\code{methylTFRann})
#' @param gcfreqs      - GC bin frequency tables (matrices for multiple motif) from (\code{methylTFRann})
#' @param gcdist       - Genome wide GC distribution from (\code{methylTFRann})
#' @param enhancer     - Specific regions like distal motif
#' @return deviation score for a given motif
#' @export 
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
compute_deviation <- function(motif, msites, tf_bindsites, gcfreqs,
                                gcdist, enhancer = NULL){
    tfbs <- tf_bindsites[[motif]]
    gcfreq <- gcfreqs[[motif]]
    w <- width(tfbs)[1]
    # This method is to increase the width as 500 bp around motif center
    tfbs <- resize(tfbs, w + 101, fix = "center")

    if (!is.null(enhancer)){
        d_hits <- findOverlaps(tfbs, enhancer)
        tfbs <- tfbs[d_hits@from]
    }
    # expected methylation
    exp_meth <- calculate_expmeth(msites, gcdist, gcfreq)

    hits <- findOverlaps(msites, tfbs, type = "within")

    mcols(tfbs)$mid_point <- round(end(tfbs) + ((start(tfbs) - end(tfbs))/2))
    x = start(msites[hits@from]) - tfbs[hits@to]$mid_point
    plot.data <- data.table(x = x, 
                            y1 = (msites[hits@from]$score/1000), 
                            y2 = msites[hits@from]$coverage)

    plot.data <- plot.data %>% group_by(x) %>% 
                    summarise(n=n(), avg_methyl=mean(y1), avg_cov = mean(y2),  
                              motif=motif)
    # GC bias correcton 
    dt_join <- left_join(plot.data, exp_meth, by='x')
    dt_join$diff <- abs(dt_join$avg_methyl - dt_join$exp_avg_methyl)

    interval_mean <- dt_join %>% mutate(cuts = cut(x, c(-200, -100, -10, 10, 100, 200))) %>% 
                        group_by(cuts) %>% 
                        summarize(n=n(), mean = mean(diff))
    
    var <- interval_mean$mean[3]/((interval_mean$mean[1] + interval_mean$mean[5])/2)
    return(var)  
}


#' run_methyltfr
#' 
#'      This function is a wrapper function to calculate the deviation in transcription factor
#'  footprint base for all given motifs using parallel package
#'
#' @param sample_ann   - a tab seperated file contains sample annotations
#' @param sample_dir   - directory where all bed file and annotation file stored
#' @param genome       - human genome version default: hg38
#' @param threads      - thread count for parallel processing
#' @param enhancer     - Specific regions like distal motif
#' @return deviation score matrix for all samples 
#' @export 
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
#' @importFrom parallel mclapply
#' @importFrom logger log_info log_error
run_methyltfr <- function(sample_ann, sample_dir, genome="hg38",
                            threads = 4, enhancer = NULL){
    
    if (!require("methylTFRann")) {
        log_error("methylTFRann package is not installed in your environment !!")
    }

    annfile = file.path(sample_dir, sample_ann)
    if (!file.exists(annfile)){
        log_error(" %s does not exist, please check the file path !!", annfile)
    }

    samples <- read.table(annfile, sep = "\t", header=TRUE)

    log_info("Loading the samples and annotation package ")
    files_list <- file.path(sample_dir, samples[, "bedFile"])
    tf_bindsites <- getTFbindsites()
    gc_dist <- getGenomeGC()
    gcfreqs <- getGCfreq()
    
    motifs <- names(gcfreqs)
    deviation <- data.table()
    for ( bedfile in files_list) {
        sample_name = basename(bedfile)
        log_info("Processing %s ", bedfile)
        msites <- read_methylome(bedfile)
        assign(sample_name, mclapply(motifs, compute_deviation, 
                                        msites = msites, 
                                        tf_bindsites = tf_bindsites, 
                                        gcfreqs = gcfreqs, 
                                        gcdist = gc_dist, mc.cores = threads))
        log_info('Done %s ', bedfile)
        deviation[sample_name] <- unlist(get(sample_name))
    }
    rownames(deviation) <- motifs
    
    return(deviation)
}