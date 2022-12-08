
#' compute_fp
#' 
#'  Compute footprint returns a data.table required to create a transcription
#' factor footprint with label.
#'
#' @param motif_name   - motif name eg: GATA string   
#' @param tf_bindsites - a GenomicRange object contains tf binding sites positions from (\code{methylTFRann})
#' @param msites       - methylation data processed from \code{RnBeads}
#' @param enhancer     - Specific regions such as distal motif, proximal motif
#' @return a \code{data.table} object to plot tf footprint
#' @export 
#' @importFrom GenomicRanges GRanges findOverlaps width resize start end
#' @importFrom data.table data.table
#' @importFrom dplyr %>% mutate
compute_fp <- function(motif_name, tf_bindsites, msites,
                                enhancer = NULL){
    tfbs <- tf_bindsites[[motif_name]]
    w <- width(tfbs)[1]
    tfbs <- resize(tfbs, w + 101, fix = "center")

    if (!is.null(enhancer)){
        d_hits <- findOverlaps(tfbs, enhancer)
        tfbs <- tfbs[d_hits@from]
    }
    
    hits <- findOverlaps(msites, tfbs, type = "within")
    
    mcols(tfbs)$mid_point <- round(end(tfbs) + (start(tfbs) - end(tfbs))/2)
    x = start(msites[hits@from]) - tfbs[hits@to]$mid_point
    plot.data <- data.table(x = x, 
                            y1 = msites[hits@from]$score/1000, 
                            y2 = msites[hits@from]$coverage)

    plot.data <- plot.data %>% group_by(x) %>% 
                    summarise(n=n(), avg_methyl=mean(y1), avg_cov = mean(y2),  
                              motif=motif_name) %>%
                    mutate(label = if_else(x == max(x), as.character(motif), NA_character_))
    
    return(plot.data)
}


#' plot_footprint
#' 
#'     Creates a footprint plot for given motifs and methylation 
#' site. It will create png files in the specified location.
#' 
#' @param motifs       - motif names as character vector
#' @param tf_bindsites - Transcript Factor binbing sites from (\code{methylTFRann})
#' @param samples       - Import methylation data 
#' @param img        - output plot png filename
#' @return image will be generated in the specified path
#' @export 
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot ggsave geom_point geom_line ggtitle
#' @importFrom ggrepel geom_label_repel
plot_footprint <- function(motifs, tf_bindsites, samples, 
                           img ="TF_footprint.png"){
    
    msites_list <- GRangesList()
    msites_list <- lapply(samples, read_methylome)
    sample_names <- unlist(lapply(samples, basename))
    names(msites_list) <- sample_names    

    plot_data <- data.table()
    
    for (i in 1:length(msites_list)){
        msites = msites_list[[i]]
        current_plot <- do.call("rbind", lapply(motifs, compute_fp, tf_bindsites, msites))
        current_plot$sample_name <- sample_names[i]
        plot_data <- rbind(plot_data, current_plot)
    }
    
    p1 <- plot_data %>%
        ggplot(aes(x = x, y = avg_methyl, group=interaction(motif, sample_name))) +
                geom_line(aes(color=sample_name)) +
                geom_point(aes(color=sample_name)) +
                geom_label_repel(aes(label = label),
                    nudge_x = 1,
                    na.rm = TRUE) +
                xlab("Distance from motif center") +
                ylab("Methylation level") +
                ggtitle("TF footprint")

    ggsave(filename = img,
            plot = p1,
            width = 11, height = 8.5)
}