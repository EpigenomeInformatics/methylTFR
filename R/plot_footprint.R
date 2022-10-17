
#' compute_fp
#' 
#'  Compute footprint returns a data.table required to create a transcription
#' factor footprint with label.
#'
#' @param motif_name  - motif_name eg GATA  
#' @param tf_bindsites - Transcription Factor binding sites as GenomicRange List object
#' @param msites  - Imported methylation data 
#' @param enhancer - - GenomicRange object containes specific regions like distal motif
#' @return tf footprint plot data (\code{data.table} object)
#' @export 
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
#' @param motifs -   motif names as character vector
#' @param tf_bindsites - Transcript Factor binbing sites from (\code{methylTFRann})   
#' @param msites - Import methylation data 
#' @param fname - output plot png filename
#' @return 
#' @export 
plot_footprint <- function(motifs, tf_bindsites, msites, 
                           fname="TF_footprint.png"){
    plot.data <- do.call("rbind", lapply(motifs, compute_fp, tf_bindsites, msites))
    p1 <- plot_data %>%
        ggplot(aes(x = x, y = avg_methyl, group=motif)) +
                geom_line(aes(color=motif)) +
                geom_point(aes(color=motif)) +
                geom_label_repel(aes(label = label),
                    nudge_x = 1,
                    na.rm = TRUE) +
                xlab("Distance from motif center") +
                ylab("Methylation level") +
                ggtitle("TF footprint")

    ggsave(filename = fname,
            plot = p1,
            width = 11, height = 8.5)
}