

# wilcoxon test for two groups
wilcoxon_helper <- function(x, groups) {
    splitx <- split(t(x), groups)
    return(wilcox.test(splitx[[1]], splitx[[2]],
                        paired = FALSE)$p.value)
}

