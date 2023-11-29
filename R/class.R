#' @title methylTFRdeviations
#' @description Class for storing results from \code{\link{run_methyltfr}} function.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @details This class inherits from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and most methods
#' for that class should work for objects of this class as well. Additionally,
#' two accessor functions are defined for extracting bias corrected deviations
#' (\code{\link{deviations}}) and deviation Z-scores
#' (\code{\link{deviationZScores}})
setClass("methylTFRdeviations", contains = "SummarizedExperiment")

# Set validity function
setValidity(
  "methylTFRdeviations",
  function(object) {
    if (!"deviations" %in% assayNames(object) || !"z" %in%
      assayNames(object)) {
      return("The assays slot must contain 'deviations' and 'z'")
    }
    return(TRUE)
  }
)

##########################################################################
########### Generic functions for methylTFRdeviations class ##############
##########################################################################

# Define generic function for deviations
setGeneric("deviations", function(x) standardGeneric("deviations"))

# Define method for methylTFRdeviations class
setMethod("deviations", "methylTFRdeviations", function(x) {
  assay(x, "deviations")
})

# Define generic function for deviationZScores
setGeneric("deviationZScores", function(x) standardGeneric("deviationZScores"))

# Define method for methylTFRdeviations class
setMethod("deviationZScores", "methylTFRdeviations", function(x) {
  assay(x, "z")
})

# Define generic function for rbind
setGeneric("rbind", function(...) standardGeneric("rbind"))

# Define generic function for cbind
setGeneric("cbind", function(...) standardGeneric("cbind"))

# Define rbind and cbind methods
setMethod("rbind", "methylTFRdeviations", function(..., deparse.level = 1) {
  inputs <- list(...)

  all_rowdata_colnames <- unique(do.call(c, lapply(inputs, function(x) colnames(rowData(x)))))
  in_common <- vapply(all_rowdata_colnames, function(x) {
    all(vapply(inputs, function(y) {
      x %in% colnames(rowData(y))
    }, rep(TRUE, length(x))))
  }, TRUE)
  common_colnames <- all_rowdata_colnames[in_common]
  inputs <- lapply(inputs, function(x) {
    rowData(x) <- rowData(x)[, common_colnames, drop = FALSE]
    x
  })

  deviations <- do.call("rbind", lapply(inputs, function(x) assay(x, "deviations")))
  z_scores <- do.call("rbind", lapply(inputs, function(x) assay(x, "z")))
  rd <- do.call("rbind", lapply(inputs, rowData))

  se <- SummarizedExperiment(
    assays = SimpleList(deviations = deviations, z = z_scores),
    colData = colData(inputs[[1]]),
    rowData = rd
  )

  as(se, "methylTFRdeviations")
})
