#' @title methylTFRdeviations
#' @description Class for storing results from \code{\link{run_methyltfr}} function.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @details This class inherits from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and most methods
#' for that class should work for objects of this class as well. Additionally,
#' two accessor functions are defined for extracting bias corrected deviations
#' (\code{\link[=deviations,methylTFRdeviations-method]{deviations}}) and deviation Z-scores
#' (\code{\link[=deviationZScores,methylTFRdeviations-method]{deviationZScores}})
#' @name methylTFRdeviations-class
#' @rdname methylTFRdeviations-class
#' @exportClass methylTFRdeviations
#' @importFrom SummarizedExperiment SummarizedExperiment
setClass("methylTFRdeviations", contains = "SummarizedExperiment")

# Set validity function
setValidity(
  "methylTFRdeviations",
  function(object) {
    if (!"deviations" %in% SummarizedExperiment::assayNames(object) || !"z" %in%
      SummarizedExperiment::assayNames(object)) {
      return("The assays slot must contain 'deviations' and 'z'")
    }
    return(TRUE)
  }
)

##########################################################################
########### Generic functions for methylTFRdeviations class ##############
##########################################################################

#' @title deviations
#' @description Function to get deviations from a methylTFRdeviations object.
#' @param x A methylTFRdeviations object.
#' @return A matrix of deviations.
#' @export
setGeneric("deviations", function(x) standardGeneric("deviations"))

#' @title deviations
#' @description Extract bias corrected deviations from methylTFRdeviations object.
#' @param x methylTFRdeviations object.
#' @return A matrix of bias corrected deviations.
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom methods setMethod
# Define method for methylTFRdeviations class
setMethod("deviations", "methylTFRdeviations", function(x) {
  SummarizedExperiment::assay(x, "deviations")
})

#' @title deviationZScores
#' @description Function to get deviation Z-scores from a methylTFRdeviations object.
#' @param x A methylTFRdeviations object.
#' @return A matrix of deviation Z-scores.
#' @export
setGeneric("deviationZScores", function(x) standardGeneric("deviationZScores"))

#' @title deviationZScores
#' @description Extract deviation Z-scores from methylTFRdeviations object.
#' @param x methylTFRdeviations object.
#' @return A matrix of deviation Z-scores.
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom methods setMethod
setMethod("deviationZScores", "methylTFRdeviations", function(x) {
  SummarizedExperiment::assay(x, "z")
})

#' @title rbind
#' @description Combine methylTFRdeviations objects by row.
#' @param ... methylTFRdeviations objects.
#' @return A methylTFRdeviations object.
#' @export
setGeneric("rbind", function(...) standardGeneric("rbind"))

#' @title cbind
#' @description Combine methylTFRdeviations objects by column.
#' @param ... methylTFRdeviations objects.
#' @return A methylTFRdeviations object.
#' @export
setGeneric("cbind", function(...) standardGeneric("cbind"))

#' @title rbind
#' @description Combine methylTFRdeviations objects by row.
#' @param ... methylTFRdeviations objects.
#' @param deparse.level An integer specifying the level of deparse to use.
#' @return A methylTFRdeviations object.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList
#' @importFrom methods setMethod as
#'
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

  se <- SummarizedExperiment::aSummarizedExperiment(
    assays = SimpleList(deviations = deviations, z = z_scores),
    colData = colData(inputs[[1]]),
    rowData = rd
  )

  return(as(se, "methylTFRdeviations"))
})
