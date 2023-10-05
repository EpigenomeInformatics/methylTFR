#' methylTFRDeviations
#'
#' Class for storing results from \code{\link{computeDeviations}} function.
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @details This class inherits from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and methods
#' for that SummarizedExperiment work for methylTFRDeviations objects as well. Additionally,
#' two accessor functions are defined for extracting bias-corrected deviations
#' (\code{\link{deviations}}) and deviation Z-scores
#' (\code{\link{deviationScores}}).
#' @slot deviations: A slot in the SummarizedExperiment containing bias-corrected deviations.
#' @slot z: A slot in the SummarizedExperiment containing deviation Z-scores.
#' @rdname methylTFRDeviations-class
#' @exportClass methylTFRDeviations
setClass("methylTFRDeviations", contains = "SummarizedExperiment")

setValidity(
  "methylTFRDeviations",
  function(object) {
    if (!all(c("deviations", "z") %in% assayNames(object))) {
      return("The assays slot must contain 'deviations' and 'z'")
    }
    return(TRUE)
  }
)
