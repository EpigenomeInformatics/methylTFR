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
        if (any(!c("deviations", "z") %in% SummarizedExperiment::assayNames(object))) {
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

#' @title cbind
#' @description Combine two methylTFRdeviations objects by column.
#' @param x methylTFRdeviations object
#' @param y methylTFRdeviations object
#' @return A methylTFRdeviations object.
#' @export
setGeneric("cbind", function(x, y) standardGeneric("cbind"))

#' @title rbind
#' @description Combine two methylTFRdeviations objects by row.
#' @param x methylTFRdeviations object
#' @param y methylTFRdeviations object
#' @return A methylTFRdeviations object.
#' @export
setGeneric("rbind", function(x, y) standardGeneric("rbind"))

#' @title cbind
#' @description Combine two methylTFRdeviations objects by column.
#' @param x methylTFRdeviations object
#' @param y methylTFRdeviations object
#' @return A methylTFRdeviations object.
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom methods setMethod
setMethod(
    "cbind", signature(x = "methylTFRdeviations", y = "methylTFRdeviations"),
    function(x, y) {
        BiocGenerics::cbind(x, y)
    }
)

#' @title rbind
#' @description Combine two methylTFRdeviations objects by row.
#' @param x methylTFRdeviations object
#' @param y methylTFRdeviations object
#' @return A methylTFRdeviations object.
#' @export
#' @importFrom BiocGenerics rbind
#' @importFrom methods setMethod
setMethod(
    "rbind", signature(x = "methylTFRdeviations", y = "methylTFRdeviations"),
    function(x, y) {
        BiocGenerics::rbind(x, y)
    }
)
