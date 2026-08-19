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
#' @examples
#' # Load example data
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' # Get deviations
#' deviations(tc_mem)
setGeneric("deviations", function(x) standardGeneric("deviations"))

#' @title deviations
#' @description Extract bias corrected deviations from methylTFRdeviations object.
#' @param x methylTFRdeviations object.
#' @return A matrix of bias corrected deviations.
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom methods setMethod
#' @examples
#' # Load the data
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' # Get deviations
#' deviations(tc_mem)
setMethod("deviations", "methylTFRdeviations", function(x) {
    SummarizedExperiment::assay(x, "deviations")
})

#' @title deviationZScores
#' @description Function to get deviation Z-scores from a methylTFRdeviations object.
#' @param x A methylTFRdeviations object.
#' @return A matrix of deviation Z-scores.
#' @export
#' @examples
#' # Load example data
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' # Get deviation Z-scores
#' deviationZScores(tc_mem)
setGeneric("deviationZScores", function(x) standardGeneric("deviationZScores"))

#' @title deviationZScores
#' @description Extract deviation Z-scores from methylTFRdeviations object.
#' @param x methylTFRdeviations object.
#' @return A matrix of deviation Z-scores.
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom methods setMethod
#' @examples
#' # Load the data
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' # Get deviation Z-scores
#' deviationZScores(tc_mem)
setMethod("deviationZScores", "methylTFRdeviations", function(x) {
    SummarizedExperiment::assay(x, "z")
})

# cbind() and rbind() are deliberately not defined here.
#
# methylTFRdeviations extends SummarizedExperiment, which already
# provides both, and those methods return an object of the subclass with
# all assays intact. Defining setGeneric("cbind", function(x, y)) shadowed
# the base generic, which made R CMD INSTALL emit
#   Creating a new generic function for 'cbind' in package 'methylTFR'
# and restricted the call to exactly two arguments. Relying on the
# inherited methods removes the warning and allows cbind(a, b, c).
#
# This requires SummarizedExperiment to be attached rather than merely
# imported: the bindCOLS and bindROWS methods it defines are only
# reachable for a subclass when its namespace is on the search path.
# SummarizedExperiment is therefore in Depends, as it is for other
# packages whose central class extends it.
