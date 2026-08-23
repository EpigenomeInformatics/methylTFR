#' @title methylTFRdeviations
#' @description Class for storing results from
#' \code{\link{run_methyltfr}} function.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @details This class inherits from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and most methods
#' for that class should work for objects of this class as well. Additionally,
#' two accessor functions are defined for extracting bias corrected
#' deviations (\code{\link[=deviations,methylTFRdeviations-method]{deviations}})
#' and deviation Z-scores
#' (\code{\link[=deviationZScores,methylTFRdeviations-method]{deviationZScores}}).
#' @name methylTFRdeviations-class
#' @rdname methylTFRdeviations-class
#' @exportClass methylTFRdeviations
#' @importFrom SummarizedExperiment SummarizedExperiment
setClass("methylTFRdeviations", contains = "SummarizedExperiment")

# Set validity function
setValidity(
  "methylTFRdeviations",
  function(object) {
    req_assays <- c("deviations", "z")
    se_assays <- SummarizedExperiment::assayNames(object)
    if (any(!req_assays %in% se_assays)) {
      return("The assays slot must contain 'deviations' and 'z'")
    }
    return(TRUE)
  }
)

##########################################################################
########### Generic functions for methylTFRdeviations class ##############
##########################################################################

#' @title deviations
#' @description Function to get deviations from a
#' methylTFRdeviations object.
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
#' @description Extract bias corrected deviations from
#' methylTFRdeviations object.
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
#' @description Function to get deviation Z-scores from a
#' methylTFRdeviations object.
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
#' @description Extract deviation Z-scores from
#' methylTFRdeviations object.
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

#' @title cbind
#' @description Combine methylTFRdeviations objects by column.
#' @details
#' Defined explicitly rather than left to inheritance. Although
#' \code{methylTFRdeviations} extends \code{SummarizedExperiment},
#' whether that class's \code{cbind} method is found for a subclass
#' defined in another package depends on the search path and on the
#' Bioconductor version: dispatch can fall through to the
#' \code{S4Vectors} method, which then fails with "unable to find an
#' inherited method for function 'bindCOLS'".
#'
#' The objects are therefore coerced to \code{SummarizedExperiment}
#' explicitly, combined there, and re-wrapped. This is defined with
#' \code{setMethod} on the existing generic rather than
#' \code{setGeneric}, so it does not shadow the base function and any
#' number of objects can be combined.
#' @param ... methylTFRdeviations objects.
#' @param deparse.level passed to the underlying method.
#' @return A methylTFRdeviations object.
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom methods as new setMethod
#' @examples
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
#' devs <- cbind(tc_mem, tc_naive)
setMethod("cbind", "methylTFRdeviations", function(..., deparse.level = 1) {
  args <- list(...)
  if (length(args) == 1L) {
    return(args[[1L]])
  }
  combined <- do.call(
    BiocGenerics::cbind,
    lapply(args, function(x) as(x, "SummarizedExperiment"))
  )
  new("methylTFRdeviations", combined)
})


#' @title rbind
#' @description Combine methylTFRdeviations objects by row.
#' @details See \code{\link{cbind}} for why this is defined explicitly
#' rather than inherited.
#' @param ... methylTFRdeviations objects.
#' @param deparse.level passed to the underlying method.
#' @return A methylTFRdeviations object.
#' @export
#' @importFrom BiocGenerics rbind
#' @importFrom methods as new setMethod
#' @examples
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' devs <- rbind(tc_mem[seq_len(2), ], tc_mem[3:4, ])
setMethod("rbind", "methylTFRdeviations", function(..., deparse.level = 1) {
  args <- list(...)
  if (length(args) == 1L) {
    return(args[[1L]])
  }
  combined <- do.call(
    BiocGenerics::rbind,
    lapply(args, function(x) as(x, "SummarizedExperiment"))
  )
  new("methylTFRdeviations", combined)
})
