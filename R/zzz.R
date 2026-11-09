#' methylTFR: transcription factor activity from DNA methylation
#'
#' Scores the methylation footprint at transcription factor binding sites
#' and returns one bias-corrected score per motif and sample.
#'
#' @section Getting started:
#' 1. Import methylation calls with [readMethylome()].
#' 2. Load annotations from methylTFRAnnotationHg38 (`getTFbindsites()`,
#'    `getGCfreq()`, `getGenomeGC()`).
#' 3. Compute scores with [runMethylTFR()] or [runMethylTFRFromRnBeads()].
#' 4. Explore with [deviations()], [deviationZScores()],
#'    [computeZScoreVariability()], [differentialDeviationTest()],
#'    [plotMotifFootprint()] and [plotExpectedFootprint()].
#'
#' See `vignette("methylTFR")` for a worked example.
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
    ".", "V1", "V2", "V3", "V6", "avg", "avg_methyl", "cuts", "gcbin",
    "mscore", "type", "x", "y1", "y2"
))
