#' @title calibrateDeviations
#' @description Standardise bias-corrected deviation scores against a
#' within-sample null, so that the resulting scores are comparable across
#' motifs and approximately standard normal in the absence of TF-specific
#' signal.
#' @details
#' Each column (sample) is centred and scaled independently, using the
#' distribution of deviation scores across all motifs in that sample as the
#' null. With several hundred motifs, the great majority of which carry no
#' strong methylation signal in any single sample, this distribution is a
#' usable estimate of the sample's noise level.
#'
#' This is deliberately different from \code{computeRowZScore}, which
#' standardises each motif against its own spread across samples. Row-wise
#' Z-scores are the right transform for visualising a heatmap, but they force
#' the standard deviation of every row to one and therefore cannot be used to
#' rank motifs by variability.
#' @param devs A numeric matrix of bias-corrected deviation scores, with motifs
#' in rows and samples in columns.
#' @param method Either \code{"robust"} (median and median absolute deviation,
#' the default) or \code{"gaussian"} (mean and standard deviation).
#' @return A numeric matrix of the same dimensions as \code{devs}.
#' @importFrom matrixStats colMedians colMads colMeans2 colSds
#' @keywords internal
calibrateDeviations <- function(devs, method = c("robust", "gaussian")) {
  method <- match.arg(method)
  devs <- as.matrix(devs)
  if (method == "robust") {
    centre <- matrixStats::colMedians(devs, na.rm = TRUE)
    scale <- matrixStats::colMads(devs, na.rm = TRUE)
  } else {
    centre <- matrixStats::colMeans2(devs, na.rm = TRUE)
    scale <- matrixStats::colSds(devs, na.rm = TRUE)
  }
  bad <- !is.finite(scale) | scale <= 0
  if (any(bad)) {
    if (all(bad)) {
      stop(
        "The null scale could not be estimated for any sample. ",
        "Check that the deviation matrix is not constant."
      )
    }
    warning(
      sum(bad), " sample(s) had a degenerate null scale; ",
      "substituting the median scale of the remaining samples."
    )
    scale[bad] <- stats::median(scale[!bad])
  }
  z <- sweep(devs, 2, centre, FUN = "-")
  z <- sweep(z, 2, scale, FUN = "/")
  return(z)
}


#' @title computeZScoreVariability
#' @description Identify transcription factor motifs whose methylTFR activity
#' varies across samples more than expected by chance.
#' @details
#' Variability is defined as the standard deviation, across samples, of
#' calibrated deviation scores. Calibration is what makes the number
#' interpretable: deviation scores are first standardised against a
#' within-sample null estimated across all motifs (see
#' \code{calibrateDeviations}), so that a calibrated score of 1 corresponds to
#' one unit of sample-level noise. A motif with variability greater than 1
#' therefore varies across samples by more than the background spread, and
#' variability near 1 is what a motif carrying no signal is expected to show.
#'
#' Under the null hypothesis that a motif's calibrated scores are independent
#' draws from a standard normal distribution, \eqn{(n - 1) s^2} follows a
#' chi-squared distribution with \eqn{n - 1} degrees of freedom, where \eqn{s}
#' is the observed variability and \eqn{n} the number of samples with a
#' non-missing score for that motif. P-values are obtained from the upper tail
#' of that distribution and adjusted for multiple testing across motifs.
#'
#' Note that this function must not be applied to the \code{z} assay returned
#' by \code{\link{deviationZScores}}. Those are row-wise Z-scores computed for
#' visualisation, and their standard deviation across samples is 1 for every
#' motif by construction, which would make the test degenerate. The default
#' input is therefore the \code{deviations} assay.
#'
#' The approach is analogous to \code{chromVAR::computeVariability}, with one
#' difference: chromVAR calibrates each motif against a set of GC- and
#' width-matched background peak sets, whereas here the null is estimated
#' across motifs within each sample. The latter needs no additional deviation
#' computations and can be applied to existing results, at the cost of
#' assuming that most motifs are inactive in any given sample.
#' @param object A \code{methylTFRdeviations} object, or a numeric matrix or
#' data.frame of bias-corrected deviation scores with motifs in rows and
#' samples in columns.
#' @param method Either \code{"robust"} (median and median absolute deviation,
#' the default) or \code{"gaussian"} (mean and standard deviation), used to
#' estimate the within-sample null.
#' @param bootstrap logical, if TRUE bootstrap confidence bounds for the
#' variability estimates are computed by resampling samples with replacement.
#' @param niterations Number of bootstrap iterations, used only when
#' \code{bootstrap} is TRUE.
#' @param conf_level Confidence level for the bootstrap bounds.
#' @param padjMethod Method for p-value adjustment, passed to
#' \code{\link[stats]{p.adjust}}. Default is "BH".
#' @return A \code{data.frame} with one row per motif and the columns
#' \code{motifs}, \code{variability}, \code{p_value} and
#' \code{p_value_adjusted}, plus \code{bootstrap_lower_bound} and
#' \code{bootstrap_upper_bound} when \code{bootstrap} is TRUE. Rows are
#' returned in the order of the input.
#' @importFrom matrixStats rowSds
#' @importFrom stats pchisq p.adjust median quantile
#' @importFrom methods is
#' @seealso \code{\link{differential_deviation_test}} for testing activity
#' differences between predefined groups.
#' @examples
#' # Load example data
#' load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
#' load(system.file("extdata", "tc_naive.rda", package = "methylTFR"))
#' devs <- cbind(tc_mem, tc_naive)
#'
#' # Rank motifs by how much their activity varies across samples
#' var_res <- computeZScoreVariability(devs)
#' head(var_res[order(-var_res$variability), ])
#' @author Irem Gunduz
#' @export
computeZScoreVariability <- function(
  object,
  method = c("robust", "gaussian"),
  bootstrap = FALSE,
  niterations = 1000L,
  conf_level = 0.95,
  padjMethod = "BH"
) {
  method <- match.arg(method)
  if (is(object, "methylTFRdeviations")) {
    devs <- deviations(object)
  } else if (is.matrix(object) || is.data.frame(object)) {
    devs <- as.matrix(object)
  } else {
    stop(
      "object must be a methylTFRdeviations object, ",
      "a matrix or a data.frame"
    )
  }
  if (!is.numeric(devs)) {
    stop("The deviation scores must be numeric")
  }
  if (ncol(devs) < 3) {
    stop(
      "At least three samples are required to estimate variability; ",
      "found ", ncol(devs), "."
    )
  }
  if (!is.logical(bootstrap) || length(bootstrap) != 1) {
    stop("bootstrap must be a single logical value")
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a number between 0 and 1")
  }

  if (nrow(devs) < 50) {
    message(
      "Only ", nrow(devs), " motifs supplied. The within-sample null is ",
      "estimated across motifs, so variability estimates from small ",
      "motif sets should be treated as indicative only."
    )
  }

  motif_names <- rownames(devs)
  if (is.null(motif_names)) {
    motif_names <- paste0("motif_", seq_len(nrow(devs)))
  }

  z <- calibrateDeviations(devs, method = method)
  variability <- matrixStats::rowSds(z, na.rm = TRUE)
  n_obs <- rowSums(!is.na(z))

  p_value <- rep(NA_real_, length(variability))
  testable <- is.finite(variability) & n_obs > 1
  p_value[testable] <- stats::pchisq(
    variability[testable]^2 * (n_obs[testable] - 1),
    df = n_obs[testable] - 1,
    lower.tail = FALSE
  )
  p_adj <- stats::p.adjust(p_value, method = padjMethod)

  res <- data.frame(
    motifs = motif_names,
    variability = variability,
    p_value = p_value,
    p_value_adjusted = p_adj,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (bootstrap) {
    if (!is.numeric(niterations) || niterations < 2) {
      stop("niterations must be a number greater than 1")
    }
    niterations <- as.integer(niterations)
    boot <- vapply(seq_len(niterations), function(k) {
      idx <- sample.int(ncol(z), ncol(z), replace = TRUE)
      matrixStats::rowSds(z[, idx, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(z)))
    alpha <- (1 - conf_level) / 2
    res$bootstrap_lower_bound <- apply(boot, 1, stats::quantile,
      probs = alpha, na.rm = TRUE
    )
    res$bootstrap_upper_bound <- apply(boot, 1, stats::quantile,
      probs = 1 - alpha, na.rm = TRUE
    )
  }

  return(res)
}
