#' @title computeVariability
#' @description Function to compute overall deviation score variability across samples per motif set
#' Adapted from chromVAR::computeVariability
#' @param object methylTFRdeviations object, or deviation score matrix
#' @param computeError if TRUE, it will compute bootstrap confidence interval
#' @param bootNsamples number of bootstrap samples to use, default is 1000
#' @param quantiles quantiles for bootstrap, default is 0.025 and 0.975
#' @param padjMethod method for p-value adjustment, default is "BH"
#' @param na.rm logical, if TRUE, remove NA values
#' @return data.frame with columns for name, variability, bootstrap lower bound,
#' bootstrap upper bound, raw p value, adjust p value.
#' @export
computeVariability <- function(object,
                               computeError = TRUE,
                               bootNsamples = 1000,
                               quantiles = c(0.025, 0.975),
                               padjMethod = "BH",
                               na.rm = TRUE) {
  if (!any(class(object) %in% c("methylTFRDeviations", "matrix", "data.frame"))) {
    stop("object must be a methylTFRDeviations object or adeviation score matrix")
  }
  if (is(object, "methylTFRDeviations")) {
    object <- assays(object)$deviations
  }
  if (!is.logical(computeError)) {
    stop("computeError must be TRUE or FALSE")
  }
  if (!is.numeric(bootNsamples) || bootNsamples < 1 || length(bootNsamples) != 1) {
    stop("bootNsamples must be a numeric value")
  }
  if (bootNsamples %% 1 != 0) {
    stop("bootNsamples must be a whole number")
  }
  if (!is.numeric(quantiles) || length(quantiles) != 2) {
    stop("quantiles must be a numeric vector of length 2")
  }
  if (!all(quantiles > 0) || !all(quantiles < 1)) {
    stop("quantiles must be between 0 and 1")
  }
  if (!quantiles[2] > quantiles[1]) {
    stop("quantiles[2] must be greater than quantiles[1]")
  }
  if (!is.character(padjMethod) || length(padjMethod) != 1) {
    stop("padjMethod must be a character value")
  }
  if (!is.logical(na.rm)) {
    message("Invalid na.rm parameter, switching to default")
    na.rm <- TRUE
  }
  motif_names <- rownames(object)
  # Compute standard deviation across samples
  sd_devs <- row_sds(object, na.rm = TRUE)
  # Simulate null distribution of standard deviation
  pval_sim <- pchisq((ncol(object) - 1) * (sd_devs^2),
    df = (ncol(object) - 1),
    lower.tail = FALSE
  )
  # Adjust p-values
  padj <- p.adjust(p = pval_sim, method = padjMethod)

  if (computeError) {
    # Compute bootstrap confidence interval
    boot_sd <- t(replicate(bootNsamples, row_sds_perm(object, na.rm = na.rm)))

    sd_error <- apply(boot_sd, 2, quantile_helper,
      quantiles = quantiles,
      na.rm = na.rm
    )

    out <- data.frame(
      name = motif_names,
      variability = sd_devs,
      bootstrap_lower_bound = sd_error[1, ],
      bootstrap_upper_bound = sd_error[2, ],
      p_value = pval_sim,
      p_value_adj = padj,
      row.names = rownames(object)
    )
  } else {
    out <- data.frame(
      name = motif_names,
      variability = sd_devs,
      p_value = pval_sim,
      p_value_adj = padj,
      row.names = rownames(object)
    )
  }

  return(out)
}
