#' @title computeZScoreVariability
#' @description Function to compute overall Z-score variability across samples per motif set
#' Adapted from chromVAR::computeVariability
#' @param object methylTFRdeviations object, or z-score matrix
#' @param computeError if TRUE, it will compute bootstrap confidence interval
#' @param bootNsamples number of bootstrap samples to use, default is 1000
#' @param quantiles quantiles for bootstrap, default is 0.025 and 0.975
#' @param padjMethod method for p-value adjustment, default is "BH"
#' @param na.rm logical, if TRUE, remove NA values
#' @return data.frame with columns for name, variability, bootstrap lower bound,
#' bootstrap upper bound, raw p value, adjust p value.
#' @export
computeZScoreVariability <- function(object,
                                     computeError = TRUE,
                                     bootNsamples = 1000,
                                     quantiles = c(0.025, 0.975),
                                     padjMethod = "BH",
                                     na.rm = TRUE) {
  if (!any(class(object) %in% c("methylTFRDeviations", "matrix", "data.frame"))) {
    stop("object must be a methylTFRDeviations object or a Z-score matrix")
  }
  if (is(object, "methylTFRDeviations")) {
    object <- assays(object)$z
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

#' @title row_sds
#' @description Compute standard deviation across rows of a matrix
#' @param X matrix
#' @param na.rm logical, if TRUE, remove NA values
#' @return numeric vector of standard deviations
#' @keywords internal
row_sds <- function(X, na.rm = FALSE) {
  res <- numeric(nrow(X))
  if (na.rm) {
    for (j in 1:nrow(X)) {
      tmp <- X[j, , drop = FALSE]
      keep <- which(!is.na(tmp))

      if (length(keep) > 1) {
        res[j] <- sd(tmp[keep])
      } else {
        res[j] <- NaN
      }
    }
  } else {
    res <- apply(X, 1, sd, na.rm = TRUE)
  }

  return(res)
}

#' @title row_sds_perm
#' @description Compute standard deviation across rows of a matrix
#' @param X matrix
#' @param na.rm logical, if TRUE, remove NA values
#' @return numeric vector of standard deviations
#' @keywords internal
row_sds_perm <- function(X, na.rm = FALSE) {
  ix <- sample(seq_len(ncol(X)), ncol(X), replace = TRUE)
  shuffled <- X[, ix, drop = FALSE]
  return(row_sds(shuffled, na.rm))
}

#' @title quantile_helper
#' @description Helper function for computeZScoreVariability
#' @param values numeric vector
#' @param quantiles numeric vector of quantiles
#' @param na.rm logical, if TRUE, remove NA values
#' @return numeric vector of quantiles
#' @keywords internal
quantile_helper <- function(values, quantiles, na.rm) {
  if (na.rm) {
    out <- quantile(values, quantiles, na.rm = TRUE)
  } else {
    if (!all_false(is.na(values))) {
      out <- rep(NA, length(quantiles))
    } else {
      out <- quantile(values, quantiles)
    }
  }
  return(out)
}
