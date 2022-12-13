#' Log-likelihood for IUPM from multiple-dilution assay data
#' @name loglik_md
#' @param tau Vector of DVL-specific IUPM parameters
#' @param assay_summary Summary of assay data from multiple dilutions
#' @return A scalar value of the log-likelihood
#'
loglik_md = function(tau, assay_summary) {
  # Total number of dilutions
  D = nrow(assay_summary)

  # Compute negative log-likelihoods for each dilution
  logliks_by_dilution <- vapply(X = 1:D,
                                FUN.VALUE = numeric(1),
                                FUN = function(d) {
                                  as.numeric(loglik_sd(l = assay_summary$dilution[d] * tau,
                                                       n = assay_summary$n[d],
                                                       M = assay_summary$M[d],
                                                       MP = assay_summary$MP[d],
                                                       m = assay_summary$m[d],
                                                       Y = assay_summary[d, grepl("Y", colnames(assay_summary))]))
                                })

  # Return the sum of (negative) log-likelihoods
  ll <- sum(logliks_by_dilution)

  return(ll)
}
