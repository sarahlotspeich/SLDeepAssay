#' Log-likelihood for IUPM from multiple-dilution assay data
#' @name loglik_md
#' @param tau Parameter
#' @param assay.summary Data
#' @return A scalar
#'
loglik_md = function(tau, assay.summary) {
  # Total number of dilutions
  D = nrow(assay.summary)

  # Compute negative log-likelihoods for each dilution
  logliks_by_dilution <- vapply(X = 1:D,
                                FUN.VALUE = numeric(1),
                                FUN = function(d) {
                                  as.numeric(loglik_sd(l = assay.summary$dilution[d] * tau,
                                                       n = assay.summary$n[d],
                                                       M = assay.summary$M[d],
                                                       MP = assay.summary$MP[d],
                                                       m = assay.summary$m[d],
                                                       Y = assay.summary[d, grepl("Y", colnames(assay.summary))]))
                                })

  # Return the sum of (negative) log-likelihoods
  ll <- sum(logliks_by_dilution)

  return(ll)
}
