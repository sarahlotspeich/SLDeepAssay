#' Gradient of the log-likelihood for IUPM from multiple-dilution assay data
#' @name gloglik_md
#' @param tau Vector of DVL-specific IUPM parameters
#' @param assay_summary Summary of assay data from multiple dilutions
#' @return A vector (of length \code{n})
#'
gloglik_md = function(tau, assay_summary) {

  D = nrow(assay_summary)
  n = length(tau)

  # compute gradients for each dilution
  glogliks_byDilution = vapply(X = 1:D, FUN.VALUE = numeric(n),
                                FUN = function(d) {
                                  as.numeric(gloglik_sd(l = assay_summary$u[d] * tau,
                                                        M = assay_summary$M[d],
                                                        MP = assay_summary$MP[d],
                                                        m = assay_summary$m[d],
                                                        Y = as.numeric(assay_summary[d, grepl("Y", colnames(assay_summary))]))) * assay_summary$u[d]
                                  }
                                )

  # return the sum of (negative) log-likelihoods
  gradient = rowSums(glogliks_byDilution)

  return(gradient)
}
