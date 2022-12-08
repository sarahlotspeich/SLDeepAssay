#' Gradient of the log-likelihood for IUPM from multiple-dilution assay data
#' @name gloglik_md
#' @param tau Parameter
#' @param assay.summary Data
#' @return A vector (of length \code{n})
#'
gloglik_md = function(tau, assay.summary) {

  D = nrow(assay.summary)
  n = length(tau)

  # compute gradients for each dilution
  glogliks_byDilution <- vapply(X = 1:D, FUN.VALUE = numeric(n),
                                FUN = function(d) {
                                  as.numeric(gloglik_sd(l = assay.summary$dilution[d] * tau,
                                                        n = assay.summary$n[d],
                                                        M = assay.summary$M[d],
                                                        MP = assay.summary$MP[d],
                                                        m = assay.summary$m[d],
                                                        Y = as.numeric(assay.summary[d, grepl("Y", colnames(assay.summary))]))) * assay.summary$dilution[d]
                                  }
                                )

  # return the sum of (negative) log-likelihoods
  gradient = rowSums(glogliks_byDilution)

  return(gradient)
}
