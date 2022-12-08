#' Derivative of the Fisher information of the log-likelihood for IUPM from multiple-dilution assay data
#' @name deriv_fisher_md
#' @param tau Parameter
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @param dilutions Vector of all dilution levels (in millions of cells per well).
#' @return A matrix
#'
deriv_fisher_md = function(tau, M, q, dilutions) {
  D = length(dilutions)
  n = length(tau)

  dI.d <- vapply(X = 1:D, FUN.VALUE = matrix(0, n, n^2),
                 FUN = function(d) deriv_fisher_sd(lambda = dilutions[d] * tau,
                                                M = M[d],
                                                q = q[d]) *
                   dilutions[d] ^ 3)

  # derivative of fisher information over all dilution levels
  dI <- apply(X = dI.d, MARGIN = c(1, 2), FUN = sum)

  return(dI)
}
