#' Bias correction for MLE from multiple-dilution assay data
#' @name BC_md
#' @param tau Parameter
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @param dilutions Vector of all dilution levels (in millions of cells per well).
#' @return A vector of the same length as \code{lambda}
#'
BC_md = function(tau, M, q, dilutions) {
  # fisher information (n x n), inverse (n x n), and stacked inverse (n^2 x 1)
  I = fisher_md(tau = tau,
                M = M,
                q = q,
                dilutions = dilutions)
  I.inv = solve(I)
  vec.I.inv = as.vector(I.inv)

  # derivative of fisher information (n x n^2)
  dI = deriv_fisher_md(tau = tau,
                       M = M,
                       q = q,
                       dilutions = dilutions)

  # expectation of third derivative of log-likelihood (n x n^2)
  exp.d3 = exp_third_deriv_md(tau = tau,
                              M = M,
                              q = q,
                              dilutions = dilutions)

  # A matric (n x n^2)
  A = - dI - 0.5 * exp.d3

  # bias correction
  bc = I.inv %*% A %*% vec.I.inv

  ## Bias corrected MLE
  tau.bc = tau - bc

  return(tau.bc)
}
