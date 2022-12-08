#' Expected value of the third derivative of the Fisher information of the log-likelihood for IUPM from multiple-dilution assay data
#' @name exp_third_deriv_md
#' @param tau Parameter
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @param dilutions Vector of all dilution levels (in millions of cells per well).
#' @return A matrix
#'
#'
exp_third_deriv_md = function(tau, M, q, dilutions) {
  D = length(dilutions)
  n = length(tau)
  ed3.d = vapply(X = 1:D, FUN.VALUE = matrix(0, n, n^2),
                 FUN = function(d) exp_third_deriv_sd(lambda = dilutions[d] * tau,
                                                      M = M[d],
                                                      q = q[d]) *
                   dilutions[d] ^ 3)

  ed3 = apply(ed3.d, c(1,2), sum)
  return(ed3)
}
