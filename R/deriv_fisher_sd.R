#' Derivative of the Fisher information of the log-likelihood for IUPM from single-dilution assay data
#' @name deriv_fisher_sd
#' @param lambda Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @return A matrix
#'
deriv_fisher_sd = function(lambda, M, q) {
  n = length(lambda)
  dI = vapply(X = 1:n, FUN.VALUE = matrix(0, n, n),
              FUN = function(i) {
                numDeriv::jacobian(func = function(l) {
                  lam = lambda
                  lam[i] = l
                  fisher_sd(lambda = lam, M = M, q = q)
                },
                x = lambda[i],
                method = "simple")
              })
  dI = matrix(dI, n, n^2)
  return(dI)
}