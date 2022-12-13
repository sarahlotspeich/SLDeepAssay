#' Bias correction for MLE from single-dilution assay data
#' @name BC_sd
#' @param lambda Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @return A vector of the same length as \code{lambda}
#'
BC_sd = function(lambda, M, q) {

  # fisher information (n x n), inverse (n x n), and stacked inverse (n^2 x 1)
  I = fisher_sd(lambda = lambda, M = M, q = q)
  I_inv = solve(I)
  vec_I_inv = as.vector(I_inv)

  # derivative of fisher information (n x n^2)
  dI = deriv_fisher_sd(lambda = lambda,
                       M = M,
                       q = q)

  # expectation of third derivative of log-likelihood (n x n^2)
  exp_d3 = exp_third_deriv_sd(lambda = lambda,
                              M = M,
                              q = q)

  # A matric (n x n^2)
  A = - dI - 0.5 * exp_d3

  # bias correction
  bc = I_inv %*% A %*% vec_I_inv

  ## Bias corrected MLE
  lambda_bc = lambda - bc

  return(lambda_bc)
}
