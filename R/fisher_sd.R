#' Fisher information of the log-likelihood for IUPM from single-dilution assay data
#' @name fisher_sd
#' @param lambda Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @return A matrix
#'
fisher_sd = function(lambda, M, q) {
  # Save parameters 
  n = length(lambda)
  Lambda = sum(lambda)
  
  # Construct n x n information matrix
  I = matrix(nrow = n, ncol = n,
             data = M * (1 - q) / (exp(Lambda) - 1))
  
  if (q > 0) {
    E.Y = get_EY(M = M, q = q, lambda = lambda)
    diag.terms = E.Y * exp(lambda) / (exp(lambda) - 1) ^ 2
    if (n == 1) {I = I + diag.terms}
    else {I = I + diag(diag.terms)}
  }
  return(I)
}