#' Expected value of the third derivative of the Fisher information of the log-likelihood for IUPM from single-dilution assay data
#' @name exp_third_deriv_sd
#' @param lambda Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @return A matrix
#'
exp_third_deriv_sd = function(lambda, M, q) {
  
  n = length(lambda)
  Lambda = sum(lambda)
  
  ed3 = matrix(nrow = n, ncol = n ^ 2,
               M * (1 - q) * (exp(Lambda) + 1) /
                 (exp(Lambda) - 1) ^ 2)
  
  EY = get_EY(lambda = lambda, M = M, q = q)
  
  diag.terms = exp(lambda) * (exp(lambda) + 1) * EY /
    (exp(lambda) - 1) ^ 3
  
  for (i in 1:n) {
    ed3[i, (n+1)*i - n] = diag.terms[i]
  }
  
  return(ed3)
}