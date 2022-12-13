#' Log-likelihood for IUPM from single-dilution assay data
#' @name loglik_sd
#' @param l Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MP Number of p24-positive wells.
#' @param m Number of p24-positive wells that underwent the UDSA.
#' @param y A vector of DVL-specific infection counts.
#' @return A scalar
#'
loglik_sd = function(l, M, MP, m, Y) {
  # Save n = # DVL detected
  n = length(l)

  # Check for invalid rate parameters (< 0)
  if (sum(l < 0) > 0) {
    ll = 0
  } else { # If parameters are valid, calculate log-likelihood
    ll = 0
    for(i in 1:n)
    {
      ll = ll + (Y[i] * log(1 - exp(-l[i])) - l[i] * (M - MP + m - Y[i])) # Wells without missing data
    }
    ll = ll + (MP - m) * log(1 - exp(-sum(l))) #Wells with missing data
  }
  return(-ll)
}
