#' Gradient of the log-likelihood for IUPM from single-dilution assay data
#' @name gloglik_sd
#' @param l Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MP Number of p24-positive wells.
#' @param m Number of p24-positive wells that underwent the UDSA.
#' @param y A vector of  DVL-specific infection counts.
#' @return A vector
#'
gloglik_sd = function(l, M, MP, m, Y) {
  # Save n = # DVL detected
  n = length(l)

  gradient = NULL
  for(i in 1:n)
  {
    gradient[i] = Y[i] * (exp(-l[i]) / (1 - exp(-l[i]))) - (M - MP + m - Y[i]) + (MP - m) * (exp(-sum(l)) / (1 - exp(-sum(l))))
  }
  return(-gradient)
}
