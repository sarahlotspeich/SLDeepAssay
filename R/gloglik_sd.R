#' Gradient of the log-likelihood for IUPM from single-dilution assay data
#' @name gloglik_sd
#' @param l Vector of DVL-specific parameters (of length \code{n}).
#' @param n Number of distinct viral lineages.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MP Number of p24-positive wells.
#' @param m Number of p24-positive wells that underwent the USDA.
#' @param y A vector of  DVL-specific infection counts (of length \code{n}).
#' @return A vector (of length \code{n})
#'
gloglik_sd = function(l, n, M, MP, m, Y)
{
  gradient = NULL
  for(i in 1:n)
  {
    gradient[i] = Y[i] * (exp(-l[i]) / (1 - exp(-l[i]))) - (M - MP + m - Y[i]) + (MP - m) * (exp(-sum(l)) / (1 - exp(-sum(l))))
  }
  return(-gradient)
}