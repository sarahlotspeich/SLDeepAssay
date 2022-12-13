#' Fisher information of the log-likelihood for IUPM from multiple-dilution assay data
#' @name fisher_md
#' @param tau Parameter
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @param u Vector of all dilution levels (in millions of cells per well).
#' @return A matrix
#'
fisher_md = function(tau, M, q, u) {
  D = length(u)
  n = length(tau)

  # fisher information for each dilution level
  I.d <- vapply(X = 1:D,
                FUN.VALUE = matrix(0, n, n),
                FUN = function(d) fisher(lambda = u[d] * tau,
                                         M = M[d],
                                         q = q[d]) * u[d] ^ 2
  )

  # fisher information over all dilution levels
  I <- apply(X = I.d,
             MARGIN = c(1, 2),
             FUN = sum)

  return(I)
}
