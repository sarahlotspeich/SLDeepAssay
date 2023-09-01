#' Expected value of Y
#' @name get_EY
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param q Proportion of p24-positive wells that underwent UDSA.
#' @param lambda Vector of DVL-specific parameters.
#' @return A scalar
#'
get_EY = function(M, q, lambda) {

  Lambda = sum(lambda)

  EY_terms = sapply(X = 1:M,
                    FUN = function(k) {
                       dbinom(k, M, 1 - exp(-Lambda)) * 
                        round(q * k) /
                        #floor(0.5 + q * k)
                        (1 - exp(-Lambda))
                    }
                    )

  EY = (1 - exp(-lambda)) * sum(EY_terms)

  return(EY)
}
