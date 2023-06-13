#' Simulate multiple-dilution imperfect assays data
#' @name simulate_assay_md_imperfect
#' @param M Total number of wells at each dilution level (a vector of length D).
#' @param tau DVL-specific IUPMs (a vector of length \code{n}). (Note: All elements in \code{tau} must be > 0.)
#' @param q Proportions of p24-positive wells that underwent UDSA at each dilution level (a vector of length D).
#' @param u Vector of dilution levels (in millions of cells per well).
#' @param sens_QVOA Sensitivity (i.e., true positive rate) for the QVOA (a scalar between 0 and 1). Default is \code{sens_QVOA = 1}. 
#' @param spec_QVOA Specificity (i.e., true negative rate) for the QVOA (a scalar between 0 and 1). Default is \code{spec_QVOA = 1}.
#' @param sens_UDSA Sensitivity (i.e., true positive rate) for the UDSA (a scalar between 0 and 1). Default is \code{sens_UDSA = 1}. 
#' @param spec_UDSA Specificity (i.e., true negative rate) for the UDSA (a scalar between 0 and 1). Default is \code{spec_UDSA = 1}.
#' @param k Overdispersion parameter (a positive number). Default is \code{Inf}, which corresponds to no overdispersion.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells across all dilution levels are deleted.
#' @return Named list with the following slots:
#' \item{assay_summary}{A This summary should contain one row per dilution level and the following columns: M (total number of wells), n (number of distinct viral lineages \[DVL\]), MN (number of p24-negative wells), m (number of deep sequenced wells), Y1,..., Yn (counts of wells positive for DVL i, (i = 1,...,n), and dilutions (dilution levels, in millions of cells per well).}
#' @export
#'
simulate_assay_md_imperfect = function(M, tau, q, u, sens_QVOA = 1, spec_QVOA = 1, sens_UDSA = 1, spec_UDSA = 1, k = Inf, remove_undetected = T, seed = NULL) {
  
  D = length(u)
  n = length(tau)
  
  # Observed matrix for each dilution level
  assay = lapply(X = 1:D,
                 FUN = function(d) simulate_assay_sd(M = M[d],
                                                     tau = tau,
                                                     u = u[d],
                                                     q = q[d],
                                                     sens_QVOA = sens_QVOA, 
                                                     spec_QVOA = spec_QVOA, 
                                                     sens_UDSA = sens_UDSA, 
                                                     spec_UDSA = spec_UDSA, 
                                                     k = k,
                                                     remove_undetected = F))
  
  return(assay)
}

