#' Simulate single-dilution assay data
#' @name simulate_assay_sd
#' @param M Total number of wells (a scalar).
#' @param tau Mean counts of cells per million infected with each DVL (a vector). (Note: All elements in \code{tau} must be > 0.)
#' @param q Proportion of p24-positive wells that underwent UDSA (a scalar between 0 and 1).
#' @param u Dilution level in millions of cells per well (a positive scalar). Default is \code{dilution = 1}.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells are deleted.
#' @return Named list with the following slots:
#' \item{any_DVL}{A vector containing overall (any DVL) infection indicators of across the wells.}
#' \item{DVL_specific}{A matrix of infection indicators with rows and columns representing the DVLs and wells, respectively.}
#' @export
#'
simulate_assay_sd = function(M, tau, q, u = 1, remove_undetected = TRUE) {
  # Calculate lambda: mean number of infected cells per well with each DVL
  lambda = tau * u

  # Create n: a constant for the number of underlying DVLs
  n = length(lambda)

  # Generate true UDSA statuses, Z = I(Xij >= 1), for all i, j
  Z = rbinom(n = M * n,
             size = 1,
             prob = (1 - exp(- rep(x = lambda, each = M))))
  
  # Reshape to get rows per DVL, columns per replicate well
  Z_mat = matrix(data = Z,
                    nrow = n,
                    ncol = M,
                    byrow = TRUE)
  
  # Construct true QVOA results from true UDSA results
  W = get_any_DVL(Z_mat)
  
  # Calculate MP: number of wells with >= 1 DVL (p24+)
  MP = sum(W == 1)

  # Calculate m: ## the number of p24-positive wells to deep sequence
  m = ifelse(test = q == 1, yes = MP, no = round(q * MP, 0)) 
  
  # Make Z for any unsequenced, p24-positive wells NA
  p24_pos = which(W == 1) ## ids of p24-positive wells
  p24_neg = which(W == 0) ## ids of p24-negative wells
  # If only partially sequencing (m < MP), make Z missing for unsequenced wells 
  if (m < MP) { 
    make_miss = p24_pos[-c(1:m)]
    Z_mat[, make_miss] = NA
  }
  
  # Subset to DVLs to return, either all (if remove_undetected = FALSE)  
  # or only those that were detected in >= 1 well (if remove_undetected = TRUE) 
  ## Create n_: a constant for the number of detected DVLs
  n_ = n ### initialize as equal to the number of underlying DVLs
  if (remove_undetected) {
    detected = rowSums(Z_mat, na.rm = TRUE) >= 1
    n_ = sum(detected)
    Z_mat = Z_mat[detected, ]
    if (n_ == 0) {
      Z_mat = matrix(0, ncol = M)
    } else if (n_ == 1) {
      Z_mat = matrix(Z_mat, ncol = M, byrow = TRUE)
    }
  } 
  
  # Reorder columns of assay returned to be... 
  ## (i) p24-positive + sequenced, 
  ## (ii) p24-negative, and 
  ## (iii) p24-positive + unsequenced 
  if (m < MP) {
    Z_mat = Z_mat[, c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
    W = W[c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
  } else {
    Z_mat = Z_mat[, c(p24_pos, p24_neg)]
    W = W[c(p24_pos, p24_neg)]
  }
  if (n_ == 1) {
    Z_mat = matrix(Z_mat, ncol = M, byrow = TRUE)
  }
 
  # Return simulated data
  return(list("any_DVL" = W,
              "DVL_specific" = Z_mat))
}
