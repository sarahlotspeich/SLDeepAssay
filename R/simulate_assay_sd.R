#' Simulate single-dilution assay data
#' @name simulate_assay_sd
#' @param M Total number of wells (a scalar).
#' @param tau Mean counts of cells per million infected with each DVL (a vector). (Note: All elements in \code{tau} must be > 0.)
#' @param q Proportion of p24-positive wells that underwent UDSA (a scalar between 0 and 1).
#' @param u Dilution level in millions of cells per well (a positive scalar). Default is \code{dilution = 1}.
#' @param sens_QVOA Sensitivity (i.e., true positive rate) for the QVOA (a scalar between 0 and 1). Default is \code{sens_QVOA = 1}. 
#' @param spec_QVOA Specificity (i.e., true negative rate) for the QVOA (a scalar between 0 and 1). Default is \code{spec_QVOA = 1}.
#' @param sens_UDSA Sensitivity (i.e., true positive rate) for the UDSA (a scalar between 0 and 1). Default is \code{sens_UDSA = 1}. 
#' @param spec_UDSA Specificity (i.e., true negative rate) for the UDSA (a scalar between 0 and 1). Default is \code{spec_UDSA = 1}.
#' @param sequence_all (Optional) Logical, if \code{sequence_all = FALSE} (the default), then only p24-positive wells are considered for UDSA. If \code{sequence_all = TRUE} and \code{q = 1}, then all wells will undergo UDSA. 
#' @param k (Optional) Overdispersion parameter (a positive number). Default is \code{k = Inf}, which corresponds to no overdispersion.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells are deleted.
#' @return Named list with the following slots:
#' \item{any_DVL}{A vector containing overall (any DVL) infection indicators of across the wells.}
#' \item{DVL_specific}{A matrix of infection indicators with rows and columns representing the DVLs and wells, respectively.}
#' @export
#'
simulate_assay_sd = function(M, tau, q, u = 1, sens_QVOA = 1, spec_QVOA = 1, sens_UDSA = 1, spec_UDSA = 1, sequence_all = FALSE, k = Inf, remove_undetected = TRUE) {
  # Calculate lambda: mean number of infected cells per well with each DVL
  lambda = tau * u

  # Create n: a constant for the number of underlying DVLs
  n = length(lambda)
  
  # Calculate p: probabilitiy of success for Z based on Poisson or NegBin
  p = ifelse(test = k == Inf,
             yes = 1 - exp(- rep(x = lambda, each = M)),
             no = 1 - (k / (lambda + k)) ^ k)

  # Generate Z: true UDSA statuses for all DVL and wells
  Z = rbinom(n = M * n,
             size = 1,
             prob = p)
  
  # Reshape to get rows per DVL, columns per replicate well
  Z_mat = matrix(data = Z,
                 nrow = n,
                 ncol = M,
                 byrow = TRUE)
  
  # Construct W: true QVOA results from true UDSA results
  W = get_any_DVL(Z_mat)
  
  # Calculate MP: number of wells with >= 1 DVL (p24+)
  MP = sum(W == 1)

  # Calculate m: the number of p24-positive wells to deep sequence
  m = ifelse(test = q == 1, 
             yes = MP, 
             no = round(q * MP, 0)) 
  
  # Make Z for any unsequenced, p24-positive wells NA
  p24_pos = which(W == 1) ## ids of p24-positive wells
  p24_neg = which(W == 0) ## ids of p24-negative wells
  # If only partially sequencing (m < MP), make Z missing for unsequenced wells 
  if (m < MP & m > 0) { 
    make_miss = p24_pos[-c(1:m)]
    Z_mat[, make_miss] = NA
  } else if (m == 0) {
    Z_mat[, p24_pos] = NA
  }
  
  # Generate imperfect QVOA results
  if (sens_QVOA < 1 | spec_QVOA < 1) {
    # Generate W*: observed QVOA statuses for all wells
    Wstar = rbinom(n = M, 
                   size = 1, 
                   prob = sens_QVOA * W + (1 - spec_QVOA) * (1 - W))
  } else {
    # Assume perfect assay
    Wstar = W
  }
  
  # Generate imperfect UDSA results
  if (sens_UDSA < 1 | spec_UDSA < 1) {
    # Generate Z*: observed UDSA statuses for all DVL and wells
    Zstar = rbinom(n = M * n, 
                   size = 1, 
                   prob = sens_UDSA * Z + (1 - spec_UDSA) * (1 - Z))
    
    # Reshape to get rows per DVL, columns per replicate well
    Zstar_mat = matrix(data = Zstar, 
                       nrow = n, 
                       ncol = M, 
                       byrow = TRUE)
  } else {
    # Assume perfect assay
    Zstar_mat = Z_mat 
  }
  
  # Calculate MP: number of wells with >= 1 DVL (p24+)
  MP = sum(Wstar == 1)
  
  # Calculate m: the number of p24-positive wells to deep sequence
  m = ifelse(test = q == 1, 
             yes = MP, 
             no = round(q * MP, 0)) 
  
  # Make Z for any unsequenced, p24-positive wells NA
  p24_pos = which(Wstar == 1) ## ids of p24-positive wells
  p24_neg = which(Wstar == 0) ## ids of p24-negative wells
  
  # Unless sequence_all and q = 1, make some UDSA results missing 
  if (!(sequence_all & q == 1)) {
    # If only partially sequencing (m < MP), make Z missing for unsequenced positive wells 
    if (m < MP & m > 0) { 
      make_miss = p24_pos[-c(1:m)]
      Zstar_mat[, make_miss] = NA
    } else if (m == 0) {
      Zstar_mat[, p24_pos] = NA
    }
    # If QVOA was imperfect, make Z missing for negative wells 
    if (sens_QVOA < 1 | spec_QVOA < 1) { 
      Zstar_mat[, p24_neg] = NA
    }    
  } 
 
  # Subset to DVLs to return, either all (if remove_undetected = FALSE)  
  # or only those that were detected in >= 1 well (if remove_undetected = TRUE) 
  ## Create n_: a constant for the number of detected DVLs
  n_ = n ### initialize as equal to the number of underlying DVLs
  if (remove_undetected) {
    detected = rowSums(Zstar_mat, na.rm = TRUE) >= 1
    n_ = sum(detected)
    Zstar_mat = Zstar_mat[detected, ]
    if (n_ == 0) {
      Zstar_mat = matrix(NA, ncol = M)
    } else if (n_ == 1) {
      Zstar_mat = matrix(Zstar_mat, ncol = M, byrow = TRUE)
    }
  } 
  
  # Reorder columns of assay returned to be... 
  ## (i) p24-positive + sequenced, 
  ## (ii) p24-negative, and 
  ## (iii) p24-positive + unsequenced 
  if (m < MP) {
    Zstar_mat = Zstar_mat[, c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
    Wstar = Wstar[c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
  } else {
    Zstar_mat = Zstar_mat[, c(p24_pos, p24_neg)]
    Wstar = Wstar[c(p24_pos, p24_neg)]
  }
  if (n_ == 1) {
    Zstar_mat = matrix(Zstar_mat, ncol = M, byrow = TRUE)
  }
  
  # Return simulated data
  return(list("any_DVL" = Wstar,
              "DVL_specific" = Zstar_mat))
}
