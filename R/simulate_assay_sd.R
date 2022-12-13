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

  # Create constant for number of existing DVLs
  n = length(lambda)

  # Generate Z = I(Xij >= 1) for all i, j
  Z = rbinom(n = M * n,
             size = 1,
             prob = (1 - exp(- rep(x = lambda, each = M))))

  data_mat = matrix(data = Z,
                    nrow = n,
                    ncol = M,
                    byrow = TRUE)

  # Calculate MP: number of wells with >= 1 DVL (p24+)
  MP = sum(colSums(data_mat) >= 1)

  # Calculate the number of p24 positive wells for UDSA
  if (q == 1) {
    m = MP
  } else {
    m = round(q * MP, 0)
  }

  # Randomly select columns to make missing
  p24_pos = which(colSums(data_mat) >= 1)
  p24_neg = which(colSums(data_mat) == 0)
  if (m < MP) {
    make_miss = p24_pos[1:(MP - m)]
    data_mat[, make_miss] = NA
  }

  # initiate n_
  n_ = n
  if (remove_undetected) {
    # Check for which DVL detected
    detected = rowSums(data_mat, na.rm = TRUE) >= 1

    # Number of DVL detected
    n_ = sum(detected)
    data_mat = data_mat[detected, ]
    if (n_ == 0) {
      return(list("any_DVL" = rep(x = 0, M),
                  "DVL_specific" = NULL))
    } else if (n_ == 1) {
      data_mat = matrix(data_mat, ncol = M, byrow = TRUE)
    }
  }

  # Reorder to look like the paper
  ## With non-missing data first
  if (m < MP) {
    data_mat = data_mat[, c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
  } else {
    data_mat = data_mat[, c(p24_pos, p24_neg)]
  }
  if (n_ == 1) {
    data_mat = matrix(data_mat, ncol = M, byrow = TRUE)
  }

  # Return simulated data
  return(list("any_DVL" = any_DVL(data_mat),
              "DVL_specific" = data_mat))
}
