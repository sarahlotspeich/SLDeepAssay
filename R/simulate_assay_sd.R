#' Simulate single-dilution assay data
#' @name simulate_assay_sd
#' @param M Total number of wells (a scalar).
#' @param n Total number of existing distinct viral lineages (DVL) (a scalar).
#' @param lambda DVL-specific rates of infection (a vector of length \code{n}). (Note: All elements in \code{lambda} must be > 0.)
#' @param q Proportion of p24-positive wells that underwent UDSA (a scalar between 0 and 1).
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells are deleted.
#' @param seed (Optional) An integer setting the random seed used to simulate the assay data. Default is \code{seed = NULL}.
#' @return Named list with the following slots:
#' \item{any_DVL}{A \code{Mx1} vector containing indicators of overall (any DVL) infection across the wells.}
#' \item{DVL_specific}{Standard error for the MLE}
#' @export
#'
simulate_assay_sd = function(M, n, lambda, q, remove_undetected = TRUE, seed = NULL) {
  # If supplied, set the random seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate Z = I(Xij >= 1) for all i, j
  Z <- rbinom(n = M * n,
              size = 1,
              prob = (1 - exp(- rep(x = lambda, each = M))))

  data_mat <- matrix(data = Z, nrow = n, ncol = M, byrow = TRUE)

  # Calculate MP: number of wells with >= 1 DVL (p24+)
  MP <- sum(colSums(data_mat) >= 1)

  # Calculate the number of p24 positive wells for UDSA
  if (q == 1) {
    m <- MP
  } else {
    m <- round(q * MP, 0)
  }

  # Randomly select columns to make missing
  p24_pos <- which(colSums(data_mat) >= 1)
  p24_neg <- which(colSums(data_mat) == 0)
  if (m < MP) {
    make_miss <- p24_pos[1:(MP - m)]
    data_mat[, make_miss] <- NA
  }

  # initiate n_
  n_ <- n

  if (remove_undetected) {
    # Check for which DVL detected
    detected <- rowSums(data_mat, na.rm = TRUE) >= 1
    # Number of DVL detected
    n_ <- sum(detected)
    data_mat <- data_mat[detected, ]
    if (n_ == 0) {
      return(NULL)
    } else if (n_ == 1) {
      data_mat <- matrix(data_mat, ncol = M, byrow = TRUE)
    }
  }

  # Reorder to look like the paper
  ## With non-missing data first
  if (m < MP) {
    data_mat <- data_mat[, c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
  } else {
    data_mat <- data_mat[, c(p24_pos, p24_neg)]
  }
  if (n_ == 1) {
    data_mat <- matrix(data_mat, ncol = M, byrow = TRUE)
  }
  return(list("any_DVL" = any_DVL(data_mat),
              "DVL_specific" = data_mat))
}
