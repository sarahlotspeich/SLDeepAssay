#' Maximum likelihood estimator (MLE) for single-dilution assay data
#' @name fit_SLDeepAssay_sd
#' @param assay Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells. Any columns with \code{NA} values will be treated as positive but missing deep sequencing information.
#' @param u Dilution level in millions of cells per well (a positive scalar). Default is \code{u = 1}.
#' @param M (Optional) Instead of supplying \code{assay}, supply the total number of wells (a scalar). Default is \code{M = NULL}.
#' @param MP (Optional) Instead of supplying \code{assay}, supply the total number of p24-positive wells (a scalar). Default is \code{MP = NULL}.
#' @param m (Optional) Instead of supplying \code{assay}, supply the total number of p24-positive wells that underwent deep sequencing (a scalar). Default is \code{m = NULL}.
#' @param Y (Optional) Instead of supplying \code{assay}, supply the numbers of wells (without missing data) that were infected with each DVL (a vector of length \code{n}). Default is \code{Y = NULL}.
#' @param corrected Logical, if \code{corrected = TRUE} the bias-corrected MLE will be returned. If there are <= 40 DVL in \code{assay}, default is \code{corrected = TRUE}; otherwise default is \code{corrected = FALSE}.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Parameter lower bound (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Parameter upper bound (passed to \code{optim}). Default is \code{ub = Inf}.
#' @return A named list with the following slots:
#' \item{mle}{MLE}
#' \item{se}{Standard error for the MLE}
#' \item{ci}{95\% confidence interval for the MLE}
#' \item{mle_bc}{Bias-corrected MLE}
#' \item{ci_bc}{95\% confidence interval for the bias-corrected MLE}
#' @export
#'
#'
fit_SLDeepAssay_sd = function(assay, u = 1, M = NULL, n = NULL, MP = NULL, m = NULL, Y = NULL, corrected = FALSE, maxit = 1E4, lb = 1E-6, ub = Inf) {
  # If the raw assay data are provided then compute the following values
  if (!is.null(assay)) {
    M = ncol(assay) # number of wells (total)
    MN = sum(colSums(assay) == 0, na.rm = TRUE) # number of wells (p24-negative)
    MP = M - MN # number of wells (p24-positive)
    m = M - MN - sum(is.na(colSums(assay))) # number of deep-sequenced p24-positive wells
    Y = rowSums(assay, na.rm = TRUE) # number of infected wells per DVL = Y
  } else { # If the raw assay data are not provided, compute just the following value
    MN = M - MP
  }
  q = m / MP # proportion of p24-positive wells deep sequenced
  n = length(Y) # number of observed DVLs

  # Indicator for whether bias correction should be computed:
  # User-specified value if provided, else yes if n <= 40
  corrected = ifelse(!corrected,
                     n <= 40,
                     corrected)

  # Fit MLE
  optimization = optim(par = - log(1 - Y / M),
                       fn = loglik_sd,
                       gr = gloglik_sd,
                       M = M,
                       MP = MP,
                       m = m,
                       Y = Y,
                       method = "L-BFGS-B",
                       control = list(maxit = maxit),
                       lower = rep(lb, n),
                       upper = rep(ub, n),
                       hessian = F)

  lambda_hat = optimization$par
  Lambda_hat = sum(lambda_hat) # MLE of the IUPM

  # Fisher information matrix
  I = fisher_sd(lambda = lambda_hat,
                M = M,
                q = q)

  # inverse of fisher information
  cov = solve(I)

  # variance estimate
  se = sqrt(sum(cov))

  # confidence interval
  ci = exp(c(log(Lambda_hat) + c(-1, 1) * (qnorm(0.975) * se / Lambda_hat)))

  # For large n, do not compute bias correction unless user overrides
  if (!corrected) {
    Lambda_hat_bc = NA
    ci_bc = NA
  } else {
    # Bias corrected MLE (BC MLE)
    lambda_hat_bc = BC_sd(lambda = lambda_hat,
                          M = M,
                          q = q)

    Lambda_hat_bc = sum(lambda_hat_bc)

    # confidence interval
    ci_bc = exp(c(log(Lambda_hat_bc) +
                    c(-1, 1) * (qnorm(0.975) * se / Lambda_hat_bc)))
  }

  results = list("mle" = Lambda_hat / u,
                 "se" = se / u,
                 "ci" = ci / u,
                 "mle_bc" = Lambda_hat_bc / u,
                 "ci_bc" = ci_bc / u)

  return(results)
}
