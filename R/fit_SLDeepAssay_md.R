#' Maximum likelihood estimator (MLE) for multiple-dilution assay data
#' @name fit_SLDeepAssay_md
#' @param assay List of data frames with assay data from each dilution level. Default is \code{NULL}.
#' @param dilutions Vector of dilutions. Default is \code{NULL}.
#' @param assay.summary (Optional) Instead of supplying \code{assay}, supply a summary of assay results. This summary should contain one row per dilution level and the following columns: M (total number of wells), n (number of distinct viral lineages \[DVL\]), MN (number of p24-negative wells), m (number of deep sequenced wells), Y1,..., Yn (counts of wells positive for DVL i, (i = 1,...,n), and dilutions (dilution levels, in millions of cells per well).
#' @param bc Logical indicator for returning bias corrected MLE and CI. By default, the bias correction is not computed if n>40.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Lower-bound on the IUPM (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Upper-bound on the IUPM (passed to \code{optim}). Default is \code{ub = Inf}.
#' @return Named list with the following slots:
#' \item{mle}{MLE}
#' \item{se}{Standard error for the MLE}
#' \item{ci}{95\% confidence interval for the MLE}
#' \item{mle.bc}{Bias-corrected MLE}
#' \item{ci.bc}{95\% confidence interval for the bias-corrected MLE}
#' @export
#'
#'
fit_SLDeepAssay_md <- function(assay = NULL, dilutions = NULL, assay.summary,
                      bc = NULL, maxit = 1E6, lb = 1E-6, ub = Inf) {

  # For each dilution level, compute summary data
  if (!is.null(assay)) {
    assay.summary = vapply(X = 1:D, FUN.VALUE = numeric(7 + n),
                           FUN = function(d) {
                             M = ncol(assay[[d]])
                             n = nrow(assay[[d]])
                             MN = sum(colSums(assay[[d]]) == 0, na.rm = TRUE)
                             MP = M - MN
                             m = MP - sum(is.na(colSums(assay[[d]])))
                             q = ifelse(MP == 0, 0, m / MP)
                             Y = rowSums(assay[[d]], na.rm = TRUE)
                             return((c("dilution" = dilutions[d], "M" = M, "n" = n,
                                       "MN" = MN, "MP" = MP, "m" = m, "q" = q, "Y" = Y)))
                           })
    assay.summary = as.data.frame(t(assay.summary))
  }

  # Indicator for whether bias correction should be computed:
  # user specified value if provided, else yes if n <= 40
  bc = ifelse(test = is.null(bc),
              yes = assay.summary$n[1] <= 40,
              no = bc)

  # Fit MLE
  optimization = optim(par = rep(0, assay.summary$n[1]),
                       fn = function(t)
                         loglik_md(tau = t, assay.summary = assay.summary),
                       gr = function(t)
                         gloglik_md(tau = t, assay.summary = assay.summary),
                       method = "L-BFGS-B",
                       control = list(maxit = maxit),
                       lower = rep(lb, assay.summary$n[1]),
                       upper = rep(ub, assay.summary$n[1]),
                       hessian = T)

  ### parameter estimate
  tau.hat = optimization$par
  Tau.hat = sum(tau.hat) # MLE of the IUPM

  # Fisher information matrix
  I <- fisher_md(tau = tau.hat, M = assay.summary$M, q = assay.summary$q,
                 dilutions = assay.summary$dilution)

  # inverse of fisher information
  cov <- solve(I)

  ### variance estimate 4
  se <- sqrt(sum(cov))

  ### confidence interval
  ci = exp(c(log(Tau.hat) + c(-1, 1) * (qnorm(0.975) * se / Tau.hat)))

  # For large n, do not compute bias correction unless user overrides
  if (bc == F) {

    Tau.hat.bc = NA
    se.bc = NA
    ci.bc = NA

  } else {

    ### bias correction
    tau.hat.bc <- BC_md(tau = tau.hat, M = assay.summary$M, q = assay.summary$q,
                        dilutions = assay.summary$dilution)

    # bias-corrected MLE for Tau
    Tau.hat.bc <- sum(tau.hat.bc)

    # bias corrected CI
    ci.bc <- exp(c(log(Tau.hat.bc) + c(-1, 1) *
                     (qnorm(0.975) * se / Tau.hat.bc)))
  }

  return(list("mle" = Tau.hat,
              "se" = se,
              "ci" = ci,
              "mle.bc" = Tau.hat.bc,
              "ci.bc" = ci.bc))
}
