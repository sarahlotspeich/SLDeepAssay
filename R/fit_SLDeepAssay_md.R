#' Maximum likelihood estimator (MLE) for multiple-dilution assay data
#' @name fit_SLDeepAssay_md
#' @param assay List of data frames with assay data from each dilution level. Default is \code{NULL}.
#' @param u Vector of dilution levels, in millions of cells per well. Default is \code{NULL}.
#' @param assay_summary (Optional) Instead of supplying \code{assay} and \code{u}, supply a summary of assay results in the form of a data frame. This summary should contain one row per dilution level and the following columns: M (total number of wells), n (number of distinct viral lineages \[DVL\]), MN (number of p24-negative wells), m (number of deep sequenced wells), Y1,..., Yn (counts of wells positive for DVL i, (i = 1,...,n), and u (dilution levels, in millions of cells per well).
#' @param corrected Logical, if \code{corrected = TRUE} the bias-corrected MLE will be returned. If \code{corrected = FALSE} the bias-corrected MLE will be not be returned. If \code{corrected = NULL}, the bias correction will be computed if here are <= 40 DVLs in \code{assay}. Default is \code{corrected = NULL}.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @return Named list with the following slots:
#' \item{mle}{MLE}
#' \item{se}{Standard error for the MLE}
#' \item{ci}{95\% confidence interval for the MLE}
#' \item{mle_bc}{Bias-corrected MLE}
#' \item{ci_bc}{95\% confidence interval for the bias-corrected MLE}
#' @export
#'
#'
fit_SLDeepAssay_md = function(assay = NULL, u = NULL, assay_summary, corrected = NULL, maxit = 1E6) {

  # For each dilution level, compute summary data
  if (!is.null(assay)) {
    assay_summary = vapply(X = 1:length(assay), 
                           FUN.VALUE = numeric(7 + n),
                           FUN = function(d) {
                             M = ncol(assay[[d]][[2]])
                             n = nrow(assay[[d]][[2]])
                             MP = sum(assay[[d]][[1]])
                             MN = M - MP
                             m = sum(!is.na(colSums(assay[[d]][[2]]))) - MN
                             q = ifelse(MP == 0, 0, m / MP)
                             Y = rowSums(assay[[d]][[2]], na.rm = TRUE)
                             return((c("u" = u[d], "M" = M, "n" = n,
                                       "MN" = MN, "MP" = MP, "m" = m, "q" = q, "Y" = Y)))
                           })
    assay_summary = as.data.frame(t(assay_summary))
  }

  # Indicator for whether bias correction should be computed:
  # user specified value if provided, else yes if n <= 40
  corrected = ifelse(test = is.null(corrected),
                     yes = assay_summary$n[1] <= 40,
                     no = corrected)

  # Fit MLE (BFGS)
  optimization = optim(par = log(rep(0.1, assay_summary$n[1])),
                       fn = function(theta, assay_summary) {
                         loglik_md(tau = exp(theta),
                                   assay_summary = assay_summary) },
                         gr = function(theta, assay_summary) {
                           gloglik_md(tau = exp(theta),
                                      assay_summary = assay_summary) *
                             exp(theta) },   # chain rule 
                       assay_summary = assay_summary,
                       method = "BFGS",
                       control = list(maxit = maxit),
                       hessian = F)
    
  ### parameter estimate
  tau_hat = exp(optimization$par)
  Tau_hat = sum(tau_hat)
  
  # Fisher information matrix
  I = fisher_md(tau = tau_hat,
                 M = assay_summary$M,
                 q = assay_summary$q,
                 u = assay_summary$u)

  # inverse of fisher information
  cov = tryCatch(expr = solve(I), 
                 error = function(e) rep(NA, nrow(I)))

  ### variance estimate 4
  se = sqrt(sum(cov))

  ### confidence interval
  ci = exp(c(log(Tau_hat) + c(-1, 1) * (qnorm(0.975) * se / Tau_hat)))

  # For large n, do not compute bias correction unless user overrides
  if (corrected == F) {
    Tau_hat_bc = NA
    ci_bc = NA
  } else {
    ### bias correction
    tau_hat_bc = BC_md(tau = tau_hat,
                        M = assay_summary$M,
                        q = assay_summary$q,
                        u = assay_summary$u)
    # bias-corrected MLE for Tau
    Tau_hat_bc = sum(tau_hat_bc)
    # bias corrected CI
    ci_bc = exp(c(log(Tau_hat_bc) + c(-1, 1) * (qnorm(0.975) * se / Tau_hat_bc)))
  }

  return(list("mle" = Tau_hat,
              "se" = se,
              "ci" = ci,
              "mle_bc" = Tau_hat_bc,
              "ci_bc" = ci_bc))
}

