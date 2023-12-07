#' Negative binomial log-likelihood for IUPM from single-dilution assay data
#' @name negbin_loglik_sd
#' @param l Vector of DVL-specific parameters.
#' @param gamma Dispersion parameter (a non-negative number)
#' @param k Alternative dispersion parameter, `k = 1/gamma` (a positive number)
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MP Number of p24-positive wells.
#' @param m Number of p24-positive wells that underwent the UDSA.
#' @param y A vector of DVL-specific infection counts.
#' @return A scalar
#'
negbin_loglik_sd = function(l, gamma, k = NULL, M, MP, m, Y) {
  
  # transform gamma to k if gamma supplied
  if (is.null(k)) {
    k = 1 / gamma
  }
  
  # if k = Inf (gamma = 0), then return Poisson loglik
  if (k == Inf) {
    
    return(loglik_sd(l = l, M = M, MP = MP, m = m, Y = Y))
    
  # Else for k < Inf, compute the NegBin loglik
  } else {
    
    # Save n = # DVL detected
    n = length(l)
    
    # Check for invalid rate parameters (< 0)
    if (sum(l < 0) > 0) {
      ll = 0
    } else { # If parameters are valid, calculate log-likelihood
      ll = 0
      for(i in 1:n)
      {
        ll = ll + Y[i] * log(1 - (k / (l[i] + k)) ^ k) +
          (M - MP + m - Y[i]) * k * log(k / (l[i] + k)) # Wells without missing data
      }
      ll = ll + (MP - m) * log(1 - prod(k / (l + k)) ^ k) # Wells with missing data
    }
    return(-ll)
  }
}


#' Negative Binomial log-likelihood for IUPM from multiple-dilution assay data
#' @name negbin_loglik_md
#' @param tau Vector of DVL-specific IUPM parameters
#' @param gamma Dispersion parameter (a non-negative number)
#' @param k Alternative dispersion parameter, `k = 1/gamma` (a positive number)
#' @param assay_summary Summary of assay data from multiple dilutions
#' @return A scalar value of the negative log-likelihood
#'
negbin_loglik_md = function(tau, gamma, k = NULL, assay_summary) {
  
  # transform gamma to k if gamma supplied
  if (is.null(k)) {
    k = 1 / gamma
  }
  
  # Total number of dilutions
  D = nrow(assay_summary)
  
  # Compute negative log-likelihoods for each dilution
  logliks_by_dilution <- vapply(
    X = 1:D,
    FUN.VALUE = numeric(1),
    FUN = function(d) {
      as.numeric(negbin_loglik_sd(
        l = assay_summary$u[d] * tau,
        k = k,
        M = assay_summary$M[d],
        MP = assay_summary$MP[d],
        m = assay_summary$m[d],
        Y = assay_summary[d, grepl("Y", colnames(assay_summary))]))
    })
  
  # Return the sum of (negative) log-likelihoods
  ll <- sum(logliks_by_dilution)
  
  return(ll)
}


#' Gradient of the negative binomial log-likelihood for IUPM from single-dilution assay data (with respect to tau and gamma)
#' @name negbin_gloglik_sd
#' @param l Vector of DVL-specific parameters.
#' @param gamma Dispersion parameter (a non-negative number)
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MP Number of p24-positive wells.
#' @param m Number of p24-positive wells that underwent the UDSA.
#' @param y A vector of  DVL-specific infection counts.
#' @return A vector
#'
negbin_gloglik_sd = function(l, gamma, M, MP, m, Y) {
  
  # if gamma = 0 (k = Inf), then return Poisson gloglik with NA (?) for entry n+1
  if (gamma == 0) {

    # compute this derivative numerically when needed
    d.gamma = numDeriv::grad(
      func = function(gamma) {
        negbin_loglik_sd(l = l, gamma = gamma, M = M, MP = MP, m = m, Y = Y)
      },
      x = 0
    )
    
    return(c(gloglik_sd(l = l, M = M, MP = MP, m = m, Y = Y), d.gamma))
    
  # Else for gamma > 0, compute the NegBin gloglik
  } else {
    
    k = 1 / gamma
    
    # Save n = # DVL detected
    n = length(l)
    
    # partials w.r.t. lambda_i
    gradient = numeric(n)
    for(i in 1:n) 
    { 
      gradient[i] = (Y[i] * (1 - (k / (l[i] + k)) ^ k) ^ (-1) +
                       (MP - m) * ((1 - prod((k / (l + k))) ^ k)) ^ (-1) *
                       prod(k / (l[-i] + k)) ^ k) *
        (k / (l[i] + k)) ^ (k + 1) -
        (M - MP + m - Y[i]) * (k / (l[i] + k))
    }
    
    # partial w.r.t. k
    grad.k = sum((Y * (1 - ((l + k) / k) ^ k) ^ (-1) +
                              (M - MP + m - Y) +
                              (MP - m) * (1 - (prod((l + k) / k)) ^ k) ^ (-1)) *
                             (log(k / (l + k)) + (l / (l + k))))
    
    # partial w.r.t. gamma (using chain rule)
    gradient[n + 1] = - grad.k / gamma ^ 2
    
  return(-gradient)
  }
}


#' Gradient of the log-likelihood for IUPM from multiple-dilution assay data
#' @name negbin_gloglik_md
#' @param tau Vector of DVL-specific IUPM parameters
#' @param gamma Dispersion parameter (a non-negative number)
#' @param assay_summary Summary of assay data from multiple dilutions
#' @return A vector (of length \code{n})
#'
negbin_gloglik_md = function(tau, gamma, assay_summary) {
  
  D = nrow(assay_summary)
  n = length(tau)
  
  # compute gradients for each dilution
  glogliks_byDilution = vapply(
    X = 1:D,
    FUN.VALUE = numeric(n + 1),
    FUN = function(d) {
      as.numeric(negbin_gloglik_sd(
        l = assay_summary$u[d] * tau,
        gamma = gamma,
        M = assay_summary$M[d],
        MP = assay_summary$MP[d],
        m = assay_summary$m[d],
        Y = as.numeric(
          assay_summary[d,grepl("Y", colnames(assay_summary))])))
    }
  )
  
  # return the sum of (negative) log-likelihoods
  # chain rule applied only to partials w.r.t. tau
  gradient = c(glogliks_byDilution[1:n, ] %*% assay_summary$u,
               sum(glogliks_byDilution[n + 1, ]))
  
  return(gradient)
}


#' Maximum likelihood estimator (MLE) for multiple-dilution assay data with LRT
#' @name lrt_SLDeepAssay_md
#' @param assay List of data frames with assay data from each dilution level. Default is \code{NULL}.
#' @param u Vector of dilution levels, in millions of cells per well. Default is \code{NULL}.
#' @param assay_summary (Optional) Instead of supplying \code{assay} and \code{u}, supply a summary of assay results in the form of a data frame. This summary should contain one row per dilution level and the following columns: M (total number of wells), n (number of distinct viral lineages \[DVL\]), MN (number of p24-negative wells), m (number of deep sequenced wells), Y1,..., Yn (counts of wells positive for DVL i, (i = 1,...,n), and u (dilution levels, in millions of cells per well).
#' @param corrected Logical, if \code{corrected = TRUE} the bias-corrected MLE will be returned. If \code{corrected = FALSE} the bias-corrected MLE will be not be returned. If \code{corrected = NULL}, the bias correction will be computed if here are <= 40 DVLs in \code{assay}. Default is \code{corrected = NULL}.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Lower-bound on the IUPM (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Upper-bound on the IUPM (passed to \code{optim}). Default is \code{ub = Inf}.
#' @param k0 initial value for k in optimization procedure. Default is \code{k = 1}.
#' @return Named list with the following slots:
#' \item{mle}{MLE}
#' \item{se}{Standard error for the MLE}
#' \item{ci}{95\% confidence interval for the MLE}
#' \item{mle_bc}{Bias-corrected MLE}
#' \item{ci_bc}{95\% confidence interval for the bias-corrected MLE}
#' @export
#'
lrt_SLDeepAssay_md = function(assay = NULL,
                              u = NULL,
                              assay_summary,
                              corrected = NULL,
                              maxit = 1E6,
                              lb = 1E-6,
                              ub = Inf) {
  
  # For each dilution level, compute summary data
  if (!is.null(assay)) {
    assay_summary = vapply(X = 1:length(assay),
                           FUN.VALUE = numeric(7 + n),
                           FUN = function(d) {
                             M = ncol(assay[[d]])
                             n = nrow(assay[[d]])
                             MN = sum(colSums(assay[[d]]) == 0, na.rm = TRUE)
                             MP = M - MN
                             m = MP - sum(is.na(colSums(assay[[d]])))
                             q = ifelse(MP == 0, 0, m / MP)
                             Y = rowSums(assay[[d]], na.rm = TRUE)
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
  
  # Fit MLE under Poisson model
  opt_pois = optim(par = log(rep(0.1, assay_summary$n[1])),
                       fn = function(theta, assay_summary) {
                         loglik_md(tau = exp(theta),
                                   assay_summary = assay_summary) },
                       gr = function(theta, assay_summary) {
                         gloglik_md(tau = exp(theta),
                                    assay_summary = assay_summary) }, 
                       assay_summary = assay_summary,
                       method = "BFGS",
                       control = list(maxit = maxit),
                       hessian = F)
  
  # Poisson model parameter estimates
  tau_hat = exp(opt_pois$par)
  Tau_hat = sum(tau_hat)
  
  # Fit MLE under NegBin model
  opt_negbin = optim(
    par = c(tau_hat, 0), # initiate NegBin search at Poisson MLE (gamma = 0)
    fn = function(theta, assay_summary) {
      negbin_loglik_md(tau = head(theta, assay_summary$n[1]),
                       gamma = tail(theta, 1),
                       assay_summary = assay_summary)
    },
    gr = function(theta, assay_summary) {
      negbin_gloglik_md(tau = head(theta, assay_summary$n[1]),
                        gamma = tail(theta, 1),
                        assay_summary = assay_summary)
    },
    assay_summary = assay_summary,
    method = "L-BFGS-B",
    control = list(maxit = maxit),
    lower = c(rep(lb, assay_summary$n[1]), 0),
    upper = ub,
    hessian = T)
  
  # MLE for gamma
  gamma_hat_negbin = tail(opt_negbin$par, 1)
  
  # log-likelihood values
  loglik_pois = -1 * opt_pois$value
  loglik_negbin = -1 * opt_negbin$value
  
  # negative binomial model parameter estimates
  if (gamma_hat_negbin > 0) {
    # likelihood ratio statistic
    lrt_stat = max(-2 * (loglik_pois - loglik_negbin), 0)
    tau_hat_negbin = head(opt_negbin$par, assay_summary$n[1])
    Tau_hat_negbin = sum(tau_hat_negbin)
  } else {
    lrt_stat = 0
    Tau_hat_negbin = Tau_hat
    gamma_hat_negbin = 0
  }
  
  # For large n, do not compute bias correction unless user overrides
  if (corrected == F) {
    
    Tau_hat_bc = NA
    ci_bc = NA
    
  } else {
    
    ### bias correction
    tau_hat_bc <- BC_md(tau = tau_hat,
                        M = assay_summary$M,
                        q = assay_summary$q,
                        u = assay_summary$u)
    
    # bias-corrected MLE for Tau
    Tau_hat_bc <- sum(tau_hat_bc)
    
  }
  
  return(list("mle" = Tau_hat,
              "mle_bc" = Tau_hat_bc,
              "mle_negbin" = Tau_hat_negbin,
              "mle_gamma" = gamma_hat_negbin,
              "lrt_stat" = lrt_stat))
}




