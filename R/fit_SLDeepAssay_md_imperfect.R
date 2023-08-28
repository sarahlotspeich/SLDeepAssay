#' Maximum likelihood estimator (MLE) for multiple-dilution imperfect assays data
#' @name fit_SLDeepAssay_md_imperfect
#' @param assay_md A list of assay data, with entries for each dilution level. Within each dilution level, the first entry should be a vector of QVOA results, with columns represent the replicate wells, and the second entry should be a matrix of UDSA results, with rows represent the distinct viral lineages (DVL) and columns represent the replicate wells. 
#' @param u A vector of dilution levels in millions of cells per well (all entries should be positive scalars).
#' @param sens_QVOA Sensitivity (i.e., true positive rate) for the QVOA (a scalar between 0 and 1). Default is \code{sens_QVOA = 1}. 
#' @param spec_QVOA Specificity (i.e., true negative rate) for the QVOA (a scalar between 0 and 1). Default is \code{spec_QVOA = 1}.
#' @param sens_UDSA Sensitivity (i.e., true positive rate) for the UDSA (a scalar between 0 and 1). Default is \code{sens_UDSA = 1}. 
#' @param spec_UDSA Specificity (i.e., true negative rate) for the UDSA (a scalar between 0 and 1). Default is \code{spec_UDSA = 1}.
#' @param corrected Logical, if \code{corrected = TRUE} the bias-corrected MLE will be returned. If \code{corrected = FALSE} the bias-corrected MLE will be not be returned. If \code{corrected = NULL}, the bias correction will be computed if here are <= 40 DVLs in \code{assay}. Default is \code{corrected = NULL}.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Parameter lower bound (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Parameter upper bound (passed to \code{optim}). Default is \code{ub = Inf}.
#' @return A named list with the following slots:
#' \item{mle}{MLE}
#' \item{convergence}{An integer code for the convergence status (returned by \code{optim})}
#' \item{message}{A character string giving any additional information about the convergence status (returned by \code{optim})}
#' @export
#'
#'
fit_SLDeepAssay_md_imperfect = function(assay_md, u, sens_QVOA = 1, spec_QVOA = 1, sens_UDSA = 1, spec_UDSA = 1, maxit = 1E4, lb = 1E-6, ub = Inf) {
  ########################################################################################
  # Compute constants ####################################################################
  ########################################################################################
  D = length(u) # number of dilution levels
  M = as.numeric(lapply(X = 1:D, 
                        FUN = function(d) ncol(assay_md[[d]][[2]]))) # number of wells (total)
  Y = lapply(X = 1:D, 
             FUN = function(d) rowSums(assay_md[[d]][[2]], na.rm = TRUE)) # number of infected wells per DVL = Y
  n = as.numeric(lapply(X = 1:D, 
                        FUN = function(d) nrow(assay_md[[d]][[2]])))# number of observed DVLs
  
  ########################################################################################
  # Set up for log-likelihood function ###################################################
  ########################################################################################
  # Build complete assay dataset
  cd_md = list()
  for (d in 1:D) {
    cd_md[[d]] = get_complete_data(Wstar = assay_md[[d]][[1]], 
                                   Zstar = assay_md[[d]][[2]])
    
    # Add columns of P(W*|W) to complete data for all wells
    ## < this won't change with lambda, so calculate once >
    if (nrow(cd_md[[d]]$complete_seq) > 0) {
      cd_md[[d]]$complete_seq$pWstarGivW = get_pWstarGivW(complete_data = cd_md[[d]]$complete_seq,
                                                          sens = sens_QVOA, 
                                                          spec = spec_QVOA)
      # Add column of P(Z*|Z) to complete data for sequenced wells
      ## < this won't change with lambda, so calculate once >
      cd_md[[d]]$complete_seq$pZstarGivZ = get_pZstarGivZ(complete_data = cd_md[[d]]$complete_seq,
                                                          n = n[d],
                                                          sens = sens_UDSA, 
                                                          spec = spec_UDSA)
    }
    
    if (nrow(cd_md[[d]]$complete_unseq) > 0) {
      cd_md[[d]]$complete_unseq$pWstarGivW = get_pWstarGivW(complete_data = cd_md[[d]]$complete_unseq,
                                                            sens = sens_QVOA, 
                                                            spec = spec_QVOA)
    }
  }
  
  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  optimization = optim(par = rep(0.1, n[1]),
                       fn = loglik_md_imperfect, 
                       u = u, 
                       complete_data_md = cd_md,
                       method = "L-BFGS-B",
                       control = list(maxit = maxit),
                       lower = rep(lb, n[1]),
                       upper = rep(ub, n[1]),
                       hessian = T)
  
  ### parameter estimate
  tau_hat = optimization$par
  Tau_hat = sum(tau_hat) # MLE of the IUPM
  
  # Fisher information matrix
  I = optimization$hessian
  
  # inverse of fisher information
  cov = tryCatch(expr = solve(I), 
                 error = function(e) matrix(data = NA, 
                                            nrow = length(lambda_hat), 
                                            ncol = length(lambda_hat)))
  
  # variance estimate
  se = sqrt(sum(cov))
  
  # confidence interval
  ci = exp(c(log(Tau_hat) + c(-1, 1) * (qnorm(0.975) * se / Tau_hat)))
  
  results = list("mle" = Tau_hat,
                 "se" = se,
                 "ci" = ci,
                 "convergence" = optimization$convergence,
                 "message" = optimization$message)
  return(results)
}
