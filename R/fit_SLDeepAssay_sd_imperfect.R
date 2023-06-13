#' Maximum likelihood estimator (MLE) for single-dilution imperfect assays data
#' @name fit_SLDeepAssay_sd_imperfect
#' @param assay_QVOA Assay data, with columns representing the wells. 
#' @param assay_UDSA Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells. 
#' @param u Dilution level in millions of cells per well (a positive scalar). Default is \code{u = 1}.
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
fit_SLDeepAssay_sd_imperfect = function(assay_QVOA, assay_UDSA, u = 1, sens_QVOA = 1, spec_QVOA = 1, sens_UDSA = 1, spec_UDSA = 1, maxit = 1E4, lb = 1E-6, ub = Inf) {
  ########################################################################################
  # Compute constants ####################################################################
  ########################################################################################
  M = length(assay_QVOA) # number of wells (total)
  Y = rowSums(assay_UDSA, na.rm = TRUE) # number of infected wells per DVL = Y
  n = length(Y) # number of observed DVLs
  
  ########################################################################################
  # Set up for log-likelihood function ###################################################
  ########################################################################################
  # Build complete assay dataset
  cd = get_complete_data(Wstar = assay_QVOA, 
                         Zstar = assay_UDSA)
  
  # Add columns of P(W*|W) to complete data for all wells
  ## < this won't change with lambda, so calculate once >
  cd$complete_seq$pWstarGivW = get_pWstarGivW(complete_data = cd$complete_seq,
                                              sens = sens_QVOA, 
                                              spec = spec_QVOA)
  cd$complete_unseq$pWstarGivW = get_pWstarGivW(complete_data = cd$complete_unseq,
                                                sens = sens_QVOA, 
                                                spec = spec_QVOA)
  
  # Add column of P(Z*|Z) to complete data for sequenced wells
  ## < this won't change with lambda, so calculate once >
  cd$complete_seq$pZstarGivZ = get_pZstarGivZ(complete_data = cd$complete_seq,
                                              n = n,
                                              sens = sens_UDSA, 
                                              spec = spec_UDSA)
  
  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  optimization = optim(par = - log(1 - Y / M),
                       fn = loglik_sd_imperfect,
                       complete_data = cd, 
                       method = "L-BFGS-B",
                       control = list(maxit = maxit),
                       lower = rep(lb, n),
                       upper = rep(ub, n),
                       hessian = F)
  
  lambda_hat = optimization$par
  Lambda_hat = sum(lambda_hat) # MLE of the IUPM
  
  results = list("mle" = Lambda_hat / u,
                 "convergence" = optimization$convergence,
                 "message" = optimization$message)
  
  return(results)
}
