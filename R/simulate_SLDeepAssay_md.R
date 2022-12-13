#' Simulate multiple-dilution assay data and fit new and existing methods to it.
#' @name simulate_SLDeepAssay_md
#' @param M Total number of wells at each dilution level (a vector of length D).
#' @param tau DVL-specific IUPMs (a vector of length \code{n}). (Note: All elements in \code{tau} must be > 0.)
#' @param q Proportions of p24-positive wells that underwent UDSA at each dilution level (a vector of length D).
#' @param u Vector of dilution levels in millions of cells per well (a vector of length D).
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells across all dilution levels are deleted.
#' @return Named list with the following slots:
#' \item{assay_summary}{Simulated multiple-dilution assay data (deep sequenced) in summary form.}
#' \item{MLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (without deep sequencing information).}
#' \item{BCMLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (without deep sequencing information).}
#' \item{MLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (with deep sequencing information).}
#' \item{BCMLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (with deep sequencing information).}
#' \item{Message}{(Optional) Message describing whether the \code{Assay} needed to be re-simulated due to either none of the DVL being detected (\code{Message = 1}) or at least one DVL being detected in all deep-sequenced wells (\code{Message = 2}). Otherwise, \code{Message = 0}.}
#' @importFrom SLDAssay get.mle 
#' @export
#'
simulate_SLDeepAssay_md <- function(M, tau, q, u, remove_undetected = TRUE, seed = NULL) {
  # If supplied, set the random seed
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  # Create indicators of whether data need to be re-simulated
  num_redo_all = 0 # Due to any DVL being detected in all wells or
  num_redo_none = 0 # No DVL being detected in any wells
  
  # Simulate single-dilution assay data
  assay_summary = simulate_assay_md(M = M,
                            tau = tau, 
                            q = q,
                            u = u,
                            remove_undetected = remove_undetected)
  
  # Check for need to re-simulate
  max_p24_neg = tryCatch(expr = max(assay_summary$M - assay_summary$MP), 
                         error = function(e) 0)
  
  min_DVL_neg = tryCatch(expr = min(colSums(assay_summary$m - assay_summary[, grepl("Y", colnames(assay_summary))])), 
                               error = function(e) 0) 
  
  while ((max_p24_neg == 0 & min_DVL_neg == 0) | is.null(assay_summary)) { 
    if (is.null(assay)) {
      num_redo_none = num_redo_none + 1
    } else {
      num_redo_all = num_redo_all + 1
    }
    # Re-simulate single-dilution assay data
    assay_summary = simulate_assay_md(M = M,
                                      tau = tau, 
                                      q = q,
                                      u = u,
                                      remove_undetected = remove_undetected)
    
    max_p24_neg = tryCatch(expr = max(assay_summary$M - assay_summary$MP), 
                           error = function(e) 0)
    
    min_DVL_neg = tryCatch(expr = min(colSums(assay_summary$m - assay_summary[, grepl("Y", colnames(assay_summary))])), 
                           error = function(e) 0) 
  }
  
  # Methods without UDSA
  woUDSA_res = get.mle(pos = assay_summary$MP,
                       replicates = assay_summary$M,
                       u = assay_summary$u * 10^6)
  se = (log(woUDSA_res$MLE) - log(woUDSA_res$Asymp_CI[1])) / qnorm(p = 1 - (0.05 / 2))
  MLE_woUDSA = data.frame(Est = woUDSA_res$MLE, 
                          SE = se, 
                          LB = woUDSA_res$Asymp_CI[1], 
                          UB = woUDSA_res$Asymp_CI[2])
  BCMLE_woUDSA = data.frame(Est = woUDSA_res$BC_MLE, 
                            SE = se, 
                            LB = exp(log(woUDSA_res$BC_MLE) - qnorm(1 - 0.05 / 2) * se), 
                            UB = exp(log(woUDSA_res$BC_MLE) + qnorm(1 - 0.05 / 2) * se))

  # Methods with UDSA
  wUDSA_res = fit_SLDeepAssay_md(assay_summary = assay_summary)
  
  MLE_wUDSA = data.frame(Est = wUDSA_res$mle, 
                         SE = wUDSA_res$se, 
                         LB = wUDSA_res$ci[1], 
                         UB = wUDSA_res$ci[2])
  BCMLE_wUDSA = data.frame(Est = wUDSA_res$mle.bc,
                           SE = wUDSA_res$se, 
                           LB = wUDSA_res$ci.bc[1], 
                           UB = wUDSA_res$ci.bc[2])
  
  # Construct Message
  Message = ifelse(test = num_redo_all > 0, 
                   yes = 2, 
                   no = ifelse(test = num_redo_none > 0, 
                               yes = 1, 
                               no = 0))

  # Construct Message_Detailed
  if (Message > 0) {
    Message_Detailed = paste0("Assay was resimulated ",
                              num_redo_none + num_redo_all,
                              ifelse((num_redo_none + num_redo_all) == 1, "time (", "times ("),
                              num_redo_none,
                              " due to no DVL being detected and ",
                              num_redo_all,
                              " due to no well being negative for DVL).")
  } else {
    Message_Detailed = NULL
  }
  
  # Construct list to return
  return(list(Assay = assay$DVL_specific,
              MLE_woUDSA = MLE_woUDSA,
              BCMLE_woUDSA = BCMLE_woUDSA,
              MLE_wUDSA = MLE_wUDSA,
              BCMLE_wUDSA = BCMLE_wUDSA,
              Message = Message,
              Message_Detailed = Message_Detailed
  ))
}