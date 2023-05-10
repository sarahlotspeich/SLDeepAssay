#' Simulate single-dilution assay data and fit new and existing methods to it.
#' @name simulate_SLDeepAssay_sd
#' @param M Total number of wells (a scalar).
#' @param tau Mean counts of cells per million infected with each DVL (a vector). (Note: All elements in \code{tau} must be > 0.)
#' @param q Fixed proportion of p24-positive wells to be deep sequenced (a scalar between 0 and 1).
#' @param u Dilution level in millions of cells per well (a positive scalar). Default is \code{u = 1}.
#' @param tpr_UDSA True positive rate for the UDSA (a scalar between 0 and 1). Default is \code{tpr_UDSA = 1}, corresponding to a "perfect" assay.
#' @param fpr_UDSA False positive rate for the UDSA (a scalar between 0 and 1). Default is \code{fpr_UDSA = 0}, corresponding to a "perfect" assay.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells are deleted.
#' @return Named list with the following slots:
#' \item{Assay}{Simulated single-dilution assay data (deep sequenced).}
#' \item{MLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (without deep sequencing information).}
#' \item{BCMLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (without deep sequencing information).}
#' \item{MLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (with deep sequencing information).}
#' \item{BCMLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (with deep sequencing information).}
#' \item{Message}{(Optional) Message describing whether the \code{Assay} needed to be re-simulated due to either none of the DVL being detected (\code{Message = 1}) or at least one DVL being detected in all deep-sequenced wells (\code{Message = 2}). Otherwise, \code{Message = 0}.}
#' \item{Message_Detailed}{(Optional) Count of times the \code{Assay} needed to be re-simulated due to none of the DVL being detected and/or at least one DVL being detected in all deep-sequenced wells. Otherwise, \code{Message_Detailed = NULL}.}
#' @importFrom SLDAssay get.mle
#' @export
#'
simulate_SLDeepAssay_sd <- function(M, tau, q, u = 1, tpr_UDSA = 1, fpr_UDSA = 0, remove_undetected = TRUE) {
  # Create indicators of whether data need to be re-simulated
  num_redo_all = 0 # Due to any DVL being detected in all wells or
  num_redo_none = 0 # No DVL being detected in any wells

  # Simulate single-dilution assay data
  assay = simulate_assay_sd(M = M,
                            tau = tau,
                            q = q,
                            u = u,
                            tpr_UDSA = tpr_UDSA, 
                            fpr_UDSA = fpr_UDSA,
                            remove_undetected = remove_undetected)

  # Checks for need to re-simulate
  ## Any DVL is detected in all deep-sequenced and p24 positive wells?
  prop_DVL_detected = tryCatch(expr = rowMeans(assay$DVL_specific, na.rm = TRUE),
                               error = function(e) 0)
  ## If so, continue resimulating until no longer true
  while (any(prop_DVL_detected == 1) | is.null(assay$DVL_specific)) {
    if (is.null(assay)) {
      num_redo_none = num_redo_none + 1
    } else {
      num_redo_all = num_redo_all + 1
    }
    # Re-simulate single-dilution assay data
    assay = simulate_assay_sd(M = M,
                              tau = tau,
                              q = q,
                              u = u,
                              remove_undetected = remove_undetected)
    # Check for need to continue re-simulating
    prop_DVL_detected = tryCatch(expr = rowMeans(assay$DVL_specific, na.rm = TRUE),
                                 error = function(e) 0)
  }

  # Methods without UDSA
  woUDSA_res = get.mle(pos = sum(assay$any_DVL),
                       replicates = length(assay$any_DVL),
                       dilutions = u * 1E6)
  ## Solve for the standard error estimator
  Z = qnorm(p = 1 - (0.05 / 2)) ## Critical value for 95% confidence interval
  se = woUDSA_res$MLE * (log(woUDSA_res$MLE) - log(woUDSA_res$Asymp_CI[1])) / Z
  MLE_woUDSA = data.frame(Est = woUDSA_res$MLE,
                          SE = se,
                          LB = woUDSA_res$Asymp_CI[1],
                          UB = woUDSA_res$Asymp_CI[2])
  BCMLE_woUDSA = data.frame(Est = woUDSA_res$BC_MLE,
                            SE = se,
                            LB = exp(log(woUDSA_res$BC_MLE) - Z * se),
                            UB = exp(log(woUDSA_res$BC_MLE) + Z * se))

  # Methods with UDSA
  wUDSA_res = fit_SLDeepAssay_sd(assay = assay$DVL_specific,
                                 u = u)

  MLE_wUDSA = data.frame(Est = wUDSA_res$mle,
                         SE = wUDSA_res$se,
                         LB = wUDSA_res$ci[1],
                         UB = wUDSA_res$ci[2])
  BCMLE_wUDSA = data.frame(Est = wUDSA_res$mle_bc,
                           SE = wUDSA_res$se,
                           LB = wUDSA_res$ci_bc[1],
                           UB = wUDSA_res$ci_bc[2])

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
