#' Simulate single-dilution assay data and fit new and existing methods to it.
#' @name simulate_SLDeepAssay_sd
#' @param M Total number of wells (a scalar).
#' @param n Total number of existing distinct viral lineages (DVL) (a scalar).
#' @param lambda DVL-specific rates of infection (a vector of length \code{n}). (Note: All elements in \code{lambda} must be > 0.)
#' @param q Fixed proportion of p24-positive wells to be deep sequenced (a scalar between 0 and 1).
#' @param dilution Number of cells per well (a scalar, in millions). Default is \code{dilution = 1}.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells are deleted.
#' @param seed (Optional) An integer setting the random seed used to simulate the assay data. Default is \code{seed = NULL}.
#' @return Named list with the following slots:
#' \item{Assay}{Simulated single-dilution assay data (deep sequenced).}
#' \item{MLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (without deep sequencing information).}
#' \item{BCMLE_woUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (without deep sequencing information).}
#' \item{MLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the MLE (with deep sequencing information).}
#' \item{BCMLE_wUDSA}{Point estimate, standard error estimate, and 95\% confidence interval for the bias-corrected MLE (with deep sequencing information).}
#' \item{Message}{(Optional) Message describing whether the \code{Assay} needed to be re-simulated due to either none of the DVL being detected (\code{Message = 1}) or at least one DVL being detected in all deep-sequenced wells (\code{Message = 2}). Otherwise, \code{Message = 0}.}
#' @importFrom SLDAssay get.mle
#' @export
#'
simulate_SLDeepAssay_sd <- function(M, n, lambda, q, dilution = 1, remove_undetected = TRUE, seed = NULL) {
  # If supplied, set the random seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create indicators of whether data need to be re-simulated
  num_redo_all = 0 # Due to any DVL being detected in all wells or
  num_redo_none = 0 # No DVL being detected in any wells

  # Simulate single-dilution assay data
  assay = simulate_assay_sd(M = M,
                            n = n,
                            lambda = lambda,
                            q = q,
                            remove_undetected = remove_undetected)

  # Check for need to re-simulate
  prop_DVL_detected = tryCatch(expr = rowMeans(assay$DVL_specific, na.rm = TRUE),
                               error = function(e) 0)
  while (any(prop_DVL_detected == 1) | is.null(assay$DVL_specific)) {
    if (is.null(assay)) {
      num_redo_none = num_redo_none + 1
    } else {
      num_redo_all = num_redo_all + 1
    }
    # Re-simulate single-dilution assay data
    assay = simulate_assay_sd(M = M,
                              n = n,
                              lambda = lambda,
                              q = q,
                              remove_undetected = remove_undetected)
    # Check for need to continue re-simulating
    prop_DVL_detected = tryCatch(expr = rowMeans(assay$DVL_specific, na.rm = TRUE),
                                 error = function(e) 0)
  }

  # Methods without UDSA
  woUDSA_res = get.mle(pos = sum(assay$any_DVL),
                       replicates = length(assay$any_DVL),
                       dilutions = dilution * 1E6)
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
  wUDSA_res = fit_SLDeepAssay_sd(assay$DVL_specific,
                                 dilution = dilution)

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

  return(list(Assay = assay$DVL_specific,
              MLE_woUDSA = MLE_woUDSA,
              BCMLE_woUDSA = BCMLE_woUDSA,
              MLE_wUDSA = MLE_wUDSA,
              BCMLE_wUDSA = BCMLE_wUDSA,
              Message = Message
              ))
}
