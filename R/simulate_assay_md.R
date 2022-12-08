#' Simulate multiple-dilution assay data
#' @name simulate_assay_md
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
simulate_assay_md = function(M, tau, q, dilutions, remove_undetected = T) {

  D = length(dilutions)
  n = length(tau)

  # Observed matrix for each dilution level
  assay = lapply(1:D,
                 X = function(d) simulate_assay_sd(M = M[d],
                                                   n = n,
                                                   lambda = tau * dilutions[d],
                                                   q = q[d],
                                                   remove_undetected = remove_undetected,
                                                   seed = seed))
  assay.summary = vapply(X = 1:D,
                         FUN.VALUE = numeric(7 + n),
                         FUN = function(d) {
                           M = ncol(assay[[d]]) # number of wells
                           n = nrow(assay[[d]]) # number of DVLs detected
                           MN = sum(colSums(assay[[d]]) == 0, na.rm = TRUE) # number of p24-negative wells
                           MP = M - MN # number of p24-positive wells = MP
                           m = MP - sum(is.na(colSums(assay[[d]]))) # number of deep-sequenced wells = m
                           q = ifelse(MP == 0, 0, m / MP) # proportion of p24-positive wells deep sequenced
                           Y = rowSums(assay[[d]], na.rm = TRUE) # number of infected wells per DVL = Y
                           return((c("dilution" = dilutions[d], "M"=M, "n"=n,
                                     "MN"=MN, "MP"=MP, "m"=m, "q"=q, "Y"=Y)))
                         })

  assay.summary = as.data.frame(t(assay.summary))

  ### remove undetected DVL and adjust n
  if (remove_undetected) {
    assay.summary = assay.summary[, !grepl("Y", colnames(assay.summary)) | colSums(assay.summary) > 0]
    assay.summary$n = sum(grepl("Y", colnames(assay.summary)))
  }

  return(assay.summary)
}
