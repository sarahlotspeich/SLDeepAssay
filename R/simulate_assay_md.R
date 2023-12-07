#' Simulate multiple-dilution assay data
#' @name simulate_assay_md
#' @param M Total number of wells at each dilution level (a vector of length D).
#' @param tau DVL-specific IUPMs (a vector of length \code{n}). (Note: All elements in \code{tau} must be > 0.)
#' @param q Proportions of p24-positive wells that underwent UDSA at each dilution level (a vector of length D).
#' @param u Vector of dilution levels (in millions of cells per well).
#' @param k Overdispersion parameter (a positive number). Default is \code{Inf}, which corresponds to no overdispersion.
#' @param remove_undetected Logical, if \code{remove_undetected = TRUE} (the default), then DVL which were not detected in any of the deep sequenced wells across all dilution levels are deleted.
#' @return A data frame with summary data of the simulated results. This summary contains one row per dilution level and the following columns: M (total number of wells), n (number of distinct viral lineages \[DVL\]), MN (number of p24-negative wells), m (number of deep sequenced wells), Y1,..., Yn (counts of wells positive for DVL i, (i = 1,...,n), and dilutions (dilution levels, in millions of cells per well).
#' @export
#'
simulate_assay_md = function(M, tau, q, u, k = Inf, remove_undetected = T, seed = NULL) {

  D = length(u)
  n = length(tau)

  # Observed matrix for each dilution level
  assay = lapply(X = 1:D,
                 FUN = function(d) simulate_assay_sd(M = M[d],
                                                     tau = tau,
                                                     u = u[d],
                                                     q = q[d],
                                                     k = k,
                                                     remove_undetected = F))
  assay_summary = vapply(X = 1:D,
                         FUN.VALUE = numeric(7 + n),
                         FUN = function(d) {
                           M = ncol(assay[[d]]$DVL_specific) # number of wells
                           n = nrow(assay[[d]]$DVL_specific) # number of DVLs detected
                           MN = sum(colSums(assay[[d]]$DVL_specific) == 0, na.rm = TRUE) # number of p24-negative wells
                           MP = M - MN # number of p24-positive wells = MP
                           m = MP - sum(is.na(colSums(assay[[d]]$DVL_specific))) # number of deep-sequenced wells = m
                           q = ifelse(MP == 0, 0, m / MP) # proportion of p24-positive wells deep sequenced
                           Y = rowSums(assay[[d]]$DVL_specific, na.rm = TRUE) # number of infected wells per DVL = Y
                           return((c("u" = u[d], "M"=M, "n"=n,
                                     "MN"=MN, "MP"=MP, "m"=m, "q"=q, "Y"=Y)))
                         })

  assay_summary = as.data.frame(t(assay_summary))

  ### remove undetected DVL and adjust n
  if (remove_undetected) {
    assay_summary = assay_summary[, !grepl("Y", colnames(assay_summary)) | colSums(assay_summary) > 0]
    assay_summary$n = sum(grepl("Y", colnames(assay_summary)))
  }

  return(assay_summary)
}

