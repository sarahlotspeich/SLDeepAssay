loglik_md_imperfect = function(tau, u, complete_data_md) {
  # Total number of dilutions
  D = length(complete_data_md)
  
  # Compute negative log-likelihoods for each dilution
  logliks_by_dilution = vapply(X = 1:D,
                               FUN.VALUE = numeric(1),
                               FUN = function(d) {
                                as.numeric(
                                  loglik_sd_imperfect(
                                    l = u[d] * tau,
                                    complete_data = complete_data_md[[d]]
                                    )
                                  )
                                 }
                               )
  
  # Return the sum of (negative) log-likelihoods
  ll = sum(logliks_by_dilution)
  
  return(ll)
}
