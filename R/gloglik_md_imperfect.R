gloglik_md_imperfect = function(tau, u, complete_data_md) {
  # Total number of dilutions
  D = length(complete_data_md)
  
  # Compute gradients of log-likelihoods for each dilution
  glogliks_by_dilution = vapply(X = 1:D,
                               FUN.VALUE = numeric(1),
                               FUN = function(d) {
                                 as.numeric(
                                   gloglik_sd_imperfect(
                                     l = u[d] * tau,
                                     complete_data = complete_data_md[[d]]
                                   )
                                 )
                               }
  )
  
  # Return the sum of gradients 
  gradient = rowSums(glogliks_byDilution)
  
  return(gradient)
}

