loglik_sd_imperfect2 = function(params, complete_data) {
  #print(paste0("alphas = ", paste0(params[1:4], collapse = ",")))
  #print(paste0("lambdas = ", paste0(params[-c(1:4)], collapse = ",")))
  
  # Calculate P(W*|W) for all rows in complete data 
  ## < this won't change with lambda, so calculate once >
  pWstarGivW = get_pWstarGivW(complete_data = complete_data,
                              sens = params[1], 
                              spec = params[2])
  
  # Calculate P(Z*|Z) for all rows in complete data 
  ## < this won't change with lambda, so calculate once >
  pZstarGivZ = get_pZstarGivZ(complete_data = complete_data,
                              n = length(params[-c(1:4)]),
                              sens = params[3], 
                              spec = params[4])
  
  # Calculate P(Z) for all rows in complete data 
  ## < this will change with lambda, so re-calculate within each iteration >
  pZ = get_pZ(complete_data = complete_data, 
              l = params[-c(1:4)]) 
  
  # Sum P(W*|W)P(Z*|Z)P(Z) over Z to get P(W*, Z*)
  complete_data$sum_to_pObs = pWstarGivW * pZstarGivZ * pZ
  pObs = rowsum(x = complete_data[, "sum_to_pObs"], group = complete_data[, "j"])
  
  # Replace pObs = 0 with 1 to avoid infinite log-likelihood 
  pObs[pObs == 0] = 1
  
  # Re-calculate log-likelihood 
  ll = sum(log(pObs))
  
  #print(paste0("loglik = ", ll))
  
  return(-ll)
}

loglik_sd_imperfect = function(l, complete_data, pWstarGivW, pZstarGivZ) {
  # Calculate P(Z) for all rows in complete data 
  ## < this will change with lambda, so re-calculate within each iteration >
  pZ = get_pZ(complete_data = complete_data, 
              l = l) 
  
  # Sum P(W*|W)P(Z*|Z)P(Z) over Z to get P(W*, Z*)
  complete_data$sum_to_pObs = pWstarGivW * pZstarGivZ * pZ
  pObs = rowsum(x = complete_data[, "sum_to_pObs"], group = complete_data[, "j"])
  
  # Replace pObs = 0 with 1 to avoid infinite log-likelihood 
  pObs[pObs == 0] = 1
  
  # Re-calculate log-likelihood 
  ll = sum(log(pObs))
  
  return(-ll)
}