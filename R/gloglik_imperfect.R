gloglik_sd_imperfect = function(l, complete_data) {
  # Derivatives of the sum over unsequenced wells 
  if (!any(is.na(complete_data$complete_unseq$j)) & nrow(complete_data$complete_unseq) > 0) {
    ## Save Lambda = sum of l
    L = sum(l)
    
    ## Compute P(W) based on vector of lambdas l
    complete_data$complete_unseq$pW = get_pW(complete_data = complete_data$complete_unseq, 
                                             l = l)
    
    ## Compute P(W*, W) = P(W*|W)P(W), to be summed over W to get P(W*)
    complete_data$complete_unseq$sum_to_pWstar = complete_data$complete_unseq$pWstarGivW * complete_data$complete_unseq$pW
    pWstar = rowsum(x = complete_data$complete_unseq$sum_to_pWstar, 
                    group = complete_data$complete_unseq$j, 
                    reorder = TRUE)
    
    ## Compute derivative of P(W*|W)P(W), to be summed over W to get derivative P'(W*)
    complete_data$complete_unseq$sum_to_d_pWstar = complete_data$complete_unseq$pWstarGivW * (- 1) ^ (1 - complete_data$complete_unseq$w) * exp(- L) 
    d_pWstar = rowsum(x = complete_data$complete_unseq$sum_to_d_pWstar, 
                      group = complete_data$complete_unseq$j, 
                      reorder = TRUE)
    
    ## To avoid NaN errors when dividing d_pWstar / pWstar, force replace 0 / 0 = 0
    pWstar[pWstar == 0] = 1
    
    ## Initialize gradient vector with P'(W*)/P(W*)
    ### (Note: This derivative is the same for all i)
    gradient = rep(x = sum(d_pWstar / pWstar), 
                   times = length(l)) 
  } else {
    gradient = rep(0, length(l))
  }
  
  if (!any(is.na(complete_data$complete_seq$j)) & nrow(complete_data$complete_seq) > 0) {
    # Derivatives of the sum over unsequenced wells 
    ## Save n = # DVL detected
    n = length(l)
    
    ## Compute P(Z) based on vector of lambdas l
    complete_data$complete_seq$pZ = get_pZ(complete_data = complete_data$complete_seq, 
                                           l = l)
    
    ## Multiply P(W*|W)P(Z*|Z) outside, since it's inside the sum for all derivatives
    complete_data$complete_seq$prod_probs = complete_data$complete_seq$pWstarGivW * complete_data$complete_seq$pZstarGivZ 
    
    # Add derivatives of the sum over sequenced wells 
    for(i in 1:n)
    {
      ## Compute product over P(Zi) for all lambda except i
      complete_data$complete_seq$pZ_exclDVLi = get_pZ(complete_data = complete_data$complete_seq, 
                                                      l = l, 
                                                      exclude_DVL = i)
      
      ## Compute derivative of P(Zi)
      complete_data$complete_seq$d_pZi = unlist((- 1) ^ (1 - complete_data$complete_seq[paste0("z", i)]) * exp(- l[i]))
      
      ## Compute P(W*, Z*, W, Z), to be summed over W and Z to get P(W*, Z*)
      complete_data$complete_seq$sum_to_pWstarZstar = complete_data$complete_seq$pWstarGivW * complete_data$complete_seq$pZstarGivZ * complete_data$complete_seq$pZ
      pWstarZstar = rowsum(x = complete_data$complete_seq$sum_to_pWstarZstar, 
                           group = complete_data$complete_seq$j, 
                           reorder = TRUE)
      
      ## Compute derivative of P(W*, Z*, W, Z), to be summed over W and Z to get derivative P'(W*, Z*)
      complete_data$complete_seq$sum_to_d_pWstarZstar = complete_data$complete_seq$prod_probs * complete_data$complete_seq$pZ_exclDVLi * complete_data$complete_seq$d_pZi
      d_pWstarZstar = rowsum(x = complete_data$complete_seq$sum_to_d_pWstarZstar, 
                             group = complete_data$complete_seq$j, 
                             reorder = TRUE)
      
      ## To avoid NaN errors when dividing d_pWstarZstar / pWstarZstar, force replace 0 / 0 = 0
      pWstarZstar[pWstarZstar == 0] = 1
      
      ## Add derivatives of log-likelihood for sequenced wells to gradient
      gradient[i] = gradient[i] + sum(d_pWstarZstar / pWstarZstar)
    }  
  }
  
  return(-gradient)
}

gloglik_md_imperfect = function(tau, u, complete_data_md) {
  # Number of DVLs
  n = length(tau)
  
  # Total number of dilutions
  D = length(complete_data_md)
  
  # Compute gradients of log-likelihoods for each dilution
  glogliks_by_dilution = vapply(X = 1:D,
                                FUN.VALUE = numeric(n),
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
  gradient = rowSums(glogliks_by_dilution)
  
  return(gradient)
}

