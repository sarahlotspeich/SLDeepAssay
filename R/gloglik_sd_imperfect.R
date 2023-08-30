#' Gradient of the log-likelihood for IUPM from single-dilution imperfect assay data
#' @name gloglik_sd_imperfect
#' @param l Vector of DVL-specific parameters.
#' @param M Total number of wells originally sequenced with the QVOA.
#' @param MPstar Number of p24-positive wells.
#' @param mstar Number of p24-positive wells that underwent the UDSA.
#' @param y A vector of  DVL-specific infection counts.
#' @return A vector
#'
gloglik_sd_imperfect = function(l, cd) {
  # Derivatives of the sum over unsequenced wells 
  ## Save Lambda = sum of l
  L = sum(l)
  
  ## Compute P(W) based on vector of lambdas l
  cd$complete_unseq$pW = get_pW(complete_data = cd$complete_unseq, 
                                l = l)
  
  ## Compute P(W*, W) = P(W*|W)P(W), to be summed over W to get P(W*)
  cd$complete_unseq$sum_to_pWstar = cd$complete_unseq$pWstarGivW * cd$complete_unseq$pW
  pWstar = rowsum(x = cd$complete_unseq$sum_to_pWstar, 
                  group = cd$complete_unseq$j, 
                  reorder = TRUE)
  
  ## Compute derivative of P(W*|W)P(W), to be summed over W to get derivative P'(W*)
  cd$complete_unseq$sum_to_d_pWstar = cd$complete_unseq$pWstarGivW * (- 1) ^ (1 - cd$complete_unseq$w) * exp(- L) 
  d_pWstar = rowsum(x = cd$complete_unseq$sum_to_d_pWstar, 
                    group = cd$complete_unseq$j, 
                    reorder = TRUE)
  
  ## To avoid NaN errors when dividing d_pWstar / pWstar, force replace 0 / 0 = 0
  pWstar[pWstar == 0] = 1
  
  ## Initialize gradient vector with P'(W*)/P(W*)
  ### (Note: This derivative is the same for all i)
  gradient = rep(x = sum(d_pWstar / pWstar), 
                 times = length(l)) 
  
  # Derivatives of the sum over unsequenced wells 
  ## Save n = # DVL detected
  n = length(l)
  
  ## Compute P(Z) based on vector of lambdas l
  cd$complete_seq$pZ = get_pZ(complete_data = cd$complete_seq, 
                              l = l)
  
  ## Multiply P(W*|W)P(Z*|Z)P(Z) outside, since it's inside the sum for all derivatives
  cd$complete_seq$prod_probs = cd$complete_seq$pWstarGivW * cd$complete_seq$pZstarGivZ * cd$complete_seq$pZ
  
  # Add derivatives of the sum over sequenced wells 
  for(i in 1:n)
  {
    ## Compute product over P(Zi) for all lambda except i
    cd$complete_seq$pZ_exclDVLi = get_pZ(complete_data = cd$complete_seq, 
                                         l = l, 
                                         exclude_DVL = i)
    
    ## Compute derivative of P(Zi)
    cd$complete_seq$d_pZi = unlist((- 1) ^ (1 - cd$complete_seq[paste0("z", i)]) * exp(- l[i]))
    
    ## Compute P(W*, Z*, W, Z), to be summed over W and Z to get P(W*, Z*)
    cd$complete_seq$sum_to_pWstarZstar = cd$complete_seq$pWstarGivW * cd$complete_seq$pZstarGivZ * cd$complete_seq$pZ
    pWstarZstar = rowsum(x = cd$complete_seq$sum_to_pWstarZstar, 
                         group = cd$complete_seq$j, 
                         reorder = TRUE)
    
    ## Compute derivative of P(W*, Z*, W, Z), to be summed over W and Z to get derivative P'(W*, Z*)
    cd$complete_seq$sum_to_d_pWstarZstar = cd$complete_seq$prod_probs * cd$complete_seq$pZ_exclDVLi * cd$complete_seq$d_pZi
    d_pWstarZstar = rowsum(x = cd$complete_seq$sum_to_d_pWstarZstar, 
                           group = cd$complete_seq$j, 
                           reorder = TRUE)
    
    ## To avoid NaN errors when dividing d_pWstarZstar / pWstarZstar, force replace 0 / 0 = 0
    pWstarZstar[pWstarZstar == 0] = 1
    
    ## Add derivatives of log-likelihood for sequenced wells to gradient
    gradient[i] = gradient[i] + sum(d_pWstarZstar / pWstarZstar)
  }
  return(-gradient)
}
