loglik_sd_imperfect = function(l, complete_data) {
  # Log-likelihood contributions
  ## Sequenced wells 
  if (!any(is.na(complete_data$complete_seq$j)) & nrow(complete_data$complete_seq) > 0) {
    ### Calculate P(Z) for all rows in sequenced wells' complete data 
    complete_data$complete_seq$pZ = get_pZ(complete_data = complete_data$complete_seq, 
                                           l = l) 
    ### Sum P(W*|W)P(Z*|Z)P(Z) over Z to get P(W*, Z*)
    complete_data$complete_seq$sum_to_pObs = with(complete_data$complete_seq, 
                                                  pWstarGivW * pZstarGivZ * pZ)
    pObs = rowsum(x = complete_data$complete_seq[, "sum_to_pObs"], 
                  group = complete_data$complete_seq[, "j"])
    ### Replace pObs = 0 with 1 to avoid infinite log-likelihood 
    pObs[pObs == 0] = 1
    ### Re-calculate log-likelihood 
    ll = sum(log(pObs))
  } else {
    ll = 0
  }
  
  ## Unsequenced wells
  if (!any(is.na(complete_data$complete_unseq$j)) & nrow(complete_data$complete_unseq) > 0) {
    ### Calculate P(W) for all rows in unsequenced wells' complete data 
    complete_data$complete_unseq$pW = get_pW(complete_data = complete_data$complete_unseq, 
                                             l = l) 
    ### Sum P(W*|W)P(W) over W to get P(W*)
    complete_data$complete_unseq$sum_to_pObs = with(complete_data$complete_unseq, 
                                                    pWstarGivW * pW)
    pObs = rowsum(x = complete_data$complete_unseq[, "sum_to_pObs"], 
                  group = complete_data$complete_unseq[, "j"])
    # Replace pObs = 0 with 1 to avoid infinite log-likelihood 
    pObs[pObs == 0] = 1
    # Re-calculate log-likelihood 
    ll = ll + sum(log(pObs))
  } 
  
  return(-ll)
}