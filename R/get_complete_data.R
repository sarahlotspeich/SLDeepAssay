get_complete_data = function(Wstar, Zstar) {
  # Save useful constants 
  M = length(Wstar) ## number of replicate wells 
  MP = sum(Wstar) ## number p24-positive replicate wells 
  n = nrow(Zstar) ## number of DVLs detected
  m = sum(!is.na(colSums(Zstar))) ## number of p24-positive, sequenced replicate wells 
  
  # Transform Z* in preparation to augment with cd
  Zstar_seq = t(Zstar) ## transpose to get rows per well, cols per DVL
  ## then keep only sequenced wells (rows) 
  if (n == 1) { ### if only one DVL detected, need to force Zstar_seq to a column
    Zstar_seq = matrix(data = Zstar_seq[complete.cases(Zstar_seq), drop = FALSE], 
                       ncol = 1) 
  } else if (m == 1) { ### if only one well sequenced, need to force Zstar_seq to a row
    Zstar_seq = matrix(data = Zstar_seq[complete.cases(Zstar_seq), drop = FALSE], 
                       nrow = 1) 
  } else {
    Zstar_seq = Zstar_seq[!is.na(colSums(Zstar)), ] ## keep only sequenced wells (rows) 
  }
  colnames(Zstar_seq) = paste0("zstar", 1:n) ## name cols {z*1, ..., z*n}

  # Create matrix with "complete data" for sequenced wells
  ## Includes all possible combinations of true Z
  cd_seq = expand.grid(replicate(n = n, expr = c(0, 1), simplify = F))
  colnames(cd_seq) = paste0("z", 1:n)  ## name cols {z1, ..., zn}
  ## Compute true W from true Z
  cd_seq$w = as.numeric(rowSums(cd_seq) > 0)
  
  # Create matrix with "complete data" for unsequenced wells
  ## Includes all possible combinations of true Z and imperfect Z*
  cd_unseq = expand.grid(replicate(n = 2 * n, expr = c(0, 1), simplify = F))
  colnames(cd_unseq) = c(paste0("z", 1:n), paste0("zstar", 1:n)) ## name cols {z1, ..., zn, z*1, ..., z*n}
  ## Compute true W from true Z
  if (n > 1) {
    cd_unseq$w = as.numeric(rowSums(cd_unseq[, paste0("z", 1:n)]) > 0)  
  } else {
    cd_unseq$w = as.numeric(cd_unseq[, paste0("z", 1:n)] > 0)
  }
  
  # Augment complete data with observed data 
  ## For sequenced wells 
  cd_seq_long = cbind(j = rep(x = 1:m, each = nrow(cd_seq)),
                      cd_seq[rep(x = 1:nrow(cd_seq), times = m), ], 
                      Zstar_seq[rep(x = 1:m, each = nrow(cd_seq)), , drop = F], 
                      wstar = rep(x = Wstar[1:m], each = nrow(cd_seq)))
  ## For unsequenced wells 
  cd_unseq_long = cbind(j = rep(x = (m + 1):M, each = nrow(cd_unseq)), 
                        cd_unseq[rep(x = 1:nrow(cd_unseq), times = (M - m)), ],
                        wstar = rep(x = Wstar[(m + 1):M], each = nrow(cd_unseq)))
  ## Stack them for all wells 
  cd_long = rbind(cd_seq_long, 
                  cd_unseq_long)
  
  return(cd_long)
}