get_complete_data = function(Wstar, Zstar) {
  # Transform Z* in preparation to augment with cd
  Zstar_seq = t(Zstar) ## transpose to get rows per well, cols per DVL
  if (ncol(Zstar_seq) == 1) {
    Zstar_seq = matrix(data = Zstar_seq[complete.cases(Zstar_seq), drop = FALSE], nrow = 1) ## keep only sequenced wells (rows) 
  } else {
    Zstar_seq = Zstar_seq[complete.cases(Zstar_seq), drop = FALSE] ## keep only sequenced wells (rows) 
  }
  colnames(Zstar_seq) = paste0("zstar", 1:ncol(Zstar_seq)) ## name cols {z*1, ..., z*n}
  
  # Save useful constants 
  M = length(Wstar) ## number of replicate wells 
  MP = sum(Wstar) ## number p24-positive replicate wells 
  m = nrow(Zstar_seq) ## number of p24-positive, sequenced replicate wells 
  n = ncol(Zstar_seq) ## number of DVLs detected
  
  # Create matrix with "complete data" for sequenced wells
  ## Includes all possible combinations of true W and Z
  cd_seq = expand.grid(replicate(n = n + 1, expr = c(0, 1), simplify = F))
  colnames(cd_seq) = c("w", paste0("z", 1:n))  ## name cols {w, z1, ..., zn}
  
  # Create matrix with "complete data" for unsequenced wells
  ## Includes all possible combinations of true W and Z and imperfect Z*
  cd_unseq = expand.grid(replicate(n = 2 * n + 1, expr = c(0, 1), simplify = F))
  colnames(cd_unseq) = c("w", paste0("z", 1:n), paste0("zstar", 1:n)) ## name cols {w, z1, ..., zn, z*1, ..., z*n}
  
  # Augment complete data with observed data 
  ## For sequenced wells 
  cd_seq_long = cbind(j = rep(x = 1:m, each = nrow(cd_seq)),
                      wstar = rep(x = Wstar[1:m], each = nrow(cd_seq)),
                      Zstar_seq[rep(x = 1:m, each = nrow(cd_seq)), ], 
                      cd_seq[rep(x = 1:nrow(cd_seq), times = m), ])
  ## For unsequenced wells 
  cd_unseq_long = cbind(j = rep(x = (m + 1):M, each = nrow(cd_unseq)), 
                        wstar = rep(x = Wstar[(m + 1):M], each = nrow(cd_unseq)),
                        cd_unseq[rep(x = 1:nrow(cd_unseq), times = (M - m)), ])
  ## Stack them for all wells 
  cd_long = rbind(cd_seq_long, 
                  cd_unseq_long)
  
  return(cd_long)
}