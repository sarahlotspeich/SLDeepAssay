# Build artifically "complete" dataset 
get_complete_data = function(Wstar, Zstar) {
  # Save useful constants 
  M = length(Wstar) ## number of replicate wells 
  MP = sum(Wstar) ## number p24-positive replicate wells 
  n = nrow(Zstar) ## number of DVLs detected
  m = sum(!is.na(colSums(Zstar))) ## number of p24-positive, sequenced replicate wells 
  
  if (m > 0) {
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
    
    # Augment complete data with observed data 
    ## For sequenced wells 
    cd_seq_long = cbind(j = rep(x = 1:m, each = nrow(cd_seq)),
                        cd_seq[rep(x = 1:nrow(cd_seq), times = m), ], 
                        Zstar_seq[rep(x = 1:m, each = nrow(cd_seq)), , drop = F], 
                        wstar = rep(x = Wstar[1:m], each = nrow(cd_seq)))
  } else {
    cd_seq_long = data.frame()
  }
  
  if ((M - m) > 0) {
    ## For unsequenced wells 
    cd_unseq_long = data.frame(cbind(j = rep(x = (m + 1):M, each = 2), 
                                     w = rep(x = c(0, 1), times = (M - m)),
                                     wstar = rep(x = Wstar[(m + 1):M], each = 2)))
  } else {
    cd_unseq_long = data.frame(j = NA, w = NA, wstar = NA)
  }
  
  return(list(complete_seq = cd_seq_long, 
              complete_unseq = cd_unseq_long))
}

# Get marginal probabilities of W
## (i) P(W = 0) = exp(-Lambda) (HIV-negative)
## (ii) P(W = 1) = 1 - exp(-Lambda) (HIV-positive)
get_pW = function(complete_data, l) {
  Lambda = sum(l) ## sum over DVL-specific rate parameters 
  (1 - exp(-Lambda)) ^ complete_data[, "w"] * exp(-Lambda) ^ (1 - complete_data[, "w"])  
}

# Get conditional probabilities of W* given W
## (i) P(W* = 0|W = 0) = spec (true negative)
## (ii) P(W* = 1|W = 0) = 1 - spec (false positive)
## (iii) P(W* = 1|W = 1) = sens (true positive)
## (iv) P(W* = 0|W = 1) = 1 - sens (false negative)
get_pWstarGivW = function(complete_data, sens, spec) {
  (1 - complete_data[, "wstar"]) * spec ^ (1 - complete_data[, "w"]) * (1 - sens) ^ complete_data[, "w"] + ### If W* = 0
    complete_data[, "wstar"] * sens ^ complete_data[, "w"] * (1 - spec) ^ (1 - complete_data[, "w"]) ### If W* = 1
  #(1 - complete_data[, "wstar"]) * spec ^ (1 - complete_data[, "w"]) * (1 - spec) ^ complete_data[, "w"] + 
  #  complete_data[, "wstar"] * sens ^ complete_data[, "w"] * (1 - sens) ^ (1 - complete_data[, "w"])
}

# Get probabilities of Zi* given Zi
## (i) P(Zi* = 0|Zi = 0) = spec (true negative)
## (ii) P(Zi* = 1|Zi = 0) = 1 - spec (false positive)
## (iii) P(Zi* = 1|Zi = 1) = sens (true positive)
## (iv) P(Zi* = 0|Zi = 1) = 1 - sens (false negative)
get_pZstar_iGivZ_i = function(z_i, zstar_i, sens, spec) {
  (1 - zstar_i) * spec ^ (1 - z_i) * (1 - sens) ^ z_i + 
    zstar_i * sens ^ z_i * (1 - spec) ^ (1 - z_i)
  # (1 - zstar_i) * spec ^ (1 - z_i) * (1 - spec) ^ z_i + 
  #   zstar_i * sens ^ z_i * (1 - sens) ^ (1 - z_i)
}

# Get joint probability of (Z1*,...,Zn*) given (Z1,...,Zn)
## P(Z1*,...,Zn*|Z1, ..., Zn) = P(Z1*|Z1)...P(Zn*|Zn)
get_pZstarGivZ = function(complete_data, n, sens, spec) {
  pZstarGivZ = 1
  for (i in 1:n) {
    pZstarGivZ = pZstarGivZ * 
      get_pZstar_iGivZ_i(z_i = complete_data[, paste0("z", i)], 
                         zstar_i = complete_data[, paste0("zstar", i)],
                         sens = sens, 
                         spec = spec)
  }
  return(pZstarGivZ)
}

# Get probabilities of Zi
## (i) P(Z_i = 1) = 1 - exp(-lambda_i)
## (ii) P(Z_i = 0) = exp(-lambda_i)
get_pZi = function(z_i, lambda_i) {
  (1 - exp(-lambda_i)) ^ z_i * exp(-lambda_i) ^ (1 - z_i)  
}

# Get joint probability of Z1,...,Zn
## P(Z_1, ..., Z_n) = P(Z_1)...P(Z_n)
get_pZ = function(complete_data, l) {
  pZ = 1
  for (i in 1:length(l)) {
    pZ = pZ * 
      get_pZi(z_i = complete_data[, paste0("z", i)], 
              lambda_i = l[i])
  }
  return(pZ)
}