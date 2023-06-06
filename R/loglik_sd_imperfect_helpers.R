# Get conditional probabilities of W* given W
## (i) P(W* = 0|W = 0) = spec (true negative)
## (ii) P(W* = 1|W = 0) = 1 - spec (false positive)
## (iii) P(W* = 1|W = 1) = sens (true positive)
## (iv) P(W* = 0|W = 1) = 1 - sens (false negative)
get_pWstarGivW = function(complete_data, sens, spec) {
  (1 - complete_data[, "wstar"]) * spec ^ (1 - z) * (1 - sens) ^ z + ### If Z* = 0
    complete_data[, "wstar"] * sens ^ z * (1 - spec) ^ (1 - z)
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