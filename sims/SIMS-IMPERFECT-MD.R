library(SLDeepAssay)

# Number of replicates per simulation setting
num_reps = 50 
## Note: This code was run in parallel a cluster instead of locally, as it can be slow. 

# Define parameters that vary over simulation settings (same as Section 3.1)
spec = 0.9 # Sensitivity of assays
sens = seq(0.8, 1, by = 0.1) # Specificity of assays
M = c(6, 12, 18) # Total number of wells
n = 6 # Number of existing DVLs
Tau = 1 # Overall IUPM (split among n DLVs)
u = c(0.5, 1, 2) # Dilutions in millions of cells per well
q = c(0, 0.5, 1) # Proportion of p24-positive wells to undergo UDSA
constant_Tau = TRUE # Indicator of constant IUPM across n DVLs

# Number of simulation settings
num_sett = length(sens) ^ 2 * length(spec) ^ 2 * length(n) * length(Tau) * length(constant_Tau) 

# Create dataframe of different simulation settings
Settings = apply(X = expand.grid("sensQVOA" = sens,
                                 "specQVOA" = spec,
                                 "sensUDSA" = sens,
                                 "specUDSA" = spec,
                                 "n" = n,
                                 "Tau" = Tau,
                                 "constant_Tau" = constant_Tau),
                 MARGIN = 2,
                 FUN = function(x) {
                   rep(x, each = num_reps)
                 }
) |>  as.data.frame() |>
  dplyr::mutate(sim_id = rep(x = seq(1, num_reps), times = num_sett),
                Lambda = NA, 
                conv = NA, 
                msg = NA,
                Lambda_naive = NA, 
                conv_naive = NA, 
                msg_naive = NA,
                assay_resampled = NA)

# Function one_sim() simulates one assay and returns the output (reshaped)
one_sim = function(setting_row) {
  # Save parameter values from setting_row
  n = as.numeric(setting_row["n"])
  Tau = as.numeric(setting_row["Tau"])
  constant_Tau = setting_row["constant_Tau"] == 1
  sensQVOA = as.numeric(setting_row["sensQVOA"])
  specQVOA = as.numeric(setting_row["specQVOA"])
  sensUDSA = as.numeric(setting_row["sensUDSA"])
  specUDSA = as.numeric(setting_row["specUDSA"])
  
  if (constant_Tau) {
    tau = rep(x = Tau / n, times = n)
  } else {
    tau = c(rep(Tau / (2 * n), n / 2), rep(3 * Tau / (2 * n), n / 2))
  }
  temp = simulate_assay_md_imperfect(M = M,
                                     tau = tau,
                                     q = q,
                                     u = u,
                                     sens_QVOA = sensQVOA, 
                                     spec_QVOA = specQVOA, 
                                     sens_UDSA = sensUDSA, 
                                     spec_UDSA = specUDSA) 

  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  # New likelihood (corrected IUPM estimator)
  fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                      u = u, 
                                      sens_QVOA = sensQVOA, 
                                      spec_QVOA = specQVOA, 
                                      sens_UDSA = sensUDSA, 
                                      spec_UDSA = specUDSA)
  setting_row[c("Lambda", "conv", "msg")] = with(fit1, c(mle, convergence, message))
  
  # Original likelihood (naive IUPM estimator)
  assay_summary = vapply(X = 1:length(u),
                         FUN.VALUE = numeric(7 + n),
                         FUN = function(d) {
                           M = ncol(temp[[d]]$DVL_specific) # number of wells
                           n = nrow(temp[[d]]$DVL_specific) # number of DVLs detected
                           MN = sum(colSums(temp[[d]]$DVL_specific) == 0, na.rm = TRUE) # number of p24-negative wells
                           MP = M - MN # number of p24-positive wells = MP
                           m = MP - sum(is.na(colSums(temp[[d]]$DVL_specific))) # number of deep-sequenced wells = m
                           q = ifelse(MP == 0, 0, m / MP) # proportion of p24-positive wells deep sequenced
                           Y = rowSums(temp[[d]]$DVL_specific, na.rm = TRUE) # number of infected wells per DVL = Y
                           return((c("u" = u[d], "M"=M, "n"=n,
                                     "MN"=MN, "MP"=MP, "m"=m, "q"=q, "Y"=Y)))
                         })
  assay_summary = as.data.frame(t(assay_summary))
  fit2 = fit_SLDeepAssay_md(assay_summary = assay_summary, 
                            corrected = FALSE)
  setting_row["Lambda_naive"] = sum(fit2$mle)
  return(setting_row)
}

# Be reproducible
sim_seed = 11422 
set.seed(sim_seed)

Results = data.frame()
for (i in 1:nrow(Settings)) {
  Results = rbind(Results, one_sim(setting_row = Settings[i, ]))
  write.csv(Results, "md_imperfect.csv", row.names = F)
  if (i %% 10 == 0) print(paste0("Sim ", i, " complete (", round(100 * i/nrow(Settings)), "%)"))
}
