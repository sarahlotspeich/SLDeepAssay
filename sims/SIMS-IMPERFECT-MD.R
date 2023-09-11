# Load package 
library(magrittr) ## for pipes 
library(SLDeepAssay, ## for IUPM estimators
        lib.loc = "/home/lotspes/R/x86_64-pc-linux-gnu-library/4.0") 

# Be reproducible
args = commandArgs(TRUE)
sim_seed = 11422 + as.integer(args)
set.seed(sim_seed)

# Number of replicates per simulation setting
num_reps = 100
## Note: You may want to run this code was run in parallel a cluster instead of locally, as it can be slow.

# Define parameters that vary over simulation settings (same as Section 3.1)
spec = 0.9 ## Sensitivity of assays
sens = seq(0.8, 1, by = 0.1) ## Specificity of assays
M = c(6, 12, 18) ## Total number of wells
n = 6 ## Number of existing DVLs
Tau = 1 ## Overall IUPM (split among n DLVs)
u = c(0.5, 1, 2) ## Dilutions in millions of cells per well
q = c(0, 0.5, 1) ## Proportion of p24-positive wells to undergo UDSA
constant_Tau = TRUE ## Indicator of constant IUPM across n DVLs

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
                 }) %>%  
  as.data.frame() %>%
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
  if (setting_row["constant_Tau"] == 1) {
    tau = rep(x = as.numeric(setting_row["Tau"]) / as.numeric(setting_row["n"]), 
              times = as.numeric(setting_row["n"]))
  } else {
    tau = c(rep(as.numeric(setting_row["Tau"]) / (2 * as.numeric(setting_row["n"])), 
                as.numeric(setting_row["n"]) / 2), 
            rep(3 * as.numeric(setting_row["Tau"]) / (2 * as.numeric(setting_row["n"])), 
                as.numeric(setting_row["n"]) / 2))
  }
  temp = simulate_assay_md_imperfect(M = M,
                                     tau = tau,
                                     q = q,
                                     u = u,
                                     sens_QVOA = as.numeric(setting_row["sensQVOA"]), 
                                     spec_QVOA = as.numeric(setting_row["specQVOA"]), 
                                     sens_UDSA = as.numeric(setting_row["sensUDSA"]), 
                                     spec_UDSA = as.numeric(setting_row["specUDSA"])) 

  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  # New likelihood (corrected IUPM estimator)
  fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                      u = u, 
                                      sens_QVOA = as.numeric(setting_row["sensQVOA"]), 
                                      spec_QVOA = as.numeric(setting_row["specQVOA"]), 
                                      sens_UDSA = as.numeric(setting_row["sensUDSA"]), 
                                      spec_UDSA = as.numeric(setting_row["specUDSA"]),
                                      lb = 1E-6)
  setting_row["Lambda"] = fit1$mle
  setting_row[c("conv", "msg")] = with(fit1, c(convergence, message))
  
  # Original likelihood (naive IUPM estimator)
  fit2 = fit_SLDeepAssay_md(assay = temp,
                            u = u, 
                            corrected = FALSE)
  # fit2 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
  #                                     u = u, 
  #                                     sens_QVOA = 1, 
  #                                     spec_QVOA = 1, 
  #                                     sens_UDSA = 1, 
  #                                     spec_UDSA = 1,
  #                                     lb = 1E-6)
  setting_row["Lambda_naive"] = fit2$mle
  # setting_row[c("conv_naive", "msg_naive")] = with(fit2, c(convergence, message))
  
  return(list(assay = temp, 
              result = setting_row))
}

# Be reproducible
Results = data.frame()
Assays = list()
for (i in 1:nrow(Settings)) {
  new_sim = one_sim(setting_row = Settings[i, ])
  Assays[[length(Assays) + 1]] = new_sim$assay
  Results = rbind(Results, new_sim$result)
  saveRDS(Assays, paste0("md-imperfect/Assays-md-imperfect-seed", sim_seed))
  write.csv(Results, paste0("md-imperfect/md-imperfect-seed", sim_seed, ".csv"), row.names = F)
  if (i %% 100 == 0) print(paste0("Sim ", i, " complete (", round(100 * i/nrow(Settings)), "%)"))
}
