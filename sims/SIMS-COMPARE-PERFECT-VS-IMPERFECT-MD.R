library(SLDeepAssay)

# Number of replicates per simulation setting
num_reps = 100

# Define parameters that vary over simulation settings (same as Section 3.1)
M = c(6, 12, 18) # Total number of wells
n = 6 # Number of existing DVLs
Tau = 1 # Overall IUPM (split among n DLVs)
u = c(0.5, 1, 2) # Dilutions in millions of cells per well
q = c(0, 0.5, 1) # Proportion of p24-positive wells to undergo UDSA

# Number of simulation settings
num_sett = length(n) * length(Tau)

# Create dataframe of different simulation settings
Settings = apply(X = expand.grid("n" = n,
                                 "Tau" = Tau),
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
  tau = rep(x = Tau / n, times = n)
  temp = simulate_assay_md_imperfect(M = M,
                                     tau = tau,
                                     q = q,
                                     u = u,
                                     sens_QVOA = 1, 
                                     spec_QVOA = 1, 
                                     sens_UDSA = 1, 
                                     spec_UDSA = 1) 
  saveRDS(temp, "~/Downloads/temp")
  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  # New likelihood (corrected IUPM estimator)
  fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                      u = u, 
                                      sens_QVOA = 1, 
                                      spec_QVOA = 1, 
                                      sens_UDSA = 1, 
                                      spec_UDSA = 1,
                                      lb = 0) #lb = 1E-6)
  setting_row[c("Lambda", "conv", "msg")] = with(fit1, c(as.numeric(mle), convergence, message))
  
  # Original likelihood (naive IUPM estimator)
  fit2 = fit_SLDeepAssay_md(assay = temp, 
                            u = u, 
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
  write.csv(Results, "~/Downloads/md_compare_perfect_imperfect.csv", row.names = F)
  if (i %% 10 == 0) print(paste0("Sim ", i, " complete (", round(100 * i/nrow(Settings)), "%)"))
}

# Check average difference between perfect/imperfect MLE code
with(Results, mean(as.numeric(Lambda) - Lambda_naive)) ## 1.366463e-05

# Check whether each replicate has a unique estimate
length(unique(Results$Lambda)) ## 100
length(unique(Results$Lambda_naive)) ## 100 