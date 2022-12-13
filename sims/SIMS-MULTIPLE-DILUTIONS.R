# Install packages
## (Run once)
## install.packages("SLDAssay")
## devtools::install_github("sarahlotspeich/SLDeepAssay")

# Load packages
## Functions for MLE/bias-corrected MLE without UDSA
## from Trumble et al. (2017)
library(SLDAssay)
## Functions for MLE/bias-corrected MLE with UDSA
## from Lotspeich et al. (2022+)
library(SLDeepAssay)
## Functions for running simulations
library(tidyr)
library(pbapply)



# Number of replicates per simulation setting
num_reps = 2
# Define parameters that remain constant across all settings
M.ratio = c(1, 2, 3) # ratio of number of wells at each dilution level
q = c(0, 0.5, 1) # proportion of p24-positive wells deep-sequenced
u = c(.5, 1, 2) # dilution levels
Tau = 1
# Define parameters that vary over simulation settings
M.scale = c(6, 9, 12)   # M = M.scale * M.ratio
n = c(6, 12, 18)        # number of DVLs
constant_Tau = c(TRUE, FALSE) # indicator of constant IUPM across n DVLs
# Number of simulation settings
num_sett = length(M.scale) * length(n) * length(constant_Tau)
# Create dataframe of different simulation settings
Settings = apply(X = expand.grid("M.scale" = M.scale,
                                 "n" = n,
                                 "constant_Tau" = constant_Tau),
                 MARGIN = 2,
                 FUN = function(x) {
                   rep(x, each = num_reps)
                 }
) |>
  as.data.frame() |>
  dplyr::mutate(sim_id = rep(x = seq(1, num_reps), times = num_sett))
# Function one_sim() simulates one assay and returns the output (reshaped)
one_sim = function(setting_row) {
  # Save parameter values from row
  M.scale = as.numeric(setting_row["M.scale"])
  n = as.numeric(setting_row["n"])
  constant_Tau = setting_row["constant_Tau"] == 1
  if (constant_Tau) {
    tau = rep(x = Tau / n, times = n)
  } else {
    tau = c(rep(Tau / (2 * n), n / 2), rep(3 * Tau / (2 * n), n / 2))
  }
  res = simulate_SLDeepAssay_md(M = M.scale * M.ratio,
                                tau = tau,
                                q = q,
                                u = u)
  res_rshp = cbind(rep = r,
                   method = c("MLE_woUDSA", "BCMLE_woUDSA", "MLE_wUDSA", "BCMLE_wUDSA"),
                   do.call(what = rbind, args = res_r[grep(pattern = "MLE", x = names(res_r), value = FALSE)]),
                   Message = res_r$Message,
                   Message_Detailed = ifelse(is.null(res_r$Message_Detailed), "", res_r$Message_Detailed))
  return(res_rshp)
}

one_sim(Settings[1,])

# Function to conduct simulation study
## Simulates "reps" assays and fits all estimators to them
## Results are returned as a dataframe with rows for each estimator
simulation_study = function(M, n, lambda, q, remove_undetected = TRUE, dilution = 1, seed = 1, reps = 1000) {
  set.seed(seed)
  all_res = data.frame()
  for (r in 1:reps) {
    res_r = simulate_SLDeepAssay_sd(M = M, lambda = lambda, q = q, dilution = dilution)
    all_res = rbind(all_res,
                    cbind(rep = r,
                          method = c("MLE_woUDSA", "BCMLE_woUDSA", "MLE_wUDSA", "BCMLE_wUDSA"),
                          do.call(what = rbind, args = res_r[grep(pattern = "MLE", x = names(res_r), value = FALSE)]),
                          Message = res_r$Message,
                          Message_Detailed = ifelse(is.null(res_r$Message_Detailed), "", res_r$Message_Detailed)))
  }
  return(data.frame(all_res, row.names = NULL))
}

# Simulation set 1: Assuming uniformly distributed IUPM across DVL
path_to_sims = "~/" # Set path to folder where you would like to save results (one CSV per)
seed = 11422
for (M in c(12, 24, 32)) {
  for (n in c(6, 12, 18)) {
    for (Lambda in c(0.5, 1)) {
      lambda_vec = rep(Lambda / n, n)
      for (q in c(1, 0.9, 0.75, 0.5)) {
        seed = seed + 1
        simulation_study(M = M, n = n, lambda = lambda_vec, q = q, seed = seed) |>
          dplyr::mutate(M = M, n = n, Lambda = Lambda, q = q) |>
          write.csv(file = paste0(path_to_sims, "M", M, "_n", n, "_Lambda", 10 * Lambda, "_q", 100 * q, ".csv"),
                    row.names = FALSE)
      }
    }
  }
}

# Simulation set 2: Assuming non-uniformly distributed IUPM across DVL
path_to_sims = "~/" # Set path to folder where you would like to save results (one CSV per setting)
seed = 11422
for (M in c(12, 24, 32)) {
  for (n in c(6, 12, 18)) {
    for (Lambda in c(0.5, 1)) {
      lambda_vec = c(rep(Lambda / (2 * n), n / 2), rep(3 * Lambda / (2 * n), n / 2))
      for (q in c(1, 0.9, 0.75, 0.5)) {
        seed = seed + 1
        simulation_study(M = M, n = n, lambda = lambda_vec, q = q, seed = seed) |>
          dplyr::mutate(M = M, n = n, Lambda = Lambda, q = q) |>
          write.csv(file = paste0(path_to_sims, "M", M, "_n", n, "_Lambda", 10 * Lambda, "_q", 100 * q, ".csv"),
                    row.names = FALSE)
      }
    }
  }
}
