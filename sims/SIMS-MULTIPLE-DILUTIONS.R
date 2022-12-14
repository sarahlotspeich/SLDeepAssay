## Title: Simulations for Multiple Dilutions Setting

## Date: 2022/12/14

## Author: Brian Richardson

## Purpose: produce simulation data to assess performance of SLDeepAssay estimator in the multiple dilutions setting

## Output: md_sim_data.csv

## Note: These simulations take ~12 hours to run. The simulation results can be found in sims/md_sim_data.csv.



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
library(dplyr)
library(pbapply)


# Random number seed for reproducibility
set.seed(130502)

# Number of replicates per simulation setting
num_reps = 1000
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
  res_rshp = cbind(method = c("MLE_woUDSA", "BCMLE_woUDSA", "MLE_wUDSA", "BCMLE_wUDSA"),
                   do.call(what = rbind, args = res[grep(pattern = "MLE", x = names(res), value = FALSE)]),
                   Message = res$Message,
                   Message_Detailed = ifelse(is.null(res$Message_Detailed), "", res$Message_Detailed))
  return(res_rshp)
}

# run simulations
sim_out = do.call("rbind", pbapply(Settings, 1,
                                   function(row) one_sim(setting_row = row)))

results = cbind(Settings, sim_out) |>
  mutate(M = paste0(M.ratio[1] * M.scale, ", ",
                   M.ratio[2] * M.scale, ", ",
                   M.ratio[3] * M.scale),
         Tau = Tau) |>
  separate(col = method, into = c("bc", "assay_type"), sep = "_")


# save simulations
path_to_sims = here::here() # Set path to folder where you would like to save results

write.csv(results, file = paste0(path_to_sims, "/sims/md_sim_data.csv"),
          row.names = FALSE)
