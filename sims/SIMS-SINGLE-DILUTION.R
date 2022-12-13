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
