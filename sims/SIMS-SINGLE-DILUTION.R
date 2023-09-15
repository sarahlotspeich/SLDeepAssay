# Install packages (Run once)
## install.packages("SLDAssay")
## devtools::install_github("sarahlotspeich/SLDeepAssay")

# Load packages
## Functions for MLE/bias-corrected MLE without UDSA
## from Trumble et al. (2017)
library(SLDAssay)
## Functions for MLE/bias-corrected MLE with UDSA
## from Lotspeich et al. (2022+)
library(SLDeepAssay)

# Number of replicates per simulation setting
num_reps = 1000

# Define parameters that vary over simulation settings
M = c(12, 24, 32) # Total number of wells
n = c(6, 12, 18) # Number of existing DVLs
Tau = c(0.5, 1)  # Overall IUPM (split among n DLVs)
u = 1 # Dilution in millions of cells per well
q = c(0.5, 0.75, 1) # Proportion of p24-positive wells to undergo UDSA
constant_Tau = c(TRUE, FALSE) # Indicator of constant IUPM across n DVLs

# Number of simulation settings
num_sett = length(M) * length(n) * length(Tau) * length(constant_Tau) * length(q)

# Create dataframe of different simulation settings
Settings = apply(X = expand.grid("M" = M,
                                 "n" = n,
                                 "Tau" = Tau,
                                 "u" = u,
                                 "q" = q,
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
  # Save parameter values from setting_row
  M = as.numeric(setting_row["M"])
  n = as.numeric(setting_row["n"])
  Tau = as.numeric(setting_row["Tau"])
  u = as.numeric(setting_row["u"])
  q = as.numeric(setting_row["q"])
  constant_Tau = setting_row["constant_Tau"] == 1
  set.seed(as.numeric(setting_row["sim_id"]))
  if (constant_Tau) {
    tau = rep(x = Tau / n, times = n)
  } else {
    tau = c(rep(Tau / (2 * n), n / 2), rep(3 * Tau / (2 * n), n / 2))
  }
  res = simulate_SLDeepAssay_sd(M = M,
                                tau = tau,
                                q = q,
                                u = u)
  
  # Resimulate if extreme scenario
  resim_all_neg = resim_all_pos = 0 ## initialize counter for which extreme scenario was encountered
  if (mean(colSums(res$Assay), na.rm = TRUE) == 0 | any(rowMeans(res$Assay, na.rm = TRUE) == 1)) {
    res = simulate_SLDeepAssay_sd(M = M,
                                  tau = tau,
                                  q = q,
                                  u = u)
    resim_all_neg = resim_all_neg + as.numeric(mean(colSums(res$Assay), na.rm = TRUE))
    resim_all_pos = resim_all_pos + as.numeric(any(rowMeans(res$Assay, na.rm = TRUE) == 1))
  }
  
  res_rshp = cbind(constant_Tau, Tau, n, M, q,
                   Method = c("MLE_woUDSA", "BCMLE_woUDSA", "MLE_wUDSA", "BCMLE_wUDSA"),
                   do.call(what = rbind, args = res[grep(pattern = "MLE", x = names(res), value = FALSE)]),
                   Message = res$Message,
                   Message_Detailed = ifelse(is.null(res$Message_Detailed), "", res$Message_Detailed),
                   resim_all_neg, resim_all_pos,
                   row.names = NULL)
  return(res_rshp)
}

## Run once (if needed)
### install.packages("pbapply")
library(pbapply)
set.seed(11422)
Results = do.call(what = rbind,
                  args = pbapply(X = Settings,
                                 MARGIN = 1,
                                 FUN = one_sim
                                 )
)

# Save
Results |>
  write.csv(file = "sim_data/sd_sim_data.csv", row.names = FALSE)
