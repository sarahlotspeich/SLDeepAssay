library(SLDeepAssay)

# Number of replicates per simulation setting
num_reps = 1000

# Define parameters that vary over simulation settings (same as Section 3.1)
M = 12 # Total number of wells
n = 6 # Number of existing DVLs
Tau = 1 # Overall IUPM (split among n DLVs)
u = 1 # Dilutions in millions of cells per well
q = 0.5 # Proportion of p24-positive wells to undergo UDSA

# Number of simulation settings
num_sett = length(M) * length(n) 

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
  tau = rep(x = as.numeric(setting_row["Tau"]) / as.numeric(setting_row["n"]), 
            times = as.numeric(setting_row["n"]))
  temp = simulate_assay_sd(M = M, 
                           tau = tau, 
                           q = q, 
                           u = u, 
                           sens_QVOA = 1, 
                           spec_QVOA = 1, 
                           sens_UDSA = 1, 
                           spec_UDSA = 1)
  return(temp$DVL_specific)
}

# Be reproducible
sim_seed = 11422 
set.seed(sim_seed)

assays = list()
for (i in 1:nrow(Settings)) {
  assays[[i]] = one_sim(setting_row = Settings[i, ])
  if (i %% 10 == 0) print(paste0("Sim ", i, " complete (", round(100 * i/nrow(Settings)), "%)"))
}

length(unique(lapply(X = assays, FUN = colSums)))
