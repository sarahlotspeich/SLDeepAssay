## Title: Simulations for Testing Overdispersion

## Date: 2023/06/20

## Author: Brian Richardson

## Purpose: produce simulation data to assess negative binomial MLE and overdispersion LRT

## Output: overdispersion_sim_data.csv

## Note: These simulations take ~12 hours to run. The simulation results can be found in sims/sim_data/md_sim_data.csv.



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
#library(tidyr)
library(dplyr)
library(pbapply)
library(here)

# Number of replicates per simulation setting
num_reps = 500
# Define parameters that remain constant across all settings
u. <- c(0.5, 1, 2)               # dilution levels
D. <- length(u.)                 # number of dilution levels
q. <- c(0, 0.5, 1)               # proportion deep-sequenced
n. <- 12                         # number of DVLs
Tau. <- 1                        # IUPM
tau. <- rep(1/n., n.)            # DVL-specific IUPMs
# Define parameters that vary over simulation settings
M. <- matrix(nrow = 2, byrow = TRUE,         # number of replicate wells
             data = c(6, 12, 18,             # small sample size
                      30, 60, 90))           # large sample

gamma. <- c(0, 1/5, 1, 4)                    # dispersion parameters
# Number of simulation settings
num_sett = length(M.) * length(gamma.)
# Create dataframe of different simulation settings
Settings = expand.grid("gamma" = gamma.,
                        "M" = 1:2,
                        "seed" = 1:num_reps) |> 
  as.data.frame() |> 
  mutate(M.label = factor(paste0("(", ifelse(M == 1,
                                             paste0(unique(M.[1,]), collapse = ", "),
                                             paste0(unique(M.[2,]), collapse = ", ")),
                                 ")"),
                          levels = paste0("(",
                                          c(paste0(unique(M.[1,]), collapse = ", "),
                                            paste0(unique(M.[2,]), collapse = ", ")),
                                          ")")))
# Function lrt_one_sim() simulates one assay and returns the output (reshaped)
lrt_one_sim = function(M, tau, q, u, gamma, seed = NULL, remove_undetected = TRUE) {
  
  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create indicators of whether data need to be re-simulated
  num_redo_all = 0 # Due to any DVL being detected in all wells or
  num_redo_none = 0 # No DVL being detected in any wells
  
  # Simulate single-dilution assay data
  assay_summary = simulate_assay_md(
    M = M,
    tau = tau,
    q = q,
    u = u,
    k = 1 / gamma,
    remove_undetected = remove_undetected)
  
  # Check for need to re-simulate
  max_p24_neg = tryCatch(expr = max(assay_summary$M - assay_summary$MP), 
                         error = function(e) 0)
  
  min_DVL_neg = tryCatch(expr = min(colSums(assay_summary$m - assay_summary[, grepl("Y", colnames(assay_summary))])), 
                         error = function(e) 0) 
  
  while ((max_p24_neg == 0 & min_DVL_neg == 0) | is.null(assay_summary)) { 
    if (is.null(assay)) {
      num_redo_none = num_redo_none + 1
    } else {
      num_redo_all = num_redo_all + 1
    }
    # Re-simulate single-dilution assay data
    # Simulate single-dilution assay data
    assay_summary = simulate_assay_md(
      M = M,
      tau = tau,
      q = q,
      u = u,
      gamma = gamma,
      remove_undetected = remove_undetected)
    
    max_p24_neg = tryCatch(expr = max(assay_summary$M - assay_summary$MP), 
                           error = function(e) 0)
    
    min_DVL_neg = tryCatch(expr = min(colSums(assay_summary$m - assay_summary[, grepl("Y", colnames(assay_summary))])), 
                           error = function(e) 0) 
  }
  
  # analyze assay
  res <- lrt_SLDeepAssay_md(assay_summary = assay_summary)
  
  # Construct list to return
  return(c("num_redo_all" = num_redo_all,
           "num_redo_none" = num_redo_none,
           "mle_pois" = res$mle,
           "mle_pois_bc" = res$mle_bc,
           "mle_negbin" = res$mle_negbin,
           "mle_gamma" = res$mle_gamma,
           "lrt_stat" = res$lrt_stat))
}

# run simulations
sim.out <- pbvapply(
  X = 1:nrow(Settings),
  FUN = function(i) {
    tryCatch(expr = c("error" = 0,
                      lrt_one_sim(tau = tau., q = q., u = u.,
                                  M = M.[Settings$M[i], ],
                                  gamma = Settings$gamma[i],
                                  seed = Settings$sim_id[i])),
             error = function(e) c(1, rep(NA, 7)))
  },
  FUN.VALUE = numeric(8)
) |> 
  t()

results <- cbind(Settings, sim.out)

# save simulations
path_to_sims = here::here() # Set path to folder to save results

write.csv(results, file = paste0(path_to_sims, "/sim_data/overdispersion_sim_data.csv"),
          row.names = FALSE)
