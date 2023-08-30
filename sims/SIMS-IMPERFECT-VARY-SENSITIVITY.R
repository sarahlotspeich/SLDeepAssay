library(SLDeepAssay)

# Number of replicates per simulation setting
num_reps = 1000
## Note: This code was run in parallel a cluster instead of locally, as it can be slow. 

# Define parameters that vary over simulation settings (same as Section 3.1)
spec = 0.9 # Sensitivity of assays
sens = seq(0.8, 1, by = 0.1) # Specificity of assays
M = c(12, 24, 32) # Total number of wells
n = 6 # Number of existing DVLs
Tau = 1 ## Overall IUPM (split among n DLVs)
u = 1 # Dilution in millions of cells per well
q = 1 # Proportion of p24-positive wells to undergo UDSA
constant_Tau = TRUE # Indicator of constant IUPM across n DVLs

# Number of simulation settings
num_sett = length(sens) ^ 2 * length(spec) ^ 2 * length(M) * length(n) * length(Tau) * length(constant_Tau) * length(q)

# Create dataframe of different simulation settings
Settings = apply(X = expand.grid("M" = M,
                                 "sensQVOA" = sens,
                                 "specQVOA" = spec,
                                 "sensUDSA" = sens,
                                 "specUDSA" = spec,
                                 "n" = n,
                                 "Tau" = Tau,
                                 "u" = u,
                                 "q" = q,
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
  if (setting_row["constant_Tau"] == 1) {
    tau = rep(x = as.numeric(setting_row["Tau"]) / as.numeric(setting_row["n"]), 
              times = as.numeric(setting_row["n"]))
  } else {
    tau = c(rep(as.numeric(setting_row["Tau"]) / (2 * as.numeric(setting_row["n"])), 
                as.numeric(setting_row["n"]) / 2), 
            rep(3 * as.numeric(setting_row["Tau"]) / (2 * as.numeric(setting_row["n"])), 
                as.numeric(setting_row["n"]) / 2))
  }
  temp = simulate_assay_sd(M = as.numeric(setting_row["M"]),
                           tau = as.numeric(setting_row["Tau"]),
                           q = as.numeric(setting_row["q"]),
                           u = as.numeric(setting_row["u"]),
                           sens_QVOA = as.numeric(setting_row["sensQVOA"]), 
                           spec_QVOA = as.numeric(setting_row["specQVOA"]), 
                           sens_UDSA = as.numeric(setting_row["sensUDSA"]), 
                           spec_UDSA = as.numeric(setting_row["specUDSA"])) 
  
  if (is.null(nrow(temp$DVL_specific))) {
    prevDVL = mean(temp$DVL_specific, na.rm = TRUE)
  } else {
    prevDVL = rowMeans(temp$DVL_specific, na.rm = TRUE)
  }
  
  assay_count = 1
  while(sum(temp$any_DVL) == 0 | is.null(dim(temp$DVL_specific)) | any(prevDVL == 1)) {
    assay_count = assay_count + 1
    temp = simulate_assay_sd(M = as.numeric(setting_row["M"]),
                             tau = as.numeric(setting_row["Tau"]),
                             q = as.numeric(setting_row["q"]),
                             u = as.numeric(setting_row["u"]),
                             sens_QVOA = as.numeric(setting_row["sensQVOA"]), 
                             spec_QVOA = as.numeric(setting_row["specQVOA"]), 
                             sens_UDSA = as.numeric(setting_row["sensUDSA"]), 
                             spec_UDSA = as.numeric(setting_row["specUDSA"])) 
    if (is.null(nrow(temp$DVL_specific))) {
      prevDVL = mean(temp$DVL_specific, na.rm = TRUE)
    } else {
      prevDVL = rowMeans(temp$DVL_specific, na.rm = TRUE)
    }
  }  
  setting_row["assay_resampled"] = assay_count > 1
  
  saveRDS(temp, "~/Downloads/temp-assay")
  
  ########################################################################################
  # Find MLEs ############################################################################
  ########################################################################################
  # New likelihood (corrected IUPM estimator)
  fit1 = fit_SLDeepAssay_sd_imperfect(assay_QVOA = temp$any_DVL, 
                                      assay_UDSA = temp$DVL_specific,
                                      sens_QVOA = as.numeric(setting_row["sensQVOA"]), 
                                      spec_QVOA = as.numeric(setting_row["specQVOA"]), 
                                      sens_UDSA = as.numeric(setting_row["sensUDSA"]), 
                                      spec_UDSA = as.numeric(setting_row["specUDSA"]))
  setting_row[c("Lambda", "conv", "msg")] = with(fit1, c(mle, convergence, message))
  
  # Original likelihood (naive IUPM estimator)
  fit2 = fit_SLDeepAssay_sd(assay = temp$DVL_specific,
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
  write.csv(Results, "sd_imperfect_vary_sensitivity.csv", row.names = F)
  if (i %% 10 == 0) print(paste0("Sim ", i, " complete (", round(100 * i/nrow(Settings)), "%)"))
}
