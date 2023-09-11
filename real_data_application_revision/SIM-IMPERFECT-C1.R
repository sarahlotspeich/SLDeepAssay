library(SLDeepAssay)

# Define parameters based on Subject C1
subject_params = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application_revision/data/UDSA_MD_C1.csv")

M = subject_params$M # Total number of wells
MP = subject_params$MP # Number of QVOA-positive wells
n = subject_params$n[subject_params$deepseq == 1] # Number of detected DVLs
u = subject_params$u # Dilutions in millions of cells per well
m = subject_params$m # Number of p24-positive wells to undergo UDSA
Y = subject_params[subject_params$deepseq == 1, -c(1:9)]

## Reproducibility
set.seed(123)

# Simulate more detailed well-level, per-DVL data 
temp = list() ## create empty list to hold multiple dilution assays
for (d in 1:length(u)) {
  ## Create vector of any DVL postive/negative results
  W = c(rep(1, MP[d]), 
        rep(0, (M[d] - MP[d])))
  
  ## Create matrix to hold DVL-specific positive/negative results
  Z = matrix(data = NA, 
             nrow = n, 
             ncol = M[d])
  
  ## If deep sequencing was done, replace NA values for m deep-sequenced wells
  if (m[d] > 0) {
    ## Take the number of wells infected with each DVL from the real data
    Z_rowsum = Y
    
    ## For each DVL... 
    for(i in 1:nrow(Z)) {
      ### Create a vector of positive/negative statuses with sum = Z_rowsum
      Zi = c(rep(x = 1, times = Z_rowsum[i]), 
             rep(x = 0, times = MP[d] - Z_rowsum[i]))
      
      ### And order them randomly before putting into Z matrix
      Z[i, 1:MP[d]] = Zi[sample(x = 1:MP[d], size = MP[d], replace = FALSE)]
    }
  } 
  
  ## Save these data to multiple dilutions assay
  temp[[d]] = list(any_DVL = W, 
                 DVL_specific = Z)
}

# Try a variety of sensitivity/specificity combinations
## (same sensitivity/specificity combinations as in Figure S5)
try_sens_spec = expand.grid(sens_qvoa = seq(0.8, 1, by = 0.1), 
                            spec_qvoa = 0.9, 
                            sens_udsa = seq(0.8, 1, by = 0.1), 
                            spec_udsa = 0.9,
                            mle = NA, 
                            ci_lb = NA, 
                            ci_ub = NA,
                            code = NA,
                            msg = NA
                            )

for (t in 1:nrow(try_sens_spec)) {
  # Estimate IUPM 
  ## New likelihood (corrected IUPM estimator for imperfect assays)
  fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                      u = u, 
                                      sens_QVOA = try_sens_spec$sens_qvoa[t], 
                                      spec_QVOA = try_sens_spec$spec_qvoa[t], 
                                      sens_UDSA = try_sens_spec$sens_udsa[t], 
                                      spec_UDSA = try_sens_spec$spec_udsa[t],
                                      lb = 1E-6)
  
  try_sens_spec$mle[t] = fit1$mle
  try_sens_spec$ci_lb[t] = fit1$ci[1]
  try_sens_spec$ci_ub[t] = fit1$ci[2]
  try_sens_spec$code[t] = fit1$convergence
  try_sens_spec$msg[t] = fit1$message
}

# For comparison - 
## Estimate IUPM using original likelihood ("perfect assays" MLE)
fit2 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                    u = u, 
                                    sens_QVOA = 1, 
                                    spec_QVOA = 1, 
                                    sens_UDSA = 1, 
                                    spec_UDSA = 1,
                                    lb = 1E-6)

try_sens_spec = try_sens_spec |> 
  dplyr::bind_rows(
    data.frame(sens_qvoa = 1, spec_qvoa = 1, sens_udsa = 1, spec_udsa = 1, 
               mle = fit2$mle, ci_lb = fit2$ci[1], ci_ub = fit2$ci[2],
               code = fit2$convergence, msg = fit2$message)
  )

# Save results 
try_sens_spec |> 
  write.csv("~/Documents/SLDeepAssay/sim_data/subject_c1_imperfect_vary_sens.csv", 
            row.names = F)

# Try a variety of sensitivity/specificity combinations
## (same sensitivity/specificity combinations as in Figure S5)
try_sens_spec = expand.grid(sens_qvoa = 0.9, 
                            spec_qvoa = seq(0.8, 1, by = 0.1), 
                            sens_udsa = 0.9, 
                            spec_udsa = seq(0.8, 1, by = 0.1),
                            mle = NA, 
                            ci_lb = NA, 
                            ci_ub = NA,
                            code = NA,
                            msg = NA
)

for (t in 1:nrow(try_sens_spec)) {
  # Estimate IUPM 
  ## New likelihood (corrected IUPM estimator for imperfect assays)
  fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                      u = u, 
                                      sens_QVOA = try_sens_spec$sens_qvoa[t], 
                                      spec_QVOA = try_sens_spec$spec_qvoa[t], 
                                      sens_UDSA = try_sens_spec$sens_udsa[t], 
                                      spec_UDSA = try_sens_spec$spec_udsa[t],
                                      lb = 1E-6)
  
  try_sens_spec$mle[t] = fit1$mle
  try_sens_spec$ci_lb[t] = fit1$ci[1]
  try_sens_spec$ci_ub[t] = fit1$ci[2]
  try_sens_spec$code[t] = fit1$convergence
  try_sens_spec$msg[t] = fit1$message
}

# For comparison - Estimate IUPM using original likelihood ("perfect assays" MLE)
try_sens_spec = try_sens_spec |> 
  dplyr::bind_rows(
    data.frame(sens_qvoa = 1, spec_qvoa = 1, sens_udsa = 1, spec_udsa = 1, 
               mle = fit2$mle, ci_lb = fit2$ci[1], ci_ub = fit2$ci[2],
               code = fit2$convergence, msg = fit2$message)
  )

# Save results 
try_sens_spec |> 
  write.csv("~/Documents/SLDeepAssay/sim_data/subject_c1_imperfect_vary_spec.csv", 
            row.names = F)
