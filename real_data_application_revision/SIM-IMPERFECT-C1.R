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

# Estimate IUPM 
## New likelihood (corrected IUPM estimator for imperfect assays)
fit1 = fit_SLDeepAssay_md_imperfect(assay_md = temp,
                                    u = u, 
                                    sens_QVOA = 1, 
                                    spec_QVOA = 1, 
                                    sens_UDSA = 1, 
                                    spec_UDSA = 1)
#with(fit1, c(mle, convergence, message))

# Original likelihood (naive IUPM estimator)
assay_summary = vapply(X = 1:length(u),
                       FUN.VALUE = numeric(7 + n),
                       FUN = function(d) {
                         M = ncol(temp[[d]]$DVL_specific) # number of wells
                         n = nrow(temp[[d]]$DVL_specific) # number of DVLs detected
                         MP = sum(temp[[d]]$any_DVL) # number of p24-positive wells
                         MN = M - MP # number of p24-negative wells
                         m = sum(!is.na(colSums(temp[[d]]$DVL_specific))) # number of deep-sequenced wells
                         q = ifelse(MP == 0, 0, m / MP) # proportion of p24-positive wells deep sequenced
                         Y = rowSums(temp[[d]]$DVL_specific, na.rm = TRUE) # number of infected wells per DVL
                         return((c("u" = u[d], "M"=M, "n"=n,
                                   "MN"=MN, "MP"=MP, "m"=m, "q"=q, "Y"=Y)))
                       })
assay_summary = as.data.frame(t(assay_summary))
fit2 = fit_SLDeepAssay_md(assay_summary = assay_summary, 
                          corrected = FALSE)

# Check: fit1 with perfect sens/spec should be approximately equal to fit2
fit1$mle; fit2$mle

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
                                      spec_UDSA = try_sens_spec$spec_udsa[t])
  
  try_sens_spec$mle[t] = fit1$mle
  try_sens_spec$ci_lb[t] = fit1$ci[1]
  try_sens_spec$ci_ub[t] = fit1$ci[2]
  try_sens_spec$code[t] = fit1$convergence
  try_sens_spec$msg[t] = fit1$message
}

# Save results 
try_sens_spec |> 
  write.csv("~/Documents/SLDeepAssay/sim_data/subject_c1_imperfect.csv", 
            row.names = F)
