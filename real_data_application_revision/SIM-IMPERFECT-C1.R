library(SLDeepAssay)

# Define parameters based on Subject C1
M = c(36, 6, 6) # Total number of wells
MP = c(4, 0, 0) # Number of QVOA-positive wells
n = 12 # Number of detected DVLs
u = c(2.5, 0.5, 0.1) # Dilutions in millions of cells per well
q = c(1, 0, 0) # Proportion of p24-positive wells to undergo UDSA

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
  if (q[d] > 0) {
    ## Simulate the number of wells infected with each DVL 
    Z_rowsum = sample(x = 1:MP[d], 
                      size = n, 
                      replace = TRUE)
    
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
with(fit1, c(mle, convergence, message))

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
fit2$mle
