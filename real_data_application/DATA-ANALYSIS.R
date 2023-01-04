## Title: IUPM Real Data Analysis

## Author:  Brian Richardson

## Date:  01/04/2023

## Purpose: Analyze real data for SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(SLDeepAssay)
library(SLDAssay)
library(kableExtra)
library(magick)

# load one data set for each of the 17 subjects
for (id in paste0("C", 1:17)) {
  
  file.name <- paste0("real_data_application/data/UDSA_MD_", id, ".csv")
  dat.name <- paste0("dat_md_", id)
  
  assign(dat.name, read.csv(file.name, row.names = 1))
  
}


## tore experimental setup data
exp_setup = matrix(nrow = 0, ncol = 5)
colnames(exp_setup) = c("id", "n", "M", "MP", "m")

# matrix to store analysis results
results = matrix(nrow = 0, ncol = 9)

for (id in paste0("C", 1:17)) {
  
  # get subject data
  assay_summary = eval(parse(text = paste0("dat_md_", id)))
  
  ### experimental setup
  exp_setup = rbind(exp_setup,
                    cbind(id,
                          assay_summary[1, "n"],
                          paste(as.character(assay_summary[, "M"]), sep="' '", collapse=", "),
                          paste(as.character(assay_summary[, "MP"]), sep="' '", collapse=", "),
                          paste(as.character(assay_summary[, "m"]), sep="' '", collapse=", ")))
  
  ### SLDAssay results
  # my.log = capture.output({ ... avoids the line printed by get.mle()
  my_log = capture.output({
    res_SLDAssay = get.mle(pos = as.numeric(assay_summary$MP),
                           replicates = as.numeric(assay_summary$M),
                           dilutions = as.numeric(assay_summary$u * 10^6),
                           monte=10)
  })
  
  # compute standard error
  SLDAssay_se = (log(res_SLDAssay$Asymp_CI[2]) - log(res_SLDAssay$MLE)) *
    res_SLDAssay$MLE / qnorm(0.975)
  
  SLDAssay_bc_ci = exp(c(log(res_SLDAssay$BC_MLE) + c(-1, 1) * (qnorm(0.975) *
                                                                   SLDAssay_se / res_SLDAssay$BC_MLE)))
  
  results = rbind(results,
                  c(id, "Without UDSA (Multiple Dilutions)",
                    res_SLDAssay$MLE, SLDAssay_se, res_SLDAssay$Asymp_CI,
                    res_SLDAssay$BC_MLE, SLDAssay_bc_ci))
  
  ### deepIUPM results
  # dilution index where deep sequencing is done
  d_deepseq = which(assay_summary$deepseq == 1) 
  # dilution level where deep sequencing is done
  u_deepseq = assay_summary$u[d_deepseq]
  
  res_SLDeepAssay_sd = fit_SLDeepAssay_sd(assay = NULL,
                        M = assay_summary$M[d_deepseq],
                        MP = assay_summary$MP[d_deepseq],
                        m = assay_summary$m[d_deepseq],
                        Y = as.numeric(assay_summary[d_deepseq, grepl("Y", colnames(assay_summary))]),
                        u = assay_summary$u[d_deepseq],
                        corrected = T)
  
  results = rbind(results,
                  c(id, "With UDSA (Single Dilution)", unlist(res_SLDeepAssay_sd)))
  
  ### deepIUPM MD results
  res_SLDeepAssay_md = fit_SLDeepAssay_md(assay = NULL, assay_summary = assay_summary, corrected = T)
  
  results = rbind(results,
                  c(id, "With UDSA (Multiple Dilutions)", unlist(res_SLDeepAssay_md)))
  
}

exp_setup = as.data.frame(exp_setup)
colnames(results)[1:2] = c("id", "method")

# Save results
write.csv(exp_setup, "real_data_application/tableS4_data.csv", row.names = F)
write.csv(results, "real_data_application/tableS5_data.csv", row.names = F)

