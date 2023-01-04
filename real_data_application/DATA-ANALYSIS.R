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


## store experimental setup data
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

colnames(results)[1:2] = c("id", "method")

## create plot data for Figure 2
figure2_data = rbind(results[, c("id", "method", "mle", "se", "ci1", "ci2")],
                 results[, c("id", "method", "mle_bc", "se", "ci_bc1", "ci_bc2")]) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("id", "method",
                 "mle", "se", "ci.lower", "ci.upper")) %>% 
  mutate_at(c("mle", "se", "ci.lower", "ci.upper"), as.numeric) %>% 
  mutate(bias.correction = factor(rep(c("MLE", "BC-MLE"),
                                      each = nrow(results))))

## data for Table S4
tableS4_data = as.data.frame(exp_setup)

# Create data for Table S5
tableS5_data = figure2_data %>%
  mutate(mle.ci = paste0("$", format(round(mle, 2), nsmall = 2), "$ $(",
                         format(round(ci.lower, 2), nsmall = 2), ", ",
                         format(round(ci.upper, 2), nsmall = 2), ")$"),
         id = factor(id, levels = paste0("C", 1:17))) %>% 
  select(id, method, bias.correction, mle.ci) %>% 
  pivot_wider(names_from = "method",
              values_from = "mle.ci") %>% 
  arrange(id)

# Save data
write.csv(figure2_data, "real_data_application/figure2_data.csv", row.names = F)
write.csv(tableS4_data, "real_data_application/tableS4_data.csv", row.names = F)
write.csv(tableS5_data, "real_data_application/tableS5_data.csv", row.names = F)

