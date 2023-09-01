## Title: IUPM Real Data Analysis (Revision)

## Author:  Brian Richardson

## Date:  2023/08/30

## Purpose: Analyze real data for SLDeepAssay paper with Biometrics revisions

# prepare workspace -------------------------------------------------------

rm(list=ls())
setwd(here::here())

# Install packages (Run once)
## install.packages("SLDAssay")
## devtools::install_github("sarahlotspeich/SLDeepAssay")

## load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(SLDeepAssay)
library(SLDAssay)
library(kableExtra)
library(magick)


# load data ---------------------------------------------------------------

for (id in paste0("C", 1:17)) {
  
  file.name <- paste0("real_data_application/data/UDSA_MD_", id, ".csv")
  dat.name <- paste0("dat_md_", id)
  
  assign(dat.name, read.csv(file.name, row.names = 1))
  
}

## store experimental setup data
exp_setup = matrix(nrow = 0, ncol = 5)
colnames(exp_setup) = c("id", "n", "M", "MP", "m")


# analyze data ------------------------------------------------------------

# 95% quantile of null distribution of LRT stat
cutoff <- emdbook::qchibarsq(q = 0.95, df = 1, mix = 0.5)

# matrix to store analysis results
results = matrix(nrow = 0, ncol = 13)

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
                    res_SLDAssay$BC_MLE, SLDAssay_bc_ci,
                    NA, NA, NA, NA))
  
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
                  c(id, "With UDSA (Single Dilution)", unlist(res_SLDeepAssay_sd),
                    NA, NA, NA, NA))
  
  ### SLDeepAssay MD results
  res_SLDeepAssay_md = fit_SLDeepAssay_md(assay = NULL, assay_summary = assay_summary, corrected = T)
  res_lrt_SLDeepAssay_md = lrt_SLDeepAssay_md(assay_summary = assay_summary, corrected = T)
  lrt_sig = ifelse(res_lrt_SLDeepAssay_md$lrt_stat <= cutoff, 0, 1)
  
  results = rbind(results,
                  c(id, "With UDSA (Multiple Dilutions)",
                    as.character(unlist(res_SLDeepAssay_md)),
                    as.character(unlist(res_lrt_SLDeepAssay_md[c("mle_negbin", "mle_gamma", "lrt_stat")])),
                    lrt_sig))
  
}

colnames(results)[c(1:2, 10:13)] = c("id", "method", "mle_negbin", 
                                     "mle_gamma", "lrt_stat", "lrt_sig")

## create plot data for Figure 2
figure2_data = rbind(results[, c("id", "method", "mle", "se", "ci1", "ci2")],
                     results[, c("id", "method", "mle_bc", "se", "ci_bc1", "ci_bc2")]) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("id", "method",
                 "mle", "se", "ci.lower", "ci.upper")) %>% 
  mutate_at(c("mle", "se", "ci.lower", "ci.upper"), as.numeric) %>% 
  mutate(bias.correction = factor(rep(c("MLE", "BC-MLE"),
                                      each = nrow(results))))


## data for Table S4 (real data assay summaries)
TableS4_realDataSummary_data = as.data.frame(exp_setup)

## Data for Table S5: real data main analysis
TableS5_realDataMainAnalysis_data = figure2_data %>%
  mutate(mle.ci = paste0("$", format(round(mle, 2), nsmall = 2), "$ $(",
                         format(round(ci.lower, 2), nsmall = 2), ", ",
                         format(round(ci.upper, 2), nsmall = 2), ")$"),
         id = factor(id, levels = paste0("C", 1:17))) %>% 
  select(id, method, bias.correction, mle.ci) %>% 
  pivot_wider(names_from = "method",
              values_from = "mle.ci") %>% 
  arrange(id)

# Data for Table S7: real data overdispersion results
TableS7_realDataLRT_data = results %>% 
  as.data.frame() %>% 
  filter(method == "With UDSA (Multiple Dilutions)") %>% 
  mutate(lrt_pval = 0.5 * pchisq(as.numeric(lrt_stat), df = 1, lower.tail = F) +
                    0.5 * (lrt_stat == 0)) %>% 
  select(id, lrt_stat, lrt_pval, mle, mle_bc, mle_negbin)


# LRT sensitivity analysis for C13 ----------------------------------------

# Table S8 data, subject C13 sensitivity analysis
TableS8_realDataLRTC13_data <- data.frame(MP1 = c(16, 17, 18, 16, 16),
                           MP2 = c(4,  4,  4,  4,  4),
                           MP3 = c(3,  3,  3,  2,  1),
                           MP4 = c(0,  0,  0,  0,  0),
                           lrt_stat = rep(NA, 5),
                           lrt_pval = rep(NA, 5),
                           mle_negbin = rep(NA, 5))

# change MP and record new LRT stat and negative binomial MLE
for (i in 1:nrow(TableS8_realDataLRTC13_data)) {
  
  dat_sens_C13 <- dat_md_C13 %>% 
    mutate(MP = as.numeric(TableS8_realDataLRTC13_data[i, 1:4]))
  
  sens_res <- lrt_SLDeepAssay_md(assay_summary = dat_sens_C13, corrected = T)
  
  TableS8_realDataLRTC13_data[i, 5:7] <- c(sens_res$lrt_stat,
                            0.5 * pchisq(as.numeric(sens_res$lrt_stat), df = 1, lower.tail = F),
                            sens_res$mle_negbin)
  
}


# Save data
write.csv(figure2_data,
          "real_data_application_revision/figure2_data.csv",
          row.names = F)

write.csv(TableS4_realDataSummary_data,
          "real_data_application_revision/TableS4_realDataSummary_data.csv",
          row.names = F)

write.csv(TableS5_realDataMainAnalysis_data,
          "real_data_application_revision/TableS5_realDataMainAnalysis_data.csv",
          row.names = F)

write.csv(TableS7_realDataLRT_data,
          "real_data_application_revision/TableS7_realDataLRT_data.csv",
          row.names = F)

write.csv(TableS8_realDataLRTC13_data,
          "real_data_application_revision/TableS8_realDataLRTC13_data.csv", 
        row.names = F)

# Compare CI widths (on log scale) for selected subjects
ci.widths = figure2_data %>% 
  mutate(ci.width = log(ci.upper) - log(ci.lower)) %>% 
  select(id, method, bias.correction, ci.width)

ci.widths.udsa.md = ci.widths %>% 
  filter(method == "With UDSA (Multiple Dilutions)") %>% 
  rename(ci.width.new = ci.width) %>% 
  select(!method)

ci.change = ci.widths %>% 
  filter(method != "With UDSA (Multiple Dilutions)") %>% 
  left_join(ci.widths.udsa.md,
            by = c("id", "bias.correction")) %>% 
  mutate(pct.change = 100 * (ci.width.new - ci.width) / ci.width)

# compare UDSA + QVOA at multiple dilutions to QVOA at multiple dilutions
ci.change %>% 
  filter(method == "Without UDSA (Multiple Dilutions)") %>% 
  arrange(by = pct.change)

# compare UDSA + QVOA at multiple dilutions to UDSA + QVOA at one dilution
ci.change %>% 
  filter(method == "With UDSA (Single Dilution)") %>% 
  arrange(by = pct.change)




