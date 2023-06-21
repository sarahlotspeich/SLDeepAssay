## Title: Table S1: IUPM overdispersion likelihood ratio test sim results

## Author:  Brian Richardson

## Date:  2023/06/20

## Purpose: Produce Table S1 in SLDeepAssay paper

## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)
library(emdbook)

## Load data
overdispersion_sim_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/overdispersion_sim_data.csv") %>% 
  mutate(Gamma = 1/k,
         gamma = ordered(Gamma))

## check for re-run sims
overdispersion_sim_data %>% 
  group_by(M.label, gamma) %>% 
  summarize(mean.error = mean(error),
            mean.redo.all = mean(num_redo_all, na.rm = T),
            mean.redo.none = mean(num_redo_none, na.rm = T)) %>% 
  kable(col.names = c("M", "gamma", "Prop. Error",
                      "Prop. Redo All", "Prop. Redo None")) %>% 
  kable_classic(full_width = F)

# 95% quantile of null distribution of LRT stat
cutoff <- emdbook::qchibarsq(q = 0.95, df = 1, mix = 0.5)

## create Table S1
overdispersion_sim_data %>% 
  group_by(M, M.label, gamma) %>% 
  summarise(power = mean(lrt_stat > cutoff, na.rm = T)) %>% 
  ungroup() %>% 
  select(M.label, gamma, power) %>% 
  kable(col.names = c("M", "gamma", "Power of LRT"),
        format = "latex",
        digits = 3,
        align = rep("c", 3),
        booktabs = TRUE,
        linesep = c("", "", "", "\\addlinespace"),
        escape = FALSE) %>% 
  row_spec(row = 0, bold = TRUE) -> tableS1

tableS1
