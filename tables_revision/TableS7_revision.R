## Title: Table S2: IUPM real data overdispersion results

## Author:  Brian Richardson

## Date:  2023/06/20

## Purpose: Produce Table S2 in SLDeepAssay paper

## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)

# Load data
tableS2_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application_revision/tableS7_data.csv")

## create Table S2
tableS2_data %>% 
  select(-mle_bc) %>% 
  kable(format = "latex",
        booktabs = T,
        escape = F,
        align = "ccccc",
        digits = rep(3, 5),
        col.names = c("Subject ID",
                      "LRT Statistic",
                      "LRT P-Value",
                      "Poisson MLE",
                      "Negative Binomial MLE")) %>%
  kable_classic_2(full_width = F) %>%
  kable_styling(font_size = 10) %>% 
  add_header_above(c(" " = 3, "Estimated IUPM" = 2), bold = T) %>% 
  row_spec(row = 0, bold = TRUE) -> tableS2

tableS2
