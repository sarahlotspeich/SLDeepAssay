## Title: Table S3: IUPM real data overdispersion results for subject C13

## Author:  Brian Richardson

## Date:  2023/06/20

## Purpose: Produce Table S3 in SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)

# Load data
tableS3_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application_revision/tableS2_data.csv")

## create Table S3
tableS3_data$MP <- paste0(tableS3_data$MP1,
                          tableS3_data$MP2,
                          tableS3_data$MP3,
                          tableS3_data$MP4,
                          sep = ", ")


tableS3_data %>% 
  select(MP, lrt_stat, lrt_pval, mle_negbin) %>% 
  kable(format = "latex",
        booktabs = T,
        escape = F,
        align = "ccccc",
        digits = rep(3, 4),
        col.names = c("MP",
                      "LRT Statistic",
                      "LRT P-Value",
                      "Negative Binomial MLE")) %>%
  kable_classic_2(full_width = F) %>%
  kable_styling(font_size = 10) %>% 
  row_spec(row = 0, bold = TRUE) -> tableS3

tableS3
