## Title: IUPM Real Data Table S4

## Author:  Brian Richardson

## Date:  01/04/2023

## Purpose: Produce Table S4 SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)

# Load data
tabS4_dat = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application/tableS4_data.csv")

## create Table S4 with experimental setup
tabS4_dat %>% 
  kable(format = "latex",
        booktabs = T,
        escape = F,
        align = "ccccc",
        col.names = linebreak(c("Subject\nID",
                                "DLVs\n(n)",
                                "Wells\n(M)",
                                "p24-Positive\nWells (MP)",
                                "Deep-Seqeunced\nWells (m)"),
                              align = "c")) %>%
  kable_classic_2(full_width = F) %>%
  kable_styling(font_size = 10) %>% 
  row_spec(row = 0, bold = TRUE) -> tableS4
