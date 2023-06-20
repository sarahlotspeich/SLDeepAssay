## Title: Table S7: IUPM real data assay summaries

## Author:  Brian Richardson

## Date:  2023/06/20

## Purpose: Produce Table S7 in SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)

# Load data
tabS7_dat = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application_revision/tableS7_data.csv")

## create Table S7 with experimental setup
tabS7_dat %>% 
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
  row_spec(row = 0, bold = TRUE) -> tableS7

tableS7
