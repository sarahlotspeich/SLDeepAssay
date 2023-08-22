## Title: Table S5: IUPM real data results

## Author:  Brian Richardson

## Date:  2023/06/20

## Purpose: Produce Table S5 in SLDeepAssay paper

## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(magick)

# Load data
tableS8_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application_revision/tableS5_data.csv")

## create Table S8
tableS8_data |> 
  kable(format = "latex",
        booktabs = T,
        escape = F,
        align = "ccccc",
        col.names = linebreak(c("Subject\nID",
                                "Bias\nCorrection",
                                "Without UDSA\n(Multiple Dilutions)",
                                "With UDSA\n(Single Dilution)",
                                "With UDSA\n(Multiple Dilutions)"),
                              align = "c")) |>
  kable_classic_2(full_width = F) |>
  kable_styling(font_size = 10) |> 
  add_header_above(c(" " = 2, "MLE (95% Confidence Interval)" = 3), bold = T) |> 
  row_spec(row = 0, bold = TRUE) |> 
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") -> tableS8

tableS8
