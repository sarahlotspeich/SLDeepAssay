## Title: Table S4: IUPM single dilution simulations for non-constant tau

## Author:  Sarah Lotspeich

## Date:  2023/06/20

## Purpose: Produce Table S4 in SLDeepAssay paper


# Load data
res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/single_dilution_results.csv")

# Run once
## install.packages("kableExtra)
library(kableExtra)

# Function format_nums() rounds and formats numbers for LaTex table
format_nums = function(x, digits = 2) {
  paste0("$", format(round(x, digits = digits), nsmall = digits) , "$")
}

# Create Table S4: Single-dilution assay with non-constant DVL-specific IUPMs & IUPM = 1
res |>
  dplyr::ungroup() |>
  dplyr::filter(!constant_Tau, Tau == 1) |>
  dplyr::select(-constant_Tau, -Tau) |>
  dplyr::mutate_at(dplyr::vars(-c("n", "M")), format_nums) |>
  magrittr::set_colnames(c("$\\pmb{M}$", "$\\pmb{n'}$", "$\\pmb{q}$", rep(c("Bias", "ASE", "ESE", "CP"), times = 4))) |>
  kable(format = "latex", align = c(rep("c", 3), rep("r", 16)), booktabs = TRUE, escape = FALSE) |>
  add_header_above(header = c(" " = 3, "MLE" = 4, "Bias-Corrected MLE" = 4, "MLE" = 4, "Bias-Corrected MLE" = 4), bold = TRUE) |>
  add_header_above(header = c(" " = 3, "Without UDSA" = 8, "With UDSA" = 8), bold = TRUE) |>
  row_spec(row = 0, bold = TRUE)
