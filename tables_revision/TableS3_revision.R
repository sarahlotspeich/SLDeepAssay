## Title: Table S3: Multiple Dilutions Setting Simulations (Non-constant Tau)

## Date: 2023/06/20

## Author: Brian Richardson

## Purpose: Produce Table S3 in SLDeepAssay paper

## Input: https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/md_sim_data.csv (produced by sims/SIMS-MULTIPLE-DILUTIONS.R)

## Output: Latex table

# load packages
library(kableExtra)
library(dplyr)
library(tidyr)

# load data
setwd(here::here())
md_sim_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/md_sim_data.csv")

# format numbers (default to 2 decimal places)
format_nums = function(x, digits = 2) {
  paste0("$", format(round(x, digits = digits), nsmall = digits) , "$")
}

# check number of sims removed
md_sim_data |>
  select(Message) |>
  sum()

# summarize sim results
md_sim_summ = md_sim_data |>
  mutate(M = factor(M, levels = c("6, 12, 18", "9, 18, 27", "12, 24, 36"))) |>
  group_by(constant_Tau, M, n, assay_type, bc) |>
  summarise(n_removed = sum(Message),
            rel_bias = mean(Est - Tau),
            ase = mean(SE),
            ese = sd(Est),
            cp = mean(LB <= Tau & Tau <= UB)) |>
  pivot_wider(id_cols = c("M", "n", "constant_Tau"),
              names_from = c("bc", "assay_type"),
              values_from = c("n_removed", "rel_bias", "ase", "ese", "cp"))

analysis_cols = as.vector(outer(c("rel_bias", "ase", "ese", "cp"),
                                as.vector(outer(c("_MLE", "_BCMLE"), c("_woUDSA", "_wUDSA"), paste0)), paste0))

col_order = c("M", "n", analysis_cols)

# produce table with simulation summary
md_sim_summ |> 
  filter(constant_Tau == 0) |> # Subset to columns with constant Tau
  dplyr::ungroup() |>
  dplyr::select(all_of(col_order)) |>
  dplyr::mutate_at(analysis_cols, format_nums) |>
  magrittr::set_colnames(c("M", "$\\pmb{n'}$", rep(c("Bias", "ASE", "ESE", "CP"), times = 4))) |>
  kable(format = "latex", digits = 3, align = c(rep("c", 4), rep("r", 16)),
        booktabs = TRUE, linesep = c("", "", "\\addlinespace"), escape = FALSE) |>
  add_header_above(header = c(" " = 2, "MLE" = 4, "Bias-Corrected MLE" = 4, "MLE" = 4, "Bias-Corrected MLE" = 4), bold = TRUE) |>
  add_header_above(header = c(" " = 2, "Without UDSA" = 8, "With UDSA" = 8), bold = TRUE) |>
  row_spec(row = 0, bold = TRUE)
