## Title: Table 1: Single dilution setting simulations (constant Tau)

## Date: 2023/07/18

## Author: Sarah Lotspeich

## Purpose: produce a table to summarize simulation data from the single dilution setting

## Input: https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_sim_data.csv (produced by sims/SIMS-SINGLE-DILUTION.R)

## Output: Latex table

# load packages
library(kableExtra)
library(dplyr)
library(tidyr)

# load data
setwd(here::here())
sd_sim_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_sim_data.csv")

# Function format_nums() rounds and formats numbers for LaTex table
format_nums = function(x, digits = 2) {
  paste0("$", format(round(x, digits = digits), nsmall = digits) , "$")
}

# check number of sims removed
sd_sim_data |>
  select(Message) |>
  sum()

# check for replicates MLE and BC-MLE without UDSA were infinite (these are excluded)
## these counts are noted in the Table footnote
sd_sim_data |> 
  filter(abs(Est) == Inf) |> 
  group_by(Method, Tau, constant_Tau) |> 
  summarize(reps = n())

# summarize sim results
sd_sim_summ = sd_sim_data |>
  mutate(
    Message = ifelse(abs(Est) == Inf, 1, Message), 
    Est = ifelse(abs(Est) == Inf, NA, Est)) |> 
  group_by(constant_Tau, Tau, n, M, q, Method) |>
  summarise(n_removed = sum(Message),
            rel_bias = mean(Est - Tau, na.rm = TRUE),
            ase = mean(SE, na.rm = TRUE),
            ese = sd(Est, na.rm = TRUE),
            cp = mean(LB <= Tau & Tau <= UB, na.rm = TRUE)) |>
  pivot_wider(id_cols = c("M", "n", "q", "Tau", "constant_Tau"),
              names_from = "Method",
              values_from = c("n_removed", "rel_bias", "ase", "ese", "cp"))
  
analysis_cols = as.vector(outer(c("rel_bias", "ase", "ese", "cp"),
                as.vector(outer(c("_MLE", "_BCMLE"), c("_woUDSA", "_wUDSA"), paste0)), paste0))

col_order = c("n", "M", "q", analysis_cols)

# produce table with simulation summary
sd_sim_summ |> 
  filter(constant_Tau == 1, Tau == 1) |> # Subset to columns with constant Tau
  dplyr::ungroup() |>
  dplyr::select(all_of(col_order)) |>
  dplyr::mutate_at(analysis_cols, format_nums) |>
  magrittr::set_colnames(c("$\\pmb{n'}$", "M", "q", rep(c("Bias", "ASE", "ESE", "CP"), times = 4))) |>
  kable(format = "latex", digits = 3, align = c(rep("c", 3), rep("r", 16)),
        booktabs = TRUE, linesep = c("", "", "\\addlinespace"), escape = FALSE) |>
  add_header_above(header = c(" " = 3, "MLE" = 4, "Bias-Corrected MLE" = 4, "MLE" = 4, "Bias-Corrected MLE" = 4), bold = TRUE) |>
  add_header_above(header = c(" " = 3, "Without UDSA" = 8, "With UDSA" = 8), bold = TRUE) |>
  row_spec(row = 0, bold = TRUE)
