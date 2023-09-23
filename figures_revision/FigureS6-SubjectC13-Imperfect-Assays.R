## Title: Subject C13 with Imperfect Assays (Figure S6)

## Author:  Sarah Lotspeich

## Date:  08/28/2023

## Purpose: Plot results from a sensitivity analysis of imperfect assays in the real data example for SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
figureS6_data = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/subject_c13_imperfect_vary_sens.csv") 
## create Figure S6
plot_data = figureS6_data |> 
  dplyr::mutate(sensQVOA = factor(x = sens_qvoa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0(seq(80, 100, by = 10), "%")),
                                  #labels = paste0("QVOA Sensitivity = ", seq(80, 100, by = 10), "%")),
                sensUDSA = factor(x = sens_udsa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(80, 100, by = 10), "%")),
                specQVOA_UDSA = "UDSA Specificity =\nQVOA Specificity= 90%"
    ) 
naive_iupm = plot_data |> 
  dplyr::filter(spec_qvoa == 1) |> 
  dplyr::pull(mle)
plot_vary_sens = plot_data |> 
  dplyr::filter(spec_qvoa < 1) |> 
  ggplot(aes(x = sensQVOA, y = mle, col = sensQVOA, ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = naive_iupm, linetype = 2, color = "grey") + 
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.3,
                position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(sensUDSA)) + 
  theme_minimal(base_size = 16) + 
  theme(legend.position = "none") +
  xlab("QVOA Sensitivity") +
  ylab("IUPM Estimate") + 
  ggthemes::scale_colour_colorblind(name = "")

ggsave(plot_vary_sens, 
       filename = "figures_revision/real_data_plot_imperfect_subjectC13_vary_sens.png", 
       width = 8, height = 4)
