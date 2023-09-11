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
figureS6A_data = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/subject_c13_imperfect_vary_sens.csv") 
## create Figure S6
plot_data = figureS6A_data |> 
  dplyr::mutate(sensQVOA = factor(x = sens_qvoa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("QVOA Sensitivity = ", seq(80, 100, by = 10), "%")),
                sensUDSA = factor(x = sens_udsa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(80, 100, by = 10), "%")),
                specQVOA_UDSA = "UDSA Specificity = QVOA Specificity = 90%"
    ) 
naive_iupm = plot_data |> 
  dplyr::filter(spec_qvoa == 1) |> dplyr::pull(mle)
plot_vary_sens = plot_data |> 
  dplyr::filter(spec_qvoa < 1) |> 
  ggplot(aes(x = 1, y = mle, col = sensQVOA, ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = naive_iupm, linetype = 2, color = "grey") + 
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.3,
                position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(sensUDSA), 
             rows = vars(specQVOA_UDSA)) + 
  theme_light(base_size = 24) +
  theme(legend.position = "top", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = "black")) +
  ylab("IUPM Estimate") + 
  ggthemes::scale_colour_colorblind(name = "")

ggsave(plot_vary_sens, 
       filename = "figures_revision/real_data_plot_imperfect_subjectC13_vary_sens.png", 
       width = 8, height = 4)

# Load data
figureS6B_data = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/subject_c13_imperfect_vary_spec.csv") 
## create Figure S6
plot_data = figureS6B_data |> 
  dplyr::mutate(specQVOA = factor(x = spec_qvoa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("QVOA Specificity = ", seq(80, 100, by = 10), "%")),
                specUDSA = factor(x = spec_udsa, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Specificity = ", seq(80, 100, by = 10), "%")),
                sensQVOA_UDSA = "UDSA Sensitivity = QVOA Sensitivity = 90%"
  ) 
plot_vary_spec = plot_data |> 
  dplyr::filter(sens_qvoa < 1) |> 
  ggplot(aes(x = 1, y = mle, col = specQVOA, ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = naive_iupm, linetype = 2, color = "grey") + 
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.3,
                position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(specUDSA), 
             rows = vars(sensQVOA_UDSA)) + 
  theme_light(base_size = 24) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = "black")) +
  ylab("IUPM Estimate") + 
  ggthemes::scale_colour_colorblind(name = "")

ggsave(plot_vary_spec, 
       filename = "figures_revision/real_data_plot_imperfect_subjectC13_vary_spec.png", 
       width = 8, height = 4)

ggsave(ggpubr::ggarrange(plot_vary_sens, plot_vary_spec, labels = "AUTO", nrow = 2), 
       filename = "figures_revision/real_data_plot_imperfect_subjectC13.png", 
       width = 6, height = 8)
