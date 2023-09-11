## Title: Subject C1 with Imperfect Assays (Figure S6)

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
figureS6_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/subject_c1_imperfect.csv")

## create Figure S6

dplyr::mutate(sensQVOA = factor(x = sensQVOA, 
                                  levels = seq(1, 0.8, by = -0.1), 
                                  labels = paste0("QVOA Sensitivity = ", seq(100, 80, by = -10), "%")),
                sensUDSA = factor(x = sensUDSA, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(80, 100, by = 10), "%"))
  ) |>
  dplyr::select(sensQVOA, sensUDSA, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:2)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Imperfect Assays)", "MLE (Perfect Assays)")),
                M = "(6, 12, 18)") 

plot_data |> 
  ggplot(aes(x = M, y = value, fill = estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, 
             linetype = 2, 
             color = "blue") +
  facet_grid(rows = vars(sensQVOA), 
             cols = vars(sensUDSA), 
             scales = "free") + 
  ggthemes::scale_fill_colorblind(name = "Method") +
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("IUPM Estimate") + 
  xlab("Number of Replicate Wells (M)") 


figureS6_data |>
  mutate(bias.correction = factor(bias.correction,
                                  levels = c("MLE", "BC-MLE")),
         method = factor(method,
                         levels = c("Without UDSA (Multiple Dilutions)",
                                    "With UDSA (Multiple Dilutions)",
                                    "With UDSA (Single Dilution)")),
         id.lab = factor(paste0("Subject ", id),
                         levels = paste0("Subject C", 1:17)),
         method.bc = factor(paste0(method, "_", bias.correction)),
         method.bc = factor(method.bc, levels = rev(levels(method.bc)))) |>
  ggplot(aes(x = bias.correction, y = log(mle),
             ymin = log(ci.lower), ymax = log(ci.upper),
             color = method, shape = bias.correction)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.3,
                position = position_dodge(width = 0.5)) +
  theme_light() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.95, 0.07),
        legend.justification = c("right", "bottom"),
        legend.box = "horizontal",
        legend.key.width = unit(1.5, "cm")) +
  ylab("IUPM of HIV") +
  labs(color = "Method",
       shape = "Bias Correction") +
  ggthemes::scale_colour_colorblind() +
  facet_wrap(~ id.lab, ncol = 5, scales = "free") -> real_data_plot_log_IUPM

real_data_plot_log_IUPM +
  ggtitle("IUPM MLEs and 95% Confidence Intervals by Method",
          subtitle="MLEs and Bias-Corrected MLEs")

ggsave(real_data_plot_log_IUPM, filename = "figures_revision/real_data_plot_log_IUPM.png", width = 10, height = 10)
