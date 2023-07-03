## Title: IUPM Real Data Plot (Figure 2)

## Author:  Brian Richardson

## Date:  01/04/2023

## Purpose: Plot real data results for SLDeepAssay paper


## Prepare workspace
rm(list=ls())
setwd(here::here())

## load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
figure2_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real_data_application/figure2_data.csv")

## create Figure 2
figure2_data |>
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
