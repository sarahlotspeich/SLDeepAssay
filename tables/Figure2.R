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
results = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/real-data-application/tableS5_data.csv")

## create plot data
plot.dat = rbind(results[, c("id", "method", "mle", "se", "ci1", "ci2")],
                 results[, c("id", "method", "mle_bc", "se", "ci_bc1", "ci_bc2")]) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("id", "method",
                 "mle", "se", "ci.lower", "ci.upper")) %>% 
  mutate_at(c("mle", "se", "ci.lower", "ci.upper"), as.numeric) %>% 
  mutate(bias.correction = factor(rep(c("MLE", "BC-MLE"),
                                      each = nrow(results)),
                                  levels = c("MLE", "BC-MLE")),
         method = factor(method,
                         levels = c("Without UDSA (Multiple Dilutions)",
                                    "With UDSA (Single Dilution)",
                                    "With UDSA (Multiple Dilutions)")),
         id.lab = factor(paste0("Subject ", id),
                         levels = paste0("Subject C", 1:17))) %>% 
  mutate(method.bc = factor(paste0(method, "_", bias.correction)),
         method.bc = factor(method.bc, levels = rev(levels(method.bc))))

## create plot
ggplot(plot.dat, aes(x = method.bc, y = mle, ymin = ci.lower, ymax = ci.upper,
                     color = method, shape = bias.correction)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  theme_light() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(.9, .07),
        legend.justification = c("right", "bottom"),
        legend.box = "horizontal",
        legend.key.width = unit(1.5, "cm")) +
  ylab("IUPM MLE and 95% CI") +
  labs(color = "Method",
       shape = "Bias Correction") +
  scale_color_grey() +
  facet_wrap(~ id.lab, ncol = 5, scales = "free") +
  scale_y_continuous(trans = "log10") -> real_data_plot

real_data_plot +
  ggtitle("IUPM MLEs and 95% Confidence Intervals by Method",
          subtitle="MLEs and Bias-Corrected MLEs")

ggsave(real_data_plot, filename = "real_data_application/real_data_plot.png", width = 10, height = 10)
