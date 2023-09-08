library(ggplot2)
library(dplyr)
library(tidyr)

# load data
overdispersion_sim_data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/overdispersion_sim_data.csv")

# create long data set
dat_long <- overdispersion_sim_data |> 
  filter(error == 0) |> 
  pivot_longer(cols = c(mle_pois, mle_pois_bc, mle_negbin)) |> 
  mutate(method.label = ifelse(name == "mle_negbin", "Negative Binomial",
                               ifelse(name == "mle_pois", "Poisson",
                                      "Poisson (BC)")),
         gamma = factor(gamma),
         M.label = factor(M.label,
                          levels = c("(6, 12, 18)", "(30, 60, 90)")))

# label number of points beyond range of plot
y.max <- 5

# facet labels
M.labs <- paste0("M = ", levels(dat_long$M.label))
names(M.labs) <- 1:length(M.labs)

upper.counts <- dat_long |> 
  group_by(gamma, M, method.label) |> 
  summarise(upper.count = sum(value > y.max)) |> 
  mutate(upper.count = ifelse(upper.count > 0,
                              paste0(upper.count, "*"),
                              ""))

# boxplot of estimated IUPM
negbin_mle_boxplot <- ggplot(
  filter(dat_long, value <= y.max),
  aes(y = value,
      fill = gamma,
      x = method.label)) +
  geom_boxplot(position = position_dodge2(width = 0.8)) +
  geom_text(data = upper.counts,
            y = y.max + 0.1,
            aes(label = upper.count,
                x = method.label,
                group = gamma),
            hjust = 0.5,
            size = 3,
            position = position_dodge2(width = 0.8)) +
  geom_hline(mapping = aes(yintercept = 1),
             color = "blue",
             linetype = "dashed") +
  coord_cartesian(ylim = c(0, y.max + 0.1)) +
  facet_wrap(~ M, ncol = 4,
             labeller = labeller(M = M.labs)) +
  labs(x = "Method",
       y = "Estimated IUPM",
       fill = "\u03b3") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = "black")) +
  ggthemes::scale_fill_colorblind()


negbin_mle_boxplot +
  ggtitle("Distribution of Estimated IUPM",
          subtitle = "Compared to True IUPM (Blue)")

ggsave("figures_revision/figS1-boxplot-overdispersion.png",
       plot = negbin_mle_boxplot, 
       dpi = 300,
       width = 8, height = 5)

