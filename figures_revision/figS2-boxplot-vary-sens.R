library(ggplot2)

Results = read.csv("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_imperfect_vary_sensitivity.csv") |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, ## Make any reps that didn't converge NA
                    yes = as.numeric(Lambda), 
                    no = NA)
  )

Results |> 
  dplyr::group_by(M, sensQVOA, sensUDSA) |> 
  dplyr::summarize(
    corrected_reps = sum(!is.na(Lambda)),
    uncorrected_reps = sum(!is.na(Lambda_naive))
  ) |> 
  dplyr::ungroup() |> 
  dplyr::summarize(
    min_corrected_reps = min(corrected_reps), 
    min_uncorrected_reps = min(uncorrected_reps) 
  )
## MLE (Imperfect Assays) converged in >= 996 / 1000 reps per setting
## MLE (Perfect Assays) converged in 1000 / 1000 reps per setting
table(Results$msg) ### only 19 / 270000 replicates total 

plot_data = Results |>
  dplyr::mutate(sensQVOA = factor(x = sensQVOA, 
                                  levels = seq(1, 0.8, by = -0.1), 
                                  labels = paste0("QVOA Sensitivity = ", seq(100, 80, by = -10), "%")),
                sensUDSA = factor(x = sensUDSA, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(80, 100, by = 10), "%"))
  ) |>
  dplyr::select(M, q, sensQVOA, sensUDSA, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:4)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Imperfect Assays)", "MLE (Perfect Assays)"))) 

# Exclude 2 replicates where MLE (Imperfect Assays) was > 10
table(plot_data$value > 10, plot_data$estimator)

plot_data |> 
  dplyr::filter(value <= 10) |> 
  ggplot(aes(x = factor(M), y = value, fill = estimator)) +
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
ggsave(filename = "figS2-boxplot-vary-sens.png", 
       device = "png", units = "in", width = 8, height = 8)