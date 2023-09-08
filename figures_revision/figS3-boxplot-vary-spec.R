library(ggplot2)

Results = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_imperfect_vary_specificity.csv") |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, ## Make any reps that didn't converge NA 
                    yes = as.numeric(Lambda), 
                    no = NA)
  )

Results |> 
  dplyr::group_by(M, specQVOA, specUDSA) |> 
  dplyr::summarize(
    corrected_reps = sum(!is.na(Lambda)),
    uncorrected_reps = sum(!is.na(Lambda_naive))
  ) |> 
  dplyr::ungroup() |> 
  dplyr::summarize(
    min_corrected_reps = min(corrected_reps), 
    min_uncorrected_reps = min(uncorrected_reps) 
  )
## MLE (Imperfect Assays) converged in 1000 / 1000 reps per setting
## MLE (Perfect Assays) converged in 1000 / 1000 reps per setting
table(Results$msg) ### 27000 / 270000 replicates total 

# Exclude 159 replicates where MLE (Imperfect Assays) was > 10, 
## 36 replicates where MLE (Perfect Assays) was > 10.
table(Results$value > 10, Results$estimator)

plot_data = Results |>
  dplyr::mutate(specQVOA = factor(x = specQVOA, 
                                  levels = seq(1, 0.8, by = -0.1), 
                                  labels = paste0("QVOA Specificity = ", seq(100, 80, by = -10), "%")),
                specUDSA = factor(x = specUDSA, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Specificity = ", seq(80, 100, by = 10), "%"))
  ) |>
  dplyr::select(M, q, specQVOA, specUDSA, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:4)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Imperfect Assays)", "MLE (Perfect Assays)")))

# Exclude 1 replicate where MLE (Imperfect Assays) was > 10, 
## 0 replicates where MLE (Perfect Assays) was > 10.
table(plot_data$value > 10, plot_data$estimator)

# Count number of extreme replicates per estimator + setting
plot_data |> 
  dplyr::group_by(estimator, specQVOA, specUDSA, M) |> 
  dplyr::summarize(num_extreme_vals = sum(value > 10)) |> 
  dplyr::filter(num_extreme_vals > 0) |> 
  dplyr::arrange(desc(num_extreme_vals))
## No more than 1 / 1000 = 0.1% per estimator + setting

plot_data |> 
  dplyr::filter(value <= 10) |> 
  ggplot(aes(x = factor(M), y = value, fill = estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, 
             linetype = 2, 
             color = "blue") +
  facet_grid(rows = vars(specQVOA), 
             cols = vars(specUDSA), 
             scales = "free") + 
  ggthemes::scale_fill_colorblind(name = "Method") +
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("Estimated IUPM") + 
  xlab("Number of Replicate Wells (M)")  
ggsave(filename = "figS3-boxplot-vary-spec.png", 
       device = "png", units = "in", width = 8, height = 8)
