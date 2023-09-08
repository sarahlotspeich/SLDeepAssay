library(ggplot2)

Results = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/sd_imperfect_vary_q.csv") |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, ## Make any reps that didn't converge NA 
                    yes = as.numeric(Lambda), 
                    no = NA)
  )

Results |> 
  dplyr::group_by(M, q) |> 
  dplyr::summarize(
    corrected_reps = sum(!is.na(Lambda)),
    uncorrected_reps = sum(!is.na(Lambda_naive))
  ) |> 
  dplyr::ungroup() |> 
  dplyr::summarize(
    min_corrected_reps = min(corrected_reps), 
    min_uncorrected_reps = min(uncorrected_reps) 
  )
## MLE (Imperfect Assays) converged in >= 1000 / 1000 reps per setting
## MLE (Perfect Assays) converged in 1000 / 1000 reps per setting
table(Results$msg) ### 9000 / 9000 replicates total 

plot_data = Results |>
  dplyr::mutate(q = factor(x = q, 
                           levels = c(1, 0.75, 0.5), 
                           labels = paste0("q = ", c(1, 0.75, 0.5)))
  ) |>
  dplyr::select(M, q, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:2)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Imperfect Assays)", "MLE (Perfect Assays)"))) 

# Exclude 0 replicates where MLE (Imperfect Assays) was > 10
table(plot_data$value > 10, plot_data$estimator)

plot_data |> 
  ggplot(aes(x = factor(M), y = value, fill = estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, 
             linetype = 2, 
             color = "blue") +
  facet_grid(cols = vars(q), 
             scales = "free") + 
  ggthemes::scale_fill_colorblind(name = "Method") +
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("IUPM Estimate") + 
  xlab("Number of Replicate Wells (M)") 
ggsave(filename = "figS4-boxplot-vary-q.png", 
       device = "png", units = "in", width = 8, height = 4)

