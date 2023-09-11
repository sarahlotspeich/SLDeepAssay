library(ggplot2)

Results = do.call(what = dplyr::bind_rows, 
                  args = lapply(X = paste0("https://raw.githubusercontent.com/sarahlotspeich/SLDeepAssay/main/sim_data/md-imperfect/md-imperfect-seed", 11422:11431, ".csv"), 
                                FUN = read.csv)) |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, ## Make any reps that didn't converge NA 
                    yes = as.numeric(Lambda), 
                    no = NA),
    Lambda_naive = ifelse(conv_naive == 0, ## Make any reps that didn't converge NA 
                    yes = as.numeric(Lambda_naive), 
                    no = NA)
  )

Results |> 
  dplyr::group_by(sensQVOA, sensUDSA) |> 
  dplyr::summarize(
    corrected_reps = sum(!is.na(Lambda)),
    uncorrected_reps = sum(!is.na(Lambda_naive))
  ) |> 
  dplyr::ungroup() |> 
  dplyr::summarize(
    min_corrected_reps = min(corrected_reps), 
    min_uncorrected_reps = min(uncorrected_reps) 
  )
## MLE (Imperfect Assays) converged in >= 999 / 1000 reps per setting
table(Results$msg) ### only 1 / 9000 replicates total 
## MLE (Perfect Assays) converged in >=992 / 1000 reps per setting
table(Results$msg_naive) ### only 14 / 9000 replicates total 

plot_data = Results |>
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

ggsave(filename = "~/Documents/SLDeepAssay/figures/figS5-boxplot-imperfect-md.png", 
       device = "png", units = "in", width = 8, height = 8)

