d = "~/Downloads/vary-sens/"
f = paste0(d, list.files(d))
Results = do.call(rbind, 
                  lapply(X = f, 
                         FUN = read.csv)) |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, 
                    yes = as.numeric(Lambda), 
                    no = NA),
    Lambda_naive = ifelse(conv_naive == 0, 
                          yes = as.numeric(Lambda_naive), 
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
## Corrected IUPM estimator converged in >= 997 / 1000 reps per setting
## Uncorrected IUPM estimator converged in >= 957 / 1000 reps per setting
table(Results$msg)
table(Results$msg_naive)

Results |>
  dplyr::mutate(sensQVOA = factor(x = sensQVOA, 
                                  levels = seq(1, 0.8, by = -0.1), 
                                  labels = paste0("QVOA Sensitivity = ", seq(1, 0.8, by = -0.1))),
                sensUDSA = factor(x = sensUDSA, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(0.8, 1, by = 0.1)))
  ) |>
  dplyr::select(M, q, sensQVOA, sensUDSA, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:4)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
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
  #scale_fill_manual(values = wesanderson::wes_palette("Royal2", n = 3)[c(1, 3)], 
  #                  name = "Method") + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("IUPM Estimate") + 
  xlab("Number of Replicate Wells (M)") + 
  ggtitle(label = "Results under varied sensitivities", 
          subtitle = "Assume: specificity = 0.9, n = 6 DVLs, true IUPM = 1, q = 1") +
  labs(caption = "*Two extreme replicates where the corrected IUPM was > 10 were excluded (one with QVOA/UDSA sensitivity 0.9/0.8, one with 0.8/0.8).") 
ggsave(filename = "~/Downloads/boxplot-vary-sens.png", 
       device = "png", units = "in", width = 8, height = 8)

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
## Corrected IUPM estimator converged in >= 997 / 1000 reps per setting
## Uncorrected IUPM estimator converged in >= 957 / 1000 reps per setting
table(Results$msg)
table(Results$msg_naive)

Results |>
  dplyr::mutate(sensQVOA = factor(x = sensQVOA, 
                                  levels = seq(1, 0.8, by = -0.1), 
                                  labels = paste0("QVOA Sensitivity = ", seq(1, 0.8, by = -0.1))),
                sensUDSA = factor(x = sensUDSA, 
                                  levels = seq(0.8, 1, by = 0.1), 
                                  labels = paste0("UDSA Sensitivity = ", seq(0.8, 1, by = 0.1)))
  ) |>
  dplyr::select(M, q, sensQVOA, sensUDSA, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:4)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
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
  #scale_fill_manual(values = wesanderson::wes_palette("Royal2", n = 3)[c(1, 3)], 
  #                  name = "Method") + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("IUPM Estimate") + 
  xlab("Number of Replicate Wells (M)") 
ggsave(filename = "~/Downloads/boxplot-vary-sens-plain.png", 
       device = "png", units = "in", width = 8, height = 8)
