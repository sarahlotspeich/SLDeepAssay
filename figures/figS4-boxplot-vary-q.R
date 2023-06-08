d = "~/Downloads/vary-q/"
f = paste0(d, list.files(d))
Results = do.call(rbind, 
                  lapply(X = f, 
                         FUN = read.csv)) |> 
  dplyr::mutate(
    Lambda = ifelse(conv == 0, 
                    yes = as.numeric(Lambda), 
                    no = NA)#,
    # Lambda_naive = ifelse(conv_naive == 0, 
    #                       yes = as.numeric(Lambda_naive), 
    #                       no = NA)
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
## Corrected IUPM estimator converged in >= 993 / 1000 reps per setting
## Uncorrected IUPM estimator converged in >= 1000 / 1000 reps per setting
table(Results$msg)
table(Results$msg_naive)

Results |>
  dplyr::mutate(q = factor(x = q, 
                           levels = c(1, 0.75, 0.5), 
                           labels = paste0("q = ", c(1, 0.75, 0.5)))
  ) |>
  dplyr::select(M, q, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:2)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
  ggplot(aes(x = factor(M), y = value, fill = estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, 
             linetype = 2, 
             color = "blue") +
  facet_grid(cols = vars(q), 
             scales = "free") + 
  ggthemes::scale_fill_colorblind(name = "Method") +
  # scale_fill_manual(values = wesanderson::wes_palette("Royal2", n = 3)[c(1, 3)], 
  #                   name = "Method") + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  ylab("IUPM Estimate") + 
  xlab("Number of Replicate Wells (M)") + 
  ylim(c(-1, 6)) +
  ggtitle(label = "Results under varied proportions deep-sequenced (q)", 
          subtitle = "Assume: sensitivity = 0.9, specificity = 0.9, n = 6 DVLs, true IUPM = 1") +
  labs(caption = "*Two extreme replicates where the corrected IUPM was > 10 were excluded (one with q = 1, one with q = 0.5).")
ggsave(filename = "~/Documents/SLDeepAssay/figures/boxplot-vary-q.png", 
       device = "png", units = "in", width = 8, height = 4)

Results |>
  dplyr::mutate(q = factor(x = q, 
                           levels = c(1, 0.75, 0.5), 
                           labels = paste0("q = ", c(1, 0.75, 0.5)))
  ) |>
  dplyr::select(M, q, Lambda, Lambda_naive) |>
  tidyr::gather("estimator", "value", -c(1:2)) |>
  dplyr::mutate(estimator = factor(x = estimator, 
                                   levels = c("Lambda", "Lambda_naive"),
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
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
  xlab("Number of Replicate Wells (M)") + 
  ylim(c(-1, 6))
ggsave(filename = "~/Documents/SLDeepAssay/figures/boxplot-vary-q-plain.png", 
       device = "png", units = "in", width = 8, height = 4)

