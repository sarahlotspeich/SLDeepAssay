d = "~/Downloads/vary-spec/"
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
## Corrected IUPM estimator converged in >= 992 / 1000 reps per setting
## Uncorrected IUPM estimator converged in >= 960 / 1000 reps per setting
table(Results$msg)
table(Results$msg_naive)

library(ggplot2)
Results |>
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
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
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
  xlab("Number of Replicate Wells (M)") + 
  ggtitle(label = "Results under varied specificities", 
          subtitle = "Assume: sensitivity = 0.9, n = 6 DVLs, true IUPM = 1, q = 1") 
ggsave(filename = "~/Documents/SLDeepAssay/figures/boxplot-vary-spec.png", 
       device = "png", units = "in", width = 8, height = 8)

Results |>
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
                                   labels = c("MLE (Corrected)", "MLE (Uncorrected)"))) |> 
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
ggsave(filename = "~/Documents/SLDeepAssay/figures/boxplot-vary-spec-plain.png", 
       device = "png", units = "in", width = 8, height = 8)