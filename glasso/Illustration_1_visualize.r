library(tidyverse)
library(glasso)
library(ggplot2)
library(reshape2)
library(BDgraph)
library(patchwork)
RhpcBLASctl::blas_set_num_threads(16)
RhpcBLASctl::omp_set_num_threads(1)
source("functions.R")

s <- 1
load("glasso/results_LR.RData")
# ── Assemble ROC curves ──
roc_df_sp <- do.call(rbind, lapply(seq_along(results), function(i) {
  results[[i]]$roc_sp$roc$rep <- i
  results[[i]]$roc_sp$roc
})) %>% mutate(method = "Spatial")

roc_df_nn <- do.call(rbind, lapply(seq_along(results), function(i) {
  results[[i]]$roc_nn$roc$rep <- i
  results[[i]]$roc_nn$roc
})) %>% mutate(method = "Independent")

roc_df <- bind_rows(roc_df_sp, roc_df_nn)

# ── Average ROC curves (same as before) ──
fpr_grid <- seq(0, 1, length.out = 200)

roc_avg <- roc_df %>%
  group_by(rep, method) %>%
  arrange(fpr, tpr) %>%
  summarize(tpr_interp = list(approx(fpr, tpr, xout = fpr_grid, rule = 2)$y),
            .groups = "drop") %>%
  unnest_wider(tpr_interp, names_sep = "_") %>%
  pivot_longer(starts_with("tpr_interp_"), names_to = "idx", values_to = "tpr") %>%
  mutate(idx = as.integer(gsub("tpr_interp_", "", idx))) %>%
  mutate(fpr = fpr_grid[idx]) %>%
  group_by(method, fpr) %>%
  summarize(tpr = mean(tpr), .groups = "drop")

# ── Ribbon data: difference between spatial and independent averages ──
ribbon_df <- roc_avg %>%
  pivot_wider(names_from = method, values_from = tpr) %>%
  mutate(method = "Spatial")  # ribbon only in the Spatial facet

# ── Independent average line, to appear in the Spatial facet ──
indep_avg_in_spatial <- roc_avg %>%
  filter(method == "Independent") %>%
  mutate(method = "Spatial")  # assign to Spatial facet

# ── Panel A: ROC ──
p_ROC <- ggplot(roc_df, aes(fpr, tpr, group = rep)) +
  geom_path(alpha = 0.05, linewidth = 0.3) +
  geom_ribbon(data = ribbon_df,
              aes(x = fpr, ymin = Independent, ymax = Spatial, group = 1),
              fill = "#0072B2", alpha = 0.1, inherit.aes = FALSE) +
  geom_path(data = roc_avg, aes(fpr, tpr, group = 1),
            linewidth = 1, color = "black") +
  #geom_path(data = indep_avg_in_spatial, aes(fpr, tpr, group = 1),
  #          linewidth = 0.7, color = "black", lty = 3) +
  geom_abline(lty = 2, color = "grey50") +
  labs(x = "False positive rate", y = "True positive rate") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 11) +
  facet_grid(~method) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        plot.margin=margin(0,0,0,0))


# Extract number of true edges per replicate from results
q <- 10
total_pairs <- q * (q - 1) / 2

edge_counts <- sapply(seq_along(results), function(i) {
  # adjust accessor to match your results structure
  sum(results[[i]]$graph[upper.tri(results[[i]]$graph)] != 0)
})

# Add P and N to roc_df, then compute standard F1
best_f1_df <- roc_df %>%
  left_join(tibble(rep = seq_along(results),
                   P = edge_counts,
                   N = total_pairs - edge_counts),
            by = "rep") %>%
  group_by(rep, method) %>%
  arrange(fpr, tpr) %>%
  mutate(precision = (tpr * P) / (tpr * P + fpr * N),
         recall    = tpr,
         f1 = ifelse(precision + recall == 0, 0,
                     2 * precision * recall / (precision + recall))) %>%
  slice_max(f1, n = 1, with_ties = FALSE) %>%
  ungroup()


paired_df <- best_f1_df %>%
  transmute(rep, method,
            Sensitivity = tpr,
            Specificity = 1 - fpr)

# Win rates
comp_wide <- paired_df %>%
  pivot_wider(names_from = method, values_from = c(Sensitivity, Specificity))
n_reps <- nrow(comp_wide)


diff_df <- comp_wide %>%
  transmute(rep,
            dSensitivity = Sensitivity_Spatial - Sensitivity_Independent,
            dSpecificity = Specificity_Spatial - Specificity_Independent)

diff_means <- diff_df %>%
  summarize(dSensitivity = mean(dSensitivity),
            dSpecificity = mean(dSpecificity))

major_breaks <- c(-0.2, 0, 0.2, 0.4)

fmt_breaks <- function(b) {
  ifelse(b == 0, "0",
         ifelse(b %in% major_breaks, sprintf("%.1f", b), sprintf("%.2f", b)))
}

fmt_x <- fmt_breaks
fmt_y <- fmt_breaks

p_comp <- ggplot(diff_df, aes(y = dSensitivity, x = dSpecificity)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey50") +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  stat_ellipse(geom = "polygon", alpha = 0.12, fill = "grey50",
               color = "grey40", linewidth = 0.2,
               level = 0.9, type = "t") +
  geom_point(alpha = 0.3, size = 1.5, color = "grey30") +
  geom_point(data = diff_means, shape = 8, size = 4,
             stroke = 1.2, color = "grey30", fill = "white") +
  scale_y_continuous(
    breaks = sort(unique(c(major_breaks, diff_means$dSensitivity))),
    labels = fmt_y,
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(major_breaks, diff_means$dSpecificity))),
    labels = fmt_x,
    minor_breaks = NULL
  ) +
  labs(y = expression(Delta * " Sensitivity"),
       x = expression(Delta * " Specificity")) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        panel.grid.minor = element_blank(),
        plot.margin=margin(0,0,0,0))

fig <- p_ROC + p_comp +
  plot_layout(widths = c(2, 1.2)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = c(0.5, -0.05),
        plot.margin = margin(-5, 5, 10, 5))

ggsave("glasso/Illustration_1.pdf", width=8, height=3)


auc_df <- roc_df %>%
  group_by(rep, method) %>%
  arrange(fpr, tpr) %>%
  summarize(auc = sum(diff(fpr) * (tpr[-1] + tpr[-n()]) / 2),
            .groups = "drop")

auc_df %>%
  group_by(method) %>%
  summarize(mean_auc = mean(auc),
            sd_auc = sd(auc))

diff_df %>%
  summarize(
    pct_both = mean(dSensitivity >= 0 & dSpecificity >= 0) * 100,
    pct_some = mean(dSensitivity >= 0 | dSpecificity >= 0) * 100
  )
