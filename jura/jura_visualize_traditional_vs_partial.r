rm(list=ls())
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)
library(spiox)
library(tidyr)
library(igraph)
tail_mcmc <- 500
precompute_all <- function(iox_all_mcmc, q, mcmc, tail_mcmc = 500, n_threads = 16, max_dist = 100, max_range=0.5) {
  
  # max_dist: avoid computing all pairwise distances
  # max_range: target distances<max_range since that's where the action is
  
  R_mcmc   <- iox_all_mcmc$Sigma %>%
    apply(3, \(x) cov2cor(x)) %>% array(dim = c(q, q, mcmc))
  Qpc_mcmc <- iox_all_mcmc$Sigma %>%
    apply(3, \(x) -cov2cor(solve(x))) %>% array(dim = c(q, q, mcmc))
  
  pre <- list()
  for (i in 1:q) pre[[i]] <- list()
  
  for (i in 1:q) {
    for (j in i:q) {
      fij <- spiox:::iox_fij(iox_all_mcmc, i, j,
                             n_threads = n_threads, matern = 1, tail_mcmc = tail_mcmc, max_dist, max_range)
      
      if (i == j) {
        pre[[i]][[j]] <- list(
          dists  = fij$dists,
          f_mean = rowMeans(fij$f_ij),
          f_low  = apply(fij$f_ij, 1, quantile, 0.025),
          f_upp  = apply(fij$f_ij, 1, quantile, 0.975)
        )
      } else {
        R_tail_ij   <- tail(R_mcmc[i, j, ],   tail_mcmc)
        Qpc_tail_ij <- tail(Qpc_mcmc[i, j, ], tail_mcmc)
        C_mat_ij    <- t(t(fij$f_ij) * as.vector(R_tail_ij))
        Qpc_mat_ij  <- t(t(fij$f_ij) * as.vector(Qpc_tail_ij))
        
        pre[[i]][[j]] <- list(
          dists    = fij$dists,
          C_mean   = rowMeans(C_mat_ij),
          C_low    = apply(C_mat_ij, 1, quantile, 0.025),
          C_upp    = apply(C_mat_ij, 1, quantile, 0.975),
          Qpc_mean = rowMeans(Qpc_mat_ij),
          Qpc_low  = apply(Qpc_mat_ij, 1, quantile, 0.025),
          Qpc_upp  = apply(Qpc_mat_ij, 1, quantile, 0.975)
        )
        pre[[j]][[i]] <- pre[[i]][[j]]
      }
      cat("Done (", i, ",", j, ") \n")
    }
  }
  pre
}
cat("Precomputing all pairs...\n")
#iox_pre <- precompute_all(iox_all_mcmc, q, mcmc, tail_mcmc = tail_mcmc, n_threads = 16)
cat("Done.\n")

#save(iox_pre, file = "jura/iox_precomputed.RData")

# ── Colors ────────────────────────────────────────────────────────────────────
bg_partial  <- "#F9F6EE"
bg_marginal <- "#EEF0F3"
bg_diag     <- "#F3F2EF"

lc_gpgpm    <- "#EE7733"
lc_iox      <- "#2166AC"

# ── Load fits ─────────────────────────────────────────────────────────────────
source("functions.R")
load("jura/GpGpm_jura_fit.RData")
load("jura/spioxmcmc_jura_fit.RData")
load("jura/iox_precomputed.RData")
#load("jura/lmc_jura_fit.RData")
#load("jura/GpGpm_jura_preds.RData")

mcmc <- dim(iox_all_mcmc$Sigma)[3]
outcome_names <- c("Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn")
q <- length(outcome_names)
tail_mcmc <- 500

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_panel <- function(bg, show_x = FALSE, show_y = FALSE) {
  theme(
    panel.background  = element_rect(fill = bg, colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(linewidth = 0.3, colour = "grey50"),
    axis.ticks        = element_line(linewidth = 0.25, colour = "grey50"),
    axis.ticks.length = unit(1.2, "pt"),
    plot.title        = element_text(size = 8, hjust = 0.5, face = "plain",
                                     margin = margin(0, 0, 1, 0, "pt")),
    plot.margin       = margin(1, 2, 1, 2, "pt"),
    axis.text.x       = if (show_x) element_text(size = 7, margin = margin(1, 0, 0, 0, "pt"))
    else element_blank(),
    axis.ticks.x      = if (show_x) element_line(linewidth = 0.25)
    else element_blank(),
    axis.text.y       = if (show_y) element_text(size = 7) else element_blank(),
    axis.ticks.y      = if (show_y) element_line(linewidth = 0.25) else element_blank()
  )
}

# ── Panel builders (unchanged) ────────────────────────────────────────────────
make_est_panel <- function(i, j, type = c("partial", "marginal"),
                           iox_pre, gpgpm_fit, names_vec,
                           show_x = FALSE, show_y = FALSE) {
  type <- match.arg(type)
  hh   <- iox_pre[[i]][[j]]$dists
  dd   <- 2
  
  if (type == "partial") {
    bg      <- bg_partial
    y_iox   <- iox_pre[[i]][[j]]$Qpc_mean
    y_low   <- iox_pre[[i]][[j]]$Qpc_low
    y_upp   <- iox_pre[[i]][[j]]$Qpc_upp
    Qpc_fit <- -cov2cor(gpgpm_fit$Q); diag(Qpc_fit) <- 1
    g_fit   <- get_gamma_mat(gpgpm_fit$nu, dd)[i, j]
    M_fit   <- matern_cor(hh, (gpgpm_fit$nu[i] + gpgpm_fit$nu[j]) / 2, gpgpm_fit$phi)
    y_gpgpm <- Qpc_fit[i, j] * g_fit * M_fit
  } else {
    bg      <- bg_marginal
    y_iox   <- iox_pre[[i]][[j]]$C_mean
    y_low   <- iox_pre[[i]][[j]]$C_low
    y_upp   <- iox_pre[[i]][[j]]$C_upp
    R_fit   <- cov2cor(gpgpm_fit$Sigma)
    g_fit   <- get_gamma_mat(gpgpm_fit$nu, dd)[i, j]
    M_fit   <- matern_cor(hh, (gpgpm_fit$nu[i] + gpgpm_fit$nu[j]) / 2, gpgpm_fit$phi)
    y_gpgpm <- R_fit[i, j] * g_fit * M_fit
  }
  
  df <- data.frame(h = hh, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  ggplot(df, aes(x = h)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp), fill = lc_iox, alpha = 0.15) +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox) +
    coord_cartesian(ylim = c(-1, 1), xlim = c(0, 0.3)) +
    labs(title = paste0(names_vec[i], ", ", names_vec[j]), x = NULL, y = NULL) +
    theme_panel(bg, show_x = show_x, show_y = show_y)
}

make_diag_panel <- function(idx, iox_pre, gpgpm_fit, names_vec,
                            show_x = TRUE, show_y = TRUE) {
  hh      <- iox_pre[[idx]][[idx]]$dists
  y_iox   <- iox_pre[[idx]][[idx]]$f_mean
  y_low   <- iox_pre[[idx]][[idx]]$f_low
  y_upp   <- iox_pre[[idx]][[idx]]$f_upp
  y_gpgpm <- matern_cor(hh, gpgpm_fit$nu[idx], gpgpm_fit$phi)
  
  df <- data.frame(h = hh, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  ggplot(df, aes(x = h)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp), fill = lc_iox, alpha = 0.15) +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.3)) +
    labs(title = names_vec[idx], x = NULL, y = NULL) +
    theme_panel(bg_diag, show_x = show_x, show_y = show_y)
}

make_effrange_panel <- function(i, j, type = c("partial", "marginal"),
                                iox_pre, gpgpm_fit, names_vec,
                                h = seq(1e-4, 1, length.out = 500),
                                thr = 0.05,
                                show_x = TRUE, show_y = TRUE) {
  type <- match.arg(type)
  dd   <- 2
  g_fit <- get_gamma_mat(gpgpm_fit$nu, dd)[i, j]
  M_fit <- matern_cor(h, (gpgpm_fit$nu[i] + gpgpm_fit$nu[j]) / 2, gpgpm_fit$phi)
  ij    <- iox_pre[[i]][[j]]
  
  if (type == "partial") {
    Qpc_fit <- -cov2cor(gpgpm_fit$Q); diag(Qpc_fit) <- 1
    y_gpgpm <- Qpc_fit[i, j] * g_fit * M_fit
    y_iox   <- approx(ij$dists, ij$Qpc_mean, xout = h, rule = 2)$y
    y_low   <- approx(ij$dists, ij$Qpc_low,  xout = h, rule = 2)$y
    y_upp   <- approx(ij$dists, ij$Qpc_upp,  xout = h, rule = 2)$y
    bg      <- bg_partial
    ttl     <- paste0("Partial eff. cross-range: ", names_vec[i], ", ", names_vec[j])
  } else {
    R_fit   <- cov2cor(gpgpm_fit$Sigma)
    y_gpgpm <- R_fit[i, j] * g_fit * M_fit
    y_iox   <- approx(ij$dists, ij$C_mean, xout = h, rule = 2)$y
    y_low   <- approx(ij$dists, ij$C_low,  xout = h, rule = 2)$y
    y_upp   <- approx(ij$dists, ij$C_upp,  xout = h, rule = 2)$y
    bg      <- bg_marginal
    ttl     <- paste0("Eff. cross-range ", names_vec[i], ", ", names_vec[j])
  }
  
  find_eff_range <- function(hh, yy, thr) {
    idx <- which(abs(yy) < thr)
    if (length(idx) == 0) return(NA)
    hh[min(idx)]
  }
  
  sgn  <- sign(y_gpgpm[1])
  er_g <- find_eff_range(h, y_gpgpm, thr)
  er_i <- find_eff_range(h, y_iox,   thr)
  
  df <- data.frame(h = h, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  p <- ggplot(df, aes(x = h)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -thr, ymax = thr),
              fill = "grey90", colour = NA, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp), fill = lc_iox, alpha = 0.12) +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox)
  
  cross_pts <- data.frame(h = c(er_g, er_i), y = rep(sgn * thr, 2),
                          col = c(lc_gpgpm, lc_iox))
  cross_pts <- cross_pts[!is.na(cross_pts$h), ]
  if (nrow(cross_pts) > 0) {
    p <- p + geom_point(data = cross_pts, aes(x = h, y = y),
                        colour = cross_pts$col, fill = cross_pts$col,
                        shape = 16, size = 2.5, inherit.aes = FALSE)
  }
  
  p + coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.3)) +
    labs(title = ttl, x = "Distance", y = NULL) +
    theme_panel(bg, show_x = show_x, show_y = show_y)
}

# ── Shared legend ─────────────────────────────────────────────────────────────
make_legend <- function(label_x = "Distance") {
  leg_df <- data.frame(
    x = rep(1:2, 2), y = rep(1:2, 2),
    method = factor(rep(c("Pars. Matérn", "Inside-out cross-cov."), each = 2),
                    levels = c("Pars. Matérn", "Inside-out cross-cov."))
  )
  ggplot(leg_df, aes(x, y, colour = method)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("Pars. Matérn" = lc_gpgpm,
                                   "Inside-out cross-cov." = lc_iox)) +
    guides(colour = guide_legend(nrow = 1, title = NULL)) +
    labs(title = label_x) +
    theme_void() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0,0,0,0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),
          legend.text = element_text(size = 8), legend.key.width = unit(18, "pt"),
          plot.title = element_text(size = 9, hjust = 0.5, margin = margin(0,0,0,0,"pt")),
          plot.margin = margin(0,0,0,0,"pt")) +
    coord_cartesian(xlim = c(0,0), ylim = c(0,0))
}

# ── Heatmap builder ──────────────────────────────────────────────────────────
make_part_corr <- function(Q) {
  R <- -cov2cor(Q)
  diag(R) <- 1
  R
}

summarize_matrix <- function(arr, names) {
  m_mean <- apply(arr, 1:2, mean)
  m_low  <- apply(arr, 1:2, quantile, 0.025)
  m_upp  <- apply(arr, 1:2, quantile, 0.975)
  
  expand.grid(row = 1:nrow(m_mean), col = 1:ncol(m_mean)) %>%
    filter(row > col) %>%
    mutate(
      row_name = factor(names[row], levels = names),
      col_name = factor(names[col], levels = names),
      mean = mapply(\(i, j) m_mean[i, j], row, col),
      low  = mapply(\(i, j) m_low[i, j],  row, col),
      upp  = mapply(\(i, j) m_upp[i, j],  row, col),
      sig  = (low > 0) | (upp < 0),
      label_mean = sprintf("%.2f", mean),
      label_ci   = sprintf("[%.2f, %.2f]", low, upp)
    ) %>%
    filter(sig) %>%
    mutate(
      row_name = droplevels(row_name),
      col_name = droplevels(col_name)
    )
}

make_heatmap <- function(df, title, limits = c(-1, 1), show_legend = FALSE, ci_size = 1.8) {
  ggplot(df, aes(x = col_name, y = row_name, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = label_mean), size = 2.5, nudge_y = 0.12) +
    geom_text(aes(label = label_ci),   size = ci_size, nudge_y = -0.15) +
    scale_fill_gradient2(low = "#C03030", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = limits, name = NULL) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(size = 9),
      plot.title = element_text(size = 11, hjust = 0.5),
      legend.position = if (show_legend) "right" else "none"
    )
}

# ── Precompute distance-zero arrays ──────────────────────────────────────────
iox_fij0 <- spiox:::iox_fij0(iox_all_mcmc, n_threads = 16, matern = 1, tail_mcmc = tail_mcmc)

iox_R0_mcmc <- iox_all_mcmc$Sigma %>% tail(c(NA,NA,tail_mcmc)) %>%
  apply(3, \(x) cov2cor(x)) %>% array(dim = c(q, q, tail_mcmc))
iox_R0_mcmc <- iox_R0_mcmc * iox_fij0

iox_Qpc0_mcmc <- iox_fij0 * (iox_all_mcmc$Sigma %>% tail(c(NA, NA, tail_mcmc)) %>%
                               apply(3, \(x) make_part_corr(solve(x))) %>%
                               array(dim = c(q, q, tail_mcmc)))

df_R   <- summarize_matrix(iox_R0_mcmc,   outcome_names)
df_Qpc <- summarize_matrix(iox_Qpc0_mcmc, outcome_names)

# align factor levels across both
all_rows <- union(levels(df_R$row_name), levels(df_Qpc$row_name))
all_cols <- union(levels(df_R$col_name), levels(df_Qpc$col_name))
df_R$row_name   <- factor(df_R$row_name,   levels = outcome_names[outcome_names %in% all_rows])
df_R$col_name   <- factor(df_R$col_name,   levels = outcome_names[outcome_names %in% all_cols])
df_Qpc$row_name <- factor(df_Qpc$row_name, levels = outcome_names[outcome_names %in% all_rows])
df_Qpc$col_name <- factor(df_Qpc$col_name, levels = outcome_names[outcome_names %in% all_cols])


###############################################################################
#
#  PART A — Traditional geostatistical analysis
#
###############################################################################

row_label <- function(txt) {
  ggplot() + annotate("text", x = 0, y = 0, label = txt,
                      angle = 90, size = 2.8, fontface = "plain") +
    theme_void() + coord_cartesian(clip = "off")
}
lw <- 0.04

# ── A1: Marginal covariance parameters (boxplots) ────────────────────────────
df_theta <- expand.grid(param = c(1,3), outcome = 1:q, iter = 1:mcmc) %>%
  mutate(value = mapply(\(p, o, m) iox_all_mcmc$Theta[p, o, m], param, outcome, iter),
         outcome = factor(outcome_names[outcome], levels = rev(outcome_names))) %>%
  mutate(
    param_label = case_when(
      param == 1 ~ "phi~(spatial~decay)",
      param == 3 ~ "log~SNR"
    ),
    param_label = factor(param_label, levels = c("phi~(spatial~decay)",
                                                 "log~SNR")),
    value = ifelse(param == 3, log((1 - value) / value), value)
  )

guide_df <- data.frame(outcome = factor(rev(outcome_names), levels = rev(outcome_names)))

figA1 <- ggplot(df_theta, aes(x = value, y = outcome)) +
  geom_hline(data = guide_df, aes(yintercept = outcome),
             linewidth = 0.2, linetype = "dotted", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, fill = "#D6E8F0", colour = "grey30",
               linewidth = 0.3, width = 0.7, fatten = 1.5) +
  facet_wrap(~ param_label, scales = "free_x", nrow = 1,
             labeller = label_parsed) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid       = element_blank(),
    panel.spacing.x  = unit(0, "pt"),
    panel.border     = element_rect(colour = "grey50", fill = NA, linewidth = 0.4),
    strip.text       = element_text(size = 9),
    axis.text.x      = element_text(size = 7),
    axis.text.y      = element_text(size = 8),
    axis.ticks.x     = element_line(linewidth = 0.25, colour = "grey50"),
    plot.margin      = margin(5, 5, 5, 5, "pt")
  )

# ── A2: Cross-correlation curves + diagonals ─────────────────────────────────
pick    <- which(outcome_names %in% c("Cd", "Co", "Ni", "Zn", "Cr"))
i_focus <- 2        # Co within pick
j_effr  <- 4        # Zn within pick
i_full  <- pick[i_focus]
j_full  <- pick[j_effr]
others  <- pick[-i_focus]
n_panels <- length(others)

row_marg <- list()
for (k in seq_along(others)) {
  row_marg[[k]] <- make_est_panel(i_full, others[k], "marginal",
                                  iox_pre, gpgpm_fit, outcome_names,
                                  show_x = TRUE, show_y = (k == 1))
}

diag_panels <- list()
pick_diag <- which(outcome_names %in% c("Co", "Cd", "Cr", "Ni"))
for (k in seq_along(pick_diag)) {
  diag_panels[[k]] <- make_diag_panel(pick_diag[k], iox_pre, gpgpm_fit, outcome_names,
                                      show_x = FALSE, show_y = (k == 1))
}

p_er_marg <- make_effrange_panel(i_full, j_full, "marginal",
                                 iox_pre, gpgpm_fit, outcome_names,
                                 show_x = TRUE, show_y = TRUE)

r_marg <- wrap_plots(c(list(row_label("Cross-correlation")), row_marg),
                     nrow = 1, widths = c(lw, rep(1, n_panels)))
r_diag <- wrap_plots(c(list(row_label("Marginal correlation")), diag_panels),
                     nrow = 1, widths = c(lw, rep(1, length(pick_diag))))

figA2 <- r_diag / r_marg  / make_legend() +
  plot_layout(heights = c(1, 1, 0.06))
figA2 <- wrap_elements(patchwork::patchworkGrob(figA2)) +
  theme(plot.margin = margin(-3, 0, -10, 0, "pt"))

# ── A3: Cross-correlation heatmap at distance 0 ──────────────────────────────
figA3 <- make_heatmap(df_R, "Cross-correlation at distance 0")

# ── A4: Effective range (marginal) ───────────────────────────────────────────
figA4 <- p_er_marg

# ── Save Part A ──────────────────────────────────────────────────────────────
ggsave("jura/Jura_A1_theta.pdf",    figA1, width = 4, height = 1.8)
ggsave("jura/Jura_A2_crosscorr.pdf", figA2, width = 6, height = 3.5)
ggsave("jura/Jura_A3_corrmat0.pdf",  figA3, width = 4, height = 3.5)
ggsave("jura/Jura_A4_effrange.pdf",  figA4, width = 3, height = 2)


###############################################################################
#
#  PART B — Partial correlation analysis (new contributions)
#
###############################################################################

# ── B1: Partial cross-correlation curves ─────────────────────────────────────
row_part <- list()
for (k in seq_along(others)) {
  row_part[[k]] <- make_est_panel(i_full, others[k], "partial",
                                  iox_pre, gpgpm_fit, outcome_names,
                                  show_x = TRUE, show_y = (k == 1))
}

p_er_part <- make_effrange_panel(i_full, j_full, "partial",
                                 iox_pre, gpgpm_fit, outcome_names,
                                 show_x = TRUE, show_y = TRUE)

r_part <- wrap_plots(c(list(row_label("Partial cross-correlation")), row_part),
                     nrow = 1, widths = c(lw, rep(1, n_panels)))

figB1 <- r_part / make_legend() +
  plot_layout(heights = c(1, 0.08))
figB1 <- wrap_elements(patchwork::patchworkGrob(figB1)) +
  theme(plot.margin = margin(-3, 0, -10, 0, "pt"))

# ── B2: Partial cross-correlation heatmap at distance 0 ─────────────────────
figB2 <- make_heatmap(df_Qpc, "Partial cross-correlation at distance 0")

# ── B3: Graph ────────────────────────────────────────────────────────────────
Qpc_mcmc_arr <- iox_all_mcmc$Sigma %>%
  apply(3, \(x) -cov2cor(solve(x))) %>% array(dim = c(q, q, mcmc))

Qpc_mean <- apply(Qpc_mcmc_arr, 1:2, mean)
Qpc_low  <- apply(Qpc_mcmc_arr, 1:2, quantile, 0.025)
Qpc_upp  <- apply(Qpc_mcmc_arr, 1:2, quantile, 0.975)

sig <- (Qpc_low > 0) | (Qpc_upp < 0)
diag(sig) <- FALSE
Qpc_graph <- Qpc_mean * sig

g <- graph.empty(n = q, directed = FALSE)
V(g)$name <- outcome_names

edge_data <- list()
for (i in 1:(q - 1)) {
  for (j in (i + 1):q) {
    if (abs(Qpc_graph[i, j]) > 1e-10) {
      g <- add_edges(g, c(i, j))
      edge_data[[length(edge_data) + 1]] <- data.frame(
        from = i, to = j,
        r  = Qpc_mean[i, j],
        lo = Qpc_low[i, j],
        hi = Qpc_upp[i, j]
      )
    }
  }
}
edf <- do.call(rbind, edge_data)

lay <- matrix(
  c( 1.3, 0.2,
     1.3, 1.8,
     0, 1.45,
     -1.3, 0.2,
     0.50, 1,
     -1.3, 1.8,
     -0.50, 0.75),
  ncol = 2, byrow = TRUE,
  dimnames = list(outcome_names, c("x", "y"))
)

node_df <- data.frame(id = 1:q, label = outcome_names,
                      x = lay[, 1], y = lay[, 2])

edge_plot <- data.frame(
  x    = node_df$x[edf$from], y    = node_df$y[edf$from],
  xend = node_df$x[edf$to],   yend = node_df$y[edf$to],
  r = edf$r, lo = edf$lo, hi = edf$hi,
  sign = ifelse(edf$r > 0, "positive", "negative")
)
edge_plot$label <- with(edge_plot,
                        sprintf("%+.2f\n[%+.2f, %+.2f]", r, lo, hi))

edge_colors <- c("positive" = "#2166AC", "negative" = "#C03030")

figB3 <- ggplot() +
  geom_segment(data = edge_plot,
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = sign, linewidth = abs(r)),
               lineend = "round") +
  scale_linewidth_continuous(range = c(0.4, 2.5), guide = "none") +
  geom_label(data = edge_plot,
             aes(x = (x + xend) / 2, y = (y + yend) / 2,
                 label = label, colour = sign),
             size = 2.2, fill = "white", label.size = 0,
             label.padding = unit(2, "pt"), lineheight = 0.85) +
  geom_point(data = node_df, aes(x = x, y = y),
             shape = 21, size = 14, fill = "white", colour = "grey30",
             stroke = 0.8) +
  geom_text(data = node_df, aes(x = x, y = y, label = label),
            size = 3.5) +
  scale_colour_manual(values = edge_colors, guide = "none") +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))

# ── B4: Effective range (partial) ────────────────────────────────────────────
figB4 <- p_er_part

# ── Save Part B ──────────────────────────────────────────────────────────────
ggsave("jura/Jura_B1_partcrosscorr.pdf", figB1, width = 6, height = 2.2)
ggsave("jura/Jura_B2_partcorrmat0.pdf",  figB2, width = 4, height = 3.5)
ggsave("jura/Jura_B3_graph.pdf",         figB3, width = 4.5, height = 3)
ggsave("jura/Jura_B4_effrange_part.pdf",  figB4, width = 3, height = 2)
