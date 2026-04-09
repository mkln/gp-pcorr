rm(list=ls())
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)

# ── Colors ────────────────────────────────────────────────────────────────────


bg_partial  <- "#F9F6EE"   # pale cream
bg_marginal <- "#EEF0F3"   # pale slate
bg_diag     <- "#F3F2EF"   # pale stone

lc_zero     <- "#333333"
lc_gpgpm    <- "#EE7733"
lc_iox      <- "#2166AC"

lc_partial  <- "#996515"   # raw ochre
lc_marginal <- "#4A6274"   # steel blue-grey

source("figures/make_sparse_Q_5.r")
source("functions.R")
load("sim_pmatern/GpGpm_pmatern_fit5.RData")
load("sim_pmatern/spiox_pmatern_fit5.RData")
mcmc <- 5000
q <- 5
nu <- c(0.2, 1, 0.5, 1.4, 0.75); phi <- 10 # for matern

edges <- list(
  list(from = 3, to = 1, weight = 1.1),   # confounder
  list(from = 3, to = 2, weight = 1.1),
  list(from = 1, to = 2, weight = -1.2),   # direct negative
  list(from = 5, to = 1, weight = 0.1)
)

Q_structured <- generate_network_matrices(
  build_B(q, add_chain(edges, c(3, 4, 5), weight = 0.7)),
  noise_var = c(1, 1, 1, 1, 1)  # tweakable!
)

Q <- Q_structured$Precision
Qpc <- Q_structured$Partial_Correlation
Sigma <- Q_structured$Covariance
Rc <- Q_structured$Marginal_Correlation


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

# ── Precompute IOX ────────────────────────────────────────────────────────────
precompute_iox <- function(i, iox_full_mcmc, q, mcmc, tail_mcmc, n_threads = 16) {
  R_mcmc   <- iox_full_mcmc$Sigma %>%
    apply(3, \(x) cov2cor(x)) %>% array(dim = c(q, q, mcmc))
  Qpc_mcmc <- iox_full_mcmc$Sigma %>%
    apply(3, \(x) -cov2cor(solve(x))) %>% array(dim = c(q, q, mcmc))
  
  out <- list()
  for (j in (1:q)[-i]) {
    fij      <- spiox:::iox_fij(iox_full_mcmc, i, j, matern = 1,
                                n_threads = n_threads, tail_mcmc = tail_mcmc)
    R_tail   <- tail(R_mcmc[i, j, ],   tail_mcmc)
    Qpc_tail <- tail(Qpc_mcmc[i, j, ], tail_mcmc)
    C_mat    <- t(t(fij$f_ij) * as.vector(R_tail))
    Qpc_mat  <- t(t(fij$f_ij) * as.vector(Qpc_tail))
    
    out[[j]] <- list(
      dists    = fij$dists,
      C_mean   = rowMeans(C_mat),
      C_low    = apply(C_mat, 1, quantile, 0.025),
      C_upp    = apply(C_mat, 1, quantile, 0.975),
      Qpc_mean = rowMeans(Qpc_mat),
      Qpc_low  = apply(Qpc_mat, 1, quantile, 0.025),
      Qpc_upp  = apply(Qpc_mat, 1, quantile, 0.975)
    )
  }
  out
}

# ── Cross-correlation panel (rows 1 & 2) ─────────────────────────────────────
make_est_panel <- function(i, j, type = c("partial", "marginal"),
                           iox_pre, gpgpm_fit,
                           Q_true, Sigma_true, nu_true, phi_true,
                           show_x = FALSE, show_y = FALSE) {
  type <- match.arg(type)
  hh   <- iox_pre[[j]]$dists
  dd   <- 2
  
  g_true <- gamma_coef(nu_true[i], nu_true[j], dd)
  M_true <- matern_cor(hh, (nu_true[i] + nu_true[j]) / 2, phi_true)
  
  if (type == "partial") {
    r_true   <- -Q_true[i, j] / sqrt(Q_true[i, i] * Q_true[j, j])
    y_true   <- r_true * g_true * M_true
    lc_truth <- lc_partial
    bg       <- bg_partial
    y_iox    <- iox_pre[[j]]$Qpc_mean
    y_low    <- iox_pre[[j]]$Qpc_low
    y_upp    <- iox_pre[[j]]$Qpc_upp
    Qpc_fit  <- -cov2cor(gpgpm_fit$Q); diag(Qpc_fit) <- 1
    g_fit    <- get_gamma_mat(gpgpm_fit$nu, dd)[i, j]
    M_fit    <- matern_cor(hh, (gpgpm_fit$nu[i] + gpgpm_fit$nu[j]) / 2, gpgpm_fit$phi)
    y_gpgpm  <- Qpc_fit[i, j] * g_fit * M_fit
  } else {
    rho_true <- Sigma_true[i, j] / sqrt(Sigma_true[i, i] * Sigma_true[j, j])
    y_true   <- rho_true * g_true * M_true
    lc_truth <- lc_marginal
    bg       <- bg_marginal
    y_iox    <- iox_pre[[j]]$C_mean
    y_low    <- iox_pre[[j]]$C_low
    y_upp    <- iox_pre[[j]]$C_upp
    R_fit    <- cov2cor(gpgpm_fit$Sigma)
    g_fit    <- get_gamma_mat(gpgpm_fit$nu, dd)[i, j]
    M_fit    <- matern_cor(hh, (gpgpm_fit$nu[i] + gpgpm_fit$nu[j]) / 2, gpgpm_fit$phi)
    y_gpgpm  <- R_fit[i, j] * g_fit * M_fit
  }
  
  df <- data.frame(h = hh, truth = y_true, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  ggplot(df, aes(x = h)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp), fill = lc_iox, alpha = 0.15) +
    geom_line(aes(y = truth), linewidth = 0.55, colour = lc_truth, linetype = "dotted") +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox) +
    coord_cartesian(ylim = c(-1, 1), xlim = c(0, 0.5)) +
    labs(title = bquote(y[.(i)] * ", " * y[.(j)]), x = NULL, y = NULL) +
    theme_panel(bg, show_x = show_x, show_y = show_y)
}

# ── Marginal correlation panel (row 3, left) ─────────────────────────────────
make_diag_panel <- function(idx, iox_full_mcmc, gpgpm_fit,
                            nu_true, phi_true,
                            tail_mcmc = 500, n_threads = 16,
                            show_x = TRUE, show_y = TRUE) {
  fii     <- spiox:::iox_fij(iox_full_mcmc, idx, idx, matern = 1,
                             n_threads = n_threads, tail_mcmc = tail_mcmc)
  hh      <- fii$dists
  y_iox   <- rowMeans(fii$f_ij)
  y_low   <- apply(fii$f_ij, 1, quantile, 0.025)
  y_upp   <- apply(fii$f_ij, 1, quantile, 0.975)
  y_true  <- matern_cor(hh, nu_true[idx], phi_true)
  y_gpgpm <- matern_cor(hh, gpgpm_fit$nu[idx], gpgpm_fit$phi)
  
  df <- data.frame(h = hh, truth = y_true, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  ggplot(df, aes(x = h)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp), fill = lc_iox, alpha = 0.15) +
    geom_line(aes(y = truth), linewidth = 0.55, colour = lc_zero, linetype = "dotted") +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.5)) +
    labs(title = bquote(y[.(idx)]), x = NULL, y = NULL) +
    theme_panel(bg_diag, show_x = show_x, show_y = show_y)
}

# ── Effective range panel (row 3, right) ──────────────────────────────────────
make_effrange_panel <- function(i_pick, j_pick, type = c("partial", "marginal"),
                                iox_pre, gpgpm_fit,
                                Q_true, Sigma_true, nu_true, phi_true,
                                h = seq(1e-4, 1, length.out = 500),
                                thr = 0.05,
                                show_x = TRUE, show_y = TRUE) {
  type <- match.arg(type)
  dd   <- 2
  
  lc_gpgpm_use <- lc_gpgpm
  lc_iox_use   <- lc_iox
  lc_truth_use <- "grey30"
  
  g_true <- gamma_coef(nu_true[i_pick], nu_true[j_pick], dd)
  M_true <- matern_cor(h, (nu_true[i_pick] + nu_true[j_pick]) / 2, phi_true)
  
  g_fit <- get_gamma_mat(gpgpm_fit$nu, dd)[i_pick, j_pick]
  M_fit <- matern_cor(h, (gpgpm_fit$nu[i_pick] + gpgpm_fit$nu[j_pick]) / 2, gpgpm_fit$phi)
  
  iox_j <- iox_pre[[j_pick]]
  
  if (type == "partial") {
    r_true  <- -Q_true[i_pick, j_pick] / sqrt(Q_true[i_pick, i_pick] * Q_true[j_pick, j_pick])
    y_true  <- r_true * g_true * M_true
    Qpc_fit <- -cov2cor(gpgpm_fit$Q); diag(Qpc_fit) <- 1
    y_gpgpm <- Qpc_fit[i_pick, j_pick] * g_fit * M_fit
    y_iox   <- approx(iox_j$dists, iox_j$Qpc_mean, xout = h, rule = 2)$y
    y_low   <- approx(iox_j$dists, iox_j$Qpc_low,  xout = h, rule = 2)$y
    y_upp   <- approx(iox_j$dists, iox_j$Qpc_upp,  xout = h, rule = 2)$y
    bg      <- bg_partial
    lc_truth_panel <- lc_partial
    ttl     <- bquote("Partial cross-correl.:" ~ y[.(i_pick)] * "," ~ y[.(j_pick)])
  } else {
    rho_true <- Sigma_true[i_pick, j_pick] / sqrt(Sigma_true[i_pick, i_pick] * Sigma_true[j_pick, j_pick])
    y_true   <- rho_true * g_true * M_true
    R_fit    <- cov2cor(gpgpm_fit$Sigma)
    y_gpgpm  <- R_fit[i_pick, j_pick] * g_fit * M_fit
    y_iox    <- approx(iox_j$dists, iox_j$C_mean, xout = h, rule = 2)$y
    y_low    <- approx(iox_j$dists, iox_j$C_low,  xout = h, rule = 2)$y
    y_upp    <- approx(iox_j$dists, iox_j$C_upp,  xout = h, rule = 2)$y
    bg       <- bg_marginal
    lc_truth_panel <- lc_marginal
    ttl      <- bquote("Cross-correl.:" ~ y[.(i_pick)] * "," ~ y[.(j_pick)])
  }
  
  find_eff_range <- function(hh, yy, thr) {
    idx <- which(abs(yy) < thr)
    if (length(idx) == 0) return(NA)
    hh[min(idx)]
  }
  
  sgn    <- sign(y_true[1])
  er_t   <- find_eff_range(h, y_true,  thr)
  er_g   <- find_eff_range(h, y_gpgpm, thr)
  er_i   <- find_eff_range(h, y_iox,   thr)
  
  df <- data.frame(h = h, truth = y_true, gpgpm = y_gpgpm,
                   iox = y_iox, iox_low = y_low, iox_upp = y_upp)
  
  p <- ggplot(df, aes(x = h)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -thr, ymax = thr),
              fill = "grey90", colour = NA, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = iox_low, ymax = iox_upp),
                fill = lc_iox_use, alpha = 0.12) +
    geom_line(aes(y = truth), linewidth = 0.55, colour = lc_truth_panel, linetype = "dotted") +
    geom_line(aes(y = gpgpm), linewidth = 0.5, colour = lc_gpgpm_use) +
    geom_line(aes(y = iox),   linewidth = 0.5, colour = lc_iox_use)
  
  # truth effective range
  if (!is.na(er_t)) {
    p <- p + geom_vline(xintercept = er_t, linewidth = 0.35,
                        linetype = "dotted", colour = lc_truth_use)
  }
  # estimated effective ranges
  yend <- sgn * thr

  # find crossing points
  cross_pts <- data.frame(
    h = c(er_t, er_g, er_i),
    y = rep(sgn * thr, 3),
    col = c(lc_truth_panel, lc_gpgpm_use, lc_iox_use),
    shp = c(NA, 16, 16)  # open circle for truth, filled for estimates
  )
  cross_pts <- cross_pts[!is.na(cross_pts$h), ]

  p <- p +
    geom_point(data = cross_pts, aes(x = h, y = y),
               colour = cross_pts$col, fill = cross_pts$col,
               shape = cross_pts$shp, size = 2.5, inherit.aes = FALSE)
  
  p <- p +
    coord_cartesian(ylim = c(0,1), xlim = c(0, 0.5)) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_panel(bg, show_x = show_x, show_y = show_y)
  
  p
}

# ── Unified figure ────────────────────────────────────────────────────────────
figure_unified <- function(i, j_pick,
                           iox_pre, iox_full_mcmc, gpgpm_fit,
                           Q_true, Sigma_true, nu_true, phi_true,
                           q, mcmc, tail_mcmc = 500, n_threads = 16) {
  
  others   <- (1:q)[-i]
  n_panels <- length(others)
  
  # ── Row 1: partial cross-correlations ──
  row1_panels <- list()
  for (k in seq_along(others)) {
    j <- others[k]
    row1_panels[[k]] <- make_est_panel(i, j, "partial", iox_pre, gpgpm_fit,
                                       Q_true, Sigma_true, nu_true, phi_true,
                                       show_x = FALSE, show_y = (k == 1))
  }
  
  # ── Row 2: marginal cross-correlations ──
  row2_panels <- list()
  for (k in seq_along(others)) {
    j <- others[k]
    row2_panels[[k]] <- make_est_panel(i, j, "marginal", iox_pre, gpgpm_fit,
                                       Q_true, Sigma_true, nu_true, phi_true,
                                       show_x = FALSE, show_y = (k == 1))
  }
  
  # ── Row 3: 2 marginals + 2 effective range panels ──
  p_diag1 <- make_diag_panel(i, iox_full_mcmc, gpgpm_fit,
                             nu_true, phi_true, tail_mcmc, n_threads,
                             show_x = TRUE, show_y = TRUE)
  p_diag2 <- make_diag_panel(j_pick, iox_full_mcmc, gpgpm_fit,
                             nu_true, phi_true, tail_mcmc, n_threads,
                             show_x = TRUE, show_y = FALSE)
  p_er1   <- make_effrange_panel(i, j_pick, "marginal", iox_pre, gpgpm_fit,
                                 Q_true, Sigma_true, nu_true, phi_true,
                                 show_x = TRUE, show_y = FALSE)
  p_er2   <- make_effrange_panel(i, j_pick, "partial", iox_pre, gpgpm_fit,
                                 Q_true, Sigma_true, nu_true, phi_true,
                                 show_x = TRUE, show_y = FALSE)
  
  # ── Row labels ──
  row_label <- function(txt) {
    ggplot() +
      annotate("text", x = 0, y = 0, label = txt,
               angle = 90, size = 2.8, fontface = "plain") +
      theme_void() +
      coord_cartesian(clip = "off")
  }
  
  lw <- 0.04  # label column width
  
  row1 <- wrap_plots(c(
    list(row_label("Partial cross-correlation")),
    row1_panels
  ), nrow = 1, widths = c(lw, rep(1, n_panels)))
  
  row2 <- wrap_plots(c(
    list(row_label("Cross-correlation")),
    row2_panels
  ), nrow = 1, widths = c(lw, rep(1, n_panels)))
  
  row3 <- wrap_plots(c(
    list(row_label("Marginal & Eff. range")),
    list(p_diag1, p_diag2, p_er1, p_er2)
  ), nrow = 1, widths = c(lw, 1, 1, 1, 1))
  
  # ── Legend ──
  leg_df <- data.frame(
    x = rep(1:3, 3), y = rep(1:3, 3),
    method = factor(rep(c("Truth", "Estimated (Pars. Matern)", "Estimated (Inside-out cross-cov.)"), each = 3),
                    levels = c("Truth", "Estimated (Pars. Matern)", "Estimated (Inside-out cross-cov.)"))
  )
  
  legend_plot <- ggplot(leg_df, aes(x, y, colour = method, linetype = method)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("Truth" = "grey30",
                                   "Estimated (Pars. Matern)" = lc_gpgpm,
                                   "Estimated (Inside-out cross-cov.)" = lc_iox)) +
    scale_linetype_manual(values = c("Truth" = "dotted",
                                     "Estimated (Pars. Matern)" = "solid",
                                     "Estimated (Inside-out cross-cov.)" = "solid")) +
    guides(colour   = guide_legend(nrow = 1, title = NULL),
           linetype = guide_legend(nrow = 1, title = NULL)) +
    labs(title = "Distance") +
    theme_void() +
    theme(legend.position   = "top",
          legend.direction  = "horizontal",
          legend.margin     = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(0, 0, 0, 0, "pt"),
          legend.text       = element_text(size = 8),
          legend.key.width  = unit(18, "pt"),
          plot.title        = element_text(size = 9, hjust = 0.5,
                                           margin = margin(0, 0, 0, 0, "pt")),
          plot.margin       = margin(0, 0, 0, 0, "pt")) +
    coord_cartesian(xlim = c(0, 0), ylim = c(0, 0))
  
  # ── Assemble ──
  fig <- row1 / row2 / row3 / legend_plot +
    plot_layout(heights = c(1, 1, 1, 0.06))
  
  fig_grob <- patchwork::patchworkGrob(fig)
  fig_trimmed <- wrap_elements(fig_grob) +
    theme(plot.margin = margin(-3, 0, -10, 0, "pt"))
  
  fig_trimmed
}

# ── Usage ─────────────────────────────────────────────────────────────────────
i <- 3
j_pick <- 4
iox_pre <- precompute_iox(i, iox_full_mcmc, q, mcmc, tail_mcmc = 500, n_threads = 16)

fig <- figure_unified(
  i = i, j_pick = j_pick,
  iox_pre = iox_pre, iox_full_mcmc = iox_full_mcmc,
  gpgpm_fit = gpgpm_fit,
  Q_true = Q, Sigma_true = Sigma, nu_true = nu, phi_true = phi,
  q = 5, mcmc = mcmc, tail_mcmc = 500, n_threads = 16
)
print(fig)
#
ggsave("sim_pmatern/Illustration_2.pdf", fig,
       width = 6, height = 4.5)

