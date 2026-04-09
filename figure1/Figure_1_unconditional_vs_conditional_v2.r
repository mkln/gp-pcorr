
library(ggplot2)
library(patchwork)

# setup pars
source("~/spiox-network/functions.R")
source("~/spiox-network/figures/make_sparse_Q_5.r")

edges <- list(
  list(from = 3, to = 1, weight = 1.1),   # confounder
  list(from = 3, to = 2, weight = 1.1),
  list(from = 1, to = 2, weight = -1.2)   # direct negative
)

Q_structured <- generate_network_matrices(
  build_B(q, add_chain(edges, c(3, 4, 5), weight = 0.7)),
  noise_var = c(1, 1, 1, 1, 1)  # tweakable!
)

Q <- Q_structured$Precision
Qpc <- Q_structured$Partial_Correlation
Sigma <- Q_structured$Covariance
Rc <- Q_structured$Marginal_Correlation





################# Model parameters ################ 
q   <- 5
dd  <- 2

nu <- c(0.2, 1, 0.5, 1.4, 0.75); phi <- 10 # for matern
#ell <- c(0.2, 0.3, 0.2, .1, 0.25) # for sqexp

h   <- seq(1e-4, 1, length.out = 500)

################ Precision matrix Q ################ 
Q <- Q_structured$Precision
Sigma <- solve(Q)

# ── Diagnostic ────────────────────────────────────────────────────────────────
cat("\n── Marginal corr at lag 0 ──\n")
marg_cor_mat <- matrix(NA, q, q)
for (i in 1:q) for (j in 1:q) {
  g <- gamma_coef(nu[i], nu[j], dd)
  marg_cor_mat[i, j] <- Sigma[i, j] / sqrt(Sigma[i, i] * Sigma[j, j]) * g
}
print(round(marg_cor_mat, 3))

cat("\n── Partial corr at lag 0 ──\n")
part_cor_mat <- matrix(NA, q, q)
for (i in 1:q) for (j in 1:q) {
  g <- gamma_coef(nu[i], nu[j], dd)
  if (i == j) {
    part_cor_mat[i, j] <- g
  } else {
    part_cor_mat[i, j] <- -Q[i, j] / sqrt(Q[i, i] * Q[j, j]) * g
  }
}
print(round(part_cor_mat, 3))

cat("\n── Key contrasts ──\n")
cat("  Suppression  (1,2): marginal =", round(marg_cor_mat[1, 2], 3),
    "  partial =", round(part_cor_mat[1, 2], 3), "\n")
cat("  Confounding  (4,5): marginal =", round(marg_cor_mat[4, 5], 3),
    "  partial =", round(part_cor_mat[4, 5], 3), "\n")
cat("  Confounding  (2,4): marginal =", round(marg_cor_mat[2, 4], 3),
    "  partial =", round(part_cor_mat[2, 4], 3), "\n")


bg_partial  <- "#F9F6EE"   # pale cream
bg_marginal <- "#EEF0F3"   # pale slate
bg_diag     <- "#F3F2EF"   # pale stone

lc_partial  <- "#996515"   # raw ochre
lc_marginal <- "#4A6274"   # steel blue-grey
lc_diag     <- "#333333"
lc_zero     <- "#333333"

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_panel <- function(bg, show_x = FALSE, show_y = FALSE, border = FALSE) {
  t <- theme(
    panel.background  = element_rect(fill = bg,
                                     colour = if (border) "grey40" else NA,
                                     linewidth = if (border) 0.5 else 0),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(linewidth = 0.3, colour = "grey50"),
    axis.ticks        = element_line(linewidth = 0.25, colour = "grey50"),
    axis.ticks.length = unit(1.2, "pt"),
    plot.title        = element_text(size = 8, hjust = 0.5, face = "plain",
                                     margin = margin(0, 0, 1, 0, "pt")),
    plot.margin       = margin(1, 3, 2, 3, "pt"),
    axis.text.x       = if (show_x) element_text(size = 7) else element_blank(),
    axis.ticks.x      = if (show_x) element_line(linewidth = 0.25, colour = "grey50") else element_blank(),
    axis.text.y       = if (show_y) element_text(size = 7) else element_blank(),
    axis.ticks.y      = if (show_y) element_line(linewidth = 0.25, colour = "grey50") else element_blank()
  )
  t
}

# ── Build panels ──────────────────────────────────────────────────────────────
make_panel <- function(row, col) {
  
  i <- row; j <- col
  g_ij  <- gamma_coef(nu_i = nu[i], nu[j]) #gamma_coef_se(ell[i], ell[j], dd)
  nu_ij <- (nu[i] + nu[j]) / 2
  
  if (i == j) {
    y       <- matern_cor(h, nu[i], phi)  #sq_exp_conv(h, ell[i], ell[i])
    region  <- "diag"
    bg      <- bg_diag
    lc      <- lc_diag
    ttl     <- bquote(y[.(i)])
    is_zero <- FALSE
  } else if (i < j) {
    r_ij    <- -Q[i, j] / sqrt(Q[i, i] * Q[j, j])
    y       <- r_ij * g_ij * matern_cor(h, (nu[i]+nu[j])/2, phi) #sq_exp_conv(h, ell[i], ell[j])
    region  <- "partial"
    bg      <- bg_partial
    is_zero <- abs(Q[i, j]) < 1e-10
    lc      <- if (is_zero) lc_zero else lc_partial
    ttl     <- bquote(y[.(i)] * ", " * y[.(j)])
  } else {
    rho_ij  <- Sigma[i, j] / sqrt(Sigma[i, i] * Sigma[j, j])
    y       <- rho_ij * g_ij * matern_cor(h, (nu[i]+nu[j])/2, phi) #sq_exp_conv(h, ell[i], ell[j])
    region  <- "marginal"
    bg      <- bg_marginal
    is_zero <- FALSE
    lc      <- lc_marginal
    ttl     <- bquote(y[.(i)] * ", " * y[.(j)])
  }
  
  lt <- if (is_zero) "dashed" else "solid"
  df <- data.frame(h = h, y = y)
  
  show_x  <- (row == q)
  show_y  <- (col == 1)
  border  <- (row == col)
  
  p <- ggplot(df, aes(x = h, y = y)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed",
               colour = "grey70") +
    geom_line(linewidth = 0.6, colour = lc, linetype = lt) +
    coord_cartesian(ylim = c(-1, 1), xlim=c(0, 0.5)) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_panel(bg, show_x = show_x, show_y = show_y, border = border)
  
  if (is_zero) {
    p <- p + annotate("text", x = max(h) * 0.25, y = 0.2,
                      label = expression(Q[ij] == 0),
                      size = 2.5, colour = lc_zero, fontface = "italic")
  }
  
  if ((row == q) & (col == 3)) {
    p <- p + xlab(expression(paste("||", bold(h), "||")))
  }
  if ((col == 1) & (row == 3)) {
    p <- p + ylab("Correlation")
  }
  
  p
}
# ── Assemble ──────────────────────────────────────────────────────────────────
plot_list <- list()
for (row in 1:q) {
  for (col in 1:q) {
    key <- paste0("p", row, col)
    plot_list[[key]] <- make_panel(row, col)
  }
}

combined <- wrap_plots(plot_list, ncol = q, nrow = q)

print(combined)

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("~/spiox-network/figures/figure_partial.pdf", combined,
       width = 6, height = 5.5, units = "in")




# ──────────────────────────────────────────────────────────────────────────────
# Side-by-side: marginal vs partial spatial correlation for a chosen pair
# with ±0.05 band marking effective range / partial effective range
#
# Usage: set i_pick and j_pick below (i_pick < j_pick)
# Requires: run figure_partial_matern.R first (or at least its setup)
# ──────────────────────────────────────────────────────────────────────────────

# ── Pick pair here ────────────────────────────────────────────────────────────
i_pick <- 3
j_pick <- 4

# ── Compute curves ────────────────────────────────────────────────────────────
g_ij  <- #gamma_coef_se(ell[i_pick], ell[j_pick], dd)#
  gamma_coef(nu[i_pick], nu[j_pick], dd)
nu_ij <- (nu[i_pick] + nu[j_pick]) / 2
M_h   <- #sq_exp_conv(h, ell[i], ell[j]) #
  matern_cor(h, nu_ij, phi)

# Marginal spatial correlation
rho_ij <- Sigma[i_pick, j_pick] / sqrt(Sigma[i_pick, i_pick] * Sigma[j_pick, j_pick])
y_marg <- rho_ij * g_ij * M_h

# Partial spatial correlation
r_ij   <- -Q[i_pick, j_pick] / sqrt(Q[i_pick, i_pick] * Q[j_pick, j_pick])
y_part <- r_ij * g_ij * M_h

df <- data.frame(
  h        = rep(h, 2),
  y        = c(y_marg, y_part),
  type     = rep(c("Cross-correlation", "Partial cross-correlation"), each = length(h))
)
df$type <- factor(df$type, levels = c("Cross-correlation", "Partial cross-correlation"))

# ── Effective range: first h where |y| drops below 0.05 ──────────────────────
find_eff_range <- function(hh, yy, thr = 0.05) {
  idx <- which(abs(yy) < thr)
  if (length(idx) == 0) return(NA)
  hh[min(idx)]
}

er_marg <- find_eff_range(h, y_marg)
er_part <- find_eff_range(h, y_part)
cat("Effective range (unconditional):", round(er_marg, 2), "\n")
cat("Effective range (partial):", round(er_part, 2), "\n")

# ── Plot ──────────────────────────────────────────────────────────────────────
colors <- c("Cross-correlation" = lc_marginal, "Partial cross-correlation" = lc_partial)
fills  <- c("Cross-correlation" = bg_marginal, "Partial cross-correlation" = bg_partial)
# ── Annotation data: one line per panel ───────────────────────────────────────
vline_df <- data.frame(
  type = factor(c("Cross-correlation", "Partial cross-correlation"), 
                levels = c("Cross-correlation", "Partial cross-correlation")),
  xint = c(er_marg, er_part),
  col  = c(lc_marginal, lc_partial)
)
vline_df <- vline_df[!is.na(vline_df$xint), ]

label_df <- data.frame(
  type  = vline_df$type,
  x     = vline_df$xint + 0.02,
  y     = max(df$y) * 0.37,
  label = paste0("Delta[eff] == ", round(vline_df$xint, 2)),
  col   = vline_df$col
)

# build parseable labels
lab_marg <- bquote("Cross-correlation:" ~ y[.(i_pick)] * "," ~ y[.(j_pick)])
lab_part <- bquote("Partial cross-correlation:" ~ y[.(i_pick)] * "," ~ y[.(j_pick)])

lvls <- c(deparse(lab_marg), deparse(lab_part))

df$type <- factor(
  ifelse(df$type == "Cross-correlation", lvls[1], lvls[2]),
  levels = lvls
)

# update vline_df and label_df to match
vline_df$type <- factor(
  ifelse(vline_df$type == "Cross-correlation", lvls[1], lvls[2]),
  levels = lvls
)
label_df$type <- factor(
  ifelse(label_df$type == "Cross-correlation", lvls[1], lvls[2]),
  levels = lvls
)


pair_plot <- ggplot(df, aes(x = h, y = y)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.05, ymax = 0.05),
            fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
             colour = "grey60") +
  geom_line(aes(colour = type), linewidth = 0.7) +
  geom_vline(data = vline_df, aes(xintercept = xint),
             linewidth = 0.4, linetype = "dotted", colour = vline_df$col) +
  geom_text(data = label_df, aes(x = x, y = y, label = label),
            hjust = -0.15, size = 4, 
            colour = label_df$col,parse = TRUE) +
  coord_cartesian(ylim = c(0, .5), xlim=c(0, 0.75)) +
  scale_colour_manual(values = colors) +
  facet_wrap(~ type, ncol = 2, labeller = label_parsed) +
  labs(
    x = expression(paste("|| ", bold(h), " ||")),
    y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(linewidth = 0.3, colour = "grey50"),
    axis.ticks       = element_line(linewidth = 0.25, colour = "grey50"),
    legend.position  = "none",
    plot.title       = element_text(size = 9, hjust = 0.5)
  )
print(pair_plot)

ggsave("~/spiox-network/figures/figure_pair_comparison.pdf", pair_plot,
       width = 5, height = 2.45, units = "in")



############ GRAPH

# Node positions (pentagon-ish, process 3 central)
node_df <- data.frame(
  id    = 1:q,
  label = paste0("y[", 1:q, "]"),
  x     = c(-1.0,  1.0,  0.0, -0.8,  0.8),
  y     = c( 1.0,  1.0,  0.0, -1.0, -1.0)
)

# Edges from Q
edge_list <- list()
for (i in 1:(q - 1)) {
  for (j in (i + 1):q) {
    if (abs(Q[i, j]) > 1e-10) {
      r_ij <- -Q[i, j] / sqrt(Q[i, i] * Q[j, j])
      edge_list[[length(edge_list) + 1]] <- data.frame(
        x    = node_df$x[i], y    = node_df$y[i],
        xend = node_df$x[j], yend = node_df$y[j],
        r    = r_ij,
        sign = ifelse(r_ij > 0, "positive", "negative")
      )
    }
  }
}
edge_df <- do.call(rbind, edge_list)

edge_colors <- c("positive" = "#2166AC", "negative" = "#C03030")

graph_plot <- ggplot() +
  # Edges
  geom_segment(data = edge_df,
               aes(x = x, y = y, xend = xend, yend = yend, colour = sign),
               linewidth = 1.2, lineend = "round") +
  # Edge weight labels
  geom_label(data = edge_df,
             aes(x = (x + xend) / 2, y = (y + yend) / 2,
                 label = sprintf("%+.2f", r),
                 colour = sign),
             size = 2.5, #family = "serif",
             fill = "white", label.size = 0, label.padding = unit(2, "pt")) +
  # Nodes
  geom_point(data = node_df, aes(x = x, y = y),
             shape = 21, size = 12, fill = "white", colour = "grey30",
             stroke = 0.8) +
  # Node labels
  geom_text(data = node_df, aes(x = x, y = y, label = label),
            size = 3.5, #family = "serif", 
            parse = TRUE) +
  scale_colour_manual(values = edge_colors, guide = "none") +
  coord_fixed(clip = "off") +
  theme_void()+
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))

print(graph_plot)

ggsave("~/spiox-network/figures/figure_graph.pdf", graph_plot,
       width = 2, height = 2, units = "in")


cat("\nDone.\n")