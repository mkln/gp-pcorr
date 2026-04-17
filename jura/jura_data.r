library(spiox)
library(gstat)
library(tidyverse)
library(meshed)
library(scico)

data(jura)


jura.all <- rbind(jura.pred, jura.val)
n_all <- nrow(jura.all)




coords <- jura.all %>% dplyr::select(Xloc, Yloc) %>%
  apply(2, \(x) (x - min(x)) / (max(x) - min(x))) %>% as.matrix()

X <- matrix(1, ncol = 1, nrow = n_all)
Y_all <- jura.all %>% dplyr::select(Cd, Co, Cr, Cu, Ni, Pb, Zn) %>% as.matrix() %>% log()
Y_all <- Y_all %>% apply(2, \(y) y - mean(y))
q <- ncol(Y_all)

# place 200 NAs uniformly at random
set.seed(2026)
n_na <- 200
all_cells <- expand.grid(row = 1:n_all, col = 1:q)
na_idx <- all_cells[sample(nrow(all_cells), n_na), ] %>% as.matrix()

Y_m <- Y_all
Y_m[na_idx] <- NA

df <- jura.all %>% dplyr::select(long, lat) %>% bind_cols(Y_m)

df_long <- df %>% 
  pivot_longer(cols=-c(long, lat), names_to = "Metal", values_to="Concentration")

pointsize <- 1.5
(jura_plot <- ggplot(df_long %>% filter(!is.na(Concentration)), aes(long, lat)) +
  geom_point(aes(color = Concentration), size=pointsize) +
  scale_color_scico(palette = "roma", direction=-1) +
  geom_point(data = df_long %>% filter(is.na(Concentration)),
             aes(shape = "Test set"), color = "red", size=pointsize) +
  facet_wrap(~Metal, nrow=2) +
  theme_minimal() +
  labs(color="log(Concentration)", x="Longitude", y="Latitude") +
  guides(color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)) +
  scale_shape_manual(name = NULL, values = c("Test set" = 16)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.87, 0.2),
    legend.justification = c(0.5, 0.5)
  ))

ggsave(filename = "jura/jura_plot.pdf", plot=jura_plot, width=10, height=4.7)


mcmc <- 5000

################################################################################
# IOX via spiox ~~ fixing theta -- not in the paper
################################################################################
set.seed(2026)
marginal_fit_time <- system.time({
  marginal_fix_Theta <- autostart(Y_m, X[, 1, drop = F], coords, m = 30, method = "response", nu = 0)
})["elapsed"]

set.seed(2026)
spiox_time_mcmc <- system.time({
  iox_fix_mcmc <- spiox(Y = Y_m, X = X[, 1, drop = F], coords = coords, m = 30,
                         method = "response", fit = "mcmc", iter = mcmc,
                         starting = marginal_fix_Theta,
                         opts = list(update_Theta = c(0, 0, 0), nu = 0),
                         print_every = 100)
})["elapsed"]
iox_fix_mcmc$time_elapsed <- spiox_time_mcmc + marginal_fit_time
save(iox_fix_mcmc, file = "jura/spiox_jura_fit.RData")


################################################################################
# IOX full MCMC via spiox
################################################################################
set.seed(2026)
spiox_time_mcmc <- system.time({
  iox_all_mcmc <- spiox(Y = Y_m, X = X[, 1, drop = F], coords = coords, m = 30,
                        method = "response", fit = "mcmc", iter = mcmc,
                        starting = list(Theta = rbind(rep(40, q), marginal_fix_Theta$Theta[2,], 1e-2)),
                        opts = list(update_Theta = c(1, 0, 1), nu = 0),
                        print_every = 100)
})["elapsed"]
iox_all_mcmc$time_elapsed <- spiox_time_mcmc + marginal_fit_time
save(iox_all_mcmc, file = "jura/spioxmcmc_jura_fit.RData")

load("jura/spioxmcmc_jura_fit.RData")
Yhat <- Y_m
Yhat[is.na(Y_m)] <- iox_all_mcmc$Y_missing_samples %>% apply(1, mean)

# make predictions at grid
coords_grid <- expand.grid(xg <- seq(0,1,length.out=10), xg) %>% as.matrix(); nout <- nrow(coords_grid)
X_grid <- matrix(1, nrow=nout, ncol=1)

Theta4 <- array(1, dim=c(4, q, dim(iox_all_mcmc$Sigma)[3]))
Theta4[1,,] <- iox_all_mcmc$Theta[1,,]
Theta4[3:4,,] <- iox_all_mcmc$Theta[2:3,,]
tail_mcmc <- 4501:5000
preds <- spiox:::spiox_predict(X_grid, coords_grid, Yhat, X[,1,drop=F], coords, dag = spiox:::dag_vecchia_predict(coords, coords_grid, m=30), B = iox_all_mcmc$Beta[,,tail_mcmc,drop=F], 
                               Sigma = iox_all_mcmc$Sigma[,,tail_mcmc], theta = Theta4[,,tail_mcmc], matern=1, num_threads=16)

Y_grid <- preds$Y %>% apply(1:2, mean)
colnames(Y_grid) <- colnames(Y_m)
df_grid <- data.frame(coords_grid) %>% bind_cols(Y_grid) %>% 
  pivot_longer(cols=-c(Var1, Var2), names_to = "Metal", values_to="Concentration")


(jura_pred_plot <- ggplot(df_grid, aes(Var1, Var2)) +
    geom_raster(aes(fill = Concentration)) +
    scale_fill_scico(palette = "roma", direction=-1) +
    facet_wrap(~Metal, nrow=2) +
    theme_minimal() +
    labs(fill="log(Concentration)", x="Longitude", y="Latitude") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "inside",
      legend.position.inside = c(0.87, 0.2),
      legend.justification = c(0.5, 0.5)
    ))



################################################################################
# Parsimonious Matern via GpGpm
################################################################################
source("fit_parsimonious_matern.R")
set.seed(2026)
gpgpm_fit <- parsimonious_matern(Y_m, X, coords, m = 60,
                                 filename = "jura/GpGpm_jura_fit.RData")

# build prediction inputs for NA locations
n_test <- nrow(na_idx)
Xo_list <- list()
locs_list <- list()
for (k in 1:n_test) {
  i <- na_idx[k, "row"]
  j <- na_idx[k, "col"]
  Xo_list[[k]] <- c(1, j)
  locs_list[[k]] <- c(as.numeric(coords[i, ]), j)
}
Xo_pmatern <- do.call(rbind, Xo_list) %>% as.data.frame()
colnames(Xo_pmatern) <- c("Intercept", "j")
Xo_pmatern$j <- factor(Xo_pmatern$j)
Xo_pmatern <- model.matrix(data = Xo_pmatern, ~ . - 1)[, -(q + 1)]
colnames(Xo_pmatern)[1] <- "Intercept"
locs_test <- do.call(rbind, locs_list)

set.seed(2026)
pmatern_sim_time <- system.time({
  pmatern_samps <- GpGpm::cond_sim(
    fit = gpgpm_fit$fit,
    locs_pred = as.matrix(locs_test),
    X_pred = Xo_pmatern,
    nsims = mcmc / 4
  )
})["elapsed"]
pmatern_time_elapsed <- gpgpm_fit$fit$time_elapsed + pmatern_sim_time

save(pmatern_samps, pmatern_time_elapsed, na_idx,
     file = "jura/GpGpm_jura_preds.RData")

################################################################################
# LMC via meshed
################################################################################
set.seed(2026)
lmc_time_mcmc <- system.time({
  lmc_meshed <- spmeshed(y = Y_m, x = X,
                         coords = coords, k = q, block_size = 60,
                         n_samples = mcmc, n_threads = 16, verbose = 10,
                         prior = list(phi = c(0.3, 200)))
})
lmc_meshed$v_mcmc <- NULL
lmc_meshed$w_mcmc <- NULL
lmc_meshed$lp_mcmc <- NULL
lmc_meshed$time_elapsed <- lmc_time_mcmc["elapsed"]
save(lmc_meshed, file = "jura/lmc_jura_fit.RData")