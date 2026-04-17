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

load("jura/spioxmcmc_jura_fit.RData")
Yhat <- Y_m
Yhat[is.na(Y_m)] <- iox_all_mcmc$Y_missing_samples %>% apply(1, mean)

# make predictions at grid
set.seed(2026)
coords_grid <- expand.grid(xg <- seq(0,1,length.out=50), xg) %>% as.matrix(); nout <- nrow(coords_grid)
X_grid <- matrix(1, nrow=nout, ncol=1)

Theta4 <- array(1, dim=c(4, q, dim(iox_all_mcmc$Sigma)[3]))
Theta4[1,,] <- iox_all_mcmc$Theta[1,,]
Theta4[3:4,,] <- iox_all_mcmc$Theta[2:3,,]
tail_mcmc <- (mcmc-499):mcmc
preds <- spiox:::spiox_predict(X_grid, coords_grid, Yhat, X[,1,drop=F], coords, dag = spiox:::dag_vecchia_predict(coords, coords_grid, m=30), B = iox_all_mcmc$Beta[,,tail_mcmc,drop=F], 
                               Sigma = iox_all_mcmc$Sigma[,,tail_mcmc], theta = Theta4[,,tail_mcmc], matern=1, num_threads=16)

Y_grid <- preds$Y %>% apply(1:2, mean)
colnames(Y_grid) <- colnames(Y_m)

coords_grid[,1] <- with(jura.all, (max(long)-min(long))*coords_grid[,1] +min(long))
coords_grid[,2] <- with(jura.all, (max(lat)-min(lat))*coords_grid[,2] +min(lat))

df_grid <- data.frame(coords_grid) %>% bind_cols(Y_grid) %>% 
  pivot_longer(cols=-c(Var1, Var2), names_to = "Metal", values_to="Concentration")

(jura_pred_plot <- ggplot(df_grid, aes(Var1, Var2)) +
    geom_raster(aes(fill = Concentration), interpolate = TRUE) +
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

ggsave(filename = "jura/jura_pred_plot.pdf", plot=jura_pred_plot, width=10, height=4.7)

