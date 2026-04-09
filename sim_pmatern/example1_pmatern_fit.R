rm(list=ls())
library(spiox)
library(tidyverse)
library(spBayes)
library(magrittr)

set.seed(2026)

source("functions.R")
source("figures/make_sparse_Q_5.r")

df <- read_csv("sim_pmatern/pmatern5_alt.csv")
n_all <- nrow(df)

edges <- list(
  list(from = 3, to = 1, weight = 1.1),   # confounder
  list(from = 3, to = 2, weight = 1.1),
  list(from = 1, to = 2, weight = -1.2),   # direct negative
  list(from = 5, to = 1, weight = 0.1)
)
q <- 5

Q_structured <- generate_network_matrices(
  build_B(q, add_chain(edges, c(3, 4, 5), weight = 0.7)),
  noise_var = c(1, 1, 1, 1, 1)  # tweakable!
)

Q <- Q_structured$Precision
Qpc <- Q_structured$Partial_Correlation
Sigma <- Q_structured$Covariance
Rc <- Q_structured$Marginal_Correlation

# A matrix of LMC
A <- Sigma %>% chol() %>% t() %>% solve()
A[abs(A)<1e-10] <- 0

coords <- df %>% dplyr::select(x, y) %>% as.matrix()
nu <- c(0.2, 1, 0.5, 1.4, 0.75); phi <- 10 # for matern
C_str <- pmatern_cov(coords, Sigma, nu, phi)
#########

X <- matrix(1, ncol=1, nrow=n_all) 
Y <- df %>% dplyr::select(-x, -y) %>% as.matrix() 
Y <- Y %>% apply(2, \(y) y-mean(y))
q <- ncol(Y)

Y_m <- Y
#Y_m[sample(1:(n_all*q), 200, replace=FALSE)] <- NA
#na_mask <- is.na(Y_m)


df <- data.frame(coords, Y) 
df %>% pivot_longer(cols=-c(x,y)) %>%
  ggplot(aes(x, y, fill=value)) +
  geom_raster() +
  facet_grid( ~ name) +
  scale_fill_viridis_c() +
  theme_minimal() 


################################################################################
# IOX Response via spiox
################################################################################
q <- ncol(Y)
set.seed(2026)
mcmc <- 5000
spiox_time_mcmc <- system.time({ iox_full_mcmc <- spiox(Y = Y_m, X = X[,1,drop=F], coords = coords, m = 30, 
                                     method = "response", fit = "mcmc", iter = mcmc, 
                                     starting=list(Theta=rbind(10, nu, 1e-3)),
                                     opts=list(update_Theta=c(1,1,1), nu=0),
                                     print_every = 100) })["elapsed"]

iox_full_mcmc$time_elapsed <- spiox_time_mcmc
save(iox_full_mcmc, file="sim_pmatern/spiox_pmatern_fit5_alt.RData")

################################################################################
# parsimonious Matern via GpGpm
################################################################################

source("fit_parsimonious_matern.R")

gpgpm_fit <- parsimonious_matern(Y, X, coords, m=30, 
                                 filename="sim_pmatern/GpGpm_pmatern_fit5_alt.RData")


################################################################################
# LMC via meshed
################################################################################
library(meshed)

n <- nrow(coords)
q <- ncol(Y)
mcmc <- 5000

set.seed(2026)

lmc_time_mcmc <- system.time({
  lmc_meshed <- spmeshed(y=Y, x=X, 
                         coords=coords, k=q, block_size=40, n_samples=mcmc, n_threads=16, verbose=10,
                       prior=list(phi=c(0.3, 200))) })

Omega_mcmc <- lmc_meshed$lambda_mcmc %>% meshed:::cube_correl_from_lambda()
Qpc_mcmc <- lmc_meshed$lambda_mcmc %>% meshed:::cube_tcrossprod() %>% apply(3, \(x) -cov2cor(solve(x))) %>%
  array(dim=c(q,q,mcmc))


save(lmc_meshed, file="sim_pmatern/lmc_pmatern_fit5_alt.RData")













