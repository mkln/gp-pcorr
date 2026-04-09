rm(list=ls())
library(spiox)
library(tidyverse)
library(magrittr)

g <- glue::glue

for( s in 1:200 ){
  cat(s, "\n")
  set.seed(s)
  source("functions.R")
  
  ################################################################################
  # Data generation
  ################################################################################
  coords <- expand.grid(x <- seq(0, 1, length.out=25), x)
  n <- nrow(coords)
  
  q <- 5
  source("figures/make_sparse_Q_5.r")
  edges <- list(
    list(from = 3, to = 1, weight = 1.1),
    list(from = 3, to = 2, weight = 1.1),
    list(from = 1, to = 2, weight = -1.2)
  )
  Q_structured <- generate_network_matrices(
    build_B(q, add_chain(edges, c(3, 4, 5), weight = 0.7)),
    noise_var = c(1, 1, 1, 1, 1)
  )
  Q <- Q_structured$Precision
  Qpc <- Q_structured$Partial_Correlation
  Sigma <- Q_structured$Covariance
  Rc <- Q_structured$Marginal_Correlation
  #nu <- c(0.2, 1, 0.5, 1.4, 0.75); phi <- 10
  nu <- runif(q, 0.2, 1.8); phi <- runif(1, 10, 40)
  C_str <- pmatern_cov(coords, Sigma, nu, phi)
  
  L <- t(chol(C_str$C))
  y <- L %*% rnorm(n * q)
  Y <- matrix(y, ncol=q)
  Y <- Y %>% apply(2, \(y) y - mean(y))
  
  Y_m <- Y
  Y_m[sample(1:(n*q), 200, replace=FALSE)] <- NA
  na_mask <- is.na(Y_m)
  
  X <- matrix(1, ncol=1, nrow=n)
  
  all_na <- apply(Y_m, 1, function(r) all(is.na(r)))
  Y      <- Y[!all_na, , drop=FALSE]
  Y_m    <- Y_m[!all_na, , drop=FALSE]
  coords <- coords[!all_na, , drop=FALSE]
  X      <- X[!all_na, , drop=FALSE]
  na_mask <- is.na(Y_m)
  n      <- nrow(Y)
  
  Gamma_nu <- get_gamma_mat(nu, 2)
  
  save(list=c("Y_m", "Y", "X", "coords", "na_mask", "q", "nu", "phi", "Sigma", "Q", "Gamma_nu"),
       file=g("sim_pmatern/reps/data_{s}.RData"))
}











