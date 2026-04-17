library(tidyverse)
library(glasso)
library(ggplot2)
library(reshape2)
library(BDgraph)

RhpcBLASctl::blas_set_num_threads(16)
RhpcBLASctl::omp_set_num_threads(1)

source("functions.R")

s <- 1


results <- list()


for(s in 1:100){
  
  set.seed(s)
  
  
  coords <- expand.grid(x <- seq(0,1,length.out=32), x) %>% as.matrix()
  colnames(coords) <- NULL
  nr <- nrow(coords)
  q <- 10
  
  # sample a random graph
  G <- bdgraph.sim(p = q, graph = "random", prob = .2)$G
  
  # sample from G-Wishart with that graph
  # b = degrees of freedom (b > 2 keeps things well-conditioned)
  # D = scale matrix
  Q <- rgwish(n = 1, adj = G, b = 5, D = diag(q))
  Q <- cov2cor(Q) 
  Sigma <- solve(Q)
  
  phi <- 10
  nu <- runif(q, 0.5, 2.5)
  
  nu_true <- outer(nu, nu, function(a, b) (a + b) / 2)
  diag(nu_true) <- nu
  
  parsim_matern_C <- pmatern_cov(coords, Sigma, nu, phi)
  C <- parsim_matern_C$C
  gamma_mat <- parsim_matern_C$gamma_mat
  L <- t(chol(C))
  
  y <- L %*% rnorm(nr * q)
  
  
  library(GpGpm)
  source("~/GpGp_multi_paper/R/helper.R")
  source("~/GpGp_multi_paper/R/fisher_multi.R")
  source("~/GpGp_multi_paper/R/fit_multi.R")
  source("~/GpGp_multi_paper/R/link_multi.R")
  source("~/GpGp_multi_paper/R/check.R")
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(16)
  
  # fit  GpGpm
  locs <- do.call(rbind, lapply(1:q, function(j) data.frame(coords, j = j)))
  
  X <- model.matrix(lm( y ~ -1 + as.factor(locs[,ncol(locs)])))   
  
  ncomp <- q
  neach <- ncomp*(ncomp+1)/2
  
  inds <- matrix(FALSE, ncomp, ncomp)    
  diag(inds) <- TRUE
  
  linkfuns <- list()
  
  ######## PARSIMONIOUS MATERN
  linkfuns$link <- link_pars
  linkfuns$dlink <- dlink_pars
  
  
  a <- log(1/phi) #log(mean(fit1$covparms[marginal_ran_inds]))
  start_logparms <- c(log1p(diag(diag(cov(matrix(y, ncol=q)))))[upper.tri(inds, diag = TRUE)], # start from cov(Y)
                      log(1/phi), # phi true
                      log(diag(nu_true)), # nu true
                      log(1e-5)*diag(q)[upper.tri(inds, diag = TRUE)]) # nuggets, no cross
  
  active_logparms <- rep(FALSE, length(start_logparms))
  active_logparms[1:neach] <- TRUE
  #active_logparms <- rep(TRUE, length(start_logparms))
  #active_logparms[ pars_log_cross_nug_inds ] <- FALSE 
  #active_logparms[(length(active_logparms)-neach+1):length(active_logparms)] <- FALSE
  
  print("starting parsimonious fit") 
  t1 <- proc.time()
  fit2 <- fit_multi(
    y, locs, X,
    neighbor_fun = nearest_multi_any, 
    order_fun = order_completely_random,
    m = 20,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
  )
  t2 <- proc.time()
  fit2$time_elapsed <- (t2-t1)[3]   
  fit2$valid <- T
  
  # this is equal to Sigma * gamma_mat in our parametrization
  SG_hat <- matrix(0, ncomp, ncomp)
  SG_hat[upper.tri(SG_hat, diag = TRUE)] <- fit2$covparms[1:neach]
  SG_hat[lower.tri(SG_hat)] <- t(SG_hat)[lower.tri(SG_hat)]
  
  phi_hat <- matrix(0, ncomp, ncomp)
  phi_hat[upper.tri(phi_hat, diag = TRUE)] <- fit2$covparms[neach+(1:neach)]
  phi_hat[lower.tri(phi_hat)] <- t(phi_hat)[lower.tri(phi_hat)]
  
  nu_hat <- matrix(0, ncomp, ncomp)
  nu_hat[upper.tri(nu_hat, diag = TRUE)] <- fit2$covparms[neach*2 + (1:neach)]
  nu_hat[lower.tri(nu_hat)] <- t(nu_hat)[lower.tri(nu_hat)]
  
  gamma_mat <- get_gamma_mat(diag(nu_hat), 2)
  
  
  Sigma_hat <- SG_hat/gamma_mat
  Q_hat <- solve(Sigma_hat)
  
  glasso_out <- glasso::glasso(Sigma_hat, .2)
  
  roc_nn <- glasso_roc(cov(matrix(y, ncol=q)), Q, 60)
  
  roc_sp <- glasso_roc(Sigma_hat, Q, 60)
  
  
  results[[s]] <- list(graph = G,
                       Sigma = Sigma,
                       Q = Q,
                       phi = phi,
                       nu = nu, 
                       y = y,
                       coords = coords, 
                       fit_GpGpm = fit2,
                       Sigma_hat = Sigma_hat,
                       roc_sp = roc_sp,
                       roc_nn = roc_nn)
  
  save(results, file = "glasso/results_LR.RData")

}


