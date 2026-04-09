library(GpGpm)
source("~/GpGp_multi_paper/R/helper.R")
source("~/GpGp_multi_paper/R/fisher_multi.R")
source("~/GpGp_multi_paper/R/fit_multi.R")
source("~/GpGp_multi_paper/R/link_multi.R")
source("~/GpGp_multi_paper/R/check.R")


source("~/spiox-network/functions.R")


RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(16)


parsimonious_matern <- function(Y, X, coords, m, filename){
  
  #X <- X[complete.cases(Y),]
  #coords <- coords[complete.cases(Y),]
  #Y <- Y[complete.cases(Y),]
  
  q <- ncol(Y)
  
  y <- as.vector(Y)
  
  X_list <- list()
  for( i in 1:q ){
    X_list[[i]] <- data.frame(X, j=i)
  } 
  
  X <- bind_rows(X_list) 
  X$j <- factor(X$j)
  X <- model.matrix(data=X, ~.-1)[,-(q+1)]
  
  # fit  GpGpm
  locs <- do.call(rbind, lapply(1:q, function(j) data.frame(coords, j = j)))
  locs <- as.matrix(locs)
  
  ncomp <- q
  neach <- ncomp*(ncomp+1)/2
  
  inds <- matrix(FALSE, ncomp, ncomp)    
  diag(inds) <- TRUE
  
  linkfuns <- list()
  
  ######## PARSIMONIOUS MATERN
  linkfuns$link <- link_pars
  linkfuns$dlink <- dlink_pars
  
  
  #a <- log(1/phi) #log(mean(fit1$covparms[marginal_ran_inds]))
  Y_a <- Y[complete.cases(Y),]
  start_logparms <- c(log1p(diag(diag(cov(Y_a))))[upper.tri(inds, diag = TRUE)], # start from cov(Y)
                      log(1/20), 
                      log(rep(exp(1), q)), 
                      log(1e-5)*diag(q)[upper.tri(inds, diag = TRUE)]) # nuggets, no cross
  
  active_logparms <- (start_logparms!=0) #rep(FALSE, length(start_logparms))
  active_logparms[1:neach] <- TRUE # sigma
  active_logparms[neach+1] <- TRUE # phi
  active_logparms[neach+1+(1:q)] <- TRUE # nu
  
  #active_logparms <- rep(TRUE, length(start_logparms))
  #active_logparms[ pars_log_cross_nug_inds ] <- FALSE 
  #active_logparms[(length(active_logparms)-neach+1):length(active_logparms)] <- FALSE
  
  print("starting parsimonious fit") 
  t1 <- proc.time()
  fit <- fit_multi(
    y, locs, X,
    neighbor_fun = nearest_multi_any, 
    order_fun = order_completely_random,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms, max_iter = 100
  )
  t2 <- proc.time()
  fit$time_elapsed <- (t2-t1)[3]   
  fit$valid <- T
  
  # this is equal to Sigma * gamma_mat in our parametrization
  SG_hat <- matrix(0, ncomp, ncomp)
  SG_hat[upper.tri(SG_hat, diag = TRUE)] <- fit$covparms[1:neach]
  SG_hat[lower.tri(SG_hat)] <- t(SG_hat)[lower.tri(SG_hat)]
  
  nug <- matrix(0, ncomp, ncomp)
  nug[upper.tri(nug, diag = TRUE)] <- head(fit$covparms, neach) * tail(fit$covparms, neach)
  nug[lower.tri(nug)] <- t(nug)[lower.tri(nug)]
  
  phi_hat <- matrix(0, ncomp, ncomp)
  phi_hat[upper.tri(phi_hat, diag = TRUE)] <- 1/fit$covparms[neach+(1:neach)]
  phi_hat[lower.tri(phi_hat)] <- t(phi_hat)[lower.tri(phi_hat)]
  
  nu_hat <- matrix(0, ncomp, ncomp)
  nu_hat[upper.tri(nu_hat, diag = TRUE)] <- fit$covparms[neach*2 + (1:neach)]
  nu_hat[lower.tri(nu_hat)] <- t(nu_hat)[lower.tri(nu_hat)]
  
  gamma_mat <- get_gamma_mat(diag(nu_hat), 2)
  
  
  Sigma_hat <- (SG_hat+nug)/gamma_mat
  Q_hat <- solve(Sigma_hat)
  
  gpgpm_fit <- list(
    Sigma=Sigma_hat,
    Q = Q_hat,
    gamma_nu = gamma_mat,
    phi = phi_hat[1,1],
    nu = diag(nu_hat),
    nugg = nug,
    fit = fit
  )
  
  save(gpgpm_fit, file = filename)
  
  return(gpgpm_fit)
}
