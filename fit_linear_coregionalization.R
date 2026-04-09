

fit_spbayes_lmc <- function(Y, X, coords, n.samples = 500, burn.in = NULL,
                            cov.model = "exponential", n.report = 10, 
                            filename = NULL, prev_fit = NULL) {
  q <- ncol(Y)
  n <- nrow(Y)
  nltr <- q * (q + 1) / 2
  if (is.null(burn.in)) burn.in <- floor(0.75 * n.samples)
  if (!is.list(X)) X <- replicate(q, X, simplify = FALSE)
  
  obs_rows <- which(complete.cases(Y))
  coords_obs <- coords[obs_rows, , drop = FALSE]
  
  formulas <- list()
  for (j in 1:q) {
    yname <- paste0("y.", j)
    xname <- paste0("x.", j)
    assign(yname, Y[obs_rows, j])
    assign(xname, X[[j]][obs_rows, , drop = FALSE])
    formulas[[j]] <- as.formula(paste0(yname, " ~ ", xname, " - 1"))
  }
  
  if (!is.null(prev_fit)) {
    # extract last sample from previous run
    theta_last <- prev_fit$p.theta.samples[nrow(prev_fit$p.theta.samples), ]
    
    K <- matrix(0, q, q)
    K[lower.tri(K, diag=TRUE)] <- theta_last %>% head(q*(q+1)/2)
    K[upper.tri(K)] <- t(K)[upper.tri(K)]
    A <- t(chol(K))
    Avec <- A[lower.tri(A, diag=T)]
    
    nms <- colnames(prev_fit$p.theta.samples)
    
    starting <- list(
      phi = theta_last[grep("^phi", nms)],
      A   = Avec,#theta_last[grep("^A", nms)],
      Psi = theta_last[grep("^Psi", nms)],
      nu  = theta_last[grep("^nu", nms)]
    )
    
  } else {
    starting <- list(phi = rep(6, q), A = diag(1, q)[lower.tri(diag(1, q), TRUE)],
                     Psi = rep(1, q), nu = rep(1, q))
    
  }
  tuning <- list(
    phi = rep(0.0001, q), A = rep(0.0001, nltr), Psi = rep(0.0001, q),
    nu  = rep(0.05, q)
  )
  
  print(tuning)
  
  priors <- list("beta.Flat",
                 phi.Unif = list(rep(0.1, q), rep(200, q)),
                 K.IW = list(q + 1, diag(0.1, q)),
                 Psi.ig = list(rep(2, q), rep(0.1, q)),
                 nu.Unif = list(rep(0.2, q), rep(2.1, q)))
  
  fitting_time <- 
    system.time({
      fit <- spMvLM(formulas, coords = coords_obs,
                    starting = starting, tuning = tuning, priors = priors,
                    n.samples = n.samples, cov.model = cov.model, n.report = n.report)
    })
  
  #fit <- spRecover(fit, start = burn.in)
  
  fit$obs_rows <- obs_rows
  fit$q <- q
  fit$elapsed <- fitting_time["elapsed"]
  
  if (!is.null(filename)) {
    lmc_fit <- fit
    save(lmc_fit, file = filename)
  }
  
  fit
}


predict_spbayes_lmc <- function(fit, X_pred, coords_pred, coords_obs, thin = 10) {
  
  # coords_pred: n_pred x 2
  # coords_obs:  n_obs x 2 (training coordinates)
  
  fit <- spRecover(fit, start = burn.in)
  
  q <- fit$q
  if (!is.list(X_pred)) X_pred <- replicate(q, X_pred, simplify = FALSE)
  pred.covars <- mkMvX(X_pred)
  
  out <- spPredict(fit, start = 1, thin = thin,
                   pred.coords = coords_pred, pred.covars = pred.covars)
  
  # posterior predictive summaries
  quants <- apply(out$p.y.predictive.samples, 1, 
                  function(x) quantile(x, c(0.5, 0.025, 0.975)))
  
  n_pred <- nrow(coords_pred)
  # unstack: spMvLM interleaves outcomes
  Y_hat <- list()
  for (j in 1:q) {
    idx <- seq(j, n_pred * q, by = q)
    Y_hat[[j]] <- t(quants[, idx])
    colnames(Y_hat[[j]]) <- c("median", "lo", "hi")
  }
  
  list(Y_hat = Y_hat, samples = out$p.y.predictive.samples, fit_out=out)
}






