
plot_pcor <- function(Q, tol = 1e-8) {
  p <- nrow(Q)
  D <- diag(1 / sqrt(diag(Q)))
  R <- -D %*% Q %*% D
  diag(R) <- 1
  
  df <- melt(R)
  df$value[abs(df$value) < tol & df$Var1 != df$Var2] <- NA
  diag_idx <- df$Var1 == df$Var2
  df$value[diag_idx] <- NA
  
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_y_reverse() +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), 
                         na.value = "transparent", name = "Partial\nCorrelation") +
    coord_equal() +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

glasso_roc <- function(Sigma_hat, Q_true, nlambda = 50) {
  library(glasso)
  p <- nrow(Sigma_hat)
  true_edges <- abs(Q_true[upper.tri(Q_true)]) > 1e-6
  n_pos <- sum(true_edges)
  n_neg <- sum(!true_edges)
  
  lambdas <- exp(seq(log(0.001), log(1), length.out = nlambda))
  path <- glassopath(Sigma_hat, rholist = lambdas)
  
  roc <- do.call(rbind, lapply(1:nlambda, function(k) {
    Wi <- path$wi[,,k]
    est_edges <- abs(Wi[upper.tri(Wi)]) > 1e-16
    tp <- sum(est_edges & true_edges)
    fp <- sum(est_edges & !true_edges)
    data.frame(lambda = lambdas[k], tpr = tp / n_pos, fpr = fp / n_neg)
  }))
  
  return(list(roc=roc,
              plot=ggplot(roc, aes(fpr, tpr)) +
    geom_path(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_abline(lty = 2, color = "grey50") +
    labs(x = "False positive rate", y = "True positive rate") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal()))
}


get_gamma_mat <- function(nu, d){
  gamma_mat <- matrix(0, q, q)
  for (i in 1:q) {
    for (j in 1:q) {
      nu_avg <- (nu[i] + nu[j]) / 2
      gamma_mat[i, j] <-
        sqrt(gamma(nu[i] + d/2) / gamma(nu[i])) *
        sqrt(gamma(nu[j] + d/2) / gamma(nu[j])) *
        gamma(nu_avg) / gamma(nu_avg + d/2)
    }
  }
  return(gamma_mat)
}

get_gamma_sqexp <- function(ell){
  gamma_coef_se <- function(ell_i, ell_j, d = 2) {
    (2 * ell_i * ell_j / (ell_i^2 + ell_j^2))^(d / 2)
  }
  
  gamma_mat <- matrix(0, q, q)
  for (i in 1:q) {
    for (j in 1:q) {
      gamma_mat[i, j] <- gamma_coef_se(ell[i], ell[j])
    }
  }
  return(gamma_mat)
}

pmatern_cov <- function(coords, Sigma, nu, phi) {
  # coords: n x d matrix of locations
  # Sigma:  q x q colocated covariance matrix
  # nu:     length-q vector of smoothness parameters
  # phi:    scalar range parameter
  
  n <- nrow(coords)
  q <- length(nu)
  d <- ncol(coords)
  h <- as.matrix(dist(coords))
  
  # gamma_ij coefficients
  gamma_mat <- get_gamma_mat(nu, d)
  
  # Matern correlation matrix for given smoothness
  matern_cor <- function(h, nu, phi) {
    out <- matrix(0, n, n)
    idx <- h > 0
    z <- phi * h[idx]
    out[idx] <- 2^(1 - nu) / gamma(nu) * z^nu * besselK(z, nu)
    out[!idx] <- 1
    out
  }
  
  # assemble qn x qn matrix: block (i,j) is sigma_ij * gamma_ij * M(h; (nu_i+nu_j)/2, phi)
  C <- matrix(0, q * n, q * n)
  for (i in 1:q) {
    for (j in 1:q) {
      nu_avg <- (nu[i] + nu[j]) / 2
      block <- Sigma[i, j] * gamma_mat[i, j] * matern_cor(h, nu_avg, phi)
      C[((i-1)*n + 1):(i*n), ((j-1)*n + 1):(j*n)] <- block
    }
  }
  return(list(C=C, gamma_mat=gamma_mat))
}

sqexp_cov <- function(coords, Sigma, ell) {
  # coords: n x d matrix of locations
  # Sigma:  q x q colocated covariance matrix
  # nu:     length-q vector of smoothness parameters
  # phi:    scalar range parameter
  
  n <- nrow(coords)
  q <- length(nu)
  d <- ncol(coords)
  h <- as.matrix(dist(coords))
  
  sq_exp_conv <- function(h, ell_i, ell_j) {
    exp(-h^2 / (2 * (ell_i^2 + ell_j^2)))
  }
  
  gamma_mat <- get_gamma_sqexp(ell)
  
  # assemble qn x qn matrix: block (i,j) is sigma_ij * gamma_ij * M(h; (nu_i+nu_j)/2, phi)
  C <- matrix(0, q * n, q * n)
  for (i in 1:q) {
    for (j in 1:q) {
      block <- Sigma[i, j] * gamma_mat[i,j] * sq_exp_conv(h, ell[i], ell[j])
      C[((i-1)*n + 1):(i*n), ((j-1)*n + 1):(j*n)] <- block
    }
  }
  return(list(C=C + 1e-6*diag(n*q), gamma_mat=gamma_mat))
}

gamma_coef <- function(nu_i, nu_j, d = 2) {
  sqrt(gamma(nu_i + d / 2) / gamma(nu_i)) *
    sqrt(gamma(nu_j + d / 2) / gamma(nu_j)) *
    gamma((nu_i + nu_j) / 2) /
    gamma((nu_i + nu_j) / 2 + d / 2)
}

################# Matern correlation function ################ 
matern_cor <- function(h, nu, phi) {
  out <- rep(1, length(h))
  idx <- h > 1e-10
  z <- phi * h[idx]
  out[idx] <- (2^(1 - nu) / gamma(nu)) * z^nu * besselK(z, nu)
  out
}

sq_exp_conv <- function(h, ell_i, ell_j) {
  exp(-h^2 / (2 * (ell_i^2 + ell_j^2)))
}

gamma_coef_se <- function(ell_i, ell_j, d = 2) {
  (2 * ell_i * ell_j / (ell_i^2 + ell_j^2))^(d / 2)
}
