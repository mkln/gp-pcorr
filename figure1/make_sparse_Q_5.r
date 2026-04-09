generate_network_matrices <- function(B, noise_var = NULL) {
  
  p <- nrow(B)
  I <- diag(p)
  
  # Default: homoskedastic noise
  if (is.null(noise_var)) {
    noise_var <- rep(1, p)
  }
  D_eps <- diag(noise_var)
  
  # SEM: X = B^T X + E  →  (I - B^T) X = E
  A <- I - t(B)
  
  # Precision and covariance
  Omega <- t(A) %*% solve(D_eps) %*% A
  Sigma <- solve(Omega)
  
  # Correlations
  D_inv_sqrt <- diag(1 / sqrt(diag(Sigma)))
  Cor_marginal <- D_inv_sqrt %*% Sigma %*% D_inv_sqrt
  
  D_omega_inv_sqrt <- diag(1 / sqrt(diag(Omega)))
  Cor_partial <- -D_omega_inv_sqrt %*% Omega %*% D_omega_inv_sqrt
  diag(Cor_partial) <- 1
  
  list(
    B = B,
    Precision = Omega,
    Covariance = Sigma,
    Partial_Correlation = Cor_partial,
    Marginal_Correlation = Cor_marginal
  )
}

build_B <- function(p, edges) {
  B <- matrix(0, p, p)
  
  for (e in edges) {
    # e = list(from=, to=, weight=)
    B[e$from, e$to] <- e$weight
  }
  
  B
}

make_sign_reversal_dag <- function(
    confounder = 3,
    i = 1,
    j = 2,
    direct = -0.5,
    conf_strength = 1.0,
    extra_edges = list()
) {
  
  edges <- list(
    list(from = confounder, to = i, weight = conf_strength),
    list(from = confounder, to = j, weight = conf_strength),
    list(from = i, to = j, weight = direct)
  )
  
  # allow arbitrary extensions
  edges <- c(edges, extra_edges)
  
  build_B(max(c(confounder, i, j, unlist(extra_edges))), edges)
}

add_chain <- function(edges, nodes, weight) {
  for (k in seq_len(length(nodes) - 1)) {
    edges[[length(edges) + 1]] <- list(
      from = nodes[k],
      to   = nodes[k + 1],
      weight = weight
    )
  }
  edges
}
