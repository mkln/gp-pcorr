library(BDgraph)
library(ggplot2)
library(patchwork)
library(reshape2)

set.seed(2026)

q <- 5

coords <- expand.grid(x <- seq(0, 1, length.out=25), x)
n <- nrow(coords)

h <- as.matrix(dist(coords))

source("figures/make_sparse_Q_5.r")

edges <- list(
  list(from = 3, to = 1, weight = 1.1),   # confounder
  list(from = 3, to = 2, weight = 1.1),
  list(from = 1, to = 2, weight = -1.2)
)

Q_structured <- generate_network_matrices(
  build_B(q, add_chain(edges, c(3, 4, 5), weight = 0.7)),
  noise_var = c(1, 1, 1, 1, 1)  # tweakable!
)

Q <- Q_structured$Precision
Qpc <- Q_structured$Partial_Correlation
Sigma <- Q_structured$Covariance
Rc <- Q_structured$Marginal_Correlation


source("functions.R")

nu <- c(0.2, 1, 0.5, 1.4, 0.75); phi <- 10 # for matern
C_str <- pmatern_cov(coords, Sigma, nu, phi)

L <- t(chol(C_str$C))
y <- L %*% rnorm(n * q)
Y <- matrix(y, ncol=q)


df <- data.frame(coords, Y)
colnames(df)[1:2] <- c("x", "y")


write_csv(df, file="sim_pmatern/pmatern5.csv")

