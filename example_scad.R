source("simu_funcs.R")
source("algo_scad.R")
set.seed(123)

############ Hyperparameters #############
K <- 30
nGrid <- 201
N <- 100
N_test <- 100
p <- 200
q <- 10
sigma <- 0.5
t <- seq(0, 1, length.out = nGrid)
true_P <- 1:q
true_N <- (q+1):(p)
thres_factor <- 10

### auto-correlation
rho <- 0.75

## true beta
beta <- create_Beta(p, q, nGrid, K)

## generate X and Y
X <- gen_X(K, nGrid, N, p, rho=rho)
Y <- gen_Y(X, beta, sigma)

lambda <- 0.5
s <- 5

## estimate beta using flm-scad
theta <- get_theta(X, s)
beta_hat <- scad_flm(Y, theta, N, p, s, nGrid, lambda)
# take around 10 seconds.


## visualize results
par(mfrow = c(4,5), mar = rep(2.5, 4))
for (kk in 1:20){
  plot(t, beta[,kk], type = "l", lty=2, lwd=2,
       ylab = "beta", main = paste("beta", kk))
  lines(t, beta_hat[,kk], lwd=2, col = "red")
}
par(mfrow = c(1,1))