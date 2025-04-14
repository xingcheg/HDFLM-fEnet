source("simu_funcs.R")
source("algo_enet_representer.R")
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

lambda <- 0.05
alpha <- 1
df <- 10
theta <- 1e-5

## estimate beta using fEnet
beta_hat <- enet_flm1(X, Y, t, N, p, df, 
                      lambda, alpha, rho = theta, thres_factor, 
                      M_iter = 1e3, tol = 1e-4)
# take around 10 seconds.



## visualize results
par(mfrow = c(4,5), mar = rep(2.5, 4))
for (kk in 1:20){
  plot(t, beta[,kk], type = "l", lty=2, lwd=2,
       ylab = "beta", main = paste("beta", kk))
  lines(t, beta_hat[,kk], lwd=2, col = "red")
}
par(mfrow = c(1,1))



## refined estimator using ridge regression
source("algo_refine.R")

## estimated true signal set
beta_norm <- apply(beta_hat, 2, function(x){
  sqrt(mean(x^2))
})
est_P <- which(beta_norm>0)
est_P

beta_hat_refine <- flm_refine(X=X[,est_P,], Y=Y, t=t, N=N, 
                              p=length(est_P), df=15, lambda=3e-9)

## visualize results
par(mfrow = c(4,5), mar = rep(2.5, 4))
for (kk in 1:10){
  plot(t, beta[,kk], type = "l", lty=2, lwd=2,
       ylab = "beta", main = paste("beta", kk))
  lines(t, beta_hat[,kk], lwd=2, col = "red")
  lines(t, beta_hat_refine[,kk], lwd=2, col = "blue")
}
par(mfrow = c(1,1))

