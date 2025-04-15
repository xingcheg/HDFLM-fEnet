library(parallel)
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

### auto-correlation
rho <- 0.75

## true beta
beta <- create_Beta(p, q, nGrid, K)

## number of replicates
M <- 200

## tuning parameters
Lambda <- c(5, 2, 1, 5e-1, 2e-1)
S <- c(3, 5, 10, 15)

nn1 <- length(Lambda)
nn2 <- length(S)


output <- mclapply(1:M, function(i){
  set.seed(111 + i*222)
  X <- gen_X(K, nGrid, N, p, rho=rho)
  Y <- gen_Y(X, beta, sigma)
  X_test <- gen_X(K, nGrid, N_test, p, rho=rho)
  Y_test0 <- gen_Y(X_test, beta, sigma=0)
  Y_test <- Y_test0 + sigma * rnorm(N_test)
  OUT_i <- NULL
  for (j in 1:nn1){
    for (k in 1:nn2){
      lambda <- Lambda[j]
      s <- S[k]
      theta <- get_theta(X, s)
      beta_hat <- scad_flm(Y, theta, N, p, s, nGrid, lambda)
      
      beta_norm <- apply(beta_hat, 2, function(x){
        sqrt(mean(x^2))
      })
      
      beta_diff_norm <- apply(beta_hat-beta, 2, function(x){
        sqrt(mean(x^2))
      })
      
      max_beta_diff <- max(beta_diff_norm)
      est_P <- which(beta_norm>0)
      est_N <- which(beta_norm==0)
      TPR <- sum( est_P %in% true_P ) / length(true_P)
      FPR <- sum( est_P %in% true_N ) / length(true_N)
      TNR <- sum( est_N %in% true_N ) / length(true_N)
      FNR <- sum( est_N %in% true_P ) / length(true_P)
      Y_test_est <- gen_Y(X_test, beta_hat, 0)
      RMSPE <- sqrt( mean( (Y_test - Y_test_est)^2 ) )
      ER <- mean( (Y_test0 - Y_test_est)^2 )
      TE <- mean( (Y_test0 )^2 )
      RER <- ER / TE
      out_i <- c(i, lambda, s, FPR, FNR, max_beta_diff, RER, RMSPE)
      OUT_i <- rbind(OUT_i, out_i)
      #cat(out, "\n")
    }
  }
  return(OUT_i)
}, mc.cores = 10)
  

output <- as.data.frame(do.call("rbind", output))
names(output) <- c("set_id", "lambda", "s", "FPR", "FNR", "MNE", "RER", "RMSPE")

  
result <- list(output = output, N=N, N_test=N_test, p=p, q=q, K=K, sigma=sigma)
#saveRDS(result, "cv_scad2_auto.rds")


