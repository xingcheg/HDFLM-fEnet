library(parallel)
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

## number of replicates
M <- 200

## tuning parameters
Df <- c(10, 15)
Lambda <- c(0.2, 1e-1, 5e-2, 2e-2)
Alpha <- c(1, 1-1e-6, 1-1e-5)
Theta <- c(0, 1e-5, 1e-4)

nn1 <- length(Df)
nn2 <- length(Lambda)
nn3 <- length(Alpha)
nn4 <- length(Theta)

tt1 <- Sys.time()

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
      for (l in 1:nn3){
        for (m in 1:nn4){
          
          df <- Df[j]
          lambda <- Lambda[k]
          alpha <- Alpha[l]
          theta <- Theta[m]
          error_flag <- try(beta_hat <- enet_flm1(X, Y, t, N, p, df, 
                                                  lambda, alpha, rho = theta, thres_factor, 
                                                  M_iter = 1e3, tol = 1e-4), 
                            silent = TRUE)
          
          if(!inherits(error_flag, "try-error")){
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
            out_i <- c(i, lambda, log10(1-alpha), log10(theta), df, FPR, FNR, max_beta_diff, RER, RMSPE)
          } else {
            out_i <- c(i, lambda, log10(1-alpha), log10(theta), df, rep(NA, 5))
          }
          
          OUT_i <- rbind(OUT_i, out_i)
          
        }
      }
    }
  }
  return(OUT_i)
}, mc.cores = 10)

tt2 <- Sys.time()
tt2 - tt1

output <- as.data.frame(do.call("rbind", output))
names(output) <- c("set_id", "lambda", "alpha", "theta", "df", "FPR", "FNR", "MNE", "RER", "RMSPE")


result <- list(output = output, N=N, N_test=N_test, p=p, q=q, K=K, sigma=sigma)
#saveRDS(result, "cv_enet2_auto.rds")


