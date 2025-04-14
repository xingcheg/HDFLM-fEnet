##### function to estimate refined coefficient estimates
flm_refine <- function(X, Y, t, N, p, df, lambda){
  nGrid <- length(t)

  # compute K^{1/2}:(nGrid, nGrid)
  Ksqrt <- K_sqrt(t)
  
  # compute Z:(N, p, nGrid)
  Z <- array(0, dim = c(N, p, nGrid))
  for (j in 1:p){
    Z[,j,] <- X[,j,] %*% Ksqrt / nGrid
  }
  
  # compute basis B:(p, nGrid, df) and matrix G:(p, df, df)
  B <- array(0, dim = c(p, nGrid, df))
  rfactor <- rep(0, p)
  for (j in 1:p){
    Tj <- cov(Z[,j,])
    eig_Tj <- eigen(Tj)
    phi_Tj <- eig_Tj$vectors * sqrt(nGrid)
    B[j,,] <- phi_Tj[,1:df]
    rfactor[j] <- sqrt(eig_Tj$values[1])
    Z[,j,] <- Z[,j,] / rfactor[j]
  }
  
  
  # compute Gamma:(N, df*p)
  Gamma <- matrix(0, N, df*p)
  for (j in 1:p){
    Gamma_j <- Z[,j,] %*% B[j,,] / nGrid
    Gamma[,(1:df)+df*(j-1)] <- Gamma_j
  }
  
  S1 <- crossprod(Gamma) / N
  diag(S1) <- diag(S1) + lambda
  inv_S1 <- solve(S1)
  S2 <- crossprod(Gamma, Y) / N
  vec_c <- inv_S1 %*% S2
  mat_c <- matrix(vec_c, df, p)
  
  beta_hat <- matrix(0, nGrid, p)
  for (j in 1:p){
    beta_hat[,j] <- (1/rfactor[j]) * Ksqrt %*% B[j,,] %*% mat_c[,j] / nGrid
  }
  
  return(beta_hat)
}

