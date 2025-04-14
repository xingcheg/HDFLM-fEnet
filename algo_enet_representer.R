library(splines)


mat_sqrt <- function(A){
  eigen_A <- eigen(A)
  D_A <- eigen_A$values
  P_A <- eigen_A$vectors
  Asqrt <- P_A %*% diag(sqrt(D_A), nrow = nrow(A)) %*% t(P_A)
}



##### function to compute K^{1/2}
bern_poly4 <- function(x){
  return( x^4 - 2*x^3 + x^2 - 1/30 )
}

K_sqrt <- function(t){
  nGrid <- length(t)
  Kl <- matrix(0, nGrid, nGrid)
  for (i in 2:nGrid){
    for (j in 1:(i-1)){
      Kl[i,j] <- - ( bern_poly4(abs(t[i]-t[j])/2) + bern_poly4((t[i]+t[j])/2) ) / 3
    }
  }
  K <- Kl + t(Kl)
  for (i in 1:nGrid){
    K[i,i] <- - ( bern_poly4(0) + bern_poly4(t[i]) ) / 3
  }
  Ksqrt <- mat_sqrt(K)
  return(Ksqrt)
}





##### function to estimate elastic-net functional linear regression
enet_flm1 <- function(X, Y, t, N, p, df, lambda, alpha, rho, thres_factor,
                      M_iter = 1e4, tol = 1e-4){
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
  
  
  # compute Gamma:(p, N, df)
  Gamma <- array(0, dim = c(p, N, df))
  for (j in 1:p){
    Gamma[j,,] <- Z[,j,] %*% B[j,,] / nGrid
  }
  
  # compute H: (p, df, df)
  H <- array(0, dim = c(p, df, df))
  L_inv <- array(0, dim = c(p, df, df))
  for (j in 1:p){
    H0 <- crossprod( Gamma[j,,] ) / N
    H1 <- H0 + rho * diag(df)
    H[j,,] <- H1
    L <- t(chol(H1))
    L_inv[j,,] <- solve(L)
  }

  # compute Omega:(p, df, df) and tau
  Omega <- array(0, dim = c(p, df, df))
  for (j in 1:p){
    Omega[j,,] <- diag(df) + L_inv[j,,] %*% t(L_inv[j,,]) * (lambda*(1-alpha)-rho)
  }
  
  # initialize D:(df, p)
  D <- matrix(0, df, p)
  r <- Y
  for (i in 1:M_iter){
    D0 <- D
    
    for (j in 1:p){
      if (j == 1){
        s1 <- t(L_inv[p,,]) %*% D[,p]
        s2 <- t(L_inv[1,,]) %*% D[,1]
        r <- r - as.numeric( Gamma[p,,] %*% s1) + as.numeric(Gamma[1,,] %*% s2 )
      } else {
        s1 <- t(L_inv[j-1,,]) %*% D[,j-1]
        s2 <- t(L_inv[j,,]) %*% D[,j]
        r <- r - as.numeric( Gamma[j-1,,] %*% s1) + as.numeric(Gamma[j,,] %*% s2 )
      }
      
      v <- crossprod(Gamma[j,,], r)
      eta <- L_inv[j,,] %*% v / N
      
      norm_eta <- sqrt(sum(eta^2))
      if (norm_eta <= thres_factor*lambda*alpha ){
        D[,j] <- 0
      } else {
        if (i <= 50){
          d_j0 <- D[,j] + rnorm(df, 0, 1e-1)
        } else if (i <= 200){
          d_j0 <- D[,j] + rnorm(df, 0, 1e-3)
        } else {
          d_j0 <- D[,j] + rnorm(df, 0, 1e-12)
        }
        
        for (kk in 1:100){
          norm_d_j0 <- sqrt(sum( d_j0^2 ))
          if (norm_d_j0 == 0){
            d_j <- d_j0
            break
          }
          Q <- Omega[j,,] + {(lambda * alpha)/norm_d_j0} * diag(df) 
          Q_inv <- solve(Q)
          d_j <- Q_inv %*% eta
          norm_diff <- sqrt( sum((d_j0 - d_j)^2) )
          if ( norm_diff < 1e-8 ){
            break
          } else {
            d_j0 <- d_j
          }
        }
        D[,j] <- d_j
      }
    }
    
    if (norm(D-D0, type = "O") < tol){
      break
    }
  }
  
  beta_hat <- matrix(0, nGrid, p)
  for (j in 1:p){
    beta_hat[,j] <- (1/rfactor[j]) * Ksqrt %*% B[j,,] %*% t(L_inv[j,,]) %*% D[,j] / nGrid
  }

  return(beta_hat)
}

