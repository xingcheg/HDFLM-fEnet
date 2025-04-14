library(splines)


## cosine basis function (no constant function)
cosine_basis <- function(K, nGrid){
  pts <- seq(0, 1, length.out = nGrid)
  res <- sapply(1:K, function(k) {
    sqrt(2) * cos( k * pi * pts )
  })
  colnames(res) <- paste("b", 1:K, sep="")
  return(res)
}



## Haar basis function
haar_basis <- function(K, nGrid){
  pts <- seq(0, 1, length.out = nGrid)
  res <- matrix(0, nGrid, K)
  for (k in 1:K){
    n <- floor(log(k)/log(2))
    l  <- k - 2^n + 1
    lower <- (l-1)/(2^n)
    mid <- (l-1/2)/(2^n)
    upper <- l/(2^n)
    res[pts>=lower & pts<mid, k] <- 2^(n/2)
    res[pts>=mid & pts<=upper, k] <- -2^(n/2)
  }
  colnames(res) <- paste("haar", 1:K, sep="")
  return(res)
}



## coefficient functions (beta)
create_Beta <- function(p, q, nGrid, K){
  Beta <- matrix(0, nGrid, p)
  pts <- seq(0, 1, length.out = nGrid)
  b_cos <- t( cosine_basis(K, nGrid) )
  for (i in 1:q){
    l1 <- exp(-(1:K)/4)
    l2 <- 2 * rbinom(K, 1, 0.5) - 1
    S <- 4 * b_cos * l1 * l2
    Beta[,i] <- apply(S, 2, sum)
  }
  return(Beta)
}



## autoregressive covariance
auto_cov <- function(rho, p){
  Sigma_l <- matrix(0, p, p)
  for (i in 2:p){
    for (j in 1:(i-1)){
      Sigma_l[i,j] <- rho^(abs(i-j))
    }
  }
  Sigma <- Sigma_l + t(Sigma_l)
  for (i in 1:p){
    Sigma[i,i] <- 1
  }
  return(Sigma)
} 


## generate multivariate functional data 
gen_X <- function(K, nGrid, N, p, basis="cosine", rho=0){
  if (basis=="cosine"){
    Phi <- cosine_basis(K, nGrid)
  }
  if (basis=="haar"){
    Phi <- haar_basis(K, nGrid)
  }
  X <- array(0, dim = c(N, p, nGrid))
  lambda <- exp(-(1:K)/4)
  
  if (rho==0){
    for (j in 1:p){
      Xi0 <- matrix(rnorm(K*N, 0, 1), K, N)
      Xi <- Xi0 * sqrt(lambda)
      X[, j, ] <- t(Phi %*% Xi)
      }
    } else if (rho>0 & rho<1){
    Sigma <- auto_cov(p, rho=rho)
    U <- chol(Sigma)
    for (i in 1:N){
      Xi0 <- matrix(rnorm(p*K, 0, 1), K, p) * sqrt(lambda)
      Xi <- Xi0 %*% U
      X[i,,] <- t(Phi %*% Xi)
      }
    } else {
    cat("invalid rho. \n")
      }

  return(X)
}




## generate Y
gen_Y <- function(X, beta, sigma){
  N <- dim(X)[1]
  y <- rep(0, N)
  for (i in 1:N){
    Temp <- t(X[i, , ]) * beta
    y[i] <- sum( apply(Temp, 2, mean) ) + rnorm(1) * sigma
  }
  return(y)
}



## recover the coefficients of X given an orthogonal basis
get_theta <- function(X, s){
  N <- dim(X)[1]
  p <- dim(X)[2]
  nGrid <- dim(X)[3]
  Phi <- cosine_basis(s, nGrid)
  theta <- vector("list", p)
  for (j in 1:p){
    tj <- matrix(0, N, s)
    for (k in 1:s){
      tj[,k] <- X[,j,] %*% Phi[,k] / nGrid
    }
    theta[[j]] <- tj
  }
  return(theta)
}




## generate multivariate functional data (*)
gen_X1 <- function(K, nGrid, N, p, U, rho=0){
  Phi0 <- cosine_basis(K, nGrid)
  Phi <- Phi0 %*% U
  
  X <- array(0, dim = c(N, p, nGrid))
  lambda <- exp(-(1:K)/4)
  
  if (rho==0){
    for (j in 1:p){
      Xi0 <- matrix(rnorm(K*N, 0, 1), K, N)
      Xi <- Xi0 * sqrt(lambda)
      X[, j, ] <- t(Phi %*% Xi)
    }
  } else if (rho>0 & rho<1){
    Sigma <- auto_cov(p, rho=rho)
    U <- chol(Sigma)
    for (i in 1:N){
      Xi0 <- matrix(rnorm(p*K, 0, 1), K, p) * sqrt(lambda)
      Xi <- Xi0 %*% U
      X[i,,] <- t(Phi %*% Xi)
    }
  } else {
    cat("invalid rho. \n")
  }
  
  return(X)
}

