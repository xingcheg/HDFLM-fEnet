scadderi <- function(lambda, a, t){
  if (t <= lambda){
    scadderi <- lambda
  } else if (t > lambda & t < a*lambda){
    scadderi <- (a*lambda-t)/(a-1)
  } else {
    scadderi <- 0
  }
  return(scadderi)
}


vec_norm <- function(a){
  return(sqrt(sum(a^2)))
}


scad_flm <- function(Y, theta1, N, p, s, nGrid, lambda, a=3.7){
  theta2 <- vector("list", p)
  theta3 <- vector("list", p)
  for (j in 1:p){
    tj <- theta1[[j]]
    theta2[[j]] <- tj %*% solve(crossprod(tj)) %*% t(tj)
    theta3[[j]] <- solve(crossprod(tj)) %*% t(tj)
  }
  distance <- 1
  H <- 1:p
  f <- matrix(0, N, p)
  f1 <- f
  R1 <- vector("list", p)
  P1 <- vector("list", p)
  idx_cum <- 0
  while(distance > 1e-4 & idx_cum < 1e3){
    idx_cum <- idx_cum + 1
    #cat("iter = ", idx_cum, "\n")
    for (j in 1:p){
      R1[[j]] <- Y - apply(f[,-j], 1, sum)
      P1[[j]] <- theta2[[j]] %*% R1[[j]]
      l1 <- lambda * sqrt(s)
      t1 <- vec_norm(f[,j]) / sqrt(N)
      norm_P1 <- vec_norm(P1[[j]])
      flag <- ifelse(norm_P1>0, 1, 0)
      f[,j] <- max(0, 1 - sqrt(N)*scadderi(l1,a,t1)/(norm_P1+1-flag)) * P1[[j]]
      f[,j] <- f[,j] - mean(f[,j])
    }
    distance <- mean(diag(crossprod(f-f1)))
    f1 <- f
  }
  f11 <- f1
  eta <- matrix(0, s, p)
  for (j in 1:p){
    eta[,j] <- theta3[[j]] %*% f11[,j]
  }
  Phi <- cosine_basis(s, nGrid) 
  beta <- Phi %*% eta
  return(beta)
}

