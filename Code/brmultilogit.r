calculateA <- function(x, y, b, ref=NULL) {
  # This is a function use to calculate the log determent of information matrix
  # Args:
  #   x: design matrix
  #   y: response variable
  #   B: J*p coefficients matrix
  n <- nrow(y)
  J <- ncol(y) - 1
  p <- ncol(x)
  
  
  if(is.null(ref) == TRUE){
    ref <- ncol(y)
  }
  
  X <- t(kronecker(t(x), diag(J)))
  Y <- array(t(y[,-ref]), dim = c(length(y[,-ref]),1))
  size <- rowSums(y)
  M <- array(0, dim = c(n ,J ,J))
  
  B <- array(b, dim = c(length(b),1))
  eta <- x %*% t(b)
  denom <- apply(eta, 1, function(x) sum(exp(x)) + 1)
  theta <- exp(eta)/denom
  
  for(i in 1:n){
    for(j in 1:J){
      for(k in 1:J){
        if(j==k){
          M[i,j,k] <- theta[i,j] * (1 - theta[i,j]) * size[i]
        } else {
          M[i,j,k] <- -theta[i,j] * theta[i,k] * size[i]
        }
      }
    }
  }
  
  ML <- list()
  for(i in 1:n){
    ML[[i]] <- M[i,,]
  }
  MM <- bdiag(ML)
  A <- t(X)%*%MM%*%X
  logA <- log(det(A))
  
  return(list(A=A, logA=logA))
}

brmultilogit <- function(x, y, ref=NULL, init=NULL, itermax = 100, bias = FALSE, tol = 1e-5) {
  # This is a function use to calculate bias corrected parameters in logistic regression and multinomial regression
  # Args:
  #   x: design matrix
  #   y: response variable
  #   ref: reference category
  #   init: initial value for parameters
  #   itermax: maximum iteration number
  #   bias: bias correction or not
  #   tol: tolerance for stopping
  
  
  n <- nrow(y)
  J <- ncol(y) - 1
  p <- ncol(x)
  
  if(is.null(init) == TRUE){
    b <- array(0, dim = c(J,p))
  } 
  
  if(is.null(ref) == TRUE){
    ref <- ncol(y)
  }
  
  X <- t(kronecker(t(x), diag(J)))
  Y <- array(t(y[,-ref]), dim = c(length(y[,-ref]),1))
  size <- rowSums(y)
  M <- array(0, dim = c(n ,J ,J))
  
  if(bias == TRUE){
    XX <- kronecker(X, X)
  }
  
  
  iter <- 1
  stop <- 0
  while(stop == 0 & iter <=itermax) {
    #print(iter)
    
    B <- array(b, dim = c(length(b),1))
    eta <- x %*% t(b)
    denom <- apply(eta, 1, function(x) sum(exp(x)) + 1)
    theta <- exp(eta)/denom
    
    for(i in 1:n){
      for(j in 1:J){
        for(k in 1:J){
          if(j==k){
            M[i,j,k] <- theta[i,j] * (1 - theta[i,j]) * size[i]
          } else {
            M[i,j,k] <- -theta[i,j] * theta[i,k] * size[i]
          }
        }
      }
    }
    
    ML <- list()
    for(i in 1:n){
      ML[[i]] <- M[i,,]
    }
    MM <- bdiag(ML)
    A <- t(X)%*%MM%*%X
    
    U <- t(X)  %*% (Y -   rep(size, each = J) * c(t(theta)))
    A_inv <- solve(A)
    
    
    
    
    if(bias == FALSE){
      B_new <- B + A_inv %*% U
    } else {
      Q <- array(0, dim = c(n, J, J^2))
      q <- array(0, dim = c(n, J, J, J))
      for(i in 1:n){
        for(j in 1:J){
          for(k in 1:J){
            for(l in 1:J){
              
              if(j==k & k==l){
                q[i,j,k,l] <- theta[i,j]*(1- theta[i,j])*(1 - 2*theta[i,j]) * size[i]
              }
              
              if(j!=k & k!=l & j!=l) {
                q[i,j,k,l] <- 2 * theta[i,j] * theta[i,k] * theta[i,l] * size[i]
              }
              
              if(j==k & j!=l){
                q[i,j,k,l] <- -theta[i,j]*(1 - 2* theta[i,j])*theta[i,l] * size[i]
              }
              
              if(j!=k & (j==l | k==l)) {
                q[i,j,k,l] <- -theta[i,j] * theta[i,k] * (1 - 2* theta[i,l]) * size[i]
              }
              
              z1 <- array(0, dim = c(J,1))
              z2 <- array(0, dim = c(J,1))
              z3 <- array(0, dim = c(J,1))
              z1[j] <- 1
              z2[k] <- 1
              z3[l] <- 1
              
              Q[i,,] <- Q[i,,] + q[i,j,k,l] * z1 %*% t(kronecker(z2, z3))
            }
          }
        }
      }
      
      QL <- array(0, dim = c(n*J, (n*J)^2))
      for(i in 1:n){
        e <- array(0, dim=c(n,1))
        e[i] <- 1
        E <- kronecker(e, diag(J))
        QL <- QL + E %*% Q[i,,] %*% t(kronecker(E,E))
      }
      
      biasE <- -1/2 * A_inv %*% (t(X)%*%QL%*% (kronecker(X, X)) %*% c(as.matrix(A_inv)))
      B_new <- B - biasE + A_inv %*% U
    }
     
    
    b_new <- matrix(B_new, nrow = J, ncol = p)
    b <- b_new
   
    if(sum(abs(B_new - B)<tol) == p*J) {
      stop <- 1
      if(bias == FALSE){
        Loglik <- sum(log(cbind(theta, 1-rowSums(theta)))*y)
      } else {
        Loglik <- sum(log(cbind(theta, 1-rowSums(theta)))*y) + 1/2 * log(det(A))
      }
    } else {
      iter <- iter + 1
    }
    
  }
  
  return(list(B = b_new, iter = iter, theta = theta, Loglik=Loglik, A= A))
}
