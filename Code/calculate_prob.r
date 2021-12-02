# This is a function to deal with the numerical issue of logSum
Raddlog <- function (a, b)
{
  result <- rep(0, length(a))
  idx1 <- a > b + 200  # whether a is 200 larger than b
  result[idx1] <- a[idx1] # if a is much larger than b
  idx2 <- b > a + 200  # whether b is 200 larger than b
  result[idx2] <- b[idx2] # if b is much larger than a
  idx0 <- !(idx1 | idx2) 
  result[idx0] <- a[idx0] + log1p(exp(b[idx0] - a[idx0]))
  result
}

Raddlog3 <- function (a, b, c)
{
  result <- Raddlog(Raddlog(a, b), c)
  result
}

# This is a Function used to calculate the probabilities in log scale
Calc_Prob_Log <- function(M, G, Beta.a, Beta.f0, Beta.f1, X){
  # Args:
  #   N: Number of Connections
  #   M: number of module
  #   G: number of Group
  #   Beta_a:  (q) * M * M ; include intercept
  #   Beta_f0: (2q) * M * M: first (q) for -1 and second (q) for 1  when anatomical = 0
  #   Beta_f1: (2q) * M * M: first (q)  when anatomial = 1
  #   X: group design matrix G*q
  
  q <- dim(Beta.a)[1]
  # probability
  pi.a <- array(1, dim=c(G, M, M))       #
  pi.a0 <- array(1, dim=c(G, M, M)) 
  pi.f1 <- array(0, dim=c(G, M, M, 2))   # F==1
  pi.f_1 <- array(0, dim=c(G, M, M, 2))  # F==-1
  pi.f0 <- array(0, dim=c(G, M, M, 2))   # F==0
  
  for (s in 1:M) {
    for (t in 1:M) {
      for (g in 1:G) {
        b <- sum(Beta.a[, s, t] * X[g, ])
        pi.a[g, s, t] <- b - Raddlog(0, b)
        pi.a0[g, s, t] <-  - Raddlog(0, b)
        
        pi.a[g, t, s] <- pi.a[g, s, t]
        pi.a0[g, t, s] <- pi.a0[g, s, t]
      }
    }
  }
  
  for (s in 1:M) {
    for (t in 1:M) {
      for(g in 1:G){
        b1 <- sum(Beta.f0[(q+1):(2*q),s,t]*X[g,])
        b2 <- sum(Beta.f0[1:q,s,t]*X[g,])
        
        denom <- Raddlog3(0, b1, b2)
        
        pi.f1[g, s, t, 1] <- b1 - denom 
        pi.f_1[g, s, t, 1] <- b2 - denom
        pi.f0[g, s, t, 1] <- - denom
        
        pi.f1[g, t, s, 1] <- pi.f1[g, s, t, 1]
        pi.f_1[g, t, s, 1] <- pi.f_1[g, s, t, 1]
        pi.f0[g, t, s, 1] <- pi.f0[g, s, t, 1]
        
        b1 <- sum(Beta.f1[(q+1):(2*q),s,t]*X[g,])
        b2 <- sum(Beta.f1[1:q,s,t]*X[g,])
        
        denom <- Raddlog3(0, b1, b2)
        pi.f1[g, s, t, 2] <- b1 - denom 
        pi.f_1[g, s, t, 2] <- b2 - denom
        pi.f0[g, s, t, 2] <- - denom
        
        pi.f1[g, t, s, 2] <- pi.f1[g, s, t, 2]
        pi.f_1[g, t, s, 2] <- pi.f_1[g, s, t, 2]
        pi.f0[g, t, s, 2] <- pi.f0[g, s, t, 2]
      }
    }
  }
  return(list(pi.a = pi.a, pi.f1=pi.f1, pi.f0=pi.f0, pi.f_1=pi.f_1, pi.a0 = pi.a0))
}


Calc_Prob <- function(M, G, Beta.a, Beta.f0, Beta.f1, X){
  # Function used to calculate the probablities
  # Args:
  #   N: Number of Connections
  #   M: number of module
  #   G: number of Group
  #   Beta_a:  (q) * M * M ; include intercept
  #   Beta_f0: (2q) * M * M: first (q) for -1 and second (q) for 1  when anatomical = 0
  #   Beta_f1: (2q) * M * M: first (q)  when anatomial = 1
  #   X: group design matrix G*q
  
  # dimension of beta parameter
  q <- dim(Beta.a)[1]
  
  # probability
  pi.a <- array(1, dim=c(G, M, M))       # 
  pi.f1 <- array(0, dim=c(G, M, M, 2))   # F==1
  pi.f_1 <- array(0, dim=c(G, M, M, 2))  # F==-1
  pi.f0 <- array(0, dim=c(G, M, M, 2))   # F==0
  
  for (s in 1:M) {
    for (t in 1:M) {
      for (g in 1:G) {
        if(exp(sum(Beta.a[, s, t] * X[g, ])) == Inf) {
          pi.a[g, s, t] <- 1
		      pi.a[g, t, s] <- pi.a[g, s, t]
        } else{
          pi.a[g, s, t] <- exp(sum(Beta.a[, s, t] * X[g, ]))/(1 + exp(sum(Beta.a[, s, t] * X[g, ])))
          pi.a[g, t, s] <- pi.a[g, s, t]
        }
      }
    }
  }
  
  for (s in 1:M) {
    for (t in 1:M) {
      for(g in 1:G){
        
        denom <- (1 + exp(sum(Beta.f0[(q+1):(2*q),s,t]*X[g,])) + exp(sum(Beta.f0[1:q,s,t]*X[g,])))
        if (denom == Inf) {
          if (exp(sum(Beta.f0[(q+1):(2*q),s,t]*X[g,]))==Inf) {
            pi.f1[g, s, t, 1] <- 1 
            pi.f_1[g, s, t, 1] <- 0
            pi.f0[g, s, t, 1] <- 0
          }
          if (exp(sum(Beta.f0[1:q,s,t]*X[g,]))==Inf) {
            pi.f1[g, s, t, 1] <- 0 
            pi.f_1[g, s, t, 1] <- 1
            pi.f0[g, s, t, 1] <- 0
          }
        } else{
          pi.f1[g, s, t, 1] <- exp(sum(Beta.f0[(q+1):(2*q),s,t]*X[g,])) / denom 
          pi.f_1[g, s, t, 1] <- exp(sum(Beta.f0[1:q,s,t]*X[g,])) / denom
          pi.f0[g, s, t, 1] <- 1 / denom
        }
        
        pi.f1[g, t, s, 1] <- pi.f1[g, s, t, 1]
        pi.f_1[g, t, s, 1] <- pi.f_1[g, s, t, 1]
        pi.f0[g, t, s, 1] <- pi.f0[g, s, t, 1]
      
        denom <- (1+exp(sum(Beta.f1[(q+1):(2*q),s,t]*X[g,]))+exp(sum(Beta.f1[1:q,s,t]*X[g,])))
        if (denom == Inf) {
          if (exp(sum(Beta.f1[(q+1):(2*q),s,t]*X[g,]))==Inf) {
            pi.f1[g, s, t, 2] <- 1 
            pi.f_1[g, s, t, 2] <- 0
            pi.f0[g, s, t, 2] <- 0
          }
          if (exp(sum(Beta.f1[1:q,s,t]*X[g,]))==Inf) {
            pi.f1[g, s, t, 2] <- 0 
            pi.f_1[g, s, t, 2] <- 1
            pi.f0[g, s, t, 2] <- 0
          }
        } else {
          pi.f1[g, s, t, 2] <- exp(sum(Beta.f1[(q+1):(2*q),s,t]*X[g,])) / denom 
          pi.f_1[g, s, t, 2] <- exp(sum(Beta.f1[1:q,s,t]*X[g,])) / denom
          pi.f0[g, s, t, 2] <- 1 / denom
        }
        
        pi.f1[g, t, s, 2] <- pi.f1[g, s, t, 2]
        pi.f_1[g, t, s, 2] <- pi.f_1[g, s, t, 2]
        pi.f0[g, t, s, 2] <- pi.f0[g, s, t, 2]
      }
    }
  }
  return(list(pi.a = pi.a, pi.f1=pi.f1, pi.f0=pi.f0, pi.f_1=pi.f_1))
}
