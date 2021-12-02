Initial <- function (II, I_G, epsilon.A, epsilon.F, D, B, NodeMod, X, M){
  # This is a function used to set inital value 
  # Args:
  #   I_G: 1*G vector, (I_1, I_2, ... , I_G), I_g is the number of samples in each group, g in 1 ... G
  #   II: cusme of I_G, (0, I_1, I_1 + I_2, ..., I_1 +...+ I_G = I)
  #   epsilon.A: threshold of anatomical measure (user defined and may obtained from emprical analysis)
  #   epsilon.F: threshold of functional measure (user defined and may obtained from emprical analysis)
  #   B: functional measure I*N
  #   D: anatomical measure I*N
  #   NodeMod: N*5, 4th and 5th element save the module the connection connecting
  #   X: group design matrix G*q
  #   M: number of module
  
  
  G <- length(I_G) # number of group
  N <- dim(B)[2] # number of connection
  q <- dim(X)[2] # number of covariates including intercept
  
  D.mean <- array(0, dim = c(G, N))
  B.mean <- array(0, dim = c(G, N)) # G*N
  A.vec <- array(0, dim = c(G, N))
  F.vec <- array(0, dim = c(G, N))
  
  for (g in 1:G) {
    D.mean[g, ] <- apply(D[(II[g] + 1):II[g + 1],], 2, mean)
    B.mean[g, ] <- apply(B[(II[g] + 1):II[g + 1],], 2, mean)
    A.vec[g, D.mean[g, ] > epsilon.A] <- 1
    F.vec[g, B.mean[g, ] > epsilon.F] <- 1
    F.vec[g, B.mean[g, ] < -epsilon.F] <- -1
  }
  
  # estimate mu and sigma
  FC.vec.A0F1<-NULL
  FC.vec.A1F1<-NULL
  FC.vec.A0F0<-NULL
  FC.vec.A1F0<-NULL
  FC.vec.A0Fn1<-NULL
  FC.vec.A1Fn1<-NULL
  for (g in 1:G) {
    for (i in (II[g]+1):II[g+1]) {
      FC.vec.A0F1<-c(FC.vec.A0F1, B[i,][intersect(which(A.vec[g,]==0), which(F.vec[g,]==1))])
      FC.vec.A1F1<-c(FC.vec.A1F1, B[i,][intersect(which(A.vec[g,]==1), which(F.vec[g,]==1))])
      FC.vec.A0F0<-c(FC.vec.A0F0, B[i,][intersect(which(A.vec[g,]==0), which(F.vec[g,]==0))])
      FC.vec.A1F0<-c(FC.vec.A1F0, B[i,][intersect(which(A.vec[g,]==1), which(F.vec[g,]==0))])
      FC.vec.A0Fn1<-c(FC.vec.A0Fn1, B[i,][intersect(which(A.vec[g,]==0), which(F.vec[g,]==-1))])
      FC.vec.A1Fn1<-c(FC.vec.A1Fn1, B[i,][intersect(which(A.vec[g,]==1), which(F.vec[g,]==-1))])
    }
  }
  
  # visualization of the observed data in initial latent group
  # par(mfrow=c(3,2))
  # hist(FC.vec.A0F1,breaks=30,xlab="transformed FC",main="A0F1")
  # hist(FC.vec.A1F1,breaks=30,xlab="transformed FC",main="A1F1")
  # hist(FC.vec.A0F0,breaks=30,xlab="transformed FC",main="A0F0")
  # hist(FC.vec.A1F0,breaks=30,xlab="transformed FC",main="A1F0")
  # hist(FC.vec.A0Fn1,breaks=30,xlab="transformed FC",main="A0Fn1")
  # hist(FC.vec.A1Fn1,breaks=30,xlab="transformed FC",main="A1Fn1")
  # par(mfrow=c(1,1))
  
  Mu.init <- matrix(c(mean(FC.vec.A0Fn1),mean(FC.vec.A0F0),mean(FC.vec.A0F1),
                      mean(FC.vec.A1Fn1),mean(FC.vec.A1F0),mean(FC.vec.A1F1)),byrow=TRUE,nrow=2)
  
  Sigma2.init <- matrix(c(var(FC.vec.A0Fn1),var(FC.vec.A0F0),var(FC.vec.A0F1),
                          var(FC.vec.A1Fn1),var(FC.vec.A1F0),var(FC.vec.A1F1)),byrow=TRUE,nrow=2)
  
  # estimate beta.a 
  SC.freq <- array(0, dim=c(G,M,M,2))
  for(s in 1:M){
    for(t in s:M){
      #cat("s=",s,"\n")
      connection <- NodeMod[NodeMod$Mod1==s&NodeMod$Mod2==t,]$Connection
      for (j in 1:2){
        for (g in 1:G) {
          sc <- NULL
          for (i in (II[g] + 1):II[g + 1]){
            sc <- c(sc, D[i,][intersect(which(A.vec[g, ]==j-1), connection)])
          }
          SC.freq[g, s, t, j] <- length(sc)
          SC.freq[g, t, s, j] <- length(sc)
        }
      }
    }
  }
  
  res <- list()
  z <- 1
  beta.a.init <- array(0,dim=c(q, M, M))
  for(s in 1:M){
    for(t in s:M){
      response <- array(SC.freq[ ,s,t, ], dim=c(g,2))
      res[[z]] <- response
      beta.a.init[,s,t] <- matrix(coef(multinom(response~X-1,family = binomial,trace=FALSE,restol=1e-20)))
      beta.a.init[,t,s] <- beta.a.init[,s,t]
      z <- z + 1
    }
  }
  
  # estimate beta_f 
  FC.freq <- array(0, dim=c(G,M,M,6))
  for(s in 1:M){
    for(t in s:M){
      #cat("s=",s,"\n")
      connection <- NodeMod[NodeMod$Mod1==s&NodeMod$Mod2==t,]$Connection
      for(j in 1:2){
        for(k in 1:3){
          l = (j-1)*3 + k
          for(g in 1:G){
            fc <- NULL
            for(i in (II[g]+1):(II[g+1])){
              fc<-c(fc, B[i,][intersect(which(A.vec[g, ]==(j-1)&F.vec[g, ]==(k-2)),connection)])
            }
            FC.freq[g, s, t, l] <- length(fc) 
            FC.freq[g, t, s, l] <- length(fc)
          } 
        }
      }
    }
  }
  
  res0 <- list()
  res1 <- list()
  z <- 1
  beta.f0.init <- array(0,dim=c(2*q ,M,M))
  beta.f1.init <- array(0,dim=c(2*q ,M,M))
  for(s in 1:M){
    for(t in s:M){
      response0 <- FC.freq[, s, t, c(2, 1, 3)]
      res0[[z]] <- response0
      beta.f0.init[, s, t] <- c(matrix(coef(multinom(response0 ~ X - 1, family = multinomial, trace = FALSE, restol = 1e-20, abstol = 1e-20, maxit = 10000)), byrow = TRUE, nrow = q))
      beta.f0.init[, t, s] <- beta.f0.init[, s, t]
      
      response1 <- FC.freq[, s, t, c(5, 4, 6)]
      res1[[z]] <- response1
      beta.f1.init[, s, t] <- c(matrix(coef(multinom(response1 ~ X - 1, family = multinomial, trace = FALSE, restol = 1e-20, abstol = 1e-20, maxit = 10000)), byrow = TRUE, nrow = q))
      beta.f1.init[, t, s] <- beta.f1.init[, s, t]
      
      z <- z + 1
    }
  }
  
  # estimate the mixture normal 
  SC.vec.0 <- NULL
  SC.vec.1 <- NULL
  for(g in 1:G){
    for(i in (II[g]+1):(II[g+1]))
      SC.vec.0 <- c(SC.vec.0, D[i,][which(A.vec[g, ]==0)])
      SC.vec.1 <- c(SC.vec.1, D[i,][which(A.vec[g, ]==1)])
  }
  SC.vec.1nz<-SC.vec.1[SC.vec.1>epsilon.A]
  SC.vec.0nz<-SC.vec.0[SC.vec.0>epsilon.A]
  
  num.table1<-matrix(c(length(SC.vec.0nz),length(SC.vec.1nz),length(SC.vec.0)-length(SC.vec.0nz),length(SC.vec.1)-length(SC.vec.1nz)), ncol=2, byrow=TRUE)
  
  Rho.init <- 1-num.table1[1, ]/colSums(num.table1)
  nme.11<-norMixEM(SC.vec.1nz, 3, maxiter = 10000)
  nme.01<-norMixEM(SC.vec.0nz, 3, maxiter = 10000)
  X2.init <- rbind(nme.01[, 1], nme.11[, 1])
  Xi2.init <- rbind(nme.01[, 2]^2, nme.11[, 2]^2)
  Gamma.init <- rbind(nme.01[, 3], nme.11[, 3])
  
  return(list(beta.a.init = beta.a.init, beta.f0.init = beta.f0.init, beta.f1.init = beta.f1.init, 
              Rho.init = Rho.init, X2.init=X2.init, Xi2.init=Xi2.init, Gamma.init=Gamma.init,
              Mu.init=Mu.init, Sigma2.init=Sigma2.init, A.init = t(A.vec), F.init = t(F.vec), res=res, res0=res0, res1=res1))
}

#############
# B <- B.obs
# D <- D.obs
# II <- II.obs
# I_G <- I_G.obs
# epsilon.A <- 3
# epsilon.F <- 0.05
#Initialization <- Initial(II = II.obs, I_G = I_G.obs, epsilon.A = 3, epsilon.F = 0.05, D = D.obs, B = B.obs, NodeMod = NodeMod, X = X.d, M = 2)