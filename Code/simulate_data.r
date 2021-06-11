NodeModule <- function(N, Q, M_st){
  # find two modules one connection connects
  # 
  # Args:
  #   N: number of connections
  #   Q: number of Nodes
  #   M_st: module of each Node
  # 
  # Return:
  #   A dataframe saves information of each connection
  
  Node_Mod <- data.frame(Connection=c(1:N))
  a<-NULL
  b<-NULL
  for (i in 1:Q) {
    a <- c(a, rep(i, Q-i))
    b <- c(b, c((i+1): Q))  
  }
  b <- b[1:N]
  Node_Mod$Node1 <- a
  Node_Mod$Node2 <- b
  rm(a)
  rm(b)
  for (i in 1:N) {
    node1 <- Node_Mod$Node1[i]
    node2 <- Node_Mod$Node2[i]
    Node_Mod$Mod1[i] <- M_st[node1]
    Node_Mod$Mod2[i] <- M_st[node2]
  }
  return(Node_Mod)
}

# Simulation Function
SimulateData <- function(num.module = 5, num.group = 2, per.group.size = rep(20,2), node.module = c(rep(c(1:5),rep(40,5))), num.node=200, seed=123, design.mat, beta.a, beta.f0, beta.f1, rho, chi, xi2, mu, sigma2, gamma){
  # set seed
  set.seed(seed)
  
  
  num.connection <- num.node*(num.node - 1)/2 # Number of Connection
  node.module <-  NodeModule(num.connection, num.node,  node.module)  # Node and Mode correspond to each other
  
  num.connection.module <- matrix(0, nrow=num.module, ncol=num.module)  # number of connections in each module
  for (i in 1:num.module) {
    for (j in 1:num.module) {
      num.connection.module[i,j] <- length(which(node.module$Mod1==i& node.module$Mod2==j))
    }
  }
  
  size <- sum(per.group.size)
  cum.size <-cumsum(per.group.size)
  cum.size <-c(0, cum.size)
  
  # calculate probablity 
  prob <- Calc_Prob(M = num.module, G = num.group, beta.a, beta.f0, beta.f1, design.mat)
  pi.a <- prob$pi.a
  pi.f1 <- prob$pi.f1
  pi.f0 <- prob$pi.f0
  pi.f_1 <- prob$pi.f_1
  
  # Anatomical hidden states 
  A <- matrix(0, ncol = num.group, nrow = num.connection)
  for (n in 1:num.connection) {
    for (g in 1:num.group) {
      s <- node.module[n,"Mod1"]
      t <- node.module[n,"Mod2"]
      A[n,g] <- rbinom(1,1,pi.a[g, s, t])
    }
  }
  
  # Functional hidden states 
  F <- matrix(0, ncol = num.group, nrow = num.connection)
  for(n in 1:num.connection) {
    s <- node.module[n,"Mod1"]
    t <- node.module[n,"Mod2"]
    for (g in 1:num.group) {
      if (A[n, g] == 0) {
        F[n, g] <- which(rmultinom(1, 1, prob=c(pi.f_1[g,s,t,1], pi.f0[g,s,t,1], pi.f1[g,s,t,1]))==1)-2
      } else {
        F[n, g]<- which(rmultinom(1, 1, prob=c(pi.f_1[g,s,t,2], pi.f0[g,s,t,2], pi.f1[g,s,t,2]))==1)-2
      }
    }
  }
  
  # Anatomial Measure 
  D <- array(0, dim = c(size, num.connection))
  for (n in 1 : num.connection) {
    for (g in 1 : num.group){
      for (i in (cum.size[g]+1) : (cum.size[g+1])) {
        if (A[n, g] == 0) {
          r1 <- rbinom(1, 1, rho[1])  # r1: random number 1 
          if (r1 == 1) {
            D[i, n] <- 0
          } else {
            r2 <- which(rmultinom(1, 1, gamma[1, ] ) == 1)  # r2: random number 2
            D[i, n] <- rnorm(1, chi[1, r2], sqrt(xi2[1, r2]))
          }
        } else {
          r1 <- rbinom(1, 1, rho[2])
          if(r1==1){
            D[i, n] <- 0
          } else {
            r2 <- which(rmultinom(1, 1, gamma[2, ]) == 1)
            D[i, n] <- rnorm(1, chi[2, r2], sqrt(xi2[2, r2]))
          }
        }
      }
    }
  }
  
  # Functional Measure 
  B <- array(0, dim = c(size, num.connection))
  for (n in 1 : num.connection) {
    for (j in 1:2) {
      for (k in 1:3) {
        for (g in 1:num.group) {
          if(A[n, g] == (j-1) & F[n, g] == (k-2)){
            for(i in (cum.size[g] + 1) : cum.size[g + 1]){
              B[i, n] <- rnorm(1, mu[j, k], sqrt(sigma2[j, k]))
            }
          }
        }
      }
    }
  }
  
  return(list(B = B, D = D, II = cum.size, N_st = num.connection.module, Node_Mod = node.module, A = A, F = F, I_G = per.group.size, pi.a=pi.a, pi.f0=pi.f0, pi.f1=pi.f1, pi.f_1=pi.f_1))
}



# ##################################################
# #
# #     check simulated data                       
# #
# ##################################################
# D <- Data$D
# A <- Data$A
# F <- Data$F
# B <- Data$B
# D11 = c(D[1:20,which(A[,1]==1)], D[21:40,which(A[,2]==1)])
# D00 = c(D[1:20,which(A[,1]==0)], D[21:40,which(A[,2]==0)])
# 
# hist(D00, breaks=100)
# hist(D11, breaks=100)
# 
# B11 = c(B[1:20,which(A[,1]==1 & F[,1]==1)],B[21:40,which(A[,2]==1 & F[,2]==1)])
# B10 = c(B[1:20,which(A[,1]==1 & F[,1]==0)],B[21:40,which(A[,2]==1 & F[,2]==0)])
# B1_1 = c(B[1:20,which(A[,1]==1 & F[,1]==-1)],B[21:40,which(A[,2]==1 & F[,2]==-1)])
# B01 = c(B[1:20,which(A[,1]==0 & F[,1]==1)],B[21:40,which(A[,2]==0 & F[,2]==1)])
# B00 = c(B[1:20,which(A[,1]==0 & F[,1]==0)],B[21:40,which(A[,2]==0 & F[,2]==0)])
# B0_1 = c(B[1:20,which(A[,1]==0 & F[,1]==-1)],B[21:40,which(A[,2]==0 & F[,2]==-1)])
# 
# mean(B11);var(B11);
# mean(B10);var(B10);
# mean(B1_1);var(B1_1);
# mean(B01);var(B01);
# mean(B00);var(B00);
# mean(B0_1);var(B0_1);
# 
# ### simulation of checking ###
# ### there is randomness in the data ###
# 
# # coef1 <- array(0, dim=c(100,2))
# # coef2 <- array(0, dim=c(100,2))
# # for(i in 1:100){
# # A<-matrix(0,ncol=G,nrow=N)
# # print(i)
# # for(n in 1:N){
# #  for(g in 1:G){
# #    s = Node_Mod[n,"Mod1"]
# #    t = Node_Mod[n,"Mod2"]
# #    A[n,g]<-rbinom(1,1,pi_a[g,s,t])
# #  }
# # }
# # t1<-table(A[Node_Mod$Mod1==1&Node_Mod$Mod2==1,1])
# # t2<-table(A[Node_Mod$Mod1==1&Node_Mod$Mod2==1,2])
# # coef1[i,] <- -coef(glm(rbind(t1,t2)~X-1,family = binomial))
# # t1<-table(A[Node_Mod$Mod1==1&Node_Mod$Mod2==2,1])
# # t2<-table(A[Node_Mod$Mod1==1&Node_Mod$Mod2==2,2])
# # coef2[i,] <- -coef(glm(rbind(t1,t2)~X-1,family = binomial))
# # }
# # 
# 
# 
# 
# pi_A = list()
# pi_A[[1]] = pi_f1[,,,1]
# pi_A[[2]] = pi_f_1[,,,1]
# pi_A[[3]] = pi_f0[,,,1]
# 
# 
# 
# # idx <- Node_Mod$Mod1 == 1 &Node_Mod$Mod2 == 1
# # YD_1 <- table(A[idx,1],F[idx,1])
# # YD_2 <- table(A[idx,2],F[idx,2])
# # Response1 = YD_1[2,c(2,1,3)]
# # Response2 = YD_2[2,c(2,1,3)]
# # Response = rbind(Response1, Response2)
# # (matrix(coef(multinom(Response~c(0,1),family = multinomial)),byrow=TRUE))[c(1,3,2,4)]
# # 
# # response = matrix(c(39,30.50,1550,4220,4230),byrow=TRUE,nrow=2)
# # multinom(response~X-1,family = multinomial)