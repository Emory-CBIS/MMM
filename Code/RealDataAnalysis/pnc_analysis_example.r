# load library and useful functions
source("./Code/model_prep.R")
# load data 
load("./Data/PNCnoise.RData")


# -----------------------------------------------------------------
# Preprocess Raw Data for the model
# -----------------------------------------------------------------
# set parameter for initial value
epsilon.A <- 2.5 # parameter to estimate initial latent structural states
epsilon.F <- 0.1 # parameter to estimate initial latent functional states
seed <- 123     # random seed
thres <- 10^(-5) # threshold to determine the point mass function
tt <- 10^(-6) # replacement for those values < 1e-5

# design matrix  first column is intercept, second column is gender and third column is age
X <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,1), nrow=4)

# read FC and SC data
N <- length(fmri.data)
Q <- length(mod) - sum(mod==0)
FC <- array(0, dim = c(Q, Q, N))
SC <- array(0, dim = c(Q, Q, N))
for (i in 1:N) {
  FC[, ,i] <- fmri.data[[i]]
  SC[, ,i] <- dmri.data[[i]]
}

# group patients
age.cut <- c(8, 16, 22)
demographic$Sexind <- as.numeric(demographic$Sex == "M") 
n <- length(age.cut)
demographic$Ageind <- demographic$Age
for (i in 1:(n-1)) {
  demographic$Ageind[demographic$Age>=age.cut[i]&demographic$Age<age.cut[i+1]] <- i - 1
}

# remove NA value
na.idx <- NULL
for(i in 1:N) {
  if (sum(is.na(SC[,,i])!=0)) {
    na.idx <- c(na.idx, i)
  }
}

remove.idx <- sort(c(115, na.idx))
SC <- SC[,,-na.idx]
FC <- FC[, ,-na.idx]
N <- N - length(na.idx)
demographic.new <- demographic[-na.idx,]

# calculte Group number
G <- (length(age.cut) - 1) * 2
demographic.new$Group <- demographic.new$Sexind*(length(age.cut)-1) + demographic.new$Ageind + 1

# set variables for function
logit <- function(x) {
  return(log(x/(1-x)))
}

# deal with zero probility
cat("the proportion of zero probility",length(SC[SC<thres])/length(SC),"\n") 
SC[SC<thres] <- tt 
SC_tran <- logit(SC) - logit(tt)


D <- array(0,c(N, Q*(Q-1)/2))
for(i in 1:N){
  D[i,] <- mat2vec(SC_tran[,,i])
}

FC_tran <- atanh(FC)
B <- array(0, c(N, Q*(Q-1)/2))
for(i in 1:N){
  B[i,] <- mat2vec(FC_tran[,,i])
}

B <- B[order(demographic.new$Group),]
D <- D[order(demographic.new$Group),]

M <- 9
M_st <- sort(mod)[39:264]
NodeMod_obs <- NodeMod(N=Q*(Q-1)/2,Q=Q,M_st = M_st)
I_G_obs <- as.numeric(table(demographic.new$Group))
II_obs <- c(0, cumsum(I_G_obs))


II <- II_obs; I_G <- I_G_obs; NodeMod <- NodeMod_obs
N_st <- matrix(0, nrow=M, ncol=M)  # number of connections in each module
for (i in 1:M) {
  for (j in 1:M) {
    N_st[i,j] <- length(which(NodeMod$Mod1==i& NodeMod$Mod2==j))
  }
}

### Initailzation ###
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

Mu.init <- matrix(c(mean(FC.vec.A0Fn1),mean(FC.vec.A0F0),mean(FC.vec.A0F1),
                    mean(FC.vec.A1Fn1),mean(FC.vec.A1F0),mean(FC.vec.A1F1)),byrow=TRUE,nrow=2)

Sigma2.init <- matrix(c(var(FC.vec.A0Fn1),var(FC.vec.A0F0),var(FC.vec.A0F1),
                        var(FC.vec.A1Fn1),var(FC.vec.A1F0),var(FC.vec.A1F1)),byrow=TRUE,nrow=2)

# random intial Beta
set.seed(seed)
Beta.a.init = array(rnorm(q*M*M),dim=c(q, M, M))
Beta.f0.init = array(rnorm(2*q*M*M),dim=c(2*q, M, M))
Beta.f1.init = array(rnorm(2*q*M*M),dim=c(2*q, M, M))

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
A.init <- t(A.vec)
F.init <- t(F.vec)

# fit MMM
time <- proc.time()
fit.obs <- EM_groupmod_log_Firth(I_G = I_G, II = II, B=B, D=D,
                                 Rho = Rho.init, Gamma = Gamma.init, X2 = X2.init, 
                                 Xi2 = Xi2.init, Mu = Mu.init, Sigma2 = Sigma2.init,
                                 Beta.a = Beta.a.init, Beta.f0 = Beta.f0.init, Beta.f1 = Beta.f1.init, X = X, 
                                 NodeMod = NodeMod, N_st = N_st, 
                                 stop = 1, tol = 1e-6, logtol = 1e-5, etol1 = 1e-3, etol2 = 1e-3, etol3 = 1e-3, maxiter=200,
                                 SIM=FALSE, posterior=FALSE)
proc.time() - time

# save fitted model
#save(fit.obs, file=paste("./Data/fitpnc.RData",sep=""))


