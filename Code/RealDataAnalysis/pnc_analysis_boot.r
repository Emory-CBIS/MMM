source("./Code/model_prep.R")

# load data
mod <- read.table("./Data/roi_id.txt")
mod <- as.matrix(mod)
mod[mod == 5] <- 0
mod[mod > 5] <- mod[mod > 5] - 1
path <- "./Data/"
load(paste(path, "fitpnc.RData",sep =""))

outpath <- "./Data/Boot/boot"
num.module <- fit.obs$M 
num.group <- fit.obs$G
per.group.size <- fit.obs$I_G
num.node <- 226
node.module <- sort(mod)[39:264]
design.mat <- fit.obs$X

bootnumber <- c(1:1000)
for (i in bootnumber) {
  set.seed(i)
  cat("bootstrap",i,"\n")
  
  # —-------------------------------------------------------------------
  # Generate Bootstrap Data
  # --------------------------------------------------------------------
  sim.data <- SimulateData(num.module = num.module, num.group = num.group, per.group.size = per.group.size, node.module = node.module, num.node=num.node, seed=i, design.mat = design.mat, beta.a = fit.obs$Beta.a, beta.f0=fit.obs$Beta.f0, beta.f1 = fit.obs$Beta.f1, rho = fit.obs$Rho, chi = fit.obs$X2, xi2 = fit.obs$Xi2, mu = fit.obs$Mu, sigma2 = fit.obs$Sigma2, gamma = fit.obs$Gamma)
  A.true <- sim.data$A 
  F.true <- sim.data$F
  B.boot <- sim.data$B
  D.boot <- sim.data$D 
  II.obs <- sim.data$II
  N_st.obs <- sim.data$N_st
  NodeMod <- sim.data$Node_Mod
  I_G.obs <- sim.data$I_G
  
  # —-------------------------------------------------------------------
  # Initialization
  # --------------------------------------------------------------------
  initial.setup <- InitialRandom(II = II.obs, I_G = I_G.obs, epsilon.A = 2.5, epsilon.F = 0.1, D = D.boot, B = B.boot, NodeMod = NodeMod, X = design.mat, M = num.module)
  beta.a.init = initial.setup$beta.a.init
  beta.f0.init = initial.setup$beta.f0.init
  beta.f1.init = initial.setup$beta.f1.init 
  rho.init = initial.setup$Rho.init
  chi.init = initial.setup$X2.init
  xi2.init = initial.setup$Xi2.init
  gamma.init = initial.setup$Gamma.init
  mu.init = initial.setup$Mu.init
  sigma2.init = initial.setup$Sigma2.init
  A.init = initial.setup$A.init
  F.init = initial.setup$F.init
  
  
  # —-------------------------------------------------------------------
  # Fit the Model
  # --------------------------------------------------------------------
  
  
  modelfit <- function() {
    out <- tryCatch(EM_groupmod_log_Firth(I_G = I_G.obs , II = II.obs , B=B.boot, D=D.boot,
                                                      Rho = rho.init, Gamma = gamma.init, X2 = chi.init, 
                                                      Xi2 = xi2.init, Mu = mu.init, Sigma2 = sigma2.init,
                                                      Beta.a = beta.a.init, Beta.f0 = beta.f0.init, Beta.f1 = beta.f1.init, X = design.mat, 
                                                      NodeMod = NodeMod, N_st = N_st.obs, 
                                                      stop = 1, tol = 1e-6, logtol = 1e-5, etol1 = 1e-3, etol2 = 1e-3, etol3 = 1e-3, maxiter=100, 
                                                      SIM=FALSE, posterior = FALSE),
                    error = function(e) NULL)
    return(out)
  }
  fit.boot <- modelfit()
  # save data
  # save(fit.boot, file= paste(outpath,i,".RData", sep=""))
}

