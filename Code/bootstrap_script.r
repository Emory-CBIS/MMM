# --------------------------------------------------------------------
# Load packageas and function; Setting parameters in simulation
# --------------------------------------------------------------------
source("./model_prep.r")
fit.boot.log <- list()
load("./Data/fit_1.RData")


num.group <- dim(fit.obs$X)[1]
num.var <- dim(fit.obs$Beta.a)[1]
num.sam <- 50
per.group.size <- rep(num.sam, num.group)
module.size <- c(15, 15, 19, 20, 39)
num.node <- sum(module.size)
node.module <- rep(c(1:5), module.size)

design.mat <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,1), nrow = num.group)
group.obs <- matrix(c(rep(design.mat[1, ], num.sam), rep(design.mat[2, ], num.sam), rep(design.mat[3, ], num.sam), rep(design.mat[4, ], num.sam)), nrow = 4*num.sam, byrow=TRUE)

# ---------------
# Beta parameter
# ---------------
beta.a.true <- fit.obs.log$Beta.a
beta.f0.true <- fit.obs.log$Beta.f0
beta.f1.true <- fit.obs.log$Beta.f1


# ---------------
# Anatomical Parameter
# ---------------
rho.true <- fit.obs.log$Rho
gamma.true <- fit.obs.log$Gamma
chi.true <- fit.obs.log$X2
xi2.true <- fit.obs.log$Xi2

# ---------------
# Functional Parameter
# ---------------
mu.true <- fit.obs.log$Mu
sigma2.true <- fit.obs.log$Sigma2


for (i in 1:200) { 
	set.seed(i)
	cat("Bootstrap Simulation", i, "\n")
	# ---------------------------------------------------------------------
	# Bootstrap data
	# ---------------------------------------------------------------------
	sim.data <- SimulateData(num.module = num.module, num.group = num.group, per.group.size = per.group.size, node.module = node.module, num.node=num.node, seed=i, design.mat = design.mat, beta.a = beta.a.true, beta.f0=beta.f0.true, beta.f1 = beta.f1.true, rho = rho.true, chi = chi.true, xi2 = xi2.true, mu = mu.true, sigma2 = sigma2.true, gamma = gamma.true)
	A.true <- sim.data$A 
	F.true <- sim.data$F
	B.obs <- sim.data$B
	D.obs <- sim.data$D 
	II.obs <- sim.data$II
	N_st.obs <- sim.data$N_st
	NodeMod <- sim.data$Node_Mod
	I_G.obs <- sim.data$I_G

	# ---------------------------------------------------------------------
	# Initialization
	# --------------------------------------------------------------------
	initial.setup <- InitialRandom(II = II.obs, I_G = I_G.obs, epsilon.A = 3, epsilon.F = 0.05, D = D.obs, B = B.obs, NodeMod = NodeMod, X = design.mat, M = num.module)
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

	# ----------------------------------------
	# Model Fitting
	# ----------------------------------------
	fit.boot.log[[i]] <- EM_groupmod_log_Firth(I_G = I_G.obs, II = II.obs, B=B.obs, D=D.obs,
								   Rho = rho.init, Gamma = gamma.init, X2 = chi.init, 
								   Xi2 = xi2.init, Mu = mu.init, Sigma2 = sigma2.init,
								   Beta.a = beta.a.init, Beta.f0 = beta.f0.init, Beta.f1 = beta.f1.init, X = design.mat, 
								   NodeMod = NodeMod, N_st = N_st.obs, 
								   stop = 1, tol = 1e-6, logtol = 1e-5, etol1 = 1e-3, etol2 = 1e-3, etol3 = 1e-3, maxiter=100, 
								   A=A.true, F=F.true, A.init = A.init, F.init = F.init, SIM = TRUE, posterior = FALSE)

								   
}								   
save(fit.boot.log, file=paste("../Data/fit_boot_1.RData",sep=""))

