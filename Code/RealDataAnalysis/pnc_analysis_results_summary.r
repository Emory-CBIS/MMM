# ------------------------
# This is a script used to summarize bootstrap results to make inference
# ------------------------
path <- "./Code/"
source(paste(path, "model_prep.r", sep = ""))
source(paste(path, "./RealDataAnalysis/analysis_prep.r", sep = ""))
output_folder <- "./Figure/"
path <- "./Data/"

# ------------------------
# load fit results
# ------------------------
load(paste(path, "fitpnc.RData",sep=""))
  
# ------------------------
# load bootstrap samples
# ------------------------
boot <- list()
n.boot <- 1000
for (i in 1:n.boot) {
  filename <- paste(path, "Boot/boot", i, ".RData",sep="")
    if (file.exists(filename)){
      load(filename)
      boot[[i]] <- fit.boot
    }
}

# ------------------------------
# filter out unqualified samples
# ------------------------------
idx.exist <- sapply(boot, function(x) length(x)!=0)
boot <- boot[idx.exist] 
check.convergence <- NULL
for (i in 1:length(boot)) {
  if (boot[[i]]$NUM != 0) {
    check.convergence[i] <- FALSE
  } else {
    check.convergence[i] <- TRUE
  }
}
boot <- boot[check.convergence==TRUE]
n.boot <- length(boot)
cat("##############", "\n")
cat("fit results and boostrap results loading finished", "\n")
cat("##############", "\n")
  
# ------------------------------
# obtain bootstrap estimation
# ------------------------------
group.idx2 <- 2; group.idx1 <- 1  # focus on the analysis of Old Female (group idx 2) v.s. Young Female (group idx 1)
q <- dim(fit.obs$X)[2]
M <- fit.obs$M
G <- fit.obs$G

# ---------------------
# Structural Results
# ---------------------
# observed data
prob.obs <- Calc_Prob(M=fit.obs$M, G=fit.obs$G, Beta.a=fit.obs$Beta.a, Beta.f0=fit.obs$Beta.f0, Beta.f1=fit.obs$Beta.f1, X=fit.obs$X)


# bootstrap data
prob.final.boot <- list()
mat.diff.obs <- PlotPiAlphaSplit(prob.obs$pi.a, G=group.idx2, M=fit$M, n=4, fig=FALSE, order=TRUE) - 
                PlotPiAlphaSplit(prob.obs$pi.a, G=group.idx1, M=fit$M, n=4, fig=FALSE, order=TRUE)
mat1.boot <-list()
mat2.boot <- list()
mat.diff.boot <-list()
l <- c(1: length(boot))
for(iter in l) {
  prob.final.boot[[iter]] <- Calc_Prob(M=fit.obs$M, G=fit.obs$G, Beta.a=boot[[iter]]$Beta.a, Beta.f0=boot[[iter]]$Beta.f0, Beta.f1=boot[[iter]]$Beta.f1, X=fit.obs$X)
  mat1.boot[[iter]] <- PlotPiAlphaSplit(prob.final.boot[[iter]]$pi.a, G=group.idx1, M=fit$M, n=2, fig=FALSE,order=TRUE)
  mat2.boot[[iter]] <- PlotPiAlphaSplit(prob.final.boot[[iter]]$pi.a, G=group.idx2, M=fit$M, n=2, fig=FALSE,order=TRUE)
  mat.diff.boot[[iter]] <- mat2.boot[[iter]] - mat1.boot[[iter]]
}
boot.diff.a <- array(0, dim = c(length(boot), M^2)) 
for (mod1 in 1:M) {
  for (mod2 in 1:M) {
    for (i in 1:length(boot)) {
      boot.diff.a[i, mod1 + M*(mod2 - 1)] <- mat.diff.boot[[i]][mod1, mod2]
    }
  }
}
  
# calculate p-value
alpha <- 0.05
p.diff.a <- round(2 * pmin(pnorm(matrix(mat.diff.obs, byrow = TRUE, nrow = 1)/apply(boot.diff.a, 2, function(x) sd(x))), pnorm(matrix(mat.diff.obs, byrow = TRUE, nrow = 1)/apply(boot.diff.a, 2, function(x) sd(x)), lower.tail = FALSE)),3)
p.diff.a <- matrix(p.diff.a, nrow = 9)
signif.diff.a <- p.diff.a < alpha

# visualization
# reproduce Figure 2
pdf(file = paste(output_folder, "signif_anatomical_diff_","G", group.idx2, " & G", group.idx1, ".pdf", sep=""))
mycol <-  colorRamp2( c(1 - alpha, 1), c("whitesmoke","yellow"))
Heatmap9by9.signif(mat.diff.obs, col = mycol, signif = signif.diff.a, pmat = p.diff.a)
dev.off()

# reproduce Figure 3
pdf(file = paste(output_folder, "legend_anantomical_diff.pdf"))
lgd = Legend(col_fun = colorRamp2(c(0, alpha), c("yellow","whitesmoke")), title = "p-value", at = c(0, alpha))
draw(lgd)
dev.off()
mat.diff.obs.a <- mat.diff.obs 
p.diff.a <- p.diff.a
circlegraph(mat.diff.obs.a,  p.diff.a, filename = paste(output_folder, "Anatomical Diff.pdf"), alpha  = 0.05, factor = -2, color.limit = 0.1)
  
# ---------------------
# Conditional Functional Results
# ---------------------

# observe data
m10 <-list()
m11 <- list()
m20 <- list()
m21 <- list()
m0.diff.obs <- PlotPiFSplit(prob.obs$pi.f_1, prob.obs$pi.f0, prob.obs$pi.f1, G=group.idx2, M=fit.obs$M, n=2, combine=FALSE, fig=1, figure=FALSE, order=TRUE)-PlotPiFSplit(prob.obs$pi.f_1, prob.obs$pi.f0, prob.obs$pi.f1, G=group.idx1, M=fit.obs$M, n=2, combine=FALSE, fig=1, figure=FALSE, order=TRUE)
m1.diff.obs <- PlotPiFSplit(prob.obs$pi.f_1, prob.obs$pi.f0, prob.obs$pi.f1, G=group.idx2, M=fit.obs$M, n=2, combine=FALSE, fig=2, figure=FALSE, order=TRUE)-PlotPiFSplit(prob.obs$pi.f_1, prob.obs$pi.f0, prob.obs$pi.f1, G=group.idx1, M=fit.obs$M, n=2, combine=FALSE, fig=2, figure=FALSE, order=TRUE)
for(i in l) {
  m10[[i]] <- PlotPiFSplit(prob.final.boot[[i]]$pi.f_1, prob.final.boot[[i]]$pi.f0, prob.final.boot[[i]]$pi.f1, G=group.idx1, M=fit.obs$M, n=2, combine=FALSE, fig=1, figure=FALSE, order=TRUE)
  m11[[i]] <- PlotPiFSplit(prob.final.boot[[i]]$pi.f_1, prob.final.boot[[i]]$pi.f0, prob.final.boot[[i]]$pi.f1, G=group.idx1, M=fit.obs$M, n=2, combine=FALSE, fig=2, figure=FALSE, order=TRUE)
  m20[[i]]<- PlotPiFSplit(prob.final.boot[[i]]$pi.f_1, prob.final.boot[[i]]$pi.f0, prob.final.boot[[i]]$pi.f1, G=group.idx2, M=fit.obs$M, n=2, combine=FALSE, fig=1, figure=FALSE, order=TRUE)
  m21[[i]] <- PlotPiFSplit(prob.final.boot[[i]]$pi.f_1, prob.final.boot[[i]]$pi.f0, prob.final.boot[[i]]$pi.f1, G=group.idx2, M=fit.obs$M, n=2, combine=FALSE, fig=2, figure=FALSE, order=TRUE)
}
  
# bootstrap data
m0.diff.boot <- list()
m1.diff.boot <- list()
for (i in l) {
  m0.diff.boot[[i]] <- m20[[i]] - m10[[i]]
  m1.diff.boot[[i]] <- m21[[i]] - m11[[i]]
}
boot.diff.f0 <- list()
boot.diff.f1 <- list()
for (state in 1:3) {
  boot.diff.f0[[state]] <- array(0, dim = c(length(boot), M^2)) 
  boot.diff.f1[[state]] <- array(0, dim = c(length(boot), M^2)) 
  for (mod1 in 1:fit.obs$M) {
    for (mod2 in 1:fit.obs$M) {
      for (i in 1:length(boot)) {
        boot.sub.f0 <- m0.diff.boot[[i]][, c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))]
        boot.sub.f1 <- m1.diff.boot[[i]][, c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))]
        boot.diff.f0[[state]][i, mod1 + fit.obs$M*(mod2 - 1)] <- boot.sub.f0[mod1, mod2]
        boot.diff.f1[[state]][i, mod1 + fit.obs$M*(mod2 - 1)] <- boot.sub.f1[mod1, mod2]
      }
    }
  }
}

# calculate significant level  
boot.diff.f0.sd <- sapply(boot.diff.f0, function(x) apply(x, 2, sd))
boot.diff.f1.sd <- sapply(boot.diff.f1, function(x) apply(x, 2, sd))
p.diff.f0.boot <- list()
p.diff.f1.boot <- list()
for(state in 1:3){
  p.diff.f0.boot[[state]] <- 2 * pmin(pnorm(matrix(m0.diff.obs[,c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))], nrow =1, byrow = TRUE)/boot.diff.f0.sd[, state]), pnorm(matrix(m0.diff.obs[,c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))], nrow =1, byrow = TRUE)/boot.diff.f0.sd[, state], lower.tail = FALSE))
  p.diff.f1.boot[[state]] <- 2 * pmin(pnorm(matrix(m1.diff.obs[,c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))], nrow =1, byrow = TRUE)/boot.diff.f1.sd[, state]), pnorm(matrix(m1.diff.obs[,c((fit.obs$M*(state-1)  + 1):(fit.obs$M*state))], nrow =1, byrow = TRUE)/boot.diff.f1.sd[, state], lower.tail = FALSE))
}
  
  
# reproduce Figure 6
alpha.cf <- 0.01
for(state in 1:3) {
    signif.f0 <- p.diff.f0.boot.combine[,((state-1)*M+1):(state*M)] < alpha.cf
    signif.f1 <- p.diff.f1.boot.combine[,((state-1)*M+1):(state*M)] < alpha.cf
    mycol <-  colorRamp2( c(1 - alpha.cf, 1), c("grey90","yellow"))
    
    png(filename = paste(output_folder, "sigif_conditional_diff0_G", group.idx2, " & G", group.idx1,"_state",state,"_", alpha.cf, ".png", sep=""), width = 600, height = 600)
    Heatmap9by9.signif(m0.diff.obs[,((state-1)*M+1):(state*M)], col = mycol, signif = signif.f0, pmat = p.diff.f0.boot.combine[,((state-1)*M+1):(state*M)])
    dev.off()
    circlegraph(m0.diff.obs[,((state-1)*M+1):(state*M)],  p.diff.f0.boot.combine[,((state-1)*M+1):(state*M)], 
                filename = paste(output_folder, "Conditional Functional Diff 0 state ", state - 2, ".pdf", sep = ""), alpha  = 0.01, factor = -1, color.limit = 0.1)
    
    png(filename = paste(output_folder, "sigif_conditional_diff1_G", group.idx2, " & G", group.idx1, "_state",state,"_", alpha.cf, ".png", sep=""), width = 600, height = 600)
    Heatmap9by9.signif(m1.diff.obs[,((state-1)*M+1):(state*M)], col = mycol, signif = signif.f1, pmat = p.diff.f1.boot.combine[,((state-1)*M+1):(state*M)])
    dev.off()
    circlegraph(m1.diff.obs[,((state-1)*M+1):(state*M)],  p.diff.f1.boot.combine[,((state-1)*M+1):(state*M)], 
                filename = paste(output_folder, "Conditional Functional Diff 1 state ", state - 2, ".pdf", sep = ""), alpha  = 0.01, factor = -1, color.limit = 0.1)
}
  

# -------------------------
# Marginal Functional results
# -------------------------
# observed data results
piAlpha.obs <- PlotPiAlpha(prob.obs$pi.a, G=G, M=M, figure=FALSE)
piF.obs <- PlotPiF(p_1 = prob.obs$pi.f_1,p0 = prob.obs$pi.f0, p1 = prob.obs$pi.f1, G = G, M=M, figure=FALSE)
piF0_1.obs <- piF.obs[ ,c(1:M)]
piF1_1.obs <- piF.obs[ ,c((3*M+1):(4*M))]
piF00.obs <- piF.obs[ ,c((M+1):(2*M))]
piF10.obs <- piF.obs[ ,c((4*M+1):(5*M))]
piF01.obs <- piF.obs[ ,c((2*M+1):(3*M))]
piF11.obs <- piF.obs[ ,c((5*M+1):(6*M))]
margProb_1.obs <- (piAlpha.obs*piF1_1.obs + (1-piAlpha.obs)*piF0_1.obs)
margProb0.obs <- (piAlpha.obs*piF10.obs + (1-piAlpha.obs)*piF00.obs) 
margProb1.obs <- (piAlpha.obs*piF11.obs + (1-piAlpha.obs)*piF01.obs)
m_1.diff.obs <- margProb_1.obs[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb_1.obs[((group.idx1-1)*M+1):(group.idx1*M), ]
m0.diff.obs <- margProb0.obs[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb0.obs[((group.idx1-1)*M+1):(group.idx1*M), ]
m1.diff.obs <- margProb1.obs[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb1.obs[((group.idx1-1)*M+1):(group.idx1*M), ]
  
# boostrap results
m_1.diff.boot <- list()
m0.diff.boot <- list()
m1.diff.boot <- list()
for(i in l) {
  piAlpha.boot <- PlotPiAlpha(prob.final.boot[[i]]$pi.a, G=G, M=M, figure=FALSE)
  piF.boot <- PlotPiF(p_1 = prob.final.boot[[i]]$pi.f_1, p0 = prob.final.boot[[i]]$pi.f0, p1 = prob.final.boot[[i]]$pi.f1, G = G, M=M, figure=FALSE)
  piF0_1.boot <- piF.boot[ ,c(1:M)]
  piF1_1.boot <- piF.boot[ ,c((3*M+1):(4*M))]
  piF00.boot <- piF.boot[ ,c((M+1):(2*M))]
  piF10.boot <- piF.boot[ ,c((4*M+1):(5*M))]
  piF01.boot <- piF.boot[ ,c((2*M+1):(3*M))]
  piF11.boot <- piF.boot[ ,c((5*M+1):(6*M))]
  margProb_1.boot <- (piAlpha.boot*piF1_1.boot + (1-piAlpha.boot)*piF0_1.boot)
  margProb0.boot <- (piAlpha.boot*piF10.boot + (1-piAlpha.boot)*piF00.boot) 
  margProb1.boot <- (piAlpha.boot*piF11.boot + (1-piAlpha.boot)*piF01.boot)
  m_1.diff.boot[[i]] <- margProb_1.boot[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb_1.boot[((group.idx1-1)*M+1):(group.idx1*M), ]
  m0.diff.boot[[i]]<- margProb0.boot[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb0.boot[((group.idx1-1)*M+1):(group.idx1*M), ]
  m1.diff.boot[[i]] <- margProb1.boot[((group.idx2-1)*M+1):(group.idx2*M), ] - margProb1.boot[((group.idx1-1)*M+1):(group.idx1*M), ]
}  
  
# function to plot legend
alpha.f <- 0.01
pdf(file = paste(output_folder, "legend_functional_diff.pdf"))
lgd = Legend(col_fun = colorRamp2(c(0, alpha.f), c("yellow","whitesmoke")), title = "p-value", at = c(0, alpha.f))
draw(lgd)
dev.off()

s <- c(1,2,3,5,6,4,7,8,9) #reorder of module

# reproduce Figure 4 and Figure 5
marginalFanalysis <- function(mat.obs, mat.boot, M, alpha = 0.01, state = "0", s = s){
  mat.boot.sd <- apply(matrix(unlist(mat.boot), nrow = M*M), 1 , sd)
  p.diff <- round(2 * pmin(pnorm(matrix(mat.obs, byrow = TRUE, nrow = 1)/mat.boot.sd), pnorm(matrix(mat.obs, byrow = TRUE, nrow = 1)/mat.boot.sd, lower.tail = FALSE)), 5)
  p.diff <- matrix(p.diff, nrow = M)
  signif.diff <- p.diff < alpha
  
  pdf(file = paste(output_folder, "signif_marginal_diff_", state, "_G", group.idx2, " & G", group.idx1, "_", alpha,".pdf", sep=""), width = 8.5 , height = 8.5)
  mycol <-  colorRamp2( c(1 - alpha, 1), c("whitesmoke","yellow"))
  Heatmap9by9.signif(mat.obs[s,s], col = mycol, signif = signif.diff[s,s], pmat = p.diff[s,s])
  dev.off()
  circlegraph(mat.obs[s,s],  p.diff[s,s], filename = paste(output_folder, "Marginal Functional Diff state", state, ".pdf", sep = ""), alpha  = alpha, factor = -2, color.limit = 0.1)
}
marginalFanalysis(mat.obs = m_1.diff.obs, mat.boot = m_1.diff.boot, M = 9, alpha = 0.01, state = "-1", s = s)
marginalFanalysis(mat.obs = m0.diff.obs, mat.boot = m0.diff.boot, M = 9, alpha = 0.01, state = "0", s = s)
marginalFanalysis(mat.obs = m1.diff.obs, mat.boot = m1.diff.boot, M = 9, alpha = 0.01, state = "1", s = s)
  
