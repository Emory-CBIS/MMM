EM_groupmod_log_Firth<-function(I_G, II, B, D, Rho, Gamma, X2, Xi2, Mu, Sigma2, NodeMod, N_st, Beta.a, Beta.f0, Beta.f1, X, tol=1e-6, logtol=1e-5, etol1 = 1e-3, etol2 = 1e-3, etol3 = 1e-3,  maxiter=100, stop, 
                                A=NULL, F=NULL, A.init = NULL, F.init = NULL, SIM=FALSE, posterior = FALSE){

  # EM algorithm for MMM
  #
  # Args:
  #   I_G:  1*G vector, (I_1, I_2, ... , I_G), I_g is the number of samples in each group, g in 1 ... G
  #   II: cusum of I_G, (0, I_1, I_1+I_2, ..., I_1+...+I_G = I)
  #   B:  Functional measure I*N
  #   D:  Anatomical measure I*N
  #   Rho:  zero proportion for each anatomical state, 1*2 vector
  #   Gamma: 2*3 matrix, proportion of mixture gaussians for each anatomial state
  #   X2: 2*3 matrix , mean of mixture gaussians for each anatomical state
  #   Xi2: 2*3 matrix, sd of mixture gaussians for each anatomial state
  #   Mu: 2*3 matrix, mean of gaussion for each anatomical state and functional state
  #   Sigma2: 2*3 matrix, sd of gaussian for each anatomial state and functional state
  #   NodeMod: save the information of node and module which help calculate N_st
  #   N_st:  N*5, 4th and 5th column save the module the connection information
  #   Beta.a:  (q) * M * M ; include intercept 
  #   Beta.f0: (2q) * M * M: first (q) for -1 and second (q) for 1  when anatomical = 0
  #   Beta.f1: (2q) * M * M: first (q) for -1 and second (q) for 1  when anatomial = 1
  #   X: group design matrix G*q including intercept
  #   stop: stopping criterion for the algorithm: 1: NUM (the number of paramteres not satisfied the convergence criterion); 2: diff (mean square difference in all parameters); 3: diff_log (difference in logliklihood);
  #   etol1: tolerence for the Beta parameter in stop 1
  #   etol2: tolerence for the DTI model stop 1
  #   etol3: tolerance for the FMRI mdoel stop 1
  #   tol: tolerance for diff
  #   diff_log: tolerance for logliklihood
  #   maxiter: maximum iteration number
  #   A: true latent structural states (only use when SIM = TRUE)
  #   F: true latent functional states (only use when SIM = TRUE)
  #   A.init: initial estimation of latent structural states (only use when SIM = TRUE)
  #   F.init: iniial estimation of latent functional states (only use when SIM = TRUE)
  #   posterior: indicator of outputting the posterior probability and original data or not 
  
  # Returns:
  #   Beta.a, Beta.f0, Beta.f1: parameters in the latent level
  #   Rho, Gamma, X2, Xi2, Mu, Sigma2: parameters in the observational level
  #   p.hat: G*N*L(=6), expectation of missing varaible
  #   Loglik: Loglikihood
  #   iter: number of iteraction
  #   time: time elapsed
  #   RES, RES0, RES1: intermediate matrix in latent model parameter fitting in M step  
  #   summary.prediction: check how accurate final estimation close to the true label
  #   NUM, diff_log, diff, relative.diff.log: stopping metrics
  #   Loglik.trace, Beta.a.trace, Beta.f0.trace, Beta.f1.trace: trace metrics
  #   N: Number of Subject
  #   G: Number of Group 
  #   M: Number of Module
  #   I_G: 1*G vector, (I_1, I_2, ... , I_G), I_g is the number of samples in each group, g in 1 ... G
  #   II: cusum of I_G, (0, I_1, I_1+I_2, ..., I_1+...+I_G = I)
  #   X: group design matrix G*q including intercept    
  
  ptm <- proc.time()
  
  G <- length(I_G)     #number of the groups 
  N <- dim(B)[2]       #number of combination of connection
  M <- dim(Beta.a)[3]  #number of Modules
  q <- dim(Beta.a)[1]  #number of covariates including intercept
  pi.a <- array(1, dim=c(G, M, M))       # A=1 
  pi.f1 <- array(0, dim=c(G, M, M, q))   # F==1
  pi.f_1 <- array(0, dim=c(G, M, M, q))  # F==-1
  pi.f0 <- array(0, dim=c(G, M, M, q))   # F==0
  

  III <- matrix(0,nrow = G, ncol = N) # III: G*N save the number of D==0 in each group for each connection
  for (g in 1:G) {
    III[g, ] <- apply(D[c((II[g]+1):II[g+1]), ], 2, function(x) sum(x==0))
  }
  
  IIII <- list()   
  for (g in 1:G) {
    IIII[[g]] <- apply(D[c((II[g]+1):II[g+1]), ], 2, function(x) which(x!=0) + II[g])
  }
  
  diff_log <- 1 # difference between two log likelihood 
  diff <- 1 # difference between parameters
  NUM <- 1 # number of parameters which are not satisfied with the convergence criteria
  
  L <- 6 # combination of anaotomical latent state and functional latent state:  l=1 j=0 k=-1; l=2 j=0 k=0; l=3 j =0 k =1; l=4 j=1 k=-1; l=5 j=1 k=0; l= 6 j=1 k=1; l = 3*j+k+2;  l = 1...6
  P.hat <- array(0, dim=c(G,N,L)) 
  p.hat <- array(0, dim=c(G,N,L))
  Loglik.trace <- NULL
  Beta.a.trace <- NULL
  Beta.f0.trace <- NULL
  Beta.f1.trace <- NULL
  Anatomical.trace <- NULL
  Functional.trace <- NULL
  
  ##################################
  #  Begin of iteration
  ##################################
  for (iter in 1:maxiter) {
    sum.na <- sum(is.na(Beta.a)) + sum(is.na(Beta.f0)) + sum(is.na(Beta.f1))
    if (sum.na >0) {
      cat("####################","\n")
      cat("beta estimation has NA","\n")
      cat("####################","\n")
    }
    
    Prob <- Calc_Prob_Log(M = M, G = G, Beta.a=Beta.a, Beta.f0=Beta.f0, Beta.f1=Beta.f1, X=X)
    pi.a <- Prob[[1]]
    pi.f1 <- Prob[[2]]
    pi.f0 <- Prob[[3]]
    pi.f_1 <- Prob[[4]]
    pi.a0 <- Prob[[5]]
    
    sum.na <- sum(is.na(pi.a)) + sum(is.na(pi.f1)) + sum(is.na(pi.f_1)) + sum(is.na(pi.f0))
    if (sum.na >0) {
      cat("####################","\n")
      cat("probability has NA","\n")
      cat("####################","\n")
    }
    
    Pi.f1 = list(); Pi.f0 = list(); Pi.f_1 = list();
    for(g in 1:G){
      Pi.f1[[g]] = pi.f1[g,,,];
      Pi.f0[[g]] = pi.f0[g,,,];
      Pi.f_1[[g]] = pi.f_1[g,,,];
    }
    
    ##############
    ### E step ###
    ##############
    
    PP <- Update_group_module_log(III=III, I_G=I_G, II=II, pi_a=pi.a, pi_a0 = pi.a0, pi_f1=Pi.f1, pi_f0=Pi.f0, pi_f_1=Pi.f_1, 
                                  Node_Mod=as.matrix(NodeMod), B=B, D=D, Rho=Rho, X2=X2, Xi2=Xi2, Mu=Mu, Sigma2=Sigma2, Gamma=Gamma)
    p.hat <- PP$p.hat
    P.hat <- PP$P.hat
    sum.na.p <- sum(is.na(p.hat))
    if (sum.na.p >0) {
      cat("####################","\n")
      cat("phat has NA","\n")
      cat("####################","\n")
    }
    if(iter==1){
      p0.hat <- p.hat # save the inital estimated probability
    }
    w <- PP$w
    sum.na.w <- sum(is.na(w))
    if (sum.na.w >0) {
      cat("####################","\n")
      cat("w has NA","\n")
      cat("####################","\n")
    }
    
    ##############
    ### M step ###
    ##############
    
    # Structral Measure parameters
    Rho.new <- PP$Rho.new;
    if (sum(is.na(Rho.new)) >0) {
      cat("####################","\n")
      cat("Rho has NA","\n")
      cat("####################","\n")
    }
    Mu.new <- matrix(0, nrow=2, ncol=3);
    Sigma2.new <- matrix(0, nrow=2, ncol=3);
    
    Sum.Gamma.new <- c(0, 0);
    Gamma.new <- matrix(0, nrow=2, ncol=3);
    X2.new <- matrix(0, nrow=2, ncol=3);
    Xi2.new <- matrix(0, nrow=2, ncol=3);
    
    for (j in 1:2) {
      for (l in 1:3) {
        o <- 3 * (j - 1) + l
        for (n in 1:N) {
          for (g in 1:G) {
            Gamma.new[j, l] <- Gamma.new[j, l] + sum(w[IIII[[g]][[n]], n, o])
            X2.new[j, l] <- X2.new[j, l] + sum(w[IIII[[g]][[n]], n, o] * D[IIII[[g]][[n]], n])
            Xi2.new[j, l] <- Xi2.new[j, l] + sum(w[IIII[[g]][[n]], n, o] * (D[IIII[[g]][[n]],n]-X2[j,l])^2) 
          }
        }
      }
    }
    Sum.Gamma.new <- rowSums(Gamma.new)
    for (j in 1:2) {
      for (l in 1:3) {
        X2.new[j, l] <- X2.new[j, l] / Gamma.new[j, l]
        Xi2.new[j, l] <- Xi2.new[j, l] / Gamma.new[j, l]
        Gamma.new[j, l] <- Gamma.new[j, l] / Sum.Gamma.new[j]
      }
    }
    
    # Functional Measure parameters
    for(j in 1:2){                                        
      for(k in 1:3){
        l <- 3*(j-1)+k;
        nume1 <- 0; nume2 <- 0
        domi <- 0
        for(n in 1:N){
          for(g in 1:G){
            nume1 <- nume1+sum(B[(II[g]+1) : II[g+1],n])*p.hat[g,n,l]
            nume2 <- nume2+sum((B[(II[g]+1) : II[g+1],n]-Mu[j,k])^2)*p.hat[g,n,l]
            domi <- domi+I_G[g]*p.hat[g,n,l]
          }
        }
        Mu.new[j,k] <- nume1/domi   
        Sigma2.new[j,k] <- nume2/domi
      }
    }
    
    # Latent structural model parameters
    Y_1 = array(0, dim = c(G, M, M))
    for (s in 1:M){
      for (t in 1:M){
        for (g in 1:G){
          Y_1[g, s, t] = sum(p.hat[g, NodeMod$Mod1==s&NodeMod$Mod2==t, 4:6])
        }
      }
    }
    
    RES <- list()
    zzz <- 1
    Beta.a.new <- array(0, dim = c(q,M,M))
    Response.anat <-  array(0, dim = c(G, 2))
    for (s in 1:M){
      for (t in s:M){
        for (g in 1:G){
          Response.anat[g,] <- c(N_st[s, t] - Y_1[g, s, t], Y_1[g, s, t])
        }
        RES[[zzz]] <- Response.anat
        brfit <- brmultilogit(X, Response.anat, bias = TRUE, ref = 1)
        Beta.a.new[, s, t] <- matrix(brfit$B)
        Beta.a.new[, t, s] <- Beta.a.new[, s, t]
        zzz = zzz + 1
      }
    }
    
    YD = array(0, dim = c(G, L, M, M))
    for(g in 1:G){
      for(s in 1:M){
        for(t in 1:M){
          for(l in 1:L){
            YD[g, l, s, t] <- sum(p.hat[g, NodeMod$Mod1==s&NodeMod$Mod2==t, l])
          }
        }
      }
    }
    
    # Latent Functional model parameters
    RES0 <- list()  # used for monitoring the number in each module  when A==0 
    RES1 <- list()  # used for monitoring the number in each modul  when A==1
    Beta.f0.new <- array(0, dim = c(2*q, M, M))
    Beta.f1.new <- array(0, dim = c(2*q, M, M))
    zzz <- 1
    Response.func <- array(0, dim = c(G, 3))
    for (s in 1:M){
      for (t in s:M){
        for (g in 1:G){
          Response.func[g, ] <- YD[g, c(2, 1, 3), s, t]
        }
        RES0[[zzz]] <- Response.func 
        brfit <- brmultilogit(X, Response.func, bias = TRUE, ref = 1)
        Beta.f0.new[, s, t] <- c(matrix(brfit$B, byrow = TRUE, nrow = q))
        Beta.f0.new[, t, s] <- Beta.f0.new[, s, t]
        
        for (g in 1:G){
          Response.func[g, ] <- YD[g, c(5, 4, 6), s, t]
        }
        RES1[[zzz]] <- Response.func
        brfit <- brmultilogit(X, Response.func, bias = TRUE, ref = 1)
        Beta.f1.new[, s, t] <- c(matrix(brfit$B, byrow = TRUE, nrow = q))
        Beta.f1.new[, t, s] <- Beta.f1.new[, s, t]
        zzz <- zzz + 1
      }
    }
    
    ##########################################
    ####### stop criterion check #############
    ##########################################
    
    # difference calculated by mean square error of each parameter 
    diff <- sqrt((sum((Beta.f1.new - Beta.f1)^2) + sum((Beta.f0.new - Beta.f0)^2) + sum((Beta.a.new - Beta.a)^2)  + sum((Rho.new - Rho)^2) + sum((Gamma.new - Gamma)^2) + sum((X2.new - X2)^2) + sum((Xi2.new - Xi2)^2)
                  + sum((Mu.new-Mu)^2)  + sum((Sigma2.new-Sigma2)^2))/(5 * q * M^2 + 2 + 6 * 5))
    
    # new stop criterion: consider different stopping criterion for different parameter
    NUM <- sum((Beta.f1.new - Beta.f1)^2/((Beta.f1)^2 + 1e-8) > etol1) + 
      sum((Beta.f0.new - Beta.f0)^2/((Beta.f0)^2 + 1e-8) > etol1) + 
      sum((Beta.a.new - Beta.a)^2/((Beta.a)^2 + 1e-8)> etol1) + 
      sum((Rho.new - Rho)^2/((Rho)^2 + 1e-8) > etol2) + 
      sum((Gamma.new - Gamma)^2/((Gamma)^2 + 1e-8)> etol2) +  
      sum((X2.new - X2)^2/((X2)^2 + 1e-8) > etol2) + 
      sum((Xi2.new - Xi2)^2/((Xi2)^2 + 1e-8) > etol2) + 
      sum((Mu.new - Mu)^2/((Mu)^2 + 1e-8) > etol3) + 
      sum((Sigma2.new - Sigma2)^2/((Sigma2)^2 + 1e-8) > etol3)
    
    if(iter>=2){
      Loglik_obs1<- Loglik_obs
    }
    if(iter==1){
      Loglik_obs1<- 0
    }
    
    
    Loglik_obs <- sum(log(apply(exp(P.hat), c(1, 2), sum)))
    print(Loglik_obs)
    
    if(iter==1){
      diff_log <- Loglik_obs
    }
    if(iter>=2){
      diff_log <- Loglik_obs - Loglik_obs1
    }
    relative.diff.log <- abs(diff_log)/abs(Loglik_obs)
    
    Loglik.trace <- c(Loglik.trace, Loglik_obs)
    
    Beta.a.trace <- c(Beta.a.trace, sum((Beta.a.new - Beta.a)^2))
    Beta.f0.trace <- c(Beta.f0.trace, sum((Beta.f0.new - Beta.f0)^2))
    Beta.f1.trace <- c(Beta.f1.trace, sum((Beta.f1.new - Beta.f1)^2))
    
    cat("Iter = ",iter,"\n")
    cat("Beta_a = ", abs(Beta.a.new[1, 1:2, 1:2] - Beta.a[1, 1:2, 1:2])/abs(Beta.a[1, 1:2, 1:2]),"\n")
    cat("diff Beta_a = ", sum((Beta.a.new - Beta.a)^2),"\n")
    cat("Beta_f0 = ", abs(Beta.f0.new[1, 1:2, 1:2] - Beta.f0[1, 1:2, 1:2])/abs(Beta.f0[1, 1:2, 1:2]),"\n")
    cat("diff Beta_f0", sum((Beta.f0.new - Beta.f0)^2),"\n")
    cat("Beta_f1 = ", abs(Beta.f1.new[1, 1:2, 1:2] - Beta.f1[1, 1:2, 1:2])/abs(Beta.f1[1, 1:2, 1:2]),"\n")
    cat("diff Beta_f1", sum((Beta.f1.new - Beta.f1)^2),"\n")
    cat("diff = ", diff,"\n","Loglik = ", Loglik_obs,"\n","Rho = ", abs(Rho.new-Rho)/abs(Rho), "\n", "Gamma = ", abs(Gamma.new-Gamma)/abs(Gamma),"\n","X2 = ", abs(X2.new-X2)/abs(X2),"\n","Xi2 = ", abs(Xi2.new-Xi2)/abs(X2),"\n","Mu = ", abs(Mu.new-Mu)/abs(Mu),"\n","Sigma2 = ", abs(Sigma2.new-Sigma2)/abs(Sigma2),"\n")
    cat("diff.log = ", diff_log,"\n")
    cat("relative.diff.log = ", abs(diff_log)/abs(Loglik_obs), "\n")
    cat("NUM", NUM, "\n")
    
    # Update parameter
    Beta.a <- Beta.a.new
    Beta.f0 <- Beta.f0.new
    Beta.f1 <- Beta.f1.new
    Rho <- Rho.new
    Gamma <- Gamma.new
    X2 <- X2.new
    Xi2 <- Xi2.new
    Mu <- Mu.new
    Sigma2 <- Sigma2.new 
    
    if (iter>=2) {
      if ((Loglik_obs1-Loglik_obs) > 0) {
        cat("#############################","\n")
        cat("No Increasing!!!!!!!!!!!!!!!!!!!!!!!!!","\n")
        cat("#############################","\n")
      }
      if((Loglik_obs1-Loglik_obs)/abs(Loglik_obs)>1e-4){
        cat("#############################","\n")
        cat("No Increasing!!!!!!!!!!!!!!!!!!!!!!!!!","\n")
        cat("#############################","\n")
        break
      }
    }
    #iter=iter+1;
    if (stop == 1) {
      if (NUM == 0 ){
        break
      }
    }
    
    if (stop == 2) {
      if (diff < tol){
        break
      }
    }
    
    if (stop == 3) {
      if (diff_log < logtol){
        break
      }
    }
  }
  time = proc.time() - ptm  # record time 
  # end of iteration
  ######################################################
  
  ######################################################
  # summarize prediction results
  mis.calssified.num.init <- 0
  mis.calssified.num.start <- 0
  mis.classified.num.final <- 0
  if (SIM == TRUE)  {
    mis.classified.num.init <- sum((A.init*3 + F.init + 2) != (A*3 + F + 2))
    mis.calssified.num.start <- sum(t(apply(p0.hat, c(1, 2), which.max)) != (A*3 + F + 2))
    mis.classified.num.final <- sum(t(apply(p.hat, c(1, 2), which.max)) != (A*3 + F + 2))
  }
  change.prediction <- sum(t(apply(p.hat, c(1, 2), which.max)) != t(apply(p0.hat , c(1, 2), which.max)))
  if (SIM == TRUE) {
    summary.prediction <- data.frame(mis.init = mis.classified.num.init, mis.start = mis.calssified.num.start, mis.final = mis.classified.num.final, change = change.prediction)
  } else {
    summary.prediction <- data.frame(change = change.prediction)
  }
  
  if (posterior == FALSE) {
    return(list(Beta.a = Beta.a, Beta.f0 = Beta.f0, Beta.f1 = Beta.f1, Rho = Rho, Gamma = Gamma, X2 = X2, Xi2 = Xi2, Mu = Mu, Sigma2 = Sigma2, Loglik = Loglik_obs,
                iter = iter, time = time, RES = RES, RES0 = RES0, RES1 = RES1, summary.prediction =  summary.prediction, NUM = NUM, diff_log = diff_log, diff=diff, relative.diff.log = relative.diff.log,
                Loglik.trace = Loglik.trace, Beta.a.trace = Beta.a.trace, Beta.f0.trace = Beta.f0.trace, Beta.f1.trace = Beta.f1.trace, 
                N = N, G = G, M = M, I_G = I_G, II = II, X = X))
  }
  if (posterior == TRUE) {
    return(list(Beta.a = Beta.a, Beta.f0 = Beta.f0, Beta.f1 = Beta.f1, Rho = Rho, Gamma = Gamma, X2 = X2, Xi2 = Xi2, Mu = Mu, Sigma2 = Sigma2, p.hat = p.hat, Loglik = Loglik_obs,
                iter = iter, time = time, RES = RES, RES0 = RES0, RES1 = RES1, summary.prediction =  summary.prediction, NUM = NUM, diff_log = diff_log, diff=diff, relative.diff.log = relative.diff.log,
                Loglik.trace = Loglik.trace, Beta.a.trace = Beta.a.trace, Beta.f0.trace = Beta.f0.trace, Beta.f1.trace = Beta.f1.trace,
                N = N, G = G, M = M, NodeMod = NodeMod, B = B, D = D, I_G = I_G, II = II, X = X))
  }
} # end of function
