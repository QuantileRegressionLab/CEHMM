# Clear the workspace
rm(list=ls())

# Define library path
# lib = getwd()
# repos = "http://cran.uk.r-project.org"
# .libPaths(c(.libPaths(), lib))

# Load necessary packages
library(foreach)
library(parallel)
library(doParallel)
library(ald)
library(quantreg)
library(expectreg)
library(mclust)
library(skewt)
library(cluster)
library(markovchain)
library(rqPen)
library(lqmm)
library(copula)
library(stats4)
library(sn)

### Let's choose the model #####
model_reg <- "expectile"
# model_reg <- "quantile"
if(model_reg == "expectile"){
  source("MainFunctions_cexp.R")
} else {
  source("MainFunctions_cquant.R")
}


#### Define simulations true values ######
n <- c(500,1000)
# n <- 1000
dist <- c("norm", "t", "st")

# Set the value of 'setting' (1 or 2)
setting <- 2  # Change this value to switch between settings

# Set the values of variables based on the value of 'setting'
if (setting == 1) {
  k <- 2
  nregs <- 2
  multi_dim <- 4
} else if (setting == 2) {
  k <- 3
  nregs <- 3
  multi_dim <- 4
} else {
  stop("Invalid value for 'setting'. Use 1 or 2.")
}

for (zz in 1:length(dist)) {
  print(paste("dist=",dist[zz]))
  
  for (ii in 1:length(n)) {
    
    
    R = 10 # restart
    MM = 250 #Montecarlo simulations
    # MM = 100 #Montecarlo simulations
    tau = list(list(.1,.1),list(.5,.5),list(.9,.9))
    tauvec = c(.1,.5,.9)
    # tauvec = c(.5)
    n.tau = length(tauvec)
    aRI <- matrix(NA,MM,n.tau)
    model <- list()
    multivariate <- list()
    set.seed(080693)
    
    
    ## Distribution parameters and copula choice #####
    wcop = "norm" #  norm or t
    # wcop = "t"
    if(dist[zz] == "norm" || dist[zz] == "t"){
      df = replicate(k, 5, simplify = F)
      gamma_skewt = rep(0, multi_dim)
    } else {
      df = replicate(k, 5, simplify = F)
      gamma_skewt = rep(2, multi_dim)
    }
    
    #### Initial data #####
    set.seed(080693)
    if(setting == 1){
      # 2X, 4Y, 2K
      beta = list(matrix(data = c(-2, 3, 1, 2, 1, -2, 4, 3, 3, 1, -2, -1), byrow =T, 1+nregs, multi_dim, dimnames = list(c("beta0", "beta1", "beta2"), c("d=1", "d=2", "d=3", "d=4"))),
                  matrix(data = c(3, -2, -3, -1, -2, 1, -3, 4, 1, 4, 1, 1), byrow = T, 1+nregs, multi_dim, dimnames = list(c("beta0", "beta1", "beta2"), c("d=1", "d=2", "d=3", "d=4"))))
    } else if(setting == 2){
      # 3X, 4Y, 3K
      beta = replicate(k, matrix(data = NA, byrow =T, 1+nregs, multi_dim, 
                                 dimnames = list(c("beta0", "beta1", "beta2", "beta3"), c("d=1", "d=2", "d=3", "d=4"))), simplify = F)
      beta[[1]][1,] <- round(runif(n = multi_dim, min = -3, max = -1),1)  
      beta[[2]][1,] <- round(runif(n = multi_dim, min = -1, max = 2),1)
      beta[[3]][1,] <- round(runif(n = multi_dim, min = 3, max = 5),1)
      beta[[1]][-1,] <- runif(n = (nregs)*multi_dim, min = -2, max = 2)
      beta[[2]][-1,] <- runif(n = (nregs)*multi_dim, min = -2, max = 2)
      beta[[3]][-1,] <- runif(n = (nregs)*multi_dim, min = -2, max = 2)
    }
    
    P <- A_diag <- Sigma <- list()
    set.seed(2)
    P[[1]] = matrix(c(runif(min = -1, max = 1, n = multi_dim^2)), multi_dim, multi_dim)
    P[[2]] = matrix(c(runif(min = -1, max = 1, n = multi_dim^2)), multi_dim, multi_dim)
    if(k == 3) P[[3]] = matrix(c(runif(min = -1, max = 1, n = multi_dim^2)), multi_dim, multi_dim)
    
    A_diag[[1]] = runif(min = 0, max = 1, n = multi_dim) 
    A_diag[[2]] = runif(min = 0, max = 1, n = multi_dim)
    if(k == 3) A_diag[[3]] = runif(min = 0, max = 1, n = multi_dim)
    Sigma[[1]] = P[[1]]*t(P[[1]]) + diag(A_diag[[1]])
    Sigma[[2]] = P[[2]]*t(P[[2]]) + diag(A_diag[[2]])
    if(k==3) Sigma[[3]] = P[[3]]*t(P[[3]]) + diag(A_diag[[3]])
    Sigma = lapply(Sigma, make.positive.definite)
    Sigma = lapply(Sigma, cov2cor)
    Sigma = lapply(Sigma, make.positive.definite)
    
    if(k == 2){
      delta = c(1,0)
      gamma = matrix(c(0,1,1,0),k,k,byrow=TRUE)
    } else if(k == 3){
      delta = c(1,0,0)
      gamma = matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),k,k,byrow=TRUE) #initial matrix
    }
    
    ### Initialization of final lists #####
    crit_states <- matrix(NA, n[ii], 3)
    beta.est <- replicate(k, array(NA,dim = c(1+nregs, multi_dim,MM)), simplify = F)
    aux_fin <- replicate(k, array(NA, dim = c(n[ii], multi_dim, MM)), simplify = F)
    sigma_fin <- replicate(k, matrix(NA, multi_dim, MM), simplify = F)
    diff_ <- tau_exp <- c()
    cpar_fin <- replicate(k, matrix(NA, MM, multi_dim*(multi_dim-1)/2), simplify = F)
    delta_fin <- matrix(NA,MM,k)
    gamma_fin <- array(NA, dim = c(k,k,MM))
    loglik_fin <- iters_fin <-  c()
    post_fin  <- array(NA, dim = c(k,n[ii],MM))
    delta_l <- gamma_l <- loglik_l <- post_l <- aux_l <- iters_l <- sigma_l <- crit <- betas.est <- taus_est <- list()
    if(setting == 1){
      #2x, 4Y
      Est.beta <- Est.beta.median <-  Bias <- Bias.median <- Stderr.beta <- replicate(k, array(NA, dim = c(1+nregs, multi_dim, n.tau),
                                                                                               dimnames = list(c("beta0", "beta1", "beta2"), c("d=1", "d=2", "d=3", "d=4"),
                                                                                                               as.character(tauvec))), simplify = F)
    } else if(setting == 2){
      #3X, 4Y
      Est.beta <- Est.beta.median <-  Bias <- Bias.median <- Stderr.beta <- replicate(k, array(NA, dim = c(1+nregs, multi_dim, n.tau),
                                                                                               dimnames = list(c("beta0", "beta1", "beta2", "beta3"), c("d=1", "d=2", "d=3", "d=4"),
                                                                                                               as.character(tauvec))), simplify = F)
    }
    cpar_est <- cpar_std <- replicate(k, matrix(NA, n.tau, multi_dim*(multi_dim-1)/2, 
                                                dimnames = list(as.character(tauvec), NULL)), simplify = F)
    
    ###Select sojourn distribution: #####
    #Geometric
    if(k == 2){
      M = 1e2; d.prob = c(0.1, 0.1); d = list(dgeom(0:M,d.prob[1]), dgeom(0:M,d.prob[2]))
      gamma.hmm = matrix(c(1-d.prob[1],d.prob[1],d.prob[2],1-d.prob[2]),2,2,byrow=TRUE)
    } else if(k == 3){
      M = 1e2; d.prob = c(0.05, 0.05, 0.05); d = list(dgeom(0:M,d.prob[1]), dgeom(0:M,d.prob[2]), dgeom(0:M,d.prob[3]))
      gamma.hmm = matrix(c(1-2*d.prob[1],d.prob[1],d.prob[1],d.prob[2],1-2*d.prob[2], d.prob[2], d.prob[3], d.prob[3], 1-2*d.prob[3]),3,3,byrow=TRUE)
    }
    
    start_time <- Sys.time()
    
    #### Loop start #####
    for (t in 1:n.tau) {
      print(paste("tau=",tauvec[t]))
      
      for (i in 1:MM) {
        print(paste("Simulation",i))
        set.seed(i)
        tau_d = tauvec[t]
        #select simulation error distribution: y = X%*%Beta + e
        if(dist[zz] == "norm"){
          multivariate[[i]] = expreg.hsmm.multi(ns = n[ii], m = k, dd= multi_dim, delta = delta, gamma = gamma,
                                                beta = beta, Sigma = Sigma, tau = tau[[t]], d=d) 
        } else {
          multivariate[[i]] = expreg.hsmm.multi.skewt(ns = n[ii], m = k, dd= multi_dim, delta = delta,
                                                      gamma = gamma, beta = beta, Sigma = Sigma, tau = tau[[t]], 
                                                      df=df,gamma_skewt = gamma_skewt,d=d)}
        
        
        
        # Y_t observations:
        Y = multivariate[[i]]$series
        x = multivariate[[i]]$regressor
        Y1 <- matrix(NA, n[ii], multi_dim)
        for (md in 1:multi_dim) {
          if(model_reg == "expectile"){Y1[,md] = Y[,md] - expectile(c(multivariate[[i]]$error[,md]), tau_d)
          } else if(model_reg == "quantile"){Y1[,md] = Y[,md] - quantile(c(multivariate[[i]]$error[,md]), tau_d)}
        }
        # pairs(Y1, col = multivariate[[i]]$state)
        
        #### EM algorithm #######
        
        ##### server #####
        ncores = detectCores()
        ncores
        # cl = makeCluster(12)
        cl = makeCluster(ncores-1)
        registerDoParallel(cl)
        # invisible(clusterEvalQ(cl = cl, .libPaths(getwd())))
        invisible(clusterEvalQ(cl = cl, .libPaths(.libPaths())))
        if(model_reg == "expectile"){
          invisible(clusterEvalQ(cl = cl, source("MainFunctions_cexp.R")))
        } else {
          invisible(clusterEvalQ(cl = cl, source("MainFunctions_cquant.R")))
        }
        parallel::clusterSetRNGStream(cl, 8693)
        
        tmp = foreach(r = 1:R, .packages = c("Rcpp", "MASS", "mvtnorm", "markovchain", "rqPen", "copula", "stats4", "ald", "quantreg")) %dopar% {
          # for (r in 1:R) {
          set.seed(r)
          # list of initial model parameters for the EM
          delta.s = runif(k)
          delta.s = delta.s/sum(delta.s)
          sigma.s = list()
          lm.init = list()
          cpar.init = replicate(k,  matrix(NA, nrow = 1,ncol = multi_dim*(multi_dim-1)/2), simplify = F)
          df_init = list()
          beta.init = replicate(k, matrix(data = NA, 1+nregs, multi_dim), simplify = F)
          sigma.s = replicate(k, c(), simplify = F)
          
          
          ## Clusterization #####
          kmeans.re <- list()
          kmeans.re$clustering = sample(k, n[ii], replace = T)
          hmm_init = markovchainFit(kmeans.re$clustering)
          init = rep(0, k)
          init[kmeans.re$clustering[1]] = 1
          gamma.s = hmm_init$estimate@transitionMatrix
          # 
          
          for(j in 1:k){
            for (dd in 1:multi_dim) {
              if(model_reg == "expectile"){
                lm.init[[j]] <- lm(Y1[kmeans.re$clustering==j,dd] ~ x[kmeans.re$clustering==j,-1])
                sigma.s[[j]][dd] <- sd(lm.init[[j]]$residuals)
                #beta initials for each cluster
                beta.init[[j]][,dd] <- lm.init[[j]]$coefficients
              } else if(model_reg == "quantile"){
                lm.init[[j]] <- quantreg::rq(Y1[kmeans.re$clustering==j,dd] ~ x[kmeans.re$clustering==j,-1], tau=.5)
                sigma.s[[j]][dd] <- mean(check(lm.init[[j]]$residuals, tau = .5)) #residuals std dev
                #beta initials for each cluster
                beta.init[[j]][,dd] <- lm.init[[j]]$coefficients
              }
            } ##dd lopp
            cpar.init[[j]] <- cor(Y1[kmeans.re$clustering==j,])[upper.tri(cor(Y1[kmeans.re$clustering==j,]))]
            df_init[[j]] <- 5
            # cpar.init[[j]] <- cbind(cpar.init[[j]])
          }
          #}

    tryCatch(em.hmm.cqereg(y = Y1, X_ = x, m = k, dd = multi_dim, delta=init, gamma=gamma.s,
                                 beta=beta.init,
                                 cpar = cpar.init, df_cop = df_init,
                                 sigma=sigma.s, tau = tau_d, which_cop = wcop,
                                 tol=10^-4, maxiter=200, trace=F),
                   error=function(e) {
                     NA
                   })
          
        } #r loop
        
        stopCluster(cl)  # Stop the cluster
        
        tmp.llk = matrix(-Inf, R, 1)
        for(r in 1:R) {
          if(is.na(tmp[[r]][1]) == F) {
            tmp.llk[r] = as.numeric(tmp[[r]]$loglik)
          } 
        }
        
        if(max(tmp.llk) == -Inf) {next
        } else {model[[i]] = tmp[[which.max(tmp.llk)]]}
        
        
        aRI[i,t] <- mclust::adjustedRandIndex(x = multivariate[[i]]$state,y = apply(model[[i]]$post,2,which.max))
        
        
        ##### state inversion sorting ######
        model[[i]]$betas <- model[[i]]$betas[order(sapply(model[[i]]$betas, `[[`, i=1))]
        
        for (j in 1:k) {
          beta.est[[j]][,,i] <- model[[i]]$betas[[j]]
          aux_fin[[j]][,,i] <- model[[i]]$aux[[j]]
          sigma_fin[[j]][,i] <- model[[i]]$sigma[[j]]
          cpar_fin[[j]][i,] <- model[[i]]$cpar[[j]]
        }
        
        delta_fin[i,] <- model[[i]]$delta
        gamma_fin[,,i] <- model[[i]]$gamma
        loglik_fin[i] <- model[[i]]$loglik
        iters_fin[i] <- model[[i]]$iteration
        post_fin[,,i] <- model[[i]]$post
        diff_[i] <- model[[i]]$dif
        crit_states[i,] <- model[[i]]$crit
        
        
      } #i loop
      
      
      # Remove any possible NA element from the list
      model <- model[!sapply(model, function(x) all(is.na(x)))]
      
      delta_l[[t]] <- delta_fin
      gamma_l[[t]] <- gamma_fin
      loglik_l[[t]] <- loglik_fin
      iters_l[[t]] <- iters_fin
      sigma_l[[t]] <- sigma_fin
      post_l[[t]] <- post_fin
      aux_l[[t]] <- aux_fin
      crit[[t]] <- colMeans(crit_states) 
      betas.est[[t]] <- list()
      taus_est[[t]] <- tau_exp
      #Estimated Coefficients, Bias and and Standard Errors
      for (j in 1:k) {
        betas.est[[t]][[j]] <- beta.est[[j]]
        # cpar_est[[j]][t,] <- colMeans(cpar_fin[[j]])
        cpar_est[[j]][t,] <- apply(cpar_fin[[j]], 2, mean, na.rm = T)
        cpar_std[[j]][t,] <- apply(cpar_fin[[j]], 2, fBasics::stdev, na.rm = T) 
        Est.beta[[j]][,,t] <- apply( beta.est[[j]], 1:2 , mean, na.rm = T)
        Est.beta.median[[j]][,,t] <- apply(beta.est[[j]], 1:2 , median, na.rm = T)
        # Bias
        Bias[[j]][,,t] <- (Est.beta[[j]][,,t] - beta[[j]])
        Bias.median[[j]][,,t] <- (Est.beta.median[[j]][,,t] - beta[[j]])
        Stderr.beta[[j]][,,t] <-  apply(beta.est[[j]], 1:2, fBasics::stdev, na.rm = T)
      }
      
      
    } # end loop in t
    end_time <- Sys.time()
    results <- list(n=n[ii], R = R, montecarlo = MM, beta = beta, Sigma.init = Sigma, aRI = aRI, betas.est = betas.est,
                    Est.beta = Est.beta, Est.beta.median = Est.beta.median, Stderr.beta = Stderr.beta, Bias = Bias, Bias.median = Bias.median,
                    delta = delta_l, gamma = gamma_l, loglik = loglik_l, iterations = iters_l,
                    sigma = sigma_l, post = post_l, aux = aux_l, Cop_par = cpar_est, Cop_std = cpar_std, 
                    sys_time = end_time - start_time, taus_est = taus_est)
    
    if(model_reg == "expectile"){
      save(results, file = paste0("c",wcop,"exp",n[ii],".",dist[zz],MM,"MC_4Y.RData"))
    } else if(model_reg == "quantile"){
      save(results, file = paste0("c",wcop,"quant",n[ii],".",dist[zz],MM,"MC_4Y.RData"))
    }
    
    
  } #ii loop
} ## zz distribution loop
