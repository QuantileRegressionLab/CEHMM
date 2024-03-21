rm(list=ls())
lib = getwd()
repos = "http://cran.uk.r-project.org"
.libPaths(c(.libPaths(), lib))


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

source("MainFunctions.R")

##### server #####
ncores = detectCores()
ncores
cl = makeCluster(ncores-1)
registerDoParallel(cl)
invisible(clusterEvalQ(cl = cl, .libPaths(getwd())))
invisible(clusterEvalQ(cl = cl, source("MainFunctions.R")))

#### Define simulations true values ######
model <- "expectile"
# model <- "quantile"
n <- c(500,1000)
dist <- c("norm", "t", "st")
for (zz in 1:3) {
  print(paste("dist=",dist[zz]))

for (ii in 1:length(n)) {

k = 2 #number of states
R = 20 # restart
nregs = 1
multi_dim = 2
MM = 100 #Montecarlo simulations
tau = list(list(.1,.1),list(.5,.5),list(.9,.9))
tauvec = c(.1,.5,.9)
n.tau = length(tauvec)
aRI <- matrix(NA,MM,n.tau)
model <- list()
multivariate <- list()
set.seed(080693)


## Distribution parameters #####
wcop = "norm" #  norm or t
# wcop = "t"
if(dist[zz] == "norm" || dist[zz] == "t"){
  df = list(5,5)
  gamma_skewt = list(0,0)
} else {
  df = list(5,5)
  # gamma_skewt = list(-2,2)
  gamma_skewt = list(2,2)
}

#### Initial data #####
beta = list(matrix(data = c(-2,1,3,-2), 1+nregs, multi_dim, dimnames = list(c("beta0", "beta1"), c("d=1", "d=2"))), 
            matrix(data = c(3,-2,-2,1), 1+ nregs, multi_dim, dimnames = list(c("beta0", "beta1"), c("d=1", "d=2"))))

Sigma = list(matrix(data = c(1,.2,.2,1), multi_dim, multi_dim), matrix(data = c(1,.7,.7,1), multi_dim, multi_dim))
delta = c(1,0)
gamma = matrix(c(0,1,1,0),k,k,byrow=TRUE)

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
Est.beta <- Est.beta.median <-  Bias <- Bias.median <- Stderr.beta <- replicate(k, array(NA, dim = c(1+nregs, multi_dim, n.tau), 
                                                                                         dimnames = list(c("beta0", "beta1"), c("d=1", "d=2"), 
                                                                                                         as.character(tauvec))), simplify = F)
cpar_est <- cpar_std <- replicate(k, matrix(NA, n.tau, multi_dim*(multi_dim-1)/2, 
                                            dimnames = list(as.character(tauvec), NULL)), simplify = F)

###Select sojourn distribution: #####
#Geometric
#M = 1e2; d.prob = c(0.1, 0.2); d = list(dgeom(0:M,d.prob[1]), dgeom(0:M,d.prob[2]))
M = 1e2; d.prob = c(0.1, 0.1); d = list(dgeom(0:M,d.prob[1]), dgeom(0:M,d.prob[2]))
gamma.hmm = matrix(c(1-d.prob[1],d.prob[1],d.prob[2],1-d.prob[2]),2,2,byrow=TRUE)

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
      if(model == "expectile"){Y1[,md] = Y[,md] - expectile(c(multivariate[[i]]$error[,md]), tau_d)
      } else if(model == "quantile"){Y1[,md] = Y[,md] - quantile(c(multivariate[[i]]$error[,md]), tau_d)}
    }
    
    
    #### EM algorithm #######
    
    
    tmp = foreach(r = 1:R, .packages = c("Rcpp", "MASS", "mvtnorm", "markovchain", "rqPen", "copula", "stats4", "ald", "quantreg")) %dopar% {
      #for (r in 1:R) {
      set.seed(r)
      # list of initial model parameters for the EM
      delta.s = runif(k)
      delta.s = delta.s/sum(delta.s)
      sigma.s = list()
      lm.init = list()
      cpar.init = replicate(k,  matrix(NA, nrow = 1,ncol = multi_dim*(multi_dim-1)/2), simplify = F)
      df_init = list()
      beta.init = replicate(k, matrix(data = NA, multi_dim, 1+nregs), simplify = F)
      sigma.s = replicate(k, matrix(data = NA, multi_dim, multi_dim), simplify = F)
      
      
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
          beta_clust <- x[kmeans.re$clustering==j,-1]
          if(model == "expectile"){
            lm.init[[j]] <- QRLM(Y = Y1[kmeans.re$clustering==j,dd], X = cbind(1, beta_clust), weights = rep(1, sum(kmeans.re$clustering==j)), 
                                 maxit = 1e3, err = 1e-06, tau = 0.5)
            sigma.s[[j]][dd] <- sd(lm.init[[j]]$residuals)
            #beta initials for each cluster
            beta.init[[j]][,dd] <- lm.init[[j]]$coefficients
          } else if(model == "quantile"){
            lm.init[[j]] <- quantreg::rq(Y1[kmeans.re$clustering==j,dd] ~ beta_clust, tau=.5)
            sigma.s[[j]][dd] <- mean(check(lm.init[[j]]$residuals, tau = .5)) #residuals std dev
            #beta initials for each cluster
            beta.init[[j]][,dd] <- lm.init[[j]]$coefficients
          }
        } ##dd lopp
          cpar.init[[j]] <- cor(Y1[kmeans.re$clustering==j,])[upper.tri(cor(Y1[kmeans.re$clustering==j,]))]
          df_init[[j]] <- 5
          cpar.init[[j]] <- cbind(cpar.init[[j]])
      }
      #}
      
      if(model == "expectile"){
        tryCatch(em.hmm.expreg(y = Y1, X_ = x, m = k, dd = multi_dim, delta=init, gamma=gamma.s,
                               beta=beta.init,
                               cpar = cpar.init, df_cop = df_init,
                               sigma=sigma.s, tau = tau_d, which_cop = wcop,
                               tol=10^-4, maxiter=200, trace=F),
                 error=function(e) {
                   NA
                 })
      } else if(model == "quantile"){
        tryCatch(em.hmm.quantreg(y = Y1, X_ = x, m = k, dd = multi_dim, delta=init, gamma=gamma.s,
                                 beta=beta.init,
                                 cpar = cpar.init,
                                 sigma=sigma.s, tau = tau_d, df_cop = df_init, which_cop = wcop,
                                 tol=10^-5, maxiter=200, trace=F),
                 error=function(e) {
                   NA
                 })
      }
    


    } #r loop
    
    tmp.llk = matrix(-Inf, R, 1)
    for(r in 1:R) {
      if(is.na(tmp[[r]][1]) == F) {
        tmp.llk[r] = as.numeric(tmp[[r]]$loglik)
      }
    }

    model[[i]] = tmp[[which.max(tmp.llk)]]
    
    
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
    cpar_est[[j]][t,] <- colMeans(cpar_fin[[j]])
    cpar_std[[j]][t,] <- apply(cpar_fin[[j]], 2, fBasics::stdev) 
    Est.beta[[j]][,,t] <- apply( beta.est[[j]], 1:2 , mean )
    Est.beta.median[[j]][,,t] <- apply(beta.est[[j]], 1:2 , median )
    # Bias
    Bias[[j]][,,t] <- (Est.beta[[j]][,,t] - beta[[j]])
    Bias.median[[j]][,,t] <- (Est.beta.median[[j]][,,t] - beta[[j]])
    Stderr.beta[[j]][,,t] <-  apply(beta.est[[j]], 1:2, fBasics::stdev)
  }
  
  
} # end loop in t
end_time <- Sys.time()
results <- list(n=n[ii], R = R, montecarlo = MM, beta = beta, Sigma.init = Sigma, aRI = aRI, betas.est = betas.est,
                Est.beta = Est.beta, Est.beta.median = Est.beta.median, Stderr.beta = Stderr.beta, Bias = Bias, Bias.median = Bias.median,
                delta = delta_l, gamma = gamma_l, loglik = loglik_l, iterations = iters_l,
                sigma = sigma_l, post = post_l, aux = aux_l, Cop_par = cpar_est, Cop_std = cpar_std, 
                sys_time = end_time - start_time, nas.llk = nas.llk_s, infs.llk = infs.llk_s, taus_est = taus_est)

if(model == "expectile"){
  save(results, file = paste0("c",wcop,"exp",n[ii],".",dist[zz],"1Y.RData"))
} else if(model == "quantile"){
  save(results, file = paste0("c",wcop,"quant",n[ii],".",dist[zz],"1Y.RData"))
}


} #ii loop
} ## zz distribution loop

