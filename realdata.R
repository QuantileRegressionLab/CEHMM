rm(list = ls())
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

ci_fun <- function(centre, stder){
  lx = centre - 1.96*stder
  ux = centre + 1.96*stder
  vec <- c(lx,ux)
  return(vec)
}

#### Dataset Upload ####
load("data/returns.RData")
ret_cexp <- ret.df
source("MainFunctions.R")

##### Server #####
ncores = detectCores()
ncores
#cl = makeCluster(ncores, outfile="")  #to print
cl = makeCluster(ncores-1) #without print
registerDoParallel(cl)
invisible(clusterEvalQ(cl = cl, .libPaths(getwd())))
invisible(clusterEvalQ(cl = cl, source("MainFunctions.R")))

#Dependent variables: returns_cryptos
#Indepedent variables: returns_indexes

model = "expectile"
# model = "quantile"
dims = 5 #number of dependent variables
ys = as.matrix(ret_cexp[2:dim(ret_cexp)[1], 1:dims])
xs = as.matrix(ret_cexp[1:dim(ret_cexp)[1]-1,-c(1:dims)]) #t-1
xs_1 = cbind(1,xs)
nregs_1 = ncol(xs_1)
K = 2 #number of states
#K = 3
n = dim(ys)[1] #sample size
B = 1e3 #bootstrap sample
R = 50 # restart
tauvec = c(.01, .05, .5, .95, .99)
n.tau = length(tauvec)
multivariate <- list()
set.seed(080693)

wcop = "t"

### Variable Inizialitazion ####
tau_exp <- matrix(NA, nrow = n.tau, ncol = dims)
model <- replicate(n.tau, list(), simplify = F)
em_time <- list()
betas.est <- replicate(n.tau, array(NA, dim = c(nregs_1, dims, K),
                                    dimnames = list(c("Intercept", colnames(xs)), colnames(ys),
                                                    as.character(paste("K =", 1:K)))),
                       simplify = F)


Corr_mat <- list()
sigmas.est <- replicate(n.tau, matrix(NA, nrow = dims, ncol = K, dimnames = list(NULL,as.character(paste("K =", 1:K)))), simplify = F)
cpar.est <- replicate(n.tau, matrix(NA, K, dims*(dims-1)/2,
                                    dimnames = list(as.character(paste("K =", 1:K)), NULL)), simplify = F)
df.est  <- matrix(NA, n.tau, K, dimnames = list(as.character(tauvec), NULL))
delta.est <- matrix(NA, nrow = n.tau, ncol = K)
gamma.est <- array(NA, dim = c(K,K,n.tau))
model_boot <- list()
######## Restart procedure ########
tmp.llk = matrix(NA, R, n.tau)

time.1 <- Sys.time()

for (t in 1:n.tau) {
  tau_d = tauvec[t]
  print(paste("tau=", tau_d))
  tmp = foreach(r = 1:R, .packages = c("Rcpp", "MASS", "mvtnorm", "markovchain", "rqPen", "copula", "stats4", "ald", "quantreg")) %dopar% {
    set.seed(r)
    #### List of initial model parameters for the EM #####
    delta.s = runif(K)
    delta.s = delta.s/sum(delta.s) 
    sigma.s = list()
    lm.init = list()
    cpar.init = replicate(K,  matrix(NA, nrow = 1,ncol = dims*(dims-1)/2), simplify = F)
    df_init = list()
    beta.init = replicate(K, matrix(data = NA, 1+dim(xs)[2], dims), simplify = F)
    sigma.s = replicate(K, matrix(data = NA, dims, 1), simplify = F)
    
    #### Initial parameters with k-means cluster: ####
    kmeans.re <- list()
    kmeans.re$clustering = sample(K, n, replace = T)
    hmm_init = markovchainFit(kmeans.re$clustering)
    init = rep(0, K)
    init[kmeans.re$clustering[1]] = 1
    gamma.s = hmm_init$estimate@transitionMatrix
    
    
    for(j in 1:K){
      for (dd in 1:dims) {
        lm.init[[j]] <- lm(ys[kmeans.re$clustering==j,dd] ~ xs[kmeans.re$clustering==j,] )
        sigma.s[[j]][dd] <- sd(lm.init[[j]]$residuals)
        #beta initials for each cluster
        beta.init[[j]][,dd] <- lm.init[[j]]$coefficients
      }
      if(wcop=="norm"){
        cpar.init[[j]] <- P2p(cor(ys[kmeans.re$clustering==j,]))
        df_init[[j]] <- 1
        cpar.init[[j]] <- cbind(cpar.init[[j]])
      }else{
        cpar.init[[j]] <- P2p(cor(ys[kmeans.re$clustering==j,]))
        df_init[[j]] <- 5
        cpar.init[[j]] <- cbind(cpar.init[[j]])
      }
    }
    

    
    if(model == "expectile"){
      tryCatch(em.hmm.expreg(y = ys, X_ = xs_1, m = K, dd = dims, delta=init, gamma=gamma.s,
                             beta=beta.init,
                             cpar = cpar.init, df_cop = df_init, 
                             sigma=sigma.s, tau = tau_d, which_cop = wcop,
                             tol=10^-6, maxiter=400, trace=F),
               error=function(e) {
                 NA
               })
    } else if(model == "quantile"){
      tryCatch(em.hmm.quantreg(y = ys, X_ = xs_1, m = K, dd = dims, delta=init, gamma=gamma.s,
                               beta=beta.init,
                               cpar = cpar.init, df_cop = df_init,
                               sigma=sigma.s, tau = tau_d, which_cop = wcop,
                               tol=10^-6, maxiter=400, trace=F),
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
  
  model[[t]] = tmp[[which.max(tmp.llk)]]
  
  ##### Expectile Estimation ####
  # abc = t(sapply(1:n, function(j) (xs_1[j,])%*%model[[t]]$betas[[apply(model[[t]]$post,2,which.max)[j]]]))
  # tau_exp[t,] <- colMeans(pmax(-(ys-abc), 0))/colMeans(abs(ys - abc))
  # colMeans(ys < abc)
  
  for (j in 1:K) {
    
    betas.est[[t]][,,j] <- model[[t]]$betas[[j]]
    sigmas.est[[t]][,j] <- model[[t]]$sigma[[j]] #variance of each y
    cpar.est[[t]][j,] <- model[[t]]$cpar[[j]] #correlations
    df.est[t,] <- unlist(model[[t]]$df)
  }
  gamma.est[,,t] <- model[[t]]$gamma
  delta.est[t,] <- model[[t]]$delta
  em_time[[t]] <- model[[t]]$time
  Corr_mat[[t]] <- lapply(model[[t]]$cpar, p2P)
}  # t loop

if(model == "expectile"){
  save.image(file = "realdata_cexp_2K_t.RData")
} else if(model == "quantile"){
  save.image(file = "realdata_cquant_2K_t.RData")
}



####### CIs Bootstraping #########

for (t in 1:n.tau) {
  
  tau_d = tauvec[t]
  print(paste("tau=", tau_d))
  M = 1e2; d.prob = 1-diag(gamma.est[,,t]); d = list(dgeom(0:M,d.prob[1]), dgeom(0:M,d.prob[2]))
  gamma.hmm = matrix(c(1-d.prob[1],d.prob[1],d.prob[2],1-d.prob[2]),2,2,byrow=TRUE)
  
  model_boot[[t]] = foreach(j = 1:B, .packages = c("Rcpp", "MASS", "mvtnorm", "markovchain", "rqPen", "copula", "stats4")) %dopar% {
    #for (j in 1:B) {
    set.seed(j)
    
    # We generate from a normal or a t copula
    
    boot_dt <- expreg.hsmm.multi.real(reg = xs_1, ns = n, m = K, dd = dims, delta = delta.est[t,], gamma = matrix(c(0,1,1,0),2,2,byrow=TRUE),
                                      beta = model[[t]]$betas, df = model[[t]]$df, Sigma = Corr_mat[[t]], sigma = model[[t]]$sigma,
                                      tau = tau_d, d=d, wcop = wcop, cpar = model[[t]]$cpar)
    
    Y = boot_dt$series 

    tryCatch(em.hmm.expreg(y = Y, X_ = xs_1, m = K, dd = dims, delta=delta.est[t,], gamma=gamma.est[,,t],
                           beta = model[[t]]$betas,
                           cpar = model[[t]]$cpar, df_cop = model[[t]]$df,
                           sigma=model[[t]]$sigma, tau = tau_d, which_cop = wcop, 
                           tol=10^-4, maxiter=200, trace=F),
             error=function(e) {
               NA
             })
    
  } ## j loop
  
} #t loop

time.2 <- Sys.time()
end.time <- time.2 - time.1

if(model == "expectile"){
  save.image(file = "realdata_cexp_2K_t.RData")
} else if(model == "quantile"){
  save.image(file = "realdata_cquant_2K_t.RData")
}

### Save estimates

betas.b <- replicate(K, array(NA, dim = c(nregs_1, dims, B)), simplify = F)
betas.b.est <- betas.b.std <- replicate(n.tau, array(NA, dim = c(nregs_1, dims, K), dimnames = list(c("Intercept", colnames(xs)), colnames(ys), as.character(paste("K =", 1:K)))), simplify = F)
cpar.b <- replicate(K, matrix(NA, nrow = B, ncol = length(model_boot[[1]][[1]]$cpar[[1]])), simplify = F)
cpar.b.est <- cpar.b.std <- replicate(n.tau, array(NA, dim = c(1, length(model_boot[[1]][[1]]$cpar[[1]]), K)), simplify = F)
df.b <- replicate(K, matrix(NA, B, 1), simplify = F)
df.b.est <- df.b.std <- replicate(n.tau, matrix(NA, K, 1, dimnames = list(as.character(c(1:K)), NULL)), simplify = F)
sigmas.b <- replicate(K, array(NA, dim = c(dims, 1, B)), simplify = F)
sigmas.b.est <- sigmas.b.std <- replicate(n.tau, array(NA, dim = c(dims, 1, K)), simplify = F)
test_statistic <- replicate(n.tau, array(NA, dim = c(nregs_1, dims, K), dimnames = list(c("Intercept", colnames(xs)), colnames(ys), as.character(paste("K =", 1:K)))), simplify = F)


for (t in 1:n.tau) {
  for (k in 1:K) {
    for (b in 1:B) {
      betas.b[[k]][,,b] <-  model_boot[[t]][[b]]$betas[[k]] #insieme delle B matrici beta per ogni tau
      cpar.b[[k]][b,] <- model_boot[[t]][[b]]$cpar[[k]]
      df.b[[k]][b] <- model_boot[[t]][[b]]$df[[k]]
      sigmas.b[[k]][,,b] <- model_boot[[t]][[b]]$sigma[[k]]
    } #b loop
    betas.b.est[[t]][,,k] <- apply(betas.b[[k]], c(1,2), mean)
    betas.b.std[[t]][,,k] <- apply(betas.b[[k]], c(1,2), fBasics::stdev)
    cpar.b.est[[t]][,,k] <- colMeans(cpar.b[[k]])
    cpar.b.std[[t]][,,k] <- apply(cpar.b[[k]], 2, fBasics::stdev)
    df.b.est[[t]][k,] <- apply(df.b[[k]], 2, mean)
    df.b.std[[t]][k,] <- apply(df.b[[k]], 2, fBasics::stdev)
    sigmas.b.est[[t]][,,k] <- apply(sigmas.b[[k]], c(1,2), mean)
    sigmas.b.std[[t]][,,k] <- apply(sigmas.b[[k]], c(1,2), fBasics::stdev)
    test_statistic[[t]][,,k] <- abs(betas.est[[t]][,,k]/betas.b.std[[t]][,,k]) # > 1.96
    
  }
}


if(model == "expectile"){
  save.image(file = "realdata_cexp_2K_t.RData")
} else if(model == "quantile"){
  save.image(file = "realdata_cquant_2K_t.RData")
}
