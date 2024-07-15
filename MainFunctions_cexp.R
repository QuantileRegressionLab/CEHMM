#Simulate from HMM
#############

nLL <- function(eta_par, u, post.hmm, dd, which_cop, t_df) {
  # llk <- c()
  # n = length(post.hmm)
  #t_df = df
  t.cop = tCopula(eta_par, dim = dd, dispstr = "un", df = t_df)
  #t.cop = tCopula(eta_par, dim = dd, dispstr = "un", df = t_df)
  llk = post.hmm * dCopula(u, t.cop, log = T)
  # for (t in 1:n) {
  #   llk[t] = post.hmm[t]*loglikCopula(param = eta_par, u[t,], norm.cop,
  #                                     error = c("-Inf"))
  # }
  llk[is.na(llk)] = 0
  llk[is.infinite(llk)] = 0
  #getTheta(norm.cop)
  sum_llk = sum(llk)
  return(-sum_llk)
}

nLL_skew <- function(eta_par, u, post.hmm, dd, which_cop, gamma_skewt, t_df) {
  
  scale_t = matrix(c(1, eta_par, eta_par, 1), 2, 2)
  
  pseudoskewt = sapply(1:dd, function(i) qst(p = u[,i], xi = 0, omega = 1, alpha = gamma_skewt[i], nu = exp(t_df)))
  
  # lnum = dskewt(pseudoskewt, scale = scale_t, gamma = gamma_skewt, df = t_df, log = T)
  lnum = dmst(pseudoskewt, xi = rep(0, dd), Omega = diag(dd), alpha = gamma_skewt, nu = exp(t_df), log = T)
  ldenom = rowSums(sapply(1:dd, function(i) dst(x = pseudoskewt[,i], xi = 0, omega = 1, alpha = gamma_skewt[i], nu = exp(t_df), log = T)))
  t.cop = lnum - ldenom
  
  # t.cop = dskewtcopula(u = u, scale = scale_t, gamma = gamma_skewt, df = t_df, log = T)
  # t.cop = tCopula(eta_par, dim = dd, dispstr = "un", df = t_df)
  #t.cop = tCopula(eta_par, dim = dd, dispstr = "un", df = t_df)
  llk = post.hmm * t.cop
  # for (t in 1:n) {
  #   llk[t] = post.hmm[t]*loglikCopula(param = eta_par, u[t,], norm.cop,
  #                                     error = c("-Inf"))
  # }
  llk[is.na(llk)] = 0
  llk[is.infinite(llk)] = 0
  #getTheta(norm.cop)
  sum_llk = sum(llk)
  return(-sum_llk)
}

nLL_bis <- function(eta_par, u, post.hmm, dd, which_cop) {
  # llk <- c()
  # n = length(post.hmm)
  #t_df = df
  if(which_cop == "norm"){
    norm.cop = normalCopula(eta_par, dim = dd, dispstr = "un")
    llk = post.hmm * dCopula(u, norm.cop, log = T)} else{
      t.cop = tCopula(eta_par[1:(dd*(dd-1)/2)], dim = dd, dispstr = "un", df = eta_par[-dd*(dd-1)/2])
      #t.cop = tCopula(eta_par, dim = dd, dispstr = "un", df = t_df)
      llk = post.hmm * dCopula(u, t.cop, log = T)
    }
  # for (t in 1:n) {
  #   llk[t] = post.hmm[t]*loglikCopula(param = eta_par, u[t,], norm.cop,
  #                                     error = c("-Inf"))
  # }
  nas_llk = sum(as.numeric(is.na(llk)))
  infs_llk = sum(as.numeric(is.infinite(llk)))
  #llk[is.na(llk)] = 0
  #llk[is.infinite(llk)] = 0
  #getTheta(norm.cop)
  sum_llk = sum(llk)
  values = list(nas = nas_llk, infs = infs_llk)
  return(values)
}
########
dAND = function(y, mu, sigma, p) {
  
  out = 2*sqrt(1-p)*sqrt(p)/(sigma*sqrt(pi)*(sqrt(p)+sqrt(1-p)))*
    exp(-abs(p-(y<=mu))*(y-mu)^2/sigma^2)
  return(out)
}

pAND = function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(q) == 0) stop("q must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  
  a = sqrt(1-p)/(sqrt(p)+sqrt(1-p))
  b = sqrt(p)/(sqrt(p)+sqrt(1-p))
  
  cdf_and =   ifelse(test = q<mu, yes = 2*b*pnorm(q,mu,sigma/sqrt(2*(1-p))), 
                     no = (b-a) + 2*a*pnorm(q,mu,sigma/sqrt(2*p)))  
  
  ifelse(test=lower.tail == TRUE,yes=return(cdf_and),no=return(1-(cdf_and)))
}


qAND = function(prob,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("prob must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in [0,1].")
  
  
  a = sqrt(1-p)/(sqrt(p)+sqrt(1-p))
  b = sqrt(p)/(sqrt(p)+sqrt(1-p))
  cdf_in_mode = b
  
  q = sapply(X=prob,FUN=function(prob,mu,sigma,p){ifelse(test=lower.tail == TRUE,yes=prob,no=(1-prob));
    ifelse(test=prob<cdf_in_mode,yes= qnorm(prob/(2*b), mu, sigma/sqrt(2*(1-p))),
           no= qnorm((prob+a-b)/(2*a), mu, sigma/sqrt(2*p)) )},mu=mu,sigma=sigma,p=p)
  
  
  return(q)
}

# and = AbscontDistribution(d = function(y, mu, sigma, p) 2*sqrt(1-p)*sqrt(p)/(sigma*sqrt(pi)*(sqrt(p)+sqrt(1-p)))*
#                             exp(-abs(p-(y<=mu))*(y-mu)^2/sigma^2), withStand = TRUE)
#############
# log-forward
l.forward = function(delta,gamma,m,f.dens){
  ns = nrow(f.dens)
  l.alpha = matrix(NA,ns,m)
  foo = delta*f.dens[1,]
  sumfoo = sum(foo)
  lscale = log(sumfoo)
  foo = foo/sumfoo
  l.alpha[1,] = log(foo)+lscale
  for(t in 2:ns)
  {
    foo = foo%*%gamma*f.dens[t,]
    sumfoo = sum(foo)
    lscale = lscale+log(sumfoo)
    foo = foo/sumfoo
    l.alpha[t,] = log(foo)+lscale
  }
  return(l.alpha)
}
#############
# log-backward
l.backward = function(delta,gamma,m,f.dens){
  ns = nrow(f.dens)
  l.beta = matrix(NA,ns,m)
  l.beta[ns,] = rep(0,m)
  foo = rep(1/m,m)
  lscale = log(m)
  for(t in (ns-1):1)
  {
    foo = gamma%*%(f.dens[t+1,]*foo)
    l.beta[t,] = log(foo)+lscale
    sumfoo = sum(foo)
    foo = foo/sumfoo
    lscale = lscale+log(sumfoo)
  }
  return(l.beta)
}
#############
# conditional M-quantile regression (WLS)
library(MASS)
QRLM = function(Y, X, weights, c = 1e200, maxit = 1e3, err = 1e-06, tau) {
  
  temp = lm.fit(x = X, y = Y, method = "qr")
  if(any(is.na(temp$coef))){
    temp$coef[which(is.na(temp$coef))] <- 0}
  coef.prev = temp$coef
  resid = temp$residuals
  done = F
  conv = NULL
  
  for(i in 1:maxit) {
    # scale = median(abs(resid)) / 0.6745
    scale = median(abs(resid - median(resid))) / 0.6745
    # scale = sd(resid)
    scale = 1
    w = psi.huber(resid / scale, k = c)
    w = rep(1, length(resid))
    ww = 2 * (1 - tau) * w
    ww[resid > 0] = 2 * tau * w[resid > 0]
    w = ww * weights
    temp <- lm.wfit(x = X, y = Y, w = w, method = "qr")
    if(any(is.na(temp$coef))){
      temp$coef[which(is.na(temp$coef))] <- 0}
    coef <- temp$coef
    resid <- temp$residuals
    
    done = (max(norm(coef - coef.prev, type = "2")) <= err)
    #print(paste("done:", done, ";coef:",coef))
    coef.prev = coef
    
    if (done) 
      break
  }
  if (!done) 
    warning(paste("rlm failed to converge in", maxit, "steps"))
  
  return(list(coefficients = coef, scale = scale, weights = w, fitted.values = temp$fitted.values, residuals = resid))
}
############################################ select simulation error distribution #######à
#Generate data from multivariate normal
expreg.hsmm.multi = function(ns, m, dd, delta, gamma,beta, Sigma, tau, d){
  # ns = sample size
  # dd = multivariate variable dimension
  # m = scalar number of states, 
  # Sigma = covariance matrix for errors
  # delta = vector of initial values for prior probabilities, (pi vector)
  # gamma = matrix of initial values for state transition probabilies. (Pi matrix)
  # mu = list of initial values for means, 
  # d = list of state duration probabilities.
  ld = sapply(d,length)
  mvect = 1:m
  state = numeric(ns)
  x = error = matrix(0, ns, dd)
  nreg = nrow(beta[[1]])
  state[1] = sample(mvect, 1, prob=delta)
  
  dur = sample(1:ld[state[1]],1,prob=d[[state[1]]])
  mu = replicate(m, array(NA, dim = c(ns,nreg,dd)), simplify = F)
  regressor = matrix(rep(rnorm(ns*(nreg-1))), ncol = (nreg-1))
  regressor = cbind(rep(1,ns), regressor) 
  
  for (j in 1:m) {
    mu[[j]] = regressor%*%beta[[j]]
  }
  
  for(t in 1:dur){
    state[t] = state[1]
    error[t,] = mvrnorm(1, mu = rep(0, dd), Sigma = Sigma[[state[t]]])
    x[t,] = mu[[state[t]]][t,] + error[t,]
  }
  
  total = dur
  
  while(total<ns){
    
    state[total+1] = sample(mvect,1,prob=gamma[state[total],])
    dur = sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(t in 1:dur){
      if(total+t>ns) break
      state[total+t] = state[total+1]
      error[total+t,] = mvrnorm(1, mu = rep(0, dd), Sigma = Sigma[[state[total+t]]])
      x[total+t,] = mu[[state[total+t]]][total+t,] + error[total+t,]
    }
    total = total + dur
  }
  
  return(list(series=x,state=state,regressor=regressor,error=error))
}


#Generate data from multivariate t Student and skew-t
expreg.hsmm.multi.skewt = function(ns, m, dd, delta, gamma,beta, Sigma, tau,df,gamma_skewt, d){
  # ns = sample size
  # dd = multivariate variable dimension
  # m = scalar number of states, 
  # Sigma = covariance matrix for errors
  # delta = vector of initial values for prior probabilities, (pi vector)
  # gamma = matrix of initial values for state transition probabilies. (Pi matrix)
  # mu = list of initial values for means, 
  # df = list of t degrees of freedom for each state,
  # gamma_skewt = list of t skewness parameter for each state,
  # d = list of state duration probabilities.
  gamma_skewt = unlist(gamma_skewt)
  
  ld = sapply(d,length)
  mvect = 1:m
  state = numeric(ns)
  x = error = matrix(0, ns, dd)
  nreg = nrow(beta[[1]])
  state[1] = sample(mvect, 1, prob=delta)
  
  dur = sample(1:ld[state[1]],1,prob=d[[state[1]]])
  mu = replicate(m, array(NA, dim = c(ns,nreg,dd)), simplify = F)
  regressor = matrix(rep(rnorm(ns*(nreg-1))), ncol = (nreg-1))
  regressor = cbind(rep(1,ns), regressor) 
  
  for (j in 1:m) {
    mu[[j]] = regressor%*%beta[[j]]
  }
  
  for(t in 1:dur){
    state[t] = state[1]
    error[t,] = rmst(1, xi = rep(0, dd), Omega = Sigma[[state[t]]], nu = df[[state[t]]], alpha = gamma_skewt)
    x[t,] = mu[[state[t]]][t,] + error[t,]
  }
  
  total = dur
  
  while(total<ns){
    
    state[total+1] = sample(mvect,1,prob=gamma[state[total],])
    dur = sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(t in 1:dur){
      if(total+t>ns) break
      state[total+t] = state[total+1]
      error[total+t,] = rmst(1, xi = rep(0, dd), Omega = Sigma[[state[total+t]]], nu = df[[state[total+t]]], alpha = gamma_skewt)
      x[total+t,] = mu[[state[total+t]]][total+t,] + error[total+t,]
       }
    total = total + dur
  }
  
  return(list(series=x,state=state,regressor=regressor,error=error))
}



# Generate data from a gaussian or t copula with AND marginals
hsmm.multi.real = function(reg, ns, m, dd, delta, gamma, beta, df, Sigma, sigma, tau, d, wcop, cpar){
  # ns = sample size
  # dd = multivariate variable dimension
  # m = scalar number of states, 
  # Sigma = covariance matrix for errors
  # delta = vector of initial values for prior probabilities, (pi vector)
  # gamma = matrix of initial values for state transition probabilies. (Pi matrix)
  # mu = list of initial values for means, 
  # d = list of state duration probabilities.
  ld = sapply(d,length)
  mvect = 1:m
  state = numeric(ns)
  x = error = matrix(0, ns, dd) 
  nreg = nrow(beta[[1]])
  state[1] = sample(mvect, 1, prob=delta)
  
  dur = sample(1:ld[state[1]],1,prob=d[[state[1]]])
  mu = replicate(m, array(NA, dim = c(ns,nreg,dd)), simplify = F)
  
  for (j in 1:m) {
    mu[[j]] = reg%*%beta[[j]]
  }
  
  for(t in 1:dur){
    state[t] = state[1]
    if(wcop == "norm"){
      mn = mvrnorm(1, mu = rep(0, dd), Sigma = Sigma[[state[t]]])
      pn = pnorm(mn)
      error[t,] = sapply(1:dd, function(i) qAND(pn[i], mu = 0, sigma = sigma[[state[t]]][i], p = tau))
    }else{
      t.cop = tCopula(param = cpar[[state[t]]], dim = dd, dispstr = "un", df = df[[state[t]]])
      ut = rCopula(1, t.cop)
      error[t,] = sapply(1:dd, function(i) qAND(ut[i], mu = 0, sigma = sigma[[state[t]]][i], p = tau))
    }

    x[t,] = mu[[state[t]]][t,] + error[t,]
  }
  
  total = dur
  
  while(total<ns){
    
    state[total+1] = sample(mvect,1,prob=gamma[state[total],])
    dur = sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(t in 1:dur){
      if(total+t>ns) break
      state[total+t] = state[total+1]
      if(wcop == "norm"){
        mn = mvrnorm(1, mu = rep(0, dd), Sigma = Sigma[[state[total+t]]])
        pn = pnorm(mn)
        error[total+t,] = sapply(1:dd, function(i) qAND(pn[i], mu = 0, sigma = sigma[[state[total+t]]][i], p = tau))
      }else{
        t.cop = tCopula(param = cpar[[state[total+t]]], dim = dd, dispstr = "un", df = df[[state[total+t]]])
        ut = rCopula(1, t.cop)
        error[total+t,] = sapply(1:dd, function(i) qAND(ut[i], mu = 0, sigma = sigma[[state[total+t]]][i], p = tau))
      }
      x[total+t,] = mu[[state[total+t]]][total+t,] + error[total+t,]
    }
    total = total + dur
  }
  
  return(list(series=x,state=state,regressor=reg,error=error))
}

###########################################
# EM for expectile HMM regression
###########################################
em.hmm.cqereg = function(y, X_, tau, m, dd, delta, gamma, beta, cpar, sigma, maxiter, df_cop, which_cop, tol, trace=F){
  # y = response variable
  # X = regressors
  # tau = quantile level
  # m = number of states
  # dd = dimension of the response variable
  # delta = initial values for prior probabilities
  # gamma = initial values for state transition probabilities
  # beta = initial values for regression coefficients
  # cpar = initial values for copula parameters
  # sigma = initial values for error standard deviations
  # maxiter = maximum number of iterations
  # df_cop = initial values for copula degrees of freedom
  # which_cop = copula type
  # tol = tolerance level
  # trace = print output
  
  start_time = Sys.time()
  X = X_
  n = dim(y)[1]
  C = dim(X)[2]
  
  if(is.list(delta)) delta = unlist(delta)

 
  startpars = list(delta = delta, beta = beta, gamma = gamma, sigma = sigma, cpar = cpar, df = df_cop)
  llk.pred = 0
  
  gamma.next = gamma
  beta.next = beta
  delta.next = delta
  sigma.next = sigma
  cpar.next = cpar
  df.next = df_cop
  
  for (iter in 1:maxiter){
    print(paste("em iter=",iter))
    lallprobs_multi = pdf_c = marg_prod = matrix(NA, nrow = n,ncol = m)
    lallprobs = u = replicate(m, matrix(NA, nrow = n,ncol = dd), simplify = F)
    for (j in 1:m){
      for (t in 1:n){
        for (index.d in 1:dd) {
          lallprobs[[j]][t,index.d] = dAND(y[t,index.d], mu = X[t,]%*%beta[[j]][,index.d],sigma=sigma[[j]][index.d],p=tau)
          #marginals for each dimension
          u[[j]][t,index.d] = pAND(y[t,index.d], mu = X[t,]%*%beta[[j]][,index.d],sigma=sigma[[j]][index.d],p=tau)
        }
      }
    }  

    # Compute the copula density
    ##Copula ##
    norm.cop <- t.cop <- list()

    for (j in 1:m){
      u[[j]] = u[[j]] * (n - 1)/n + .Machine$double.eps/n 
      coef_ <- cpar[[j]][1:(dd*(dd-1)/2)]
      if(which_cop == "indip"){ pdf_c[,j] <- 1
      } else if(which_cop == "norm"){
        norm.cop[[j]] <- normalCopula(coef_, dim = dd, dispstr = "un")
        pdf_c[,j] <- dCopula(u[[j]], norm.cop[[j]])} else {
          df_ <- df_cop[[j]]
          t.cop[[j]] <- tCopula(coef_, dim = dd, dispstr = "un", df = df_)
          pdf_c[,j] <- dCopula(u[[j]], t.cop[[j]])
          }
      marg_prod[,j] <- apply(lallprobs[[j]], 1, prod)
    }
    lallprobs_multi <- marg_prod * pdf_c
    
 # Calculate forward and backward probabilities
    #forward probs
    la = matrix(NA , m ,n)
    pa = matrix(NA , n ,m)
    foo_s = foo_t = matrix(NA , n ,m)
    foo = delta * lallprobs_multi[1,]
    if(sum(is.na(foo))>0){foo[is.na(foo)] <- .Machine$double.xmin}
    if(sum(foo==0)>0 ){foo[foo==0] <- .Machine$double.xmin}
    sumfoo = sum(foo)
    lscale = (log(sumfoo))
    foo = foo / sumfoo
    pa[1,] = foo
    la[,1] = lscale + log(foo)
    foo_s[1,] <- foo
    for (t in 2:n){
      foo_s[t,] <- foo
      foo = foo %*% gamma * lallprobs_multi[t,]
      if(sum(is.na(foo))>0){foo[is.na(foo)] <- .Machine$double.xmin}
      if(sum(foo==0)>0){foo[foo==0] <- .Machine$double.xmin}
      sumfoo = sum(foo)
      lscale = lscale + (log(sumfoo))
      
      
      ##***
      if(sum(as.integer(is.nan(log(sumfoo))))>1){print(paste("log sumfoo NaN at t=",t))}
      #***
      foo = foo / sumfoo
      pa[t,] = foo
      la[,t] = log(foo) + lscale
    }
    
    #backward probs
    lb = matrix(NA , m ,n)
    lscale_s <- c()
    foo = rep(1/m, m)
    foo_t[n,] <- foo
    lb[,n] = rep(0,m)
    lscale = log(m)
    lscale_s[1] <- lscale
    for (t in (n-1):1){ 
      foo_t[t,] <- foo
      foo = gamma %*% (lallprobs_multi[t+1,] * foo)
      if(sum(is.na(foo))>0){foo[is.na(foo)] <- .Machine$double.xmin}
      if(sum(foo==0)>0){foo[foo==0] <- .Machine$double.xmin}
      lb[,t] = log((foo)) + lscale
      sumfoo = sum(foo)
      foo = foo / sumfoo
      lscale = lscale + (log(sumfoo))
      if(is.infinite(lscale) || is.nan(lscale)){lscale <- lscale_s[t+1]}
      lscale_s[t] <- lscale
    }
    
    c   =  max(la[,n])                                      
    llk = c+log(sum(exp(la[,n]-c)))
    if(is.nan(llk)){
      print(paste("llk=",llk)); llk = c(.Machine$double.eps)}
    print(paste("llk=",llk,";c=",c))
    post.hmm = matrix(NA , ncol = n, nrow = m)
    for (j in 1:m){                                           
      post.hmm[j,] = exp(la[j,]+lb[j,]-llk)
    }
    
    
    #Delta estimate
    delta.next = exp(la[,1]+lb[,1]-llk)
    delta.next = delta.next/sum(delta.next)
    
    #Gamma estimate
    gamma.next = matrix(NA,nrow = m, ncol = m)
    for (j in 1:m){
      for (k in 1:m){                                                      
        tmp = exp(la[j,1:(n-1)]+   
                    log(lallprobs_multi[2:n,k])+lb[k,2:n]-llk)
        tmp[tmp==Inf] <- 1
        tmp[tmp==-Inf] <- 0
        gamma.next[j,k] = gamma[j,k]*sum(tmp)  
      }   
    }                                                       
    gamma.next = gamma.next/apply(gamma.next,1,sum)
    
    
    #Beta and Sigma estimate
    aux = replicate(m, matrix(data = NA, n, dd), simplify = F)
    beta.next = replicate(m,matrix(data = NA, C, dd), simplify = F)
    sigma.next = replicate(m, matrix(NA,nrow = dd, ncol = 1), simplify = F)
    gamma_skewt.next = replicate(m, matrix(NA, dd, 1), simplify = F)
    if(which_cop=="norm"){           #copula correlation coefficients
      cpar.next = replicate(m, matrix(NA, nrow = 1,ncol = dd*(dd-1)/2), simplify = F)}else{
        cpar.next = replicate(m, matrix(NA, nrow = 1,ncol = dd*(dd-1)/2 + 0), simplify = F) 
        df.next = list()
      } 
    model = list()
    for (j in 1:m){
      if(sum(is.infinite(post.hmm)) >0){
        post.hmm[post.hmm==Inf]<-1
        post.hmm[post.hmm== -Inf]<-0
      }
      ### ****************
      for (index.d in 1:dd) {
        model[[index.d]] = QRLM(Y = y[,index.d], X = X, weights = post.hmm[j,], tau=tau)
        beta.next[[j]][,index.d] = model[[index.d]]$coefficients
        aux[[j]][,index.d] = post.hmm[j,]*abs(tau-as.numeric(model[[index.d]]$residuals<=0))*model[[index.d]]$residuals^2
        sigma.next[[j]][index.d] = sqrt(2*sum(aux[[j]][,index.d])/sum(post.hmm[j,]))
      }
    }   
    
    nas_llk <- infs_llk <- c()
    for (j in 1:m){   
      if(which_cop == "indip"){ df.next[[j]] <- 1; cpar.next[[j]] <- rep(0, (dd*(dd-1)/2))
      } else if(which_cop == "norm"){
        pseudon = sapply(1:dd, function(i) qnorm(p = u[[j]][,i], mean = 0, sd = 1))
        cpar.next[[j]] = P2p(cov2cor(t(pseudon)%*%diag(post.hmm[j,])%*%pseudon))
        df.next[[j]] <- 1
      } else {
        pseudot = sapply(1:dd, function(i) qt(p = u[[j]][,i], df = df_cop[[j]]))
        scale_t = p2P(cpar[[j]][1:(dd*(dd-1)/2)], dd)
        wt = sapply(1:n, function(i) (df_cop[[j]] + dd) / (df_cop[[j]] + pseudot[i,]%*%solve(scale_t)%*%pseudot[i,]))
        cpar.next[[j]] = P2p(cov2cor(t(pseudot)%*%diag(wt * post.hmm[j,])%*%pseudot))
        
        fit <- optim(df_cop[[j]], fn = nLL, u=u[[j]],post.hmm=post.hmm[j,], dd=dd, eta_par=cpar[[j]],
                     which_cop = which_cop, method = "L-BFGS-B", lower = 2.0001, upper = 1e3,
                     control = list(factr = 1e7, maxit = 5e3)) # 1e7 * 1e-15 = 1e-8 by default. Smaller -> stricter convergence
        df.next[[j]] <- fit$par
      }
    }
   
    
    crit_beta <- c()
    crit_sigma <- c()
    crit_cpar <- c()
    crit_df <- c()
    #Set Convergence Criterion
  
    for (j in 1:m) {
      crit_beta[j] <- sum(abs(beta[[j]]-beta.next[[j]]))
      crit_sigma[j] <- sum(abs(sigma[[j]]-sigma.next[[j]]))
      crit_cpar[j] <- sum(abs(cpar[[j]][1:(dd*(dd-1)/2)]-cpar.next[[j]][1:(dd*(dd-1)/2)]))
      crit_df[j] <- sum(abs(df_cop[[j]]-df.next[[j]]))
    }
    crit = sum(sum(abs(delta-delta.next)) + sum(abs(crit_beta)) + sum(abs(crit_cpar)) + sum(abs(crit_df)) + sum(abs(gamma-gamma.next)) + sum(abs(crit_sigma))) # + sum(abs(llk-llk.pred))
    dif_ = abs(llk - llk.pred)
    if(trace) {cat("iteration ",iter,": loglik = ", llk, "\n") 
      cat("time of a single iteration: ",time_iter, "\n")}
    if(crit<tol){
      
      #Output
      return(list(m=m, dif=dif_, betas=beta, delta=delta, gamma=gamma, cpar = cpar, df=df_cop, nas_llk = nas_llk, infs_llk = infs_llk,
                  loglik=llk, iteration=iter, sigma=sigma, post = post.hmm, aux = aux, u = u, crit = crit1))
    }
    delta = delta.next
    beta = beta.next
    gamma = gamma.next
    sigma = sigma.next
    cpar = cpar.next
    df_cop = df.next
    llk.pred = llk
    #info.criteria
    if(which_cop == "indip"){
      n.par = dd * C * m + dd * m + m * (m - 1) + (m - 1)
    } else if(which_cop == "norm"){
      n.par = dd * C * m + dd * m + m * (m - 1) + (m - 1) + m*(dd*(dd-1)/2)
    } else {
      n.par = dd * C * m + dd * m + m * (m - 1) + (m - 1) + m*(dd*(dd-1)/2) + m
    }
    aic.crit = -2*llk + 2*n.par
    bic.crit = -2*llk + log(n)*n.par
    icl.crit = bic.crit - 2 * sum(post.hmm * ifelse(post.hmm > 0, log(post.hmm), 0))
    crit1 = c(aic.crit, bic.crit, icl.crit)
    names(crit1) = c("AIC", "BIC", "ICL")
  }                                                        
  cat(paste("\nNo convergence after",maxiter,"iterations\n\n"))  
  
  #Output
  return(list(m=m, dif=dif_, betas=beta, delta=delta, gamma=gamma, cpar = cpar, df=df_cop, nas_llk = nas_llk, infs_llk = infs_llk,
              loglik=llk, iteration=iter, sigma=sigma, post = post.hmm, aux = aux, u = u, crit = crit1))
}





