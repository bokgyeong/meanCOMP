
# fit the joint model of Zoop and whale with known sigmas
Rf_fitSCOMP = function(
    y, X, offset = NULL, M, A, t, th, loglam,
    new.run, n.iter, n.save, n.thin, n.update.basis,
    n.core, path.cpp, filename){
  
  p = ncol(X)
  q = ncol(M)
  n = rep(nrow(A), t)
  nt = c(0, cumsum(n))
  Q = diag(rowSums(A, n[1])) - A
  Qs = t(M) %*% Q %*% M
  
  # ===========================================================================-
  # MCMC ----
  # ===========================================================================-
  
  outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
  
  updateCOV = TRUE
  adaptInterval = 200
  adaptFactorExponent = 0.8 # tuning parameter for the adaptive MCMC
  updateInterval = n.update.basis
  
  betaInd = 1:p
  alphaInd = p+1
  deltaInd = (p+1+1):(p+1+q)
  IdeltaInd = (p+1+q+1):(p+1+q+q)
  tauInd = p+1+q+q+1
  
  # hyperparameters
  shape_s = 0.001; rate_s = 1000
  phi = 0.1
  
  if(new.run){
    
    # initial parameter values
    beta = rnorm(p, 0, 0.1)
    alpha = 0
    delta = rnorm(q, 0, 0.1)
    Idelta = rep(1, q)
    tau = 2
    
    # MH proposal covariance matrices
    sigma2 = rep(0.2^2, 3)
    COVbeta = diag(p)
    COValpha = 1
    COVdelta = diag(q)
    adapIter = 1
    
    start = 1; postSamples = c(); Accprob = 0; rtime = 0
    
  } else {
    load(filename)
    start = which(outers/n.thin == nrow(postSamples))
  }
  
  sourceCpp(path.cpp)
  
  for(i in (start+1):length(outers) ){
    outeri = outers[i]-outers[i-1]
    pre.n.iter = outers[i-1]
    
    if(is.null(offset)){
      ptm = proc.time()[3]
      dummy = scompExBS(
        outeri, y, X, M, Qs, nt, beta, alpha, delta, Idelta, 
        tau, phi, shape_s, rate_s, sigma2, COVbeta, COValpha, COVdelta, 
        updateCOV, adaptInterval, adaptFactorExponent, adapIter, 
        pre.n.iter, updateInterval, n.thin, th, loglam, 1:n.core) 
      rtime = rtime + proc.time()[3] - ptm  
      
    } else {
      ptm = proc.time()[3]
      dummy = scompExBSoffset(
        outeri, y, X, offset, M, Qs, nt, beta, alpha, delta, Idelta, 
        tau, phi, shape_s, rate_s, sigma2, COVbeta, COValpha, COVdelta, 
        updateCOV, adaptInterval, adaptFactorExponent, adapIter, 
        pre.n.iter, updateInterval, n.thin, th, loglam, 1:n.core)
      rtime = rtime + proc.time()[3] - ptm
    }
    
    postSamples = rbind(postSamples, dummy$Sample)
    Accprob = ( Accprob * pre.n.iter + colSums(dummy$Accprob) ) / outers[i]
    
    nSamples = nrow(dummy$Sample)
    beta = dummy$Sample[nSamples, betaInd]
    alpha = dummy$Sample[nSamples, alphaInd]
    delta = dummy$Sample[nSamples, deltaInd]
    Idelta = dummy$Sample[nSamples, IdeltaInd]
    tau = dummy$Sample[nSamples, tauInd]
    
    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COVbeta = dummy$COVbeta
    COValpha = dummy$COValpha
    COVdelta = dummy$COVdelta
    
    save(postSamples, rtime, Accprob, pre.n.iter, n.thin,
         sigma2, COVbeta, COValpha, COVdelta, 
         beta, alpha, delta, Idelta, tau,
         updateCOV, adaptInterval, adaptFactorExponent, adapIter, updateInterval,
         shape_s, rate_s, phi,
         file = filename)
  }
}

