
# fit the joint model of Zoop and whale with known sigmas
Rf_fitSZICOMP = function(
    y, X, M, A, t, th, loglam,
    new.run, n.iter, n.save, n.thin, n.update.basis,
    n.core, path.cpp, filename){
  
  p1 = ncol(X)
  p2 = ncol(X)
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
  
  beta1Ind = 1:p1
  beta2Ind = (p1+1):(p1+p2)
  alphaInd = p1+p2+1
  gammaInd = (p1+p2+1+1):(p1+p2+1+q)
  deltaInd = (p1+p2+1+q+1):(p1+p2+1+q+q)
  IgammaInd = (p1+p2+1+q+q+1):(p1+p2+1+q+q+q)
  IdeltaInd = (p1+p2+1+q+q+q+1):(p1+p2+1+q+q+q+q)
  kappaInd = p1+p2+1+q+q+q+q+1
  tauInd = p1+p2+1+q+q+q+q+1+1
  
  # hyperparameters
  shape_s = 0.001; rate_s = 1000
  phi = 0.1
  
  if(new.run){
    
    # initial parameter values
    beta1 = rnorm(p1)
    beta2 = rnorm(p2)
    alpha = rnorm(1)
    gamma = rnorm(q)
    delta = rnorm(q)
    Igamma = rep(1, q)
    Idelta = rep(1, q)
    kappa = 2
    tau = 2
    w = rep(0, N)
    w[y > 0] = 1
    
    # MH proposal covariance matrices
    sigma2 = rep(0.2^2, 5)
    COVbeta1 = diag(p1)
    COVbeta2 = diag(p2)
    COValpha = 1
    COVgamma = diag(q)
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
    
    ptm = proc.time()[3]
    dummy = szicompExBSnu(
      outeri, y, X, X, M, Qs, nt, w, beta1, beta2, alpha, 
      gamma, delta, Igamma, Idelta, kappa, tau, phi, shape_s, rate_s, 
      sigma2, COVbeta1, COVbeta2, COValpha, COVgamma, COVdelta, 
      updateCOV, adaptInterval, adaptFactorExponent, adapIter, 
      pre.n.iter, updateInterval, n.thin, th, loglam, 1:n.core) 
    rtime = rtime + proc.time()[3] - ptm  
    
    postSamples = rbind(postSamples, dummy$Sample)
    Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
    
    nSamples = nrow(dummy$Sample)
    beta1 = dummy$Sample[nSamples, beta1Ind]
    beta2 = dummy$Sample[nSamples, beta2Ind]
    alpha = dummy$Sample[nSamples, alphaInd]
    gamma = dummy$Sample[nSamples, gammaInd]
    delta = dummy$Sample[nSamples, deltaInd]
    Igamma = dummy$Sample[nSamples, IgammaInd]
    Idelta = dummy$Sample[nSamples, IdeltaInd]
    kappa = dummy$Sample[nSamples, kappaInd]
    tau = dummy$Sample[nSamples, tauInd]
    w = dummy$w
    
    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COVbeta1 = dummy$COVbeta1
    COVbeta2 = dummy$COVbeta2
    COValpha = dummy$COValpha
    COVgamma = dummy$COVgamma
    COVdelta = dummy$COVdelta
    
    save(postSamples, Accprob, rtime, pre.n.iter, n.thin,
         sigma2, COVbeta1, COVbeta2, COValpha, COVgamma, COVdelta,
         beta1, beta2, alpha, gamma, delta, Igamma, Idelta, kappa, tau, w,
         updateCOV, adaptInterval, adaptFactorExponent, adapIter, updateInterval,
         shape_s, rate_s, phi,
         file = filename)
  }
}

