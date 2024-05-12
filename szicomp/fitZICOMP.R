rm(list = ls())
require(Rcpp); require(RcppArmadillo); require(sitmo)

path.r = 'src/RFtns.R'
path.cpp = 'src/RcppFtns.cpp'

dataname = 'under'

# ==============================================================================
# fit the spatial mean-parameterized COMP regression model ----
# ==============================================================================

load(paste0('data/', dataname, '.RData'))
load('appx/set.RData')

q = 50
M = M0[,1:q]

source(path.r)

# y: response
# X: design matrix
# M: basis matrix
# A: adjacent matrix
# t: the number of time points
# th: particle matrix for spline approximation
# loglam: log of rate associated with the particle (th)
# new.run: whether is is a new run of MCMC
# n.iter: the number of iterations
# n.thin: thinning period
# n.save: save results every n.save'th iteration
# n.update.basis: update basis-vector selection every n.update.basis'th iteration
# n.core: the number of cores used for the exchange algorithm
# path.cpp: path to the c++ source file
# filename: name for the results

Rf_fitSZICOMP(
  y, X, M, A, t, th, loglam,
  new.run = T, n.iter = 1000000, n.save = 2000,
  # new.run = T, n.iter = 200, n.save = 100, # for a test run
  n.thin = 20, n.update.basis = 200, n.core = 20,
  path.cpp, filename = paste0('fit/', dataname, 'ZICOMP.RData'))


# new.run = T; n.iter = 200; n.save = 100; n.thin = 20; n.update.basis = 200; n.core = 2
# filename = paste0('fit/', dataname, 'ZICOMP.RData')



# =============================================================================-
# Summary ----
# =============================================================================-
require(coda); require(batchmeans); require(tidyverse)
bm_se = function(x){ bm(x)$se }
hpd = function(x){ paste0('(', round(HPDinterval(as.mcmc(x))[1], 2), ',', round(HPDinterval(as.mcmc(x))[2], 2), ')') }

load(paste0('data/', dataname, '.RData'))
load(paste0('fit/', dataname, 'ZICOMP.RData'))

p1 = p2 = ncol(X)
q = 50
M = M0[,1:q]

beta1Ind = 1:p1
beta2Ind = (p1+1):(p1+p2)
alphaInd = p1+p2+1
gammaInd = (p1+p2+1+1):(p1+p2+1+q)
deltaInd = (p1+p2+1+q+1):(p1+p2+1+q+q)
IgammaInd = (p1+p2+1+q+q+1):(p1+p2+1+q+q+q)
IdeltaInd = (p1+p2+1+q+q+q+1):(p1+p2+1+q+q+q+q)
kappaInd = p1+p2+1+q+q+q+q+1
tauInd = p1+p2+1+q+q+q+q+1+1

burn = 10000
posterior = postSamples[-(1:burn),]
rtime/60/60

ind.pars = c(beta1Ind, beta2Ind, alphaInd, kappaInd, tauInd)
name.pars = c(paste0('beta1_', 1:p1), paste0('beta2_', 1:p1), 'log(nu)', 'kappa', 'tau')

df.summary = data.frame(
  Mean = round(apply(posterior[,ind.pars], 2, mean), 2),
  Median = round(apply(posterior[,ind.pars], 2, median), 2),
  HPD = apply(posterior[,ind.pars], 2, hpd),
  MCSE = round(apply(posterior[,ind.pars], 2, bm_se), 4),
  ESS = round(apply(posterior[,ind.pars], 2, ess)),
  row.names = name.pars)

df.summary







