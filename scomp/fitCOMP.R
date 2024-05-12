rm(list = ls())
require(akima); require(Rcpp); require(RcppArmadillo); require(sitmo)

path.r = 'src/RFtns.R'
path.cpp = 'src/RcppFtns.cpp'

# ==============================================================================
# fit the spatial mean-parameterized COMP regression model ----
# ==============================================================================

load('data/sim.RData')
load('appx/set.RData')

q = 100
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

Rf_fitSCOMP(
  y, X, offset = NULL, M, A, t, th, loglam,
  new.run = T, n.iter = 1000000, n.save = 2000,
  # new.run = T, n.iter = 200, n.save = 100, # for a test run
  n.thin = 20, n.update.basis = 200, n.core = 20,
  path.cpp, filename = 'fit/simCOMP.RData')



# =============================================================================-
# Summary ----
# =============================================================================-
require(coda); require(batchmeans); require(tidyverse)
bm_se = function(x){ bm(x)$se }
hpd = function(x){ paste0('(', round(HPDinterval(as.mcmc(x))[1], 2), ',', round(HPDinterval(as.mcmc(x))[2], 2), ')') }

load('data/sim.RData')
load('fit/simCOMP.RData')

p = ncol(X)
q = 100
M = M0[,1:q]

betaInd = 1:p
alphaInd = p+1
deltaInd = (p+1+1):(p+1+q)
IdeltaInd = (p+1+q+1):(p+1+q+q)
tauInd = p+1+q+q+1

burn = 10000
posterior = postSamples[-(1:burn),]
rtime/60/60

ind.pars = c(betaInd, alphaInd, tauInd)
name.pars = c(paste0('beta', 1:p), 'log(nu)', 'tau')

df.summary = data.frame(
  Mean = round(apply(posterior[,ind.pars], 2, mean), 2),
  Median = round(apply(posterior[,ind.pars], 2, median), 2),
  HPD = apply(posterior[,ind.pars], 2, hpd),
  MCSE = round(apply(posterior[,ind.pars], 2, bm_se), 4),
  ESS = round(apply(posterior[,ind.pars], 2, ess)),
  row.names = name.pars)

df.summary







